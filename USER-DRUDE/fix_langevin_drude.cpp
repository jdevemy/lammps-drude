/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_langevin_drude.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixLangevinDrude::FixLangevinDrude(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 9) error->all(FLERR,"Illegal fix langevin/drude command");
  // TODO add options for tally and zero

  // Langevin thermostat should be applied every step
  nevery = 1;
  vector_flag = 1;
  global_freq = nevery;
  extvector = 0;
  size_vector = 6;
  //extscalar = 1;
  
  // core temperature
  tstr_core = NULL;
  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    tstr_core = new char[n];
    strcpy(tstr_core,&arg[3][2]);
    tstyle_core = EQUAL;
  } else {
    t_start_core = force->numeric(FLERR,arg[3]);
    t_target_core = t_start_core;
    tstyle_core = CONSTANT;
  }
  t_period_core = force->numeric(FLERR,arg[4]);
  int seed_core = force->inumeric(FLERR,arg[5]);  

  // drude temperature
  tstr_drude = NULL;
  if (strstr(arg[7],"v_") == arg[6]) {
    int n = strlen(&arg[6][2]) + 1;
    tstr_drude = new char[n];
    strcpy(tstr_drude,&arg[6][2]);
    tstyle_drude = EQUAL;
  } else {
    t_start_drude = force->numeric(FLERR,arg[6]);
    t_target_drude = t_start_drude;
    tstyle_drude = CONSTANT;
  }
  t_period_drude = force->numeric(FLERR,arg[7]);
  int seed_drude = force->inumeric(FLERR,arg[8]);
  
  // error checks
  if (t_period_core <= 0.0)
    error->all(FLERR,"Fix langevin/drude period must be > 0.0");
  if (seed_core  <= 0) error->all(FLERR,"Illegal langevin/drude seed");
  if (t_period_drude <= 0.0)
    error->all(FLERR,"Fix langevin/drude period must be > 0.0");
  if (seed_drude <= 0) error->all(FLERR,"Illegal langevin/drude seed"); 

  random_core  = new RanMars(lmp,seed_core);
  random_drude = new RanMars(lmp,seed_drude);
  
  int iarg = 9;
  zero = 0;
  while (iarg < narg) { 
    if (strcmp(arg[iarg],"zero") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix langevin/drude command");
      if (strcmp(arg[iarg+1],"no") == 0) zero = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) zero = 1;
      else error->all(FLERR,"Illegal fix langevin/drude command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix langevin/drude command");
  }

  tflag = 0; // no external compute/temp is specified yet (for bias)
  energy = 0.;
}

/* ---------------------------------------------------------------------- */

FixLangevinDrude::~FixLangevinDrude()
{
  delete random_core;
  delete [] tstr_core;
  delete random_drude;
  delete [] tstr_drude;
}

/* ---------------------------------------------------------------------- */

int FixLangevinDrude::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLangevinDrude::init()
{
  // check variable-style target core temperature
  if (tstr_core) {
    tvar_core = input->variable->find(tstr_core);
    if (tvar_core < 0)
      error->all(FLERR,"Variable name for fix langevin/drude does not exist");
    if (input->variable->equalstyle(tvar_core)) tstyle_core = EQUAL;
    else error->all(FLERR,"Variable for fix langevin/drude is invalid style");
  }

  // check variable-style target drude temperature
  if (tstr_drude) {
    tvar_drude = input->variable->find(tstr_drude);
    if (tvar_drude < 0)
      error->all(FLERR,"Variable name for fix langevin/drude does not exist");
    if (input->variable->equalstyle(tvar_drude)) tstyle_drude = EQUAL;
    else error->all(FLERR,"Variable for fix langevin/drude is invalid style");
  }

  // TODO: allow bias
  //if (temperature->tempbias) which = BIAS; // only for core
  //else which = NOBIAS;

  char typetag[] = "drudetype", idtag[] = "drudeid";
  int dummy;
  index_drudetype = atom->find_custom(typetag, dummy);
  if (index_drudetype == -1) 
    error->all(FLERR,"Unable to get DRUDETYPE atom property");
  index_drudeid = atom->find_custom(idtag, dummy);
  if (index_drudeid == -1) 
    error->all(FLERR,"Unable to get DRUDEID atom property");
}

/* ---------------------------------------------------------------------- */

void FixLangevinDrude::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    error->all(FLERR,"RESPA style not compatible with fix langevin/drude");
  }

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int fix_dof = 0;
  int *drudetype = atom->ivector[index_drudetype];
  //int *drudeid = atom->ivector[index_drudeid];
  int dim = domain->dimension;

  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  int dof_core_loc = 0, dof_drude_loc = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) { // Only the cores need to be in the group.
      if (drudetype[i] == 0) // Non-polarizable atom
          dof_core_loc++;
      else {
          if (drudetype[i] == 2) continue;
          dof_core_loc++;
          dof_drude_loc++;
      }
    }
  }
  int ncoreloc = dof_core_loc;
  dof_core_loc *= dim;
  dof_drude_loc *= dim;
  MPI_Allreduce(&dof_core_loc,  &dof_core,  1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&dof_drude_loc, &dof_drude, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&ncoreloc, &ncore, 1, MPI_INT, MPI_SUM, world);
  dof_core -= fix_dof;
  if (zero) dof_core -= dim; // The center of mass is not thermalized.

  langevin(vflag, false);
}

/* ---------------------------------------------------------------------- */

void FixLangevinDrude::post_force(int vflag)
{
  langevin(vflag, true);
}

/* ---------------------------------------------------------------------- */

void FixLangevinDrude::langevin(int /*vflag*/, bool thermalize=true)
{ 
  // Compute the kinetic energy and temperature of the reduced degrees of
  // freedom. Thermalize by adding the langevin force if thermalize=true.
  // Each core-Drude pair is thermalized only once: where the core is local.

  double **v = atom->v, **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal, nall = atom->nlocal + atom->nghost;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double ftm2v = force->ftm2v, mvv2e = force->mvv2e;
  double kb = force->boltz, dt = update->dt;
  
  int *drudetype = atom->ivector[index_drudetype];
  int *drudeid = atom->ivector[index_drudeid];
  double vdrude[3], vcore[3]; // velocities in reduced representation
  double fdrude[3], fcore[3]; // forces in reduced representation
  double Ccore, Cdrude, Gcore, Gdrude;
  double fcoresum[3], fcoreloc[3];
  int dim = domain->dimension;

  /*for (int i = 0; i < nlocal; i++){ // Check comm_style brick/drude
    if (drudetype[i] == 2) {
        int j = atom->map(drudeid[i]);
        if (j < 0) std::cout << "CORE n°" << drudeid[i] << " missing on proc " << 
        comm->me << " where is its DRUDE (" << atom->tag[i] << ").\n";
        else if (j >= nlocal) std::cout << "CORE n°" << drudeid[i] << " is a ghost on proc " << 
        comm->me << " whereas its DRUDE (" << atom->tag[i] << ") is owned.\n";
    } else if (drudeid[i]) {
        int j = atom->map(drudeid[i]);
        if (j < 0) std::cout << "DRUDE n°" << drudeid[i] << " missing on proc " << 
        comm->me << " where is its CORE (" << atom->tag[i] << ").\n";
        else if (j >= nlocal) std::cout << "DRUDE n°" << drudeid[i] << " is a ghost on proc " << 
        comm->me << " whereas its CORE (" << atom->tag[i] << ") is owned.\n";
    }
  }*/
  
  // Compute target core temperature
  if (thermalize) {
    if (tstyle_core == CONSTANT)
      t_target_core = t_start_core; // + delta * (t_stop-t_start_core);
    else {
      modify->clearstep_compute();
      t_target_core = input->variable->compute_equal(tvar_core);
      if (t_target_core < 0.0)
        error->one(FLERR, "Fix langevin/drude variable returned "
                          "negative core temperature");
      modify->addstep_compute(update->ntimestep + nevery);
    }

    // Compute target drude temperature
    if (tstyle_drude == CONSTANT)
      t_target_drude = t_start_drude; // + delta * (t_stop-t_start_core);
    else {
      modify->clearstep_compute();
      t_target_drude = input->variable->compute_equal(tvar_drude);
      if (t_target_drude < 0.0)
        error->one(FLERR, "Fix langevin/drude variable returned "
                          "negative drude temperature");
      modify->addstep_compute(update->ntimestep + nevery);
    }
 
    // Clear ghost forces
    // They have already been communicated if needed
    for (int i = nlocal; i < nall; i++) {
      for (int k = 0; k < dim; k++)
        f[i][k] = 0.;
    }
    if (zero) for (int k=0; k<dim; k++) fcoreloc[k] = 0.;
  }

  double kineng_core_loc = 0., kineng_drude_loc = 0.;
  // NB : the masses are the real masses, not the reduced ones
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) { // only the cores need to be in the group
      if (drudetype[i] == 0) { // Non-polarizable atom
        double mi;
        if (rmass)
          mi = rmass[i];
        else
          mi = mass[type[i]]; 
        if (thermalize) {
          Gcore  = mi / t_period_core  / ftm2v;
          Ccore  = sqrt(2.0 * Gcore  * kb * t_target_core  / dt / ftm2v / mvv2e);
          for(int k = 0; k < dim; k++){
            fcore[k] = Ccore  * random_core->gaussian()  - Gcore  * v[i][k];
            if (zero) fcoreloc[k] += fcore[k];
            f[i][k] += fcore[k];
          }
        }
        kineng_core_loc += mi * (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      } else {
        if (drudetype[i] == 2) continue; // Done together with the core
        
        int j = atom->map(drudeid[i]);
        double mi, mj, mtot, mu; // i is core, j is drude
        if (rmass) {
          mi = rmass[i];
          mj = rmass[j];
        } else {
          mi = mass[type[i]]; 
          mj = mass[type[j]];
        } 
        mtot = mi + mj;
        mu = mi * mj / mtot;
        mi /= mtot;
        mj /= mtot;
          
        if (thermalize) {
          Gcore  = mtot / t_period_core  / ftm2v;
          Gdrude = mu   / t_period_drude / ftm2v;
          Ccore  = sqrt(2.0 * Gcore  * kb * t_target_core  / dt / ftm2v / mvv2e);
          Cdrude = sqrt(2.0 * Gdrude * kb * t_target_drude / dt / ftm2v / mvv2e);
        }

        for (int k=0; k<dim; k++) {
          // TODO check whether a fix_modify temp can subtract a bias velocity
          vcore[k] = mi * v[i][k] + mj * v[j][k]; 
          vdrude[k] = v[j][k] - v[i][k];

          kineng_core_loc += mtot * vcore[k] * vcore[k];
          kineng_drude_loc += mu * vdrude[k] * vdrude[k];
            
          if (thermalize) {
            fcore[k]  = Ccore  * random_core->gaussian()  - Gcore  * vcore[k];
            fdrude[k] = Cdrude * random_drude->gaussian() - Gdrude * vdrude[k];
           
            if (zero) fcoreloc[k]  += fcore[k];
            
            f[i][k] += mi * fcore[k] - fdrude[k];
            f[j][k] += mj * fcore[k] + fdrude[k];

            // TODO tally energy if asked
          }
        }
      }
    }
  }
  
  if(zero && thermalize) { // Remove the drift
    MPI_Allreduce(fcoreloc,  fcoresum,  dim, MPI_DOUBLE, MPI_SUM, world);
    for (int k=0; k<dim; k++) fcoresum[k] /= ncore;
    for (int i=0; i<nlocal; i++) {
      if (mask[i] & groupbit) { // only the cores need to be in the group
        if (drudetype[i] == 0) {
          for (int k=0; k<dim; k++) f[i][k] -= fcoresum[k];
        } else {
          if (drudetype[i] == 2) continue; // Done together with the core
          int j = atom->map(drudeid[i]);
          double mi, mj, mtot; // i is core, j is drude
          if (rmass) {
            mi = rmass[i];
            mj = rmass[j];
          } else {
            mi = mass[type[i]]; 
            mj = mass[type[j]];
          } 
          mtot = mi + mj;
          mi /= mtot;
          mj /= mtot;
          for (int k=0; k<dim; k++) {
            f[i][k] -= mi * fcoresum[k];
            f[j][k] -= mj * fcoresum[k];
          }
        }
      }
    }
  }

  kineng_core_loc *= 0.5 * mvv2e;
  kineng_drude_loc *= 0.5 * mvv2e;
  MPI_Allreduce(&kineng_core_loc,&kineng_core,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&kineng_drude_loc,&kineng_drude,1,MPI_DOUBLE,MPI_SUM,world);
  temp_core = 2.0 * kineng_core / (dof_core * kb);
  temp_drude = 2.0 * kineng_drude / (dof_drude * kb);
  
  // Reverse communication of the forces on ghost Drude particles
  comm->reverse_comm();
}

/* ---------------------------------------------------------------------- */

void FixLangevinDrude::reset_target(double t_new)
{
  t_target_core = t_start_core = t_new;
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixLangevinDrude::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"t_target_core") == 0) {
    return &t_target_core;
  } else if (strcmp(str,"t_target_drude") == 0) {
    return &t_target_drude;
  } else error->all(FLERR, "Illegal extract string in fix langevin/drude");
  return NULL;
}

/* ---------------------------------------------------------------------- */

double FixLangevinDrude::compute_vector(int n)
{
    switch(n) {
        case 0: return temp_core;
        case 1: return temp_drude;
        case 2: return dof_core;
        case 3: return dof_drude;
        case 4: return kineng_core;
        case 5: return kineng_drude;
        default: return 0.;
    }
}

