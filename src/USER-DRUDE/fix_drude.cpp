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
#include "fix_drude.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "error.h"
#include "domain.h"

#include <set>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;


FixDrude *FixDrude::sptr = NULL;

/* ---------------------------------------------------------------------- */

FixDrude::FixDrude(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 2 + atom->ntypes) error->all(FLERR,"Illegal fix drude command");
  memory->create(drudetype, atom->ntypes+1, "fix_drude::drudetype");
  for (int i=2; i<narg; i++) drudetype[i-1] = force->inumeric(FLERR,arg[i]);
  memory->create(drudeid, atom->nmax, "fix_drude::drudeid");
  comm_border = 1; // drudeid
  
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);
  atom->add_callback(2);
}

/* ---------------------------------------------------------------------- */

FixDrude::~FixDrude()
{
  memory->destroy(drudetype);
  memory->destroy(drudeid);
}

/* ---------------------------------------------------------------------- */

int FixDrude::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   look in bond lists for Drude partner tags and fill atom->drudeid
------------------------------------------------------------------------- */
void FixDrude::build_drudeid(){
  int nlocal = atom->nlocal;
  int *type = atom->type;

  std::vector<tagint> drude_vec; // list of my Drudes' tags
  std::vector<tagint> core_drude_vec;
  partner_set = new std::set<tagint>[nlocal]; // Temporary sets of bond partner tags
  
  sptr = this;
  // Build list of my atoms' bond partners
  for (int i=0; i<nlocal; i++){
    if (drudetype[type[i]] == NOPOL_TYPE) continue;
    drudeid[i] = 0;
    for (int k=0; k<num_bond[i]; k++){
      core_drude_vec.push_back(atom->tag[i]);
      core_drude_vec.push_back(bond_atom[i][k]);
    }
  }
  // Loop on procs to fill my atoms' sets of bond partners
  comm->ring(core_drude_vec.size(), sizeof(tagint),
             (char *) core_drude_vec.data(),
             4, ring_build_partner, NULL, 1);

  // Build the list of my Drudes' tags
  // The only bond partners of a Drude particle is its core, 
  // so fill drudeid for my Drudes.
  for (int i=0; i<nlocal; i++){
    if (drudetype[type[i]] == DRUDE_TYPE){
      drude_vec.push_back(atom->tag[i]);
      drudeid[i] = *partner_set[i].begin(); // only one 1-2 neighbor, the core
    }
  }
  // At this point each of my Drudes knows its core.
  // Send my list of Drudes to other procs and myself
  // so that each core finds its Drude.
  comm->ring(drude_vec.size(), sizeof(tagint), 
             (char *) drude_vec.data(), 
             3, ring_search_drudeid, NULL, 1); 
  delete [] partner_set;
}

/* ----------------------------------------------------------------------
 * when receive buffer, build the set of received Drude tags.
 * Look in my cores' bond partner tags if there is a Drude tag.
 * If so fill this core's dureid.
------------------------------------------------------------------------- */
void FixDrude::ring_search_drudeid(int size, char *cbuf){
  // Search for the drude partner of my cores
  Atom *atom = sptr->atom;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  std::set<tagint> *partner_set = sptr->partner_set;

  tagint *first = (tagint *) cbuf;
  tagint *last = first + size;
  std::set<tagint> drude_set(first, last);
  std::set<tagint>::iterator it;
  
  for (int i=0; i<nlocal; i++) {
    if (drudetype[type[i]] != CORE_TYPE || drudeid[i] > 0) continue;
    for (it = partner_set[i].begin(); it != partner_set[i].end(); it++) { // Drude-core are 1-2 neighbors
      if (drude_set.count(*it) > 0){
        drudeid[i] = *it;
        break;
      }
    }
  }
}

/* ----------------------------------------------------------------------
 * buffer contains bond partners. Look for my atoms and add their partner's
 * tag in its set of bond partners.
------------------------------------------------------------------------- */
void FixDrude::ring_build_partner(int size, char *cbuf){
  // Add partners from incoming list
  Atom *atom = sptr->atom;
  int nlocal = atom->nlocal;
  std::set<tagint> *partner_set = sptr->partner_set;
  tagint *it = (tagint *) cbuf;
  tagint *last = it + size;

  while (it < last) {
    int j = atom->map(*it);
    if (j >= 0 && j < nlocal)
      partner_set[j].insert(*(it+1));
    j = atom->map(*(it+1));
    if (j >= 0 && j < nlocal)
      partner_set[j].insert(*it);
    it += 2;
  }
}


/* ----------------------------------------------------------------------
   allocate atom-based array for drudeid
------------------------------------------------------------------------- */

void FixDrude::grow_arrays(int nmax)
{
  memory->grow(drudeid,nmax,3,"fix_drude:drudeid");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixDrude::copy_arrays(int i, int j, int delflag)
{
    drudeid[j] = drudeid[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixDrude::pack_exchange(int i, double *buf)
{
    int m = 0;
    buf[m++] = ubuf(drudeid[i]).d;
    return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixDrude::unpack_exchange(int nlocal, double *buf)
{
    int m = 0;
    drudeid[nlocal] = (tagint) ubuf(buf[m++]).i;
    return m;
}

/* ----------------------------------------------------------------------
   pack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixDrude::pack_border(int n, int *list, double *buf)
{
    int m = 0;
    for (int i=0; i<n; i++){
        j = list[i];
        buf[m++] = ubuf(drudeid[j]).d
    }
    return m;
}

/* ----------------------------------------------------------------------
   unpack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixDrude::unpack_border(int n, int first, double *buf)
{
    int m = 0;
    for (int i=first; i<first+n; i++){
        drudeid[i] = (tagint) ubuf(buf[m++]).i;
    }
    return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixDrude::pack_restart(int i, double *buf)
{
    int m = 1;
    buf[m++] = ubuf(drudeid[i]).d;
    buf[0] = m; // first is the size, then the values
    return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixDrude::unpack_restart(int nlocal, int nth)
{
    double **extra = atom->extra;

    // skip to Nth set of extra values
    int m = 0;
    for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
    m++;

    drudeid[nlocal] = (tagint) ubuf(extra[nlocal][m++]).i;
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixDrude::maxsize_restart()
{
  return 2;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixDrude::size_restart(int nlocal)
{
  return 2;
}

