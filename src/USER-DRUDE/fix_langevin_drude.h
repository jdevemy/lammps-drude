/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(langevin/drude,FixLangevinDrude)

#else

#ifndef LMP_FIX_LANGEVIN_DRUDE_H
#define LMP_FIX_LANGEVIN_DRUDE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLangevinDrude : public Fix {
 public:
  FixLangevinDrude(class LAMMPS *, int, char **);
  virtual ~FixLangevinDrude();
  int setmask();
  void init();
  void setup(int vflag);
  virtual void post_force(int vflag);
  void langevin(int vflag, bool thermalize); 
  void reset_target(double);
  virtual void *extract(const char *, int &);
  double compute_vector(int);
  int pack_reverse_comm(int, int, double*);
  void unpack_reverse_comm(int, int*, double*);

 protected:
  double t_start_core,t_period_core,t_target_core;
  double t_start_drude,t_period_drude,t_target_drude;
  int tstyle_core, tstyle_drude;
  int tvar_core, tvar_drude;
  char *tstr_core, *tstr_drude;
  double energy;
  int tflag;

  class AtomVecEllipsoid *avec;

  class RanMars *random_core, *random_drude;
  int seed, seed_drude;
  int dof_core, dof_drude;
  double kineng_core, kineng_drude;
  double temp_core, temp_drude;
  int zero;
  bigint ncore;
};

}

#endif
#endif

