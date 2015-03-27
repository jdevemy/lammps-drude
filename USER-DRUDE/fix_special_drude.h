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

FixStyle(special/drude,FixSpecialDrude)

#else

#ifndef LMP_FIX_SPECIAL_DRUDE_H
#define LMP_FIX_SPECIAL_DRUDE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSpecialDrude : public Fix {
 public:
  FixSpecialDrude(class LAMMPS *, int, char **);
  int setmask();
  void init();
  virtual void pre_force(int);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  //void post_force_respa(int, int, int);
  
 protected:
  bool done;
  int index_drudetype, index_drudeid;
};

}

#endif
#endif

