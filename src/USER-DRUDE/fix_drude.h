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

FixStyle(drude,FixDrude)

#else

#ifndef LMP_FIX__DRUDE_H
#define LMP_FIX_DRUDE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDrude : public Fix {
 public:
  FixDrude(class LAMMPS *, int, char **);
  virtual ~FixDrude();
  int setmask();

  void grow_arrays(int nmax);
  void copy_arrays(int i, int j, int delflag);
  int pack_exchange(int i, double *buf);
  int unpack_exchange(int nlocal, double *buf);
  int pack_border(int n, int *list, double *buf);
  int unpack_border(int n, int first, double *buf);
  int pack_restart(int i, double *buf);
  void unpack_restart(int nlocal, int nth);
  int maxsize_restart();
  int size_restart(int nlocal);

  void build_drudeid();
  void ring_search_drudeid(int size, char *cbuf);
  void ring_build_partner(int size, char *cbuf);

  tagint * drudeid;
  int * drudetype;
  bool is_reduced;
};

}

#endif
#endif

