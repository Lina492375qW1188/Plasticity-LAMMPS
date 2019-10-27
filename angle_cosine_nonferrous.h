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

#ifdef ANGLE_CLASS

AngleStyle(cosine_nonferrous,AngleCosineNonferrous)

#else

#ifndef LMP_ANGLE_COSINE_NONFERROUS_H
#define LMP_ANGLE_COSINE_NONFERROUS_H

#include <stdio.h>
#include "angle.h"

namespace LAMMPS_NS {

class AngleCosineNonferrous : public Angle {
 public:
  AngleCosineNonferrous(class LAMMPS *);
  virtual ~AngleCosineNonferrous();
  virtual void compute(int, int);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, int, int, int);
  
  int *iT;   // recorder (T)

 protected:
  double *k, *mkT, *cy, *c0; // *kT = (*mkT)*(*k[type]); (*yT)=(*cy) (T)

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
