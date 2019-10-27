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

#include <math.h>
#include <stdlib.h>
#include "angle_cosine_nonferrous.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleCosineNonferrous::AngleCosineNonferrous(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleCosineNonferrous::~AngleCosineNonferrous()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(mkT);  // destroy (T)
    memory->destroy(cy);   // destroy (T)
    memory->destroy(c0);   // destroy (T)
  }
}

/* ---------------------------------------------------------------------- */

void AngleCosineNonferrous::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double rsq1,rsq2,r1,r2,c,a,a11,a12,a22;
  double vdelx1,vdely1,vdelz1,vdelx2,vdely2,vdelz2; // relative vel (T)
  double r2v1,r1v1,r1v2,r2v2;           // inner product of r and v (T)
  double kT,yT,dc;                                 // kT, yT, dc/dt (T)

  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **v = atom->v;                            // velocity (T)
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

/* ----------------------------------------------------------------------
   plasticity statement
------------------------------------------------------------------------- */
  int *iT = new int[nanglelist];                   // recorder (T)
  for (n = 0; n < nanglelist; n++) *(iT+n) = 0;    // initialize iT[n] (T)

/* -----------------------------------------------------------------------*/

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // c = cosine of angle

    c = delx1*delx2 + dely1*dely2 + delz1*delz2; // inner product (T)
    c /= r1*r2;                   // inner product of unit vector (T)
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;


    // force & energy

/* ----------------------------------------------------------------------
   calculate the tendency of angle
------------------------------------------------------------------------- */

    vdelx1 = v[i1][0] - v[i2][0];            // relative velocity (T)
    vdely1 = v[i1][1] - v[i2][1];            // relative velocity (T)
    vdelz1 = v[i1][2] - v[i2][2];            // relative velocity (T)

    vdelx2 = v[i3][0] - v[i2][0];            // relative velocity (T)
    vdely2 = v[i3][1] - v[i2][1];            // relative velocity (T)
    vdelz2 = v[i3][2] - v[i2][2];            // relative velocity (T)   
    
    r2v1 = delx2*vdelx1 + dely2*vdely1 + delz2*vdelz1; // r2*v1 (T)
    r1v1 = delx1*vdelx1 + dely1*vdely1 + delz1*vdelz1; // r1*v1 (T)
    r1v2 = delx1*vdelx2 + dely1*vdely2 + delz1*vdelz2; // r1*v2 (T)
    r2v2 = delx2*vdelx2 + dely2*vdely2 + delz2*vdelz2; // r2*v2 (T)
    
    dc = (r2v1 - c*r1v1)/r1 + (r1v2 - c*r2v2)/r2;      // dc/dt (T)

/* ----------------------------------------------------------------------
   plasticity statement
------------------------------------------------------------------------- */

    if (c > cy[type] || *(iT+n) > 0){      // *iT > 0, into plastic (T)
      if (dc >= 0) {                       // dc >= 0, angle decrease (T)
        kT = k[type]*mkT[type];            // mkT[type] (T)
        yT = c0[type];                     // c0[type] (T)
        *(iT+n) = 1;
      }
      else {                               // inn > 0, angle increase (T)
        kT = k[type];                      // mkT[type] (T)
        yT = c0[type];                     // c0[type] (T)
        *(iT+n) = 1;
      }
    }

    else {                                 // not into plastic region (T)
      kT = k[type];
      yT = 1.0;
    }

/* ---------------------------------------------------------------------- */

    if (eflag) eangle = kT*(yT+c); //
 
    a = kT;                           // (T)
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;
    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
  }
  
  delete [] iT; // release memory (T)
  
}

/* ---------------------------------------------------------------------- */

void AngleCosineNonferrous::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(k,n+1,"angle:k");
  memory->create(mkT,n+1,"angle:mkT"); // allocate memory (T)
  memory->create(cy,n+1,"angle:cy");   // allocate memory (T)
  memory->create(c0,n+1,"angle:c0");   // allocate memory (T)
  
  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void AngleCosineNonferrous::coeff(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nangletypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double mkT_one = force->numeric(FLERR,arg[2]); // set plasticity coeff(T)
  double cy_one = force->numeric(FLERR,arg[3]);  // set plasticity coeff(T)
  double c0_one = force->numeric(FLERR,arg[4]);  // set plasticity coeff(T)

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    mkT[i] = mkT_one;                            // set plasticity coeff(T)
    cy[i] = cy_one;                              // set plasticity coeff(T)
    c0[i] = c0_one;                              // set plasticity coeff(T)
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleCosineNonferrous::equilibrium_angle(int i)
{
  return MY_PI;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleCosineNonferrous::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&mkT[1],sizeof(double),atom->nangletypes,fp); // write mkT (T)
  fwrite(&cy[1],sizeof(double),atom->nangletypes,fp);  // write  cy (T)
  fwrite(&c0[1],sizeof(double),atom->nangletypes,fp);  // write  c0 (T)
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleCosineNonferrous::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nangletypes,fp);
    fread(&mkT[1],sizeof(double),atom->nangletypes,fp);    // read mkT (T)
    fread(&cy[1],sizeof(double),atom->nangletypes,fp);     // read  cy (T)
    fread(&c0[1],sizeof(double),atom->nangletypes,fp);     // read  c0 (T)
  }
  MPI_Bcast(&k[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&mkT[1],atom->nangletypes,MPI_DOUBLE,0,world); // MPI mkT (T)
  MPI_Bcast(&cy[1],atom->nangletypes,MPI_DOUBLE,0,world);  // MPI  cy (T)
  MPI_Bcast(&c0[1],atom->nangletypes,MPI_DOUBLE,0,world);  // MPI  c0 (T)

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleCosineNonferrous::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,k[i],mkT[i],cy[i],c0[i]); //
}

/* ---------------------------------------------------------------------- */

double AngleCosineNonferrous::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;
  double **v = atom->v;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double vdelx1,vdely1,vdelz1,vdelx2,vdely2,vdelz2; // relative vel (T)
  double r2v1,r1v1,r1v2,r2v2;           // inner product of r and v (T)
  double kT,yT,dc;                                 // kT, yT, dc/dt (T)

  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

/* ----------------------------------------------------------------------
   calculate the tendency of angle
------------------------------------------------------------------------- */

    vdelx1 = v[i1][0] - v[i2][0];            // relative velocity (T)
    vdely1 = v[i1][1] - v[i2][1];            // relative velocity (T)
    vdelz1 = v[i1][2] - v[i2][2];            // relative velocity (T)

    vdelx2 = v[i3][0] - v[i2][0];            // relative velocity (T)
    vdely2 = v[i3][1] - v[i2][1];            // relative velocity (T)
    vdelz2 = v[i3][2] - v[i2][2];            // relative velocity (T)   

    r2v1 = delx2*vdelx1 + dely2*vdely1 + delz2*vdelz1; // r2*v1 (T)
    r1v1 = delx1*vdelx1 + dely1*vdely1 + delz1*vdelz1; // r1*v1 (T)
    r1v2 = delx1*vdelx2 + dely1*vdely2 + delz1*vdelz2; // r1*v2 (T)
    r2v2 = delx2*vdelx2 + dely2*vdely2 + delz2*vdelz2; // r2*v2 (T)
    
    dc = (r2v1 - c*r1v1)/r1 + (r1v2 - c*r2v2)/r2;      // dc/dt (T)
    
/* ----------------------------------------------------------------------
   plasticity statement
------------------------------------------------------------------------- */
  int *iT = new int(0);         // iT: recorder(T)
  
  if (c > cy[type] || *iT > 0){ // *iT > 0, into plastic (T)
    if (dc >= 0) {              // dc >= 0, angle decrease (T)
      kT = k[type]*mkT[type];  // mkT[type] (T)
      yT = c0[type];           // cy (T)
      *iT = 1;
    }
    else {                      // dc < 0, angle increase (T)
      kT = k[type];            // mkT[type] (T)
      yT = c0[type];           // cy (T)
      *iT = 1;
    }
  }
  else {                        // not into plastic region (T)        
      kT = k[type];
      yT = 1.0;
  }

  return kT*(yT+c);

  delete iT;

/* ---------------------------------------------------------------------- */

}
