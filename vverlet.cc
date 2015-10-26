#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "md.h"

/************** Generic velocity verlet first stage ***********/
void vverlet_1(numset ns, partcl *q, ionpnt &ion, parmtr param, double dt)
{
  int i, k, id;
  double xm;
  
  //
  // UO2 lattice Verlet update
  for (i=0;i<ns.npt;i++)
    {
      id = q[i].id;
      xm = param.xm[id];
      for (k=0;k<3;k++)
        {
          q[i].xv[k] += 0.5*dt*q[i].ff[k]/xm;
          q[i].xx[k] += dt*q[i].xv[k];
        }
    }
  //
  // Ion Verlet update
  xm = param.xm[0];
  for (k=0;k<3;k++) {
    ion.xv[k] += 0.5*dt*ion.ff[k]/xm;
    ion.xx[k] += dt*ion.xv[k]; }
}

/************** Generic velocity verlet second stage ***********/
void vverlet_2(numset ns, partcl *q, ionpnt &ion, parmtr param, double dt, 
               status &sys)
{
  int i, k, id;
  double xm;
  
  for (k=0;k<NKIND;k++)
    sys.mv2[k] = 0.0;
  
  //
  // UO2 Verlet update
  for (i=0;i<ns.npt;i++)
    {
      id = q[i].id;
      xm = param.xm[id];
      for (k=0;k<3;k++)
        {
          q[i].xv[k]  += 0.5*dt*q[i].ff[k]/xm;
	  sys.mv2[id] += xm*q[i].xv[k]*q[i].xv[k]; }}
  //
  // Ion Verlet update
  xm = param.xm[0];
  for (k=0;k<3;k++) {
    ion.xv[k]  += 0.5*dt*ion.ff[k]/xm;
    sys.mv2[0] += xm*ion.xv[k]*ion.xv[k]; }
  
}

void trace_print(FILE *file, period ttime, status sys,  ionpnt ion)
{
  int k;
  double r, dx[3];
  for (k=0;k<3;k++)
    dx[k] = ion.xx[k] - sys.ix[k];
  
  r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
  
  fprintf(file, "%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n", 
          ttime.tnow*TPS, r, ion.xx[0], ion.xx[1], ion.xx[2], sys.mv2[0]*0.5);
}


