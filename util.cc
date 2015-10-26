#include <iostream>
#include <fstream>
#include <cmath>
#include "md.h"

//
// Input parameter parsing
void parameter_read(numset &ns, parmtr &param)
{
  int i, j;
  double x;
  FILE *fp;
  char dummy[255], rest[255];
  
  if ((fp = fopen("config.prm", "r")) == NULL)
    {
      printf("Error opening file - config.prm \n"); exit(1);
    }
  
  fgets(dummy, 255, fp);
  fgets(dummy, 255, fp);
  fscanf(fp, "%lf\n", &x);
  param.a = x;
  param.a_half = param.a*0.5;
  fgets(dummy, 255, fp);
  fscanf(fp, "%lf\n", &x);
  param.rc = x;
  param.rc2 = param.rc*param.rc;
  fgets(dummy, 255, fp);
  fscanf(fp, "%lf\n", &x);
  param.rs = x;
  fgets(dummy, 255, fp);
  for (i = 0; i< NKIND ; i++) {
    fscanf(fp, "%lf\n", &x); param.z[i] = x; }
  for (i = 0; i< NKIND; i ++) {
    if (param.z[i] > param.z[0]) {
      param.zz[i][0] = param.z[i];
      param.zz[i][1] = param.z[0]; 
      param.alpha[i] = 1./(1. + pow(param.z[0]/param.z[i], 1./6.)); }
    else {
      param.zz[i][0] = param.z[0];
      param.zz[i][1] = param.z[i]; 
      param.alpha[i] = 1./(1. + pow(param.z[i]/param.z[0], 1./6.)); }
    for (j = 0; j<NKIND; j++) {
      param.au[i][j] = 0.8854*aB/
        (pow(param.z[i],0.23) + pow(param.z[j],0.23)); }}
  fgets(dummy, 255, fp);
  fscanf(fp, "%lf\n", &x); param.KE = x; // initial KE of the ion in eV
  fgets(dummy, 255, fp);
  for (i=1;i <NKIND; i++) {
    fscanf(fp, "%lf\n", &x); param.t[i] = x*kb; } // Temperature in K -> eV
  fgets(dummy, 255, fp);
  for (i=0;i <NKIND; i++) {
    fscanf(fp, "%lf\n", &x); param.xm[i] = x; }
  fgets(dummy, 255, fp);
  for (i=1;i <NKIND; i++) {
    fscanf(fp, "%lf\n", &x); param.muffin[i] = x; }
  fgets(dummy, 255, fp);
  for (i=1;i <NKIND; i++) {
    fscanf(fp, "%lf\n", &x); 
    param.s2avg[i] = x/kb; // Unit = A^2/K => A^2/eV
    }
    // For Debye model
    // <s^2> = 3hbar^2 /m kB *T/Theta^2
    // Here, x is unit of K while param.t has eV - careful!
    //param.s2avg[i] = 3.*0.064654148215055*0.064654148215055*
    //	     param.t[i]/param.xm[i]/8.617343e-5/8.617343e-5/x/x; }
  fgets(dummy, 255, fp);
  fscanf(fp, "%d\n", &j); ns.nevent_pre = j;
  fscanf(fp, "%d\n", &j); ns.nevent_reed = j;
  fgets(dummy, 255, fp);
  fscanf(fp, "%s\n", rest);
  if (rest[0] == 'Y' || rest[0] == 'y') {
    ns.restart = true;}
  else {
    ns.restart = false;}
  fgets(dummy, 255, fp);
  fscanf(fp, "%d\n", &j); ns.nbin = j; ns.dn = ns.nbin;
  int fclose(FILE *fp);
	//std::cout << param.rs << std::endl;
}


//
// Allocation of internal particles of each cell
void cell_initialize(numset &ns, nbrcl cell[])
{
  int i, j, k, l, m, n;
  
  n = 0;
  m = 0;
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      for (k=0;k<3;k++) {
	
	for (l=0;l<N_u; l++) {
	  cell[n].nid[l] = m;
	  m++; }
	
	for (l=0;l<N_o; l++) {
	  cell[n].nid[l+N_u] = m;
	  m++; }
	cell[n].npt = N_u + N_o;
	n++; }}}
  ns.npt = 27*(N_u+N_o);
}

//
// Assign each q-particle into a single cell
void q_initialize(numset ns, partcl *q)
{
  int i, n;
  n = 0;
  do {
    for (i=0;i<N_u;i++) {
      q[n].id = 1; 
      n++; }
    for (i=0;i<N_o;i++) {
      q[n].id = 2;
      n++; }
  } while (n< ns.npt);
}

//
// Add dop results from pre-REED simulations
void dop_update_pre(numset &ns, ststcl *dop_all, double *depth)
{
  dop_all[ns.ndop_all].depth  = depth[0];
  dop_all[ns.ndop_all].weight = 1.0;
  ns.ndop_all ++;
}

//
// Add dop results frmo REED simulations
void dop_update_reed(numset &ns, ststcl *dop_all, int ndop, 
                     double *depth, double *weight)
{
  for (int i=0;i<ndop;i++) {
    dop_all[ns.ndop_all].depth  =  depth[i];
    dop_all[ns.ndop_all].weight = weight[i];
    ns.ndop_all ++;
  }
}
