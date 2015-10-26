#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "md.h"

int nint(double x) // Nearest integer function
{
  int y;
  x > 0.0 ? y=int(floor(x+0.5)): y=int(ceil(x-0.5));
  //
  // if x > 0 then return floor or ceil.
  return y;
}

bintree *add_tree(numset  ns, ionpnt ion, partcl    *q, nbrcl cell[], 
		  bintree *p, status sys, period ttime, double dt);

bintree *expand_tree(bintree *p, int cnt)
{
  int i, k;
  if (p == NULL) {
    p = new bintree;
    if (p == NULL) std::cout << "out of memory for pointer allocation" 
			     << std::endl;
    p->TERM = 1; p->TACT = 0; p->left = p->right = NULL; p->level = cnt;   
  }
  if (p->left  == NULL && cnt <= Nrfn) p->left  = expand_tree(p->left,  cnt+1);
  if (p->right == NULL && cnt <= Nrfn) p->right = expand_tree(p->right, cnt+1);
  
  return p;
}

//
// If the ion breaches the regional limit, clone the ion (add one more TREE)
void cloning_check(numset  ns, ionpnt  ion, partcl    *q, nbrcl cell[], 
                   bintree *p, status &sys, period ttime, double dt,
                   double range[])
{
  int i, k, id;
  double dx[3], r2;
  
  r2 = 0.0;
  for (k=0;k<3;k++) {
    dx[k] = ion.xx[k] - sys.ix[k];
    r2 += dx[k]*dx[k]; }
  
  id = sys.clone_level;
  if ( r2 > range[id] && id < Nrfn) {
    sys.clone_level++;
    p = add_tree(ns, ion, q, cell, p, sys, ttime, dt);
    
  }	
}	

//
// Add a new ion cloning pointer
bintree *add_tree(numset  ns, ionpnt ion, partcl    *q, nbrcl cell[], 
                  bintree *p, status sys, period ttime, double dt)
{
  int i, k;
  if (p->TACT == 0) {
    p->TACT    = 1; p->left->TACT = p->right->TACT = 0;
    p->TERM    = 1;
    p->level   = sys.clone_level;
    p->dt      = dt;
    p->tnow    = ttime.tnow;
    p->mv2     = sys.mv2[0];
    p->r_min   = sys.r_min;
    p->n_count = ion.n_count;
    p->pka_id  = ion.pka_id;
    p->pka[0]  = ion.pka[0]; p->pka[1] = ion.pka[1];

    for (k=0;k<3;k++) { 
      p->ccntr[k] = sys.center[k];
      p->ix[k] = ion.xx[k];
      p->iv[k] = ion.xv[k];
      for (i=0;i<ns.npt;i++) {
	p->xx[i][k] = q[i].xx[k];
	p->xv[i][k] = q[i].xv[k];
      }
      for (i=0;i<Ncomp;i++) p->fsum[i][k] = ion.fsum[i][k];
      p->xo[k] = ion.xo[k]; p->told = ion.told;
      for (i=0;i<27;i++) p->vctr[i][k] = cell[i].id[k];
    }    
  }
  else if (p->left->TACT == 0) {
    p->left = add_tree(ns, ion, q, cell, p->left, sys, ttime, dt);}
  else if (p->left->TERM != 0) {
    p->left = add_tree(ns, ion, q, cell, p->left, sys, ttime, dt);}
  else {
    p->right = add_tree(ns, ion, q, cell, p->right, sys, ttime, dt); }
  return p;
}

//
// Terminate the current ion cloning pointer
void terminate_tree(bintree *p)
{
  if (p->left->TACT != 0) {
    if (p->left->TERM != 0) terminate_tree(p->left);
    else if (p->left->TERM == 0 && p->right->TACT == 0) {
      p->right->TACT = 1; p->right->TERM = 0; }
  }
  if (p->right->TACT != 0) {
    if (p->right->TERM != 0) terminate_tree(p->right);}
  if(p->left->TACT == 0 && p->right->TACT == 0) {
    p->left->TACT = 1; p->left->TERM = 0; }
}

//
// Search and terminate intermediate branches which have all terminated ends
void terminate_check(bintree *p)
{
  if(p->left->TACT != 0) {
    if (p->left->TERM != 0) terminate_check(p->left); }
  if(p->right->TACT != 0) {
    if (p->right->TERM != 0) terminate_check(p->right);}
  if(p->left->TACT != 0 && p->right->TACT != 0) {
    if (p->left->TERM == 0 && p->right->TERM == 0) p->TERM = 0; }
}

//
// Delete all of the pointers
void delete_tree(bintree *p)
{
  if (p->left != NULL) delete_tree(p->left);
  if (p->right!= NULL) delete_tree(p->right);
  delete p;
}

void clean_tree(bintree *p)
{
  if (p->left != NULL) clean_tree(p->left);
  if (p->right!= NULL) clean_tree(p->right);
  p->TACT = 0; p->TERM = 1;
}


//
// After cloning event is done, check whether all of the single event is done.
void stopping_check_cloning(bintree *p, bool &event)
{
  if (p->TACT == 0) event = false;
  else if ((p->left->TACT != 0) && (p->right->TACT != 0)) { 
    if ((p->left->TERM == 0) && (p->right->TERM == 0)) event = false; 
  }
  else if ((p->left->TACT != 0) && (p->right->TACT == 0)) {
    if (p->left->TERM == 0 && p->right->TERM == 0) event  = false; 
  }
}

//
// After the cloned ion is done, retrace the effective TREE and replay
// the simulation
void replay(numset ns,    ionpnt &ion, partcl *q, nbrcl cell[], bintree *p, 
            period &ttime, status &sys, double &dt)
{
  int i, k;
  if (p->left->TACT != 0) {
    if (p->left->TERM != 0) {
      replay(ns, ion, q, cell, p->left, ttime, sys, dt); }}
  if (p->right->TACT != 0) {
     if (p->right->TERM!= 0) {
       replay(ns, ion, q, cell, p->right, ttime, sys, dt); }}
  if (p->left->TERM == 0 && p->right->TACT == 0) {
    dt              = p->dt;       // Retrieve the corresponding time step
    ttime.tnow      = p->tnow;     // Retrieve the current time
    sys.clone_level = p->level;    // Retrieve the cloned level
    sys.mv2[0]      = p->mv2;
    sys.r_min       = p->r_min;
    ion.n_count     = p->n_count;
    ion.pka_id      = p->pka_id;
    ion.pka[0]      = p->pka[0]; ion.pka[1] = p->pka[1];
    for (k=0;k<3;k++) {
      sys.center[k] = p->ccntr[k]; // Retrieve the center of main cell
      ion.xx[k] = p->ix[k];        // Retrieve the ion position
      ion.xv[k] = p->iv[k];        // Retrieve the ion velocity
      for (i=0;i<ns.npt;i++) {     // Retrive lattice position/velocity
	q[i].xx[k] = p->xx[i][k];
	q[i].xv[k] = p->xv[i][k];
      }
      for (i=0;i<Ncomp;i++) ion.fsum[i][k] = p->fsum[i][k];
      ion.xo[k] = p->xo[k]; ion.told = p->told;
      for(i=0;i<27;i++) {           // Retrieve cell id vector
	cell[i].id[k] = p->vctr[i][k];
      }
    }
  }
}
//
// Regional definition - estimating distribution of dop depth
// Equi-partition version. We split the region into same number of particles.
// Still weights are varied along the cloning as 1, 1/2, 1/4, 1/8 ....
void region_def_equi(numset ns, ststcl *dop_all, double range_spawn[])
{
  
  int i, j, k, id;
  double dl, lmax, dist[Nbin], sum, depth, dN, area, area_old;
  bool tag;

  lmax = 0.0;
  for (i=0;i<ns.ndop_all;i++) {
    lmax = std::max(lmax, dop_all[i].depth);
  }

  dl = lmax / ((double) (Nbin-1));
  for (k=0;k<Nbin;k++) dist[k] = 0;
  sum = 0.0;
  
  for (i=0;i<ns.ndop_all;i++) {
    id = nint (dop_all[i].depth/dl);
    dist[id] += 1.0;  
    sum      += 1.0; }
  
  area = 0.0; area_old = 0.0;
  int n = 0;
  
  std::cout << "# cloning criterion decision " << std::endl;
  
  for (i=0;i<Nrfn;i++) {
    dN = sum*((double) (i +1)) / ((double) (Nrfn+1));
    tag = true;
    while (tag) {
      area += dist[n];
      if (area > dN ) {
	tag = false;
	depth = dl*((double) n); 
	range_spawn[i] = depth*depth;  
        std::cout << i << "th at " << depth << " \\AA with " << n << "th grid"
                  << " with assigned ions of " << area - area_old << std::endl;
        area_old = area;
      }
      n++;
    }
  } 
}

//
// Regional definition - estimating distribution of dop depth
// Equi-partition version. We split the region into same number of particles.
// Still weights are varied along the cloning as 1, 1/2, 1/4, 1/8 ....
void region_def_equi_log(numset ns, ststcl *dop_all, double range_spawn[])
{
  
  int i, j, k, id;
  double dl, lmax, dist[Nbin], sum, depth, dN, area, area_old;
  bool tag;
  
  lmax = 0.0;
  for (i=0;i<ns.ndop_all;i++) {
    lmax = std::max(lmax, dop_all[i].depth);
  }
  
  dl = log10(lmax) / ((double) (Nbin-1));
  for (k=0;k<Nbin;k++) dist[k] = 0;
  sum = 0.0;
  
  for (i=0;i<ns.ndop_all;i++) {
    id = nint (log10(dop_all[i].depth)/dl);
    dist[id] += 1.0;  
    sum      += 1.0; }
  
  area = 0.0; area_old = 0.0;
  int n = 0;
  
  std::cout << "# cloning criterion decision " << std::endl;
  
  for (i=0;i<Nrfn;i++) {
    dN = sum*((double) (i +1)) / ((double) (Nrfn+1));
    tag = true;
    while (tag) {
      area += dist[n];
      if (area > dN ) {
	tag = false;
	depth = dl*((double) n); depth = pow(10, depth);
	range_spawn[i] = depth*depth;  
        std::cout << i << "th at " << depth << " \\AA with " << n << "th grid"
                  << " with assigned ions of " << area - area_old << std::endl;
        area_old = area;
      }
      n++;
    }
  }
}


//
// Read initial distribution data - when restart option is "Yes"
void distribution_read(numset &ns, ststcl *dop_all)
{
  int i, j, n;
  double x, y;
  FILE *file;
  char dummy[255];
  
  if ((file = fopen("init_dist.dat", "r")) == NULL)
    {
      printf("Error opening file - init_dist.dat \n"); exit(1);
    }
  
  fgets(dummy, 255, file);
  i = 0;
  while (!feof(file)) {
    fscanf(file, "%d %lf %lf\n", &n, &x, &y);
    dop_all[i].depth  = x;
    dop_all[i].weight = y;
    i++;
    if (i> ns.nevent_pre) { printf("Error on index of dop distribution\n");
    exit(1); }
  }
  fclose(file);
  ns.ndop_all = i;  
}

//
// Depth of Penetration calculation
// DOP and corresponding weight are stored in dop[]
void dop_calc(double *depth, double *weight, int &ndop_local, 
              ionpnt ion, status sys)
{
  
  double dx[3];
  
  for(int k=0;k<3;k++) dx[k] = ion.xx[k] - sys.ix[k];
  depth [ndop_local] = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
  weight[ndop_local] = 1./pow(2., ((double) sys.clone_level) );
  ndop_local++;
}

//
// 1. Write each DOP result at dop_result_all.dat
// 2. Write histogram results at dop_distribution_all.dat
// When restart is required, copy dop_result_all.dat into init_dist.dat
// Even though dop_distribution_all.dat is available, I recommend you to
// postprocess separately using dop_result_all.dat. Postprocessor is
// developed with python script language
void reed_result(numset ns, ststcl *dop_all)
{
  int i, j, k, id;
  double dl, lmax, dist[Nbin2];
  std::ofstream output, result;
  
  output.open("dop_distribution_all.dat", std::ios::out);
  result.open("dop_result_all.dat", std::ios::out);
  
  result << "# index   depth  weight" << std::endl;
  
  lmax = 0.0;
  
  for (i=0;i<ns.ndop_all;i++) {
    lmax = std::max(lmax, dop_all[i].depth);
    result << setiosflags(std::ios::scientific) << i << "  " 
	   << dop_all[i].depth << " " << dop_all[i].weight << std::endl;
  }
  
  dl = 1.00 * lmax / ((double) (Nbin2-1));
  
  for (k=0;k<Nbin2;k++) dist[k] = 0;
  for (i=0;i<ns.ndop_all;i++) {
    id = nint (dop_all[i].depth/dl);
    dist[id] = dist[id] + dop_all[i].weight; 
  }
  
  output << "# distance,   number of bins" << std::endl;
  for (i=0;i<Nbin2;i++)
    output <<  setiosflags(std::ios::scientific) 
           << dl*((double) i) << "  " << dist[i] << std::endl;
 
  result.close();
  output.close();
}
