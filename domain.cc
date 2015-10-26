#define E_crit 10.E0   // Potential energy criterion for ion stopping
#define c_diff 1.E-2
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <mpi.h>
#include "md.h"

typedef std::vector<double> vec1d;
typedef std::vector<std::vector<std::vector<double> > > vec3d;

void renew      (nbrcl   cell[], partcl      *q, int   n, parmtr   param);
void lattice_uo2(nbrcl cell[], parmtr param, int n, double ux[][3], 
                  double uv[][3], double ox[][3], double ov[][3]);
void vector_check(status    &sys, parmtr   param, ionpnt &ion, int vctr[]);
void neighbor_remap(partcl *q,    nbrcl   cell[], parmtr param, status &sys, 
                    int vctr[]);

double rand_double() // random number generator
{
  double x = ((double) rand())/((double) RAND_MAX);
  return x;
}

double normal_dist() // normal distribution number generator
{
  double v1, v2;
  double r = 1.0;
  do {
    v1 = 2.*rand_double() - 1.0;
    v2 = 2.*rand_double() - 1.0;
    r = v1*v1 + v2*v2;
  } while (r >= 1.);
  double x = v1*sqrt(-2.*log(r)/r);
  return x;
}
void change_size(numset &ns, int id, vec1d &nhist, vec3d &zhist)
{
  int i, j;
  i = id - ns.nbin;
  j = i/ns.dn + 1;
  ns.nbin += j*ns.dn;
  nhist.resize(ns.nbin,0.0);
  for (i=0;i<Ncomp;i++) {
    for (j=0;j<2;j++) {
      zhist[i][j].resize(ns.nbin,0.0);
    }
  }     
}


/**********************************************************************/
void reed_initialize(partcl *q, ionpnt &ion, nbrcl cell[], 
		     parmtr param, status &sys, double &dt)
{
  
  int i, j, k, n;
  double xm, lambda, mv2;
  
  // cell vctr(id) initialization
  n = 0;
  for (i=0;i<3;i++) {
    sys.center[i] = 0.0;
    for (j=0;j<3;j++) {
      for (k=0;k<3;k++) {
	cell[n].id[0] = i-1;
	cell[n].id[1] = j-1;
	cell[n].id[2] = k-1;
	cell[n].center[0] = ((double) (i-1)) * param.a;
	cell[n].center[1] = ((double) (j-1)) * param.a;
	cell[n].center[2] = ((double) (k-1)) * param.a;
	renew(cell, q, n, param);
	n++; }}}
  //
  // allocate initial ion position as one of the U of center cell (14th/4th)
  n = cell[13].nid[0];
  xm = param.xm[0];
  mv2 = 0.0;
  ion.pka_id = n;
  ion.pka[0] = ion.pka[1] = 0;
  for(k=0;k<3;k++) {
    ion.xx[k] = q[n].xx[k];
    sys.ix[k] = ion.xx[k];
    q[n].xx[k] = param.a*100.; // artificially big number - tricky!!!
    ion.xv[k] = rand_double() - 0.5; 
    mv2 += xm*ion.xv[k]*ion.xv[k]; 
  }
  //std::cout << "Xe " << ion.xx[0] << " " <<  ion.xx[1] << " " <<  ion.xx[2] << " " << std::endl;
  lambda = sqrt(2.0*param.KE/mv2);
  for (k=0;k<3;k++) {
    ion.xv[k] = ion.xv[k]*lambda;      
    ion.ff[k] = 0.0; 
    ion.xo[k] = ion.xx[k];
    for (i=0;i<Ncomp;i++) {
      ion.fsum[i][k] = 0.0;
    }
  }
  ion.told = 0.0; ion.n_count = 0;
  // Initial time step
  dt = c_diff/sqrt(param.KE*2./xm);
  sys.r_min = param.rc; sys.mv2[0] = param.KE*2.;
}

//
// Allocate new position and velocity of each particle per cell
void renew(nbrcl cell[], partcl *q, int n, parmtr param)
{
  int i, k, id;
  double ux[N_u][3], uv[N_u][3], ox[N_o][3], ov[N_o][3];
  
  lattice_uo2(cell, param, n, ux, uv, ox, ov);
  for (i=0;i<N_u;i++)
    {
      id = cell[n].nid[i];
      for (k=0;k<3;k++) {
	q[id].xx[k] = ux[i][k];
	q[id].xv[k] = uv[i][k]; 
        q[id].ff[k] = 0.0; }
      //      std::cout << "U " << q[id].xx[0] << " " << q[id].xx[1] << " " << q[id].xx[2] << std::endl;
    }
  for (i=0;i< N_o ; i++)
    {
      id = cell[n].nid[i+N_u];
      for (k=0;k<3;k++) {
	q[id].xx[k] = ox[i][k];
	q[id].xv[k] = ov[i][k]; 
        q[id].ff[k] = 0.0; }
      //      std::cout << "O " << q[id].xx[0] << " " << q[id].xx[1] << " " << q[id].xx[2] << std::endl;

    }
  
}

//
// UO2 lattice generation
// 1. Position follows the normal distribution by Debye temperature
// 2. Velocity is coupled with initial temperature by isokinetic coupling
void lattice_uo2(nbrcl cell[], parmtr param, int n, double ux[][3], 
		       double uv[][3], double ox[][3], double ov[][3])
{
  int i, k;
  double xm, sumv[3], mv2, lambda, alpha;


  // Uranium lattice generation
  // standard deviation of s2 by Debye temperature
  alpha = sqrt(param.s2avg[1]*param.t[1]); 
  //double uvect[4][3] = {-0.5,-0.5,-0.5, 0.0,0.0,-0.5, -0.5,0.0,0.0, 
  //                     0.0,-0.5,0.0};
  double uvect[4][3] = {{0.0,0.0,0.0},  {0.5,0.5,0.0}, 
			{0.0,0.5,0.5},  {0.5,0.0,0.5}};
  //double uvect[4][3] = {-0.4,-0.4,-0.4, 0.1,0.1,-0.4, -0.4,0.1,0.1, 
  //                      0.1,-0.4,0.1};  
  xm = param.xm[1];
  for (k=0;k<3;k++) sumv[k] = 0.0;
  for (i=0;i<N_u;i++) {
    for (k=0;k<3;k++) {
      ux[i][k] = uvect[i][k] * param.a + cell[n].center[k] 
                                       + alpha*normal_dist(); 
      uv[i][k] = rand_double() - 0.5;
      sumv[k] += uv[i][k]; }}
  for (k=0;k<3;k++) sumv[k] = sumv[k]/ ((double) (N_u));
  for (i=0;i<N_u;i++) {
    for (k=0;k<3;k++) {
      uv[i][k] -= sumv[k];
      mv2 += xm*uv[i][k]*uv[i][k]; }}
  lambda = sqrt(3.*((double) N_u)*param.t[1]/mv2);
  for (i=0;i<N_u;i++) {
    for (k=0;k<3;k++) {
      uv[i][k] = uv[i][k]*lambda; }}
  

  // oxygen lattice generation
  // standard deviation of s2 by Debye temperature
  alpha = sqrt(param.s2avg[2]*param.t[2]); 
  double ovect[8][3] = {{-0.25,-0.25,-0.25}, {-0.25,-0.25,0.25}, 
			{-0.25,0.25,-0.25}, {0.25,-0.25,-0.25}, 
			{-0.25,0.25,0.25}, {0.25,0.25,-0.25},
                        {0.25,-0.25,0.25}, {0.25,0.25,0.25}};
  //double ovect[8][3] = {-0.15,-0.15,-0.15, -0.15,-0.15,0.35,-0.15,0.35,-0.15,
  //                    0.35,-0.15,-0.15, -0.15,0.35,0.35,   0.35,0.35,-0.15,
  //                     0.35,-0.15,0.35,   0.35,0.35,0.35};
  xm = param.xm[2];
  for (k=0;k<3;k++) sumv[k] = 0.0;
  for (i=0;i<N_o;i++) {
    for (k=0;k<3;k++) {
      ox[i][k] = ovect[i][k] * param.a + cell[n].center[k] 
                                       + alpha*normal_dist();
      ov[i][k] = rand_double() - 0.5;
      sumv[k] += ov[i][k]; }}
  for (k=0;k<3;k++) sumv[k] = sumv[k]/ ((double) (N_o));
  for (i=0;i<N_o;i++) {
    for (k=0;k<3;k++) {
      ov[i][k] -= sumv[k];
      mv2 += xm*ov[i][k]*ov[i][k]; }}
  lambda = sqrt(3.*((double) N_o)*param.t[2]/mv2);
  for (i=0;i<N_o;i++) {
    for (k=0;k<3;k++) {
      ov[i][k] = ov[i][k]*lambda; }}
}

//
//
void energy_dump(ionpnt &ion, status    sys, period ttime, numset &ns,
                 vec1d &nhist, vec3d &zhist)
{
  int i, k, id;
  double dx[3], dtime, dN, r2, weight, dE[Ncomp-1];
  
  for (i=0;i<Ncomp-1;i++) dE[i] = 0.0;

  dtime = ttime.tnow - ion.told;
  dN = ((double) ion.n_count);

  //  std::cout << dN << " " << dtime << std::endl;

  for (k=0;k<3;k++) {
    dx[k] = ion.xx[k] - ion.xo[k];
    for (i=0;i<Ncomp-1;i++) {
      dE[i] -= ion.fsum[i][k]*dx[k]/dN;
      ion.fsum[i][k] = 0.0;
    }
    ion.xo[k] = ion.xx[k];
  }

  ion.told = ttime.tnow;
  ion.n_count = 0;

  r2 = 0.0;
  for (k=0;k<3;k++) {
    dx[k] = ion.xx[k] - sys.ix[k];
    r2 += dx[k]*dx[k];
  }
  id = int(sqrt(r2)/ns.dl);
  weight = 1./pow(2., ((double) sys.clone_level));
  //  std::cout << id << " " << weight << " " << ns.dl <<std::endl;
  if (id >= ns.nbin) change_size(ns, id, nhist, zhist);
  nhist[id] += weight;
  for (i=0;i<Ncomp-1;i++) {
    zhist[i][0][id] += dE[i]*weight;
    zhist[i][1][id] += dE[i]*weight/dtime;        
  }
  zhist[Ncomp-1][0][id] += ((double) ion.pka[0]) * weight;
  zhist[Ncomp-1][1][id] += ((double) ion.pka[1]) * weight;
  ion.pka[0] = ion.pka[1] = 0;
}

//
// 1. Time step determination
// 2. Check whether remapping is required
void dt_dop_check_pre(partcl   *q, ionpnt &ion, nbrcl cell[], parmtr param, 
                      status &sys, double  &dt, numset   &ns, period ttime)
{
  int i, vctr[3];
  double dt_2, dt_3, med[NKIND];

  
  sys.epot[1] = sys.epot[1] / ((double) N_u) / 27.;
  sys.epot[2] = sys.epot[2] / ((double) N_o) / 27.;
  sys.mv2[1]  = sys.mv2[1]  / ((double) N_u) / 27.;
  sys.mv2[2]  = sys.mv2[2]  / ((double) N_o) / 27.;
  
  dt_2 = 0.0;
  for (i=0;i<NKIND;i++) {
    med[i] = (sys.mv2[i] + 2.*std::max(0.0, sys.epot[i]))/param.xm[i];
    med[i] = sqrt(med[i]);
    dt_2 = std::max(dt_2, med[i]);
  }
  dt_2 =  c_diff/dt_2;
  dt = std::min(dt*1.05, dt_2*0.25 + dt*0.75);
  
  dt_3 = 0.5*sys.r_min/sqrt(sys.mv2[0]/param.xm[0]);
  dt = std::min(dt_3, dt);
  ttime.tnow += dt;
  vector_check(sys, param, ion, vctr);
  if (vctr[0] == 0 && vctr[1] == 0 && vctr[2] == 0) 
    {}
  else {
    neighbor_remap(q, cell, param, sys, vctr);
    for (int k=0;k<3;k++) {
      for (i=0;i<Ncomp;i++) ion.fsum[i][k] = 0.0; 
    }
  }
}
//
// 1. Time step determination
// 2. Check whether remapping is required
void dt_dop_check(partcl   *q, ionpnt &ion, nbrcl cell[], parmtr param, 
		  status &sys, double  &dt, numset   &ns, period ttime,
                  vec1d &nhist, vec3d &zhist)
{
  int i, vctr[3];
  double dt_2, dt_3, med[NKIND];

  
  sys.epot[1] = sys.epot[1] / ((double) N_u) / 27.;
  sys.epot[2] = sys.epot[2] / ((double) N_o) / 27.;
  sys.mv2[1]  = sys.mv2[1]  / ((double) N_u) / 27.;
  sys.mv2[2]  = sys.mv2[2]  / ((double) N_o) / 27.;
  
  dt_2 = 0.0;
  for (i=0;i<NKIND;i++) {
    med[i] = (sys.mv2[i] + 2.*std::max(0.0, sys.epot[i]))/param.xm[i];
    med[i] = sqrt(med[i]);
    dt_2 = std::max(dt_2, med[i]);
  }
  dt_2 =  c_diff/dt_2;
  dt = std::min(dt*1.05, dt_2*0.25 + dt*0.75);
  
  dt_3 = 0.5*sys.r_min/sqrt(sys.mv2[0]/param.xm[0]);
  dt = std::min(dt_3, dt);
  ttime.tnow += dt;
  
  vector_check(sys, param, ion, vctr);
  if (vctr[0] == 0 && vctr[1] == 0 && vctr[2] == 0) 
    {}
  else {
    //    std::cout << ion.xx[0] <<ion.xx[1] << ion.xx[2] << std::endl;
    neighbor_remap(q, cell, param, sys, vctr);
    energy_dump(ion, sys, ttime, ns, nhist, zhist);
  }
}

//
// if the ion breaches the limit of the main cell, it assigns direction vectors
// (0,0,0) => (+-1,+-1,+-1)
void vector_check(status &sys, parmtr param, ionpnt &ion, int vctr[])
{
  int k;
  double x;
  
  for (k=0;k<3;k++) {
    x = ion.xx[k] - sys.center[k];
    if ( x > param.a_half)
      vctr[k] = 1;
    else if ( -x > param.a_half)
      vctr[k] = -1;
    else
      vctr[k] = 0;
  }
  
}

// 
// Remapping of 27 cells
// If vctr has +1 or -1 component, all corresponding cells are remapped
// while 0 vctr keeps the current states
void neighbor_remap(partcl *q, nbrcl cell[], parmtr param, status &sys, 
		    int vctr[])
{
  int i, k;
  int id[3];
  bool tag;
  
  for (k=0;k<3;k++) {
    sys.center[k] += ((double) vctr[k]) * param.a; }
  
  for (i=0;i<27;i++) {
    for (k=0;k<3;k++) 
      id[k] = cell[i].id[k] - vctr[k] + 1;
    
    tag = false;
    for (k=0;k<3;k++) {
      if ( id[k] != (id[k]+3)%3)
	tag = true;
      cell[i].id[k] = (id[k]+3)%3 - 1;
      cell[i].center[k] = sys.center[k] + ((double) cell[i].id[k])*param.a;
    }
    if (tag) renew(cell, q, i, param);
  }
}

//
// Determine stopping simulations
// 1. if KE of the ion is lower than PE
// 2. or PE is lower than certain limit
// 1. only is not enough to decide - when the ion becomes stationary,
// neighboring lattice particles run away and PE becomes extremely small,
// yielding infinite loop for stopping check
void stopping_check(status sys, bool &tag)
{
  if (0.5*sys.mv2[0] <  E_crit)
    tag = true;
}

void hist_result(numset ns, int mpi_rank, vec1d &nhist, vec3d &zhist)
{
  int i, j, k, n_bin, n_max;
  double dE, dEdt;
  char fname[100];

  n_bin = ns.nbin;
  MPI_Allreduce(&n_bin, &n_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (n_max > n_bin) {
    nhist.resize(n_max,0.0);
    for (i=0;i<Ncomp;i++) {
      for (j=0;j<2;j++) {
        zhist[i][j].resize(n_max,0.0);
      }
    }
  }
  double z_send[Ncomp][2][n_max], z_hist[Ncomp][2][n_max];
  double n_send[n_max], n_hist[n_max];

  for (i=0;i<Ncomp;i++) {
    for (j=0;j<2;j++) {
      for (k=0;k<n_max;k++) {
        z_send[i][j][k] = zhist[i][j][k];
	z_hist[i][j][k] = 0.0;
      }
    }
  }

  for (k=0;k<n_max;k++) {
    n_send[k] = nhist[k]; n_hist[k] = 0.0;
  }


  //
  /*
  sprintf(fname,"hist_by_%3.3d.dat", mpi_rank);
  std::ofstream tentativ;
  tentativ.open(fname, std::ios::out);
  tentativ << "# histogram with nbin = " << n_max << " and dx = " 
	   << ns.dl << std::endl;
  for (k=0;k<n_max;k++) {
    tentativ << k << " " << setiosflags(std::ios::scientific) 
	     << ns.dl*(double) k << " " ;
    for (j=0;j<2;j++) {
      for (i=0;i<Ncomp;i++) {
	tentativ << setiosflags(std::ios::scientific) 
		 << zhist[i][j][k] << " ";
      }
    }
    tentativ << nhist[k] << std::endl;
  }  
  tentativ.close();
  */
  // 


  MPI_Reduce(n_send, n_hist, n_max, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  n_bin = Ncomp*2*n_max;
  MPI_Reduce(z_send, z_hist, n_bin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (mpi_rank == 0)
    {      
      std::ofstream result;      
      result.open("hist_E_dump.dat", std::ios::out);
      result << "# histogram with nbin = " << n_max << " and dx = " 
             << ns.dl << std::endl;
      result << "# n --- n*dl --- SUM of dE and dE/dt  --- dE and dE/dt by ZBL_U --- dE and dE/dt by ZBL_O --- dE and dE/dt by Fir_U --- dE and dE/dt by Fir_O --- dE and dE/dt by estopping --- Number of PKA of U/O --- acc. of histogram" << std::endl;
      for (k=0;k<n_max;k++) {
        result << k << " " << setiosflags(std::ios::scientific) 
               << ns.dl*((double) k) << " " ;
	dE = dEdt = 0.0;
	for (i=0;i<Ncomp-1;i++) {
	  dE += z_hist[i][0][k]; dEdt += z_hist[i][1][k];
	}
	if (n_hist[k] > 1.E-8) {
	  result << setiosflags(std::ios::scientific) << dE/n_hist[k] << " "
		 << dEdt/n_hist[k] << " " ;
	}
	else {
	  result << "0.0  0.0  ";
	}
	for (i=0;i<Ncomp;i++) {
	  for (j=0;j<2;j++) {
            if (n_hist[k] > 1.E-8) {
              result << setiosflags(std::ios::scientific) 
                     << z_hist[i][j][k]/n_hist[k] << " ";
            }
            else {
              result << "0.0  ";
            }
          }
        }
        result << n_hist[k] << std::endl;
      }  
      result.close();
    }
}
