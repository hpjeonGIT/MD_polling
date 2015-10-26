/*  REED_2007C_MPI

Parallel REED code for radiation damage analysis of UO2
Coding started: Jan. 28, 2008
By Department of Applied Science, University of California, Davis
Byoungseon Jeon, Graduate student

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Basically, the REED_2007C follows the scheme of Beardmore and Gr{\o}bech-Jensen
of PRE vol 57, pp.7278--7287, 1998

Basic operation is same as a serial version. The code is implemented with
MPI library for simulations with distributed memory systems.

*/
#define Ntrace 200
#define Rate_trace 0.95
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <vector>
#include "md.h"

typedef std::vector<double> vec1d;
typedef std::vector<std::vector<std::vector<double> > > vec3d;

int ipow(int x, int y)
{
  int z=1;
  if (y > 1)
    z = ipow(x, y-1);
  return x*z;
}

void parameter_read (numset &ns, parmtr &param);
void cell_initialize(numset &ns, nbrcl cell[]);
void q_initialize   (numset ns,  partcl *q);
void zbl_e_density  (double rad[][2], double rho[][2]);
void reed_initialize(partcl    *q, ionpnt &ion, nbrcl cell[], parmtr param, 
                     status  &sys, double &dt);
void dt_dop_check_pre(partcl    *q, ionpnt &ion, nbrcl cell[], parmtr param, 
                      status  &sys, double  &dt, numset   &ns, period ttime);
void dt_dop_check   (partcl    *q, ionpnt &ion, nbrcl cell[], parmtr param, 
		     status  &sys, double  &dt, numset   &ns, period ttime,
                     vec1d &nhist, vec3d &zhist);
void force          (numset    ns, partcl  *q,  ionpnt &ion,  parmtr param, 
                     status &sys);
void e_stopping_cai (ionpnt &ion,  partcl *q,   parmtr param, status &sys,
                     double rad[][2],           double rho[][2]);
void vverlet_1(numset ns, partcl *q, ionpnt &ion, parmtr param, double dt);
void vverlet_2(numset ns, partcl *q, ionpnt &ion, parmtr param, double dt, 
	       status &sys);
void stopping_check(status sys, bool  &tag);
void trace_print   (FILE *file, period ttime, status sys, ionpnt ion);
bintree *expand_tree (bintree *p, int n);
void  terminate_tree (bintree *p);
void  terminate_check(bintree *p);
void  delete_tree    (bintree *p);
void  clean_tree     (bintree *p);

void cloning_check(numset  ns, ionpnt   ion, partcl    *q, nbrcl cell[], 
                   bintree *p, status  &sys, period ttime, double dt, 
                   double range[]);
void replay       (numset  ns, ionpnt  &ion, partcl    *q, nbrcl cell[], 
                   bintree *p, period &ttime, status  &sys, double &dt);
void stopping_check_cloning(bintree *p, bool &event);
void region_def_equi  (numset  ns, ststcl *dop, double range[]);
void distribution_read(numset &ns, ststcl *dop);
void dop_calc         (double *depth, double *weight, int &ndop, 
                       ionpnt ion, status sys);
void reed_result      (numset  ns, ststcl *dop);

void dop_update_pre (numset &ns, ststcl *dop,           double *depth);
void dop_update_reed(numset &ns, ststcl *dop, int ndop, double *depth, 
                     double *weight);
void hist_result(numset ns, int mpi_rank, vec1d &nhist, vec3d &zhist); 

int main (int argc, char *argv[])
{
  //
  // MPI variables
  int mpi_rank, mpi_size;
  MPI_Status statuses, *stat; MPI_Request *req;

  /* Variable declaration */  
  partcl      *q;   // particles of lattices
  ionpnt     ion;   // ion itself
  nbrcl cell[27];   // neighboring 27 cells
  numset      ns;   // storing number variables
  parmtr   param;   // storing data of physical parameters
  period   ttime;   // time variables
  status     sys;   // system variables
  status    sys2;   // tentative system variables
  bintree   *isp;   // ion spawning pointer
  time_t   t0,t1;   // check wall time
  double      dt;   // time step (variable)
  FILE     *file;   // file stream to record traces of particles

  ststcl *dop_all;  // parameters about Depth Of Penetration
  int flag_collect, flag_finish, flag_temp, flag_renew;
  bool tag_break, tag_run, tag_renew, single_event;
  double rad[Nrho_max][2], rho[Nrho_max][2], 
         range_spawn[Nrfn], range_temp[Nrfn];
  double *depth_local, *weight_local, *depth_extnl, *weight_extnl;

  int    i, j, k, ndop_extnl, ndop_local, nion, ncrit, node_source;  
  char fname[25];
  double a_time, i_time, c_time, b_time, ref_ke, lmax;
  t0 = time(NULL);


  //
  // MPI_initialization and random seeding
  MPI_Init(&argc, &argv); MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); i = (int) t0;  
  MPI_Bcast(&i, 1, MPI_INT, 0, MPI_COMM_WORLD); srand(i+mpi_rank);
  a_time = MPI_Wtime(); b_time = a_time; i_time = 0.0; c_time = 0.0;
  req = new MPI_Request [mpi_size-1]; stat = new MPI_Status [mpi_size-1];

  //
  // Parameter parsing and cell/particle initialization
  parameter_read(ns, param); zbl_e_density(rad, rho); 
  cell_initialize(ns, cell); q = new partcl [ns.npt];
  ns.Nmax_all   = ns.nevent_pre + ns.nevent_reed*ipow(2, Nrfn+1) + 1000; 
  ns.Nmax_local = ipow(2, Nrfn+1); dop_all   = new ststcl [ns.Nmax_all]; 
  depth_local  = new double [ns.Nmax_local];
  weight_local = new double [ns.Nmax_local];
  depth_extnl  = new double [ns.Nmax_local];
  weight_extnl = new double [ns.Nmax_local];
  q_initialize(ns, q); sprintf(fname,"dop_traject_pe%3.3d.dat", mpi_rank);
  file = fopen(fname, "w"); 
  fputs("# time   r     x      y     z   kinetic-E   mag. of ZBL    Firsov    e-stopping\n", file);

  vec3d zhist(Ncomp, std::vector <vec1d> (2, vec1d (ns.nbin, 0.0)));
  vec1d nhist(ns.nbin, 0.0);
    
  //  ########################################################################
  //  First stage - no rare event enhancing ##################################
  //  This stage can be skipped by input option
  ns.ndop_all = 0; ns.ndop = 0; ndop_local = 0; tag_run = true; flag_temp = 0;

  if (!ns.restart) {
    do { 
      //
      // Pre-REED - preparation of initial condition %%%%%%%%%%%%%%%%%%%%%%%%%%
      a_time = MPI_Wtime();
      reed_initialize(q, ion, cell, param, sys, dt); 
      ref_ke = sys.mv2[0];
      j = 0; ttime.tnow = 0.0; sys.clone_level = 0;
      flag_collect = 0; flag_finish = 0; tag_break = false;
      fprintf(file, "\n# %d th pre-REED with weight = 1.0\n", ns.ndop+1);
      c_time += MPI_Wtime() - a_time;
      do {
        //
        // MD routine of each ion simulation *********************************
        a_time = MPI_Wtime();
        dt_dop_check_pre(q, ion, cell, param, sys, dt, ns, ttime);
        vverlet_1(ns, q, ion, param, dt); force(ns, q, ion, param, sys);
        e_stopping_cai(ion, q, param, sys, rad, rho); 
        vverlet_2(ns, q, ion, param, dt, sys); stopping_check(sys, tag_break);
        c_time += MPI_Wtime() - a_time;
        //if (sys.mv2[0] < ref_ke) {
	if (j%10 == 0) {
	   // std::cout << sqrt(ion.xx[0]*ion.xx[0] + ion.xx[1]*ion.xx[1] + 
	//	    ion.xx[2]*ion.xx[2]) << " "  << sys.mv2[0]*0.5 << " "  
	 //   << dt*TPS << std::endl;


          trace_print(file, ttime, sys, ion); 
          ref_ke = sys.mv2[0]*Rate_trace; }
	j++;

        //
        // Polling at high frequency
        a_time = MPI_Wtime();
        if ( mpi_rank == 0) {  // Collect DOP results - master - MPI_tag = 0
          MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag_collect, 
                     &statuses); 
          if (flag_collect) {
            MPI_Recv(depth_extnl, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, 
                     MPI_COMM_WORLD, &statuses);
            flag_collect = 0; dop_update_pre(ns, dop_all, depth_extnl); 
          }
        }
        else { // Check finish flag - slaves - MPI_tag = 1
          MPI_Iprobe(0, 1, MPI_COMM_WORLD, &flag_finish, &statuses);
          if (flag_finish) {
            MPI_Recv(&flag_temp,1, MPI_INT, 0, 1, MPI_COMM_WORLD, &statuses);
            tag_break = true; tag_run = false; flag_finish = 0;
          }
        }
        i_time += MPI_Wtime() - a_time;
      } 
      //
      // end of loop for each ion simulation **********************************
      while (!tag_break);       

      // Calculation of DOP of each ion simulation 
      if (tag_run) {
        a_time = MPI_Wtime();
        dop_calc(depth_local, weight_local, ndop_local, ion, sys);         
        ndop_local = 0; ns.ndop++;
        c_time += MPI_Wtime() - a_time;
      }

      // Root processor updates the results of own 
      if (mpi_rank == 0) {
        a_time = MPI_Wtime();
        dop_update_pre(ns, dop_all, depth_local); 
        c_time += MPI_Wtime() - a_time;
      }

      //
      // Collecting of results from other processors $$$$$$$$$$$$$$$$$$$$$$$$$$
      // Sending at low frequency
      a_time = MPI_Wtime();
      if ( mpi_rank == 0) { // Check and send finish flag - master - MP_tag = 1
        if (ns.ndop_all >= ns.nevent_pre) {
           
          for (k=1;k<mpi_size;k++) 
            MPI_Isend(&flag_temp, 1, MPI_INT, k, 1, MPI_COMM_WORLD, &req[k-1]);

          for (i=1;i<mpi_size;i++) {
            MPI_Iprobe(i, 0, MPI_COMM_WORLD, &flag_collect, &statuses); 
            if (flag_collect) {
              MPI_Recv(depth_extnl, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, 
                       &statuses);
              flag_collect = 0; dop_update_pre(ns, dop_all, depth_extnl); 
            }
          }
          MPI_Waitall(mpi_size-1, req, stat);
          std::cout << ns.ndop_all << " pre_REED data collected " << std::endl;
          tag_run = false; 
        }
      }
      else if (tag_run) { // Send DOP results - slaves - MPI_tag = 0
        MPI_Send(depth_local, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); 
      }
      i_time += MPI_Wtime() - a_time;
      // End of communication $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      
    } // end of pre-REED simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while (tag_run);
  }
  else {
    // 
    // With restart option, initial distribution of DOP parsing
    if (mpi_rank == 0) distribution_read(ns, dop_all);
  }
  
  //
  // Now first part of REED-MD is completed here ------------------------------
  // Spawning regional limit definition
  if (mpi_rank == 0) {
    region_def_equi(ns, dop_all, range_spawn); 
    lmax = 0.0;
    for (i=0;i<ns.ndop_all;i++) lmax = std::max(lmax, dop_all[i].depth);    
    std::cout << "max. range is determined as " << lmax << " in pre-REED"
              <<std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(range_spawn, Nrfn, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&lmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  ns.dl = lmax/((double) (ns.nbin - 1));
  if (ns.dl < param.a) {
    ns.nbin = int (lmax/param.a); ns.dn = ns.nbin;
    ns.dl = lmax/((double) (ns.nbin - 1));
  }
  if (ns.nbin <2 ) {
    MPI_Finalize();
    std::cout << "Histogram bin error - exit " << std::endl;
    exit(1);
  }
  //  std::cout << ns.dl << " at " << mpi_rank            <<std::endl;
  //
  //  ########################################################################
  //  Second stage - rare event enhancing ####################################
  nion = 0; flag_renew = 0; tag_run = true; tag_renew = false; ncrit = 10;
  isp = NULL; isp = new bintree; isp->left = isp->right = NULL;
  i = 0; isp = expand_tree(isp,i);
  do {
    //
    // REED loop - initialization of each single event  %%%%%%%%%%%%%%%%%%%%%%
    a_time = MPI_Wtime();
    reed_initialize(q, ion, cell, param, sys, dt);
    ref_ke = sys.mv2[0]*Rate_trace;
    single_event  = true; // true means keeping the loop
    sys.clone_level = 0; j = 0; ttime.tnow = 0.0; isp->TACT = 0;
    flag_collect = 0; flag_finish = 0; tag_break = false; ndop_local = 0;
    fprintf(file, "\n# %d th REED with weight = 1.0\n", ns.ndop+1);
    c_time += MPI_Wtime() - a_time;
    do {                  

      // Loop of single REED event *******************************************
      a_time = MPI_Wtime();
      dt_dop_check(q, ion, cell, param, sys, dt, ns, ttime, nhist, zhist);
      vverlet_1(ns, q, ion, param, dt);
      force(ns, q, ion, param, sys); 
      e_stopping_cai(ion, q, param, sys, rad, rho);
      vverlet_2(ns, q, ion, param, dt, sys);
      cloning_check(ns, ion, q, cell, isp, sys, ttime, dt, range_spawn);
      stopping_check(sys, tag_break);
      
      if (tag_break) { 
        terminate_tree(isp); terminate_check(isp); ns.ndop++; 
        dop_calc(depth_local, weight_local, ndop_local, ion, sys); 
        stopping_check_cloning(isp, single_event);           
      }	

      if (tag_break && single_event) {              
        replay(ns, ion, q, cell, isp, ttime, sys, dt); tag_break = false;
	ref_ke = sys.mv2[0]*Rate_trace;
        fprintf(file, "\n# %d th REED with weight of %f \n", 
                ns.ndop+1, 1./pow(2., ((double) sys.clone_level)) );
      }

      c_time += MPI_Wtime() - a_time;

      if (sys.mv2[0] < ref_ke)  {
	trace_print(file, ttime, sys, ion); 
	ref_ke = sys.mv2[0]*Rate_trace; }
      j++;

      //
      // Polling at high frequency $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      a_time = MPI_Wtime();
      if ( mpi_rank == 0) { // Collect DOP results - master - MPI_tag = 0
        for (i=1;i<mpi_size;i++) {
          MPI_Iprobe(i, 0, MPI_COMM_WORLD, &flag_collect, 
                     &statuses);
          if (flag_collect) {
            MPI_Recv(&ndop_extnl,  1,          MPI_INT,    i, 0, 
                     MPI_COMM_WORLD, &statuses);
            MPI_Recv(depth_extnl,  ndop_extnl, MPI_DOUBLE, i, 4, 
                     MPI_COMM_WORLD, &statuses);
            MPI_Recv(weight_extnl, ndop_extnl, MPI_DOUBLE, i, 5, 
                     MPI_COMM_WORLD, &statuses);
            flag_collect = 0; nion++;
            dop_update_reed(ns, dop_all, ndop_extnl, depth_extnl, 
                            weight_extnl); 
          }}}
      else { // Check finish flag - slaves - MPI_tag = 1
        MPI_Iprobe(0, 1, MPI_COMM_WORLD, &flag_finish, &statuses);
        if (flag_finish) {
          MPI_Recv(&flag_temp,1, MPI_INT, 0, 1, MPI_COMM_WORLD, &statuses);
          tag_break = true; tag_run = false; tag_renew = false;
        } // Check spawning range - slaves - MPI_tag = 2
        MPI_Iprobe(0, 2, MPI_COMM_WORLD, &flag_renew, &statuses);
        if (flag_renew) {
          MPI_Recv(range_temp, Nrfn, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, 
                   &statuses); tag_renew = true; }
      }
      i_time += MPI_Wtime() - a_time;
      // End of internal polling $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


    } // End of a single REED event including spawned ions *******************
    while (!tag_break);
    clean_tree(isp);
    
    //
    // Sending new spawning points/renewal $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    if (mpi_rank == 0) {
      a_time = MPI_Wtime(); nion ++;
      dop_update_reed(ns, dop_all, ndop_local, depth_local, weight_local);
      c_time += MPI_Wtime() - a_time;
      if ( nion > ncrit ) {
        a_time = MPI_Wtime();
        ncrit = nion+10; std::cout << nion <<" th REED completed"<<std::endl;
        region_def_equi(ns, dop_all, range_spawn);        
        c_time += MPI_Wtime() - a_time; 
        // Send spawning range - master - MPI_tag = 2
        a_time = MPI_Wtime();
        for (k=1;k < mpi_size; k++) 
          MPI_Isend(range_spawn, Nrfn, MPI_DOUBLE, k, 2, MPI_COMM_WORLD,
                    &req[k-1]); 
        for (i=1;i<mpi_size;i++) {
          MPI_Iprobe(i, 0, MPI_COMM_WORLD, &flag_collect, &statuses); 
          if (flag_collect) {
            MPI_Recv(&ndop_extnl,  1,          MPI_INT,    i, 0, 
                     MPI_COMM_WORLD, &statuses);
            MPI_Recv(depth_extnl,  ndop_extnl, MPI_DOUBLE, i, 4, 
                     MPI_COMM_WORLD, &statuses);
            MPI_Recv(weight_extnl, ndop_extnl, MPI_DOUBLE, i, 5, 
                     MPI_COMM_WORLD, &statuses);
          flag_collect = 0; nion++;
          dop_update_reed(ns, dop_all, ndop_extnl, depth_extnl, weight_extnl); 
          }
        }
        MPI_Waitall(mpi_size-1, req, stat);
        i_time += MPI_Wtime() - a_time;
        std::cout << nion << " ion REED data collected " << std::endl;
      }
    }
    else if (tag_renew) {
      for (k=0;k<Nrfn;k++) range_spawn[k] = range_temp[k]; tag_renew = false; }
    // End of spawning point communication $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    //
    // Sending finish flag/sending results $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    a_time = MPI_Wtime();
    if ( mpi_rank == 0) { // Check and send finish flag - master - MPI_tag = 1
      if (nion >= ns.nevent_reed) {
        for (k=1;k<mpi_size;k++) 
          MPI_Isend(&flag_temp, 1, MPI_INT, k, 1, MPI_COMM_WORLD, &req[k-1]);
        for (i=1;i<mpi_size;i++) {
          MPI_Iprobe(i, 0, MPI_COMM_WORLD, &flag_collect, &statuses); 
          if (flag_collect) {
            MPI_Recv(&ndop_extnl,  1,          MPI_INT,    i, 0, 
                     MPI_COMM_WORLD, &statuses);
            MPI_Recv(depth_extnl,  ndop_extnl, MPI_DOUBLE, i, 4, 
                     MPI_COMM_WORLD, &statuses);
            MPI_Recv(weight_extnl, ndop_extnl, MPI_DOUBLE, i, 5, 
                     MPI_COMM_WORLD, &statuses);
            flag_collect = 0; nion++;
            dop_update_reed(ns, dop_all, ndop_extnl, depth_extnl, 
                            weight_extnl); 
          }}
        MPI_Waitall(mpi_size-1, req, stat);
        std::cout << nion << " ion REED data collected " << std::endl;        
        tag_run = false; }}
    else if (tag_run) { // Send DOP results - slaves - MPI_tag = 0
      MPI_Send(&ndop_local,  1,          MPI_INT,    0, 0, MPI_COMM_WORLD);
      MPI_Send(depth_local,  ndop_local, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
      MPI_Send(weight_local, ndop_local, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
    }
    i_time += MPI_Wtime() - a_time;
    // End of finish flag/sending results $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    a_time = MPI_Wtime();
    fprintf(file, "## single ion REED closed\n\n");

    c_time += MPI_Wtime() - a_time;
  } // End of a single REED simulation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  while (tag_run);  
  
  if (mpi_rank == 0) reed_result(ns, dop_all);

  delete(q); delete(dop_all); fclose(file);
  delete(depth_local); delete(weight_local); 
  delete(depth_extnl); delete(weight_extnl);
  delete(req); delete(stat); delete_tree(isp);

  hist_result(ns, mpi_rank, nhist, zhist);

  /* Finish */
  std::cout << c_time << " " << i_time << " " << MPI_Wtime() - b_time 
            << " at " << mpi_rank << std::endl;

  MPI_Finalize();
  t1 = time(NULL);
  printf ("Elapsed wall clock time: %ld seconds\n", (long) (t1 - t0));
  return 0;
}
