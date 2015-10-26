//#define eps0 8.854187817E-12  // in SI unit with C^2/J/m - electric constant
//#define eps_0 0.005526349867722157 // in MD unit with eV/1e/A
//#define e_charge 1.602176487e-19 // in Si unit with C - elementary charge
#define hbar 0.064654148215055  // planck constant/2pi in eV*MD time unit
#define aa 0.4685024785301      // radius parameter - in Firsov potential
#define PI 3.1415926535898
#define firsov 0.051444896719640566 // Firsov constant = 0.7hbar/PI/aB^2
#define vf_cai 226.18681903945458 // hbar*(9pi/4)^(1/3)/me in Cai's formula
#define x0_cai 0.014105695605393126 // pi*aB*me/hbar in Cai's formula
#define dEdx_cai_old 0.6156910851288977 // 2*me^2*e^4/4/pi^2/eps0^2/3/hbar^3
#define dEdx_cai 0.048995139807938501 //2me^2e^4/3/4^2/pi^3/eps0^2/hbar^3
#define me  5.4857990943E-4    // electron mass in amu
#define vB 222.71802408586265   // Bohr velocity in amu/eV/A/MD time unit
#define a00 0.24005195         // Lambda constant from Brandt and Kitagawa

#include <iostream>
#include <cmath>
#include "md.h"
 
//
//  force-field
// 1. ZBL
// 2. Firsov inelastic collision
// 3. Cai's electron stopping
void force(numset ns, partcl *q, ionpnt &ion, parmtr param, status &sys)
{    
  
  /* Variable declaration */
  int i, j, k, id, jd, id_min;
  double r2, r, dx[3], dv[3], z1, z2, ff[3], med, x, phi[4], alpha, epot, 
         r_min, zf1, zf2;
  
  /* Force initialization */
  for(i=0;i<ns.npt;i++)
    {
      for(k=0;k<3;k++)
        {
	  q[i].ff[k] = 0.0;
	}
    }        
  for (k=0;k<3;k++)
    ion.ff[k] = 0.0;
  
  r_min = param.rc;
  id_min = sys.id_min; // it gaurantees that id_min will be a certain number
  for (i=0;i<NKIND;i++)
    sys.epot[i] = 0.0;
  
  //
  // ion <-> lattice (UO2) interactions
  z1 = param.z[0];
  for(i=0;i<ns.npt;i++) {
    id = q[i].id;
    z2 = param.z[id];
    r2 = 0.0;
    for(k=0;k<3;k++) {
      dx[k] = ion.xx[k] - q[i].xx[k];
      dv[k] = ion.xv[k] - q[i].xv[k]; 
      r2 += dx[k]*dx[k]; }
    if (r2 < param.rc2) {
       
      // ZBL potential
      r = sqrt(r2);
      x = r/param.au[0][id];
      phi[0] = 0.18175*exp(-3.1998*x);
      phi[1] = 0.50986*exp(-0.94229*x);
      phi[2] = 0.28022*exp(-0.4029*x);
      phi[3] = 0.028171*exp(-0.20162*x);
      med = z1*z2*( (x*3.1998+1.0)*phi[0] + (x*0.94229+1.0)*phi[1] + 
		    (x*0.4029+1.0)*phi[2] + (x*0.20162+1.0)*phi[3] 
		    )*EPS/r2/r;
      
      for (k=0;k<3;k++) {
	ff[k] = dx[k]*med;
	ion.ff[k]  += ff[k];
	q[i].ff[k] -= ff[k]; 
        if (id == 1) ion.fsum[0][k] += ff[k];
        else ion.fsum[1][k] += ff[k];
      }
      epot = 0.5*z1*z2*EPS*(phi[0]+phi[1]+phi[2]+phi[3])/r;
      sys.epot[0]  += epot;
      sys.epot[id] += epot;
      
      // Firsov collision
      alpha = param.alpha[id];
      zf1 = param.zz[id][0];
      zf2 = param.zz[id][1];
      med = (zf1*zf1/pow(1.+0.8*alpha*pow(zf1,1./3.)*r/aa,4) + 
             zf2*zf2/pow(1.+0.8*(1.-alpha)*pow(zf2,1./3.)*r/aa,4) )*firsov;
      for (k=0;k<3;k++) {
	ff[k] = dv[k]*med;
        ion.ff[k]  -= ff[k];
        q[i].ff[k] += ff[k];        
        if (id == 1) ion.fsum[2][k] -= ff[k];
        else ion.fsum[3][k] -= ff[k];
      }
      
      // Closest particle search
      if (r < r_min) {
        r_min = r;
        id_min = i; }
    }}
  sys.r_min = r_min;
  sys.id_min = id_min;
  
  
  //
  // UO2 <-> UO2 interactions
  for(i=0; i<ns.npt-1;i++) {
    id = q[i].id;
    z1 = param.z[id];
    for(j=i+1;j<ns.npt;j++) {
      jd = q[j].id;
      z2 = param.z[jd];
      r2 = 0.0;
      for(k=0;k<3;k++) {
	dx[k] = q[i].xx[k] - q[j].xx[k];
	r2 += dx[k]*dx[k]; }
      
      // ZBL
      if (r2 < param.rc2) {
	r = sqrt(r2);
	x = r/param.au[id][jd];
	phi[0] = 0.18175*exp(-3.1998*x);
	phi[1] = 0.50986*exp(-0.94229*x);
	phi[2] = 0.28022*exp(-0.4029*x);
	phi[3] = 0.028171*exp(-0.20162*x);
	med = z1*z2*( (x*3.1998+1.0)*phi[0] + (x*0.94229+1.0)*phi[1] + 
		      (x*0.4029+1.0)*phi[2] + (x*0.20162+1.0)*phi[3] 
		      )*EPS/r2/r;
	for (k=0;k<3;k++) {
	  ff[k] = med*dx[k];
	  q[i].ff[k] += ff[k];
	  q[j].ff[k] -= ff[k]; 
        }
	epot = 0.5*z1*z2*EPS*(phi[0]+phi[1]+phi[2]+phi[3])/r;
        sys.epot[id] += epot;
        sys.epot[jd] += epot; 
      }}}
}

// Electron stopping by Cai et al (PRB vol54, pp.17147--17157, 1996)
void e_stopping_cai(ionpnt &ion, partcl *q, parmtr param, status &sys, 
                    double rad[][2], double rho[][2])
{
  int i, k, n, id;
  double rs, v1, v2, v4, vf, vf2, vf4, vr, yr, qeff, lambda, gamma, C, 
         z_eff, dEdx, G, x, x2, z1, x0, Aforce;
  bool tag;
  
  // Cai's electron stopping
  rs = param.rs;
  z1 = param.z[0];
  v1 = sqrt(ion.xv[0]*ion.xv[0]+ion.xv[1]*ion.xv[1]+ion.xv[2]*ion.xv[2]);
  
  // Fermi velocity
  // vf = hbar / m * (9 pi/4)^(1/3) / rs
  vf= vf_cai/rs;
  if (v1 >= vf) 
    vr = v1*(1.+vf*vf/5./v1/v1);
  else {
    v2 = v1*v1;
    v4 = v2*v2;
    vf2 = vf*vf;
    vf4 = vf2*vf2;
    vr = 0.75*vf*(1.+(2.*v2/(3.*vf2)) - (v4/(15.*vf4))); }
  yr = vr/vB/pow(z1,2./3.);
  if (yr < 0.1 ) yr = 0.1;

  qeff = 1.-exp(-0.95*(yr-0.07));
  lambda = aB*2.*a00*pow((1.-qeff),2./3.) / pow(z1, 1./3.) / (1.-(1.-qeff)/7.);

  C = vB*vB/vf/vf/2.;
  gamma = qeff + C*(1.-qeff)*log(1.+16.*lambda*lambda/rs/rs);     
  z_eff = z1*gamma;
  x0 = x0_cai*vf;
  dEdx = dEdx_cai*v1*(log(1.+x0) - x0/(1.+x0));
  
  //
  // Now, we decide rs depending on the breach of ion to the closest atom
  
  id = q[sys.id_min].id;
  if (sys.r_min < param.muffin[id]) {
    tag = true;
    i = Nrho_max;
    if (sys.id_min != ion.pka_id) {
      //std::cout << "struck with " << sys.id_min << std::endl;
      ion.pka_id = sys.id_min;
      ion.pka[id-1] ++;
    }
    do {        
      i--;
      if (sys.r_min < rad[i][id-1] || i==0) {        
	n = i;
	tag = false; }
    } while (tag);
    
    rs = pow(3./4./PI/rho[n][id-1], 1./3.); 
    //std::cout << rs << "rs" << rho[n][id-1] << std::endl;
  }
  else ion.pka_id = -1;
  
  
  // end of breach check
  
  x = rs/aB;
  x2 = x*x;
  G = 1. + 0.717*x - 0.125*x2 - 0.0124*x2*x + 0.00212*x2*x2;
  dEdx = dEdx*z_eff*z_eff*G;
  //std::cout << z_eff << " " << gamma  << " " << yr << std::endl;
  
  for (k=0;k<3;k++) 
    {
      Aforce = dEdx*ion.xv[k]/v1;
      ion.ff[k] -= Aforce;
      ion.fsum[4][k] -= Aforce;
    }
  ion.n_count ++;
}

//
// ZBL electron denisty data loading - U and O
void zbl_e_density(double rad[][2], double rho[][2])
{
  int i;
  double radius[Nrho_max][2] = {
    { 0.00104000000283  ,  0.00233999988995 },
    { 0.00206999992952  ,  0.00467999977991 },
    { 0.00310999993235  ,  0.0070199999027 },
    { 0.00414999993518  ,  0.00935999955982 },
    { 0.00518999993801  ,  0.0116999996826 },
    { 0.0062199998647  ,  0.0140000004321 },
    { 0.00725999986753  ,  0.0164000000805 },
    { 0.00829999987036  ,  0.0186999998987 },
    { 0.00932999979705  ,  0.021099999547 },
    { 0.0104000000283  ,  0.0233999993652 },
    { 0.0124000003561  ,  0.0281000006944 },
    { 0.0144999995828  ,  0.0328000001609 },
    { 0.0165999997407  ,  0.0375000014901 },
    { 0.0186999998987  ,  0.0421000011265 },
    { 0.0207000002265  ,  0.0467999987304 },
    { 0.0228000003844  ,  0.0515000000596 },
    { 0.0249000005424  ,  0.0562000013888 },
    { 0.0270000007004  ,  0.0608999989927 },
    { 0.0289999991655  ,  0.0654999986291 },
    { 0.0310999993235  ,  0.0702000036836 },
    { 0.035300001502  ,  0.0795999988914 },
    { 0.0394000001252  ,  0.0890000015497 },
    { 0.0476999990642  ,  0.0983000025153 },
    { 0.0518999993801  ,  0.108000002801 },
    { 0.0560000017285  ,  0.116999998689 },
    { 0.0601000003517  ,  0.126000002027 },
    { 0.0643000006676  ,  0.136000007391 },
    { 0.068400003016  ,  0.144999995828 },
    { 0.0725999996066  ,  0.153999999166 },
    { 0.0808999985456  ,  0.16400000453 },
    { 0.0891999974847  ,  0.182999998331 },
    { 0.0974999964237  ,  0.201000005007 },
    { 0.105999998748  ,  0.219999998808 },
    { 0.11400000006  ,  0.238999992609 },
    { 0.122000001371  ,  0.256999999285 },
    { 0.130999997258  ,  0.270000010729 },
    { 0.138999998569  ,  0.294999986887 },
    { 0.146999999881  ,  0.31400001049 },
    { 0.156000003219  ,  0.331999987364 },
    { 0.172000005841  ,  0.351000010967 },
    { 0.188999995589  ,  0.388999998569 },
    { 0.204999998212  ,  0.425999999046 },
    { 0.222000002861  ,  0.462999999523 },
    { 0.238999992609  ,  0.500999987125 },
    { 0.254999995232  ,  0.537999987602 },
    { 0.272000014782  ,  0.575999975204 },
    { 0.287999987602  ,  0.612999975681 },
    { 0.305000007153  ,  0.651000022888 },
    { 0.354999989271  ,  0.688000023365 },
    { 0.388000011444  ,  0.726000010967 },
    { 0.421000003815  ,  0.800999999046 },
    { 0.453999996185  ,  0.875 },
    { 0.486999988556  ,  0.949999988079 },
    { 0.521000027657  ,  1.02999997139 },
    { 0.554000020027  ,  1.10000002384 },
    { 0.587000012398  ,  1.17999994755 },
    { 0.620000004768  ,  1.25 },
    { 0.652999997139  ,  1.32000005245 },
    { 0.72000002861  ,  1.39999997616 },
    { 0.786000013351  ,  1.47000002861 },
    { 0.851999998093  ,  1.49000000954 },
    { 0.919000029564  ,  1.5 },
    { 0.985000014305  ,  1.50999999046 },
    { 1.04999995232  ,  1.51999998093 },
    { 1.12000000477  ,  1.52999997139 },
    { 1.25  ,  1.53999996185 },
    { 1.32000005245  ,  1.54999995232 }       
  };
  
  double  density[Nrho_max][2] = {
    { 2650000.0  ,  2010.0 },
    { 1850000.0  ,  1880.0 },
    { 1300000.0  ,  1750.0 },
    { 912000.0  ,  1630.0 },
    { 647000.0  ,  1520.0 },
    { 464000.0  ,  1410.0 },
    { 340000.0  ,  1320.0 },
    { 254000.0  ,  1230.0 },
    { 196000.0  ,  1140.0 },
    { 156000.0  ,  1070.0 },
    { 110000.0  ,  926.0 },
    { 86200.0  ,  805.0 },
    { 72700.0  ,  699.0 },
    { 63100.0  ,  608.0 },
    { 55100.0  ,  528.0 },
    { 48000.0  ,  460.0 },
    { 41500.0  ,  400.0 },
    { 35600.0  ,  348.0 },
    { 30400.0  ,  303.0 },
    { 26000.0  ,  264.0 },
    { 19100.0  ,  201.0 },
    { 14600.0  ,  153.0 },
    { 9820.0  ,  118.0 },
    { 8570.0  ,  90.5 },
    { 7620.0  ,  70.1999969482 },
    { 6820.0  ,  54.7999992371 },
    { 6070.0  ,  43.2000007629 },
    { 5370.0  ,  34.5 },
    { 4700.0  ,  27.8999996185 },
    { 3520.0  ,  23.0 },
    { 2600.0  ,  16.3999996185 },
    { 1950.0  ,  12.6999998093 },
    { 1530.0  ,  10.5 },
    { 1270.0  ,  9.10000038147 },
    { 1110.0  ,  8.19999980927 },
    { 1010.0  ,  7.53999996185 },
    { 925.0  ,  7.0 },
    { 850.0  ,  6.53000020981 },
    { 775.0  ,  6.09000015259 },
    { 617.0  ,  5.65000009537 },
    { 465.0  ,  4.84999990463 },
    { 336.0  ,  4.09999990463 },
    { 237.0  ,  3.43000006676 },
    { 168.0  ,  2.83999991417 },
    { 122.0  ,  2.33999991417 },
    { 93.3000030518  ,  1.91999995708 },
    { 75.5999984741  ,  1.55999994278 },
    { 64.5  ,  1.26999998093 },
    { 46.7999992371  ,  1.03999996185 },
    { 38.0999984741  ,  0.841000020504 },
    { 30.0  ,  0.555000007153 },
    { 23.0  ,  0.368000000715 },
    { 17.2999992371  ,  0.246000006795 },
    { 12.8999996185  ,  0.165999993682 },
    { 9.68000030518  ,  0.112999998033 },
    { 7.40999984741  ,  0.0780000016093 },
    { 5.82999992371  ,  0.0549000017345 },
    { 4.71999979019  ,  0.0386000014842 },
    { 3.32999992371  ,  0.0284000001848 },
    { 2.51999998093  ,  0.0206000003964 },
    { 1.95000004768  ,  0.0206000003964 },
    { 1.50999999046  ,  0.0206000003964 },
    { 1.15999996662  ,  0.0206000003964 },
    { 0.897000014782  ,  0.0206000003964 },
    { 0.688000023365  ,  0.0206000003964 },
    { 0.418999999762  ,  0.0206000003964 },
    { 0.342000007629  ,  0.0206000003964 },    
  };
  
  for (i=0;i<Nrho_max;i++) {
    rad[i][0] = radius[i][0];
    rad[i][1] = radius[i][1];
    rho[i][0] = density[i][0];
    rho[i][1] = density[i][1]; 
  }
}
