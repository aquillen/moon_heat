
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h"
#include "spring.h"
#include "tools.h"
#include "output.h"
#include "kepcart.h"


// extern int NS; // number of springs
// extern int NPERT;  // number of external perturbing bodies


// make a PlutoCharon binary with two masses m1,m2 
// add those two masses 
// m1 = mp and m2 = mratio*mp
// masses are separated by distance sep, circular orbit for binary
// center of mass of PC and velocity set to be in orbit with given orbital elements
// return mean motion around PC
// r1 is radius for display of Pluto
// orbital elements are aa,ee,ii, longnode, argperi, meananom
// these refer to orbit of resolved body [0,N) about binary
// the resolved body is assumed to be alone and at origin
double add_pluto_charon_kep(struct reb_simulation* const r,double mp, 
    double mratio, double sep,
    double aa,double ee, double ii,double longnode,double argperi,double meananom, 
    double r1)
{
   // const int i1= r->N;
   // const int i2= i1+1;
   struct reb_particle pt;
   double m1 = mp; // 
   pt.ax = 0.0; pt.ay = 0.0; pt.az = 0.0;
   pt.vx = 0.0; pt.vy = 0.0; pt.vz = 0.0;
   pt.x  = 0.0; pt.y  = 0.0; pt.z  = 0.0;
   double m2 = mratio*m1;
   double mtot  = m1 + m2;
   double vc_bin = sqrt(r->G*mtot/sep);
   double x1 = sep*m2/mtot;  // initial orientation of binary is here
   double x2 = -sep*m1/mtot;
   double vy1 =  vc_bin*m2/mtot;
   double vy2 = -vc_bin*m1/mtot;
   double r2 = r1*pow(m2/m1,0.333333); // display radius
   if (r2 < 1.0) r2=1.0;
   //
   double GM = r->G*(mtot+1.0); // +1 here for extended body
   OrbitalElements orbel;
   orbel.a = aa;
   orbel.e = ee;
   orbel.i = ii;
   orbel.argperi = argperi;
   orbel.longnode= longnode;
   orbel.meananom= meananom;
   PhaseState ps;
   cartesian(GM, orbel, &ps);
   double om_orb = sqrt(GM/aa)/aa;

   // ratio of reduced mass to binary mass
   double muBratio = mratio/pow(1.0 + mratio,2.0); 
   double fac = 1.0 + (3.0/8.0)*muBratio*pow(sep/aa,2.0);
   om_orb *= fac; // correct mean motion by binary quad moment
   
   pt.x = ps.x + x2; pt.y = ps.y + 0.0; pt.z = ps.z + 0.0;
   pt.vx = ps.xd + 0.0; pt.vy = ps.yd + vy2; pt.vz = ps.zd + 0.0; 
   pt.m = m2; pt.r = r2;
   reb_add(r,pt);   // add Charon

   pt.x = ps.x + x1; pt.y = ps.y + 0.0; pt.z = ps.z + 0.0;
   pt.vx = ps.xd + 0.0; pt.vy = ps.yd + vy1; pt.vz = ps.zd + 0.0; 
   pt.m = m1; pt.r = r1;
   reb_add(r,pt);   // Pluto is last
   return om_orb;
}

// add a point mass 
// the point mass added has mass m1, radius r1 (for display)
// if ip>=0 put point mass in orbit about another point mass with index ip
// if ip==-1 
//    the extended mass has indices [il,ih)
//    add a single point mass with extended body orbiting it
//    the pt mass added is at origin
//    if extended body is rotating then it still rotates
//       however if orbit is tilted then obliquity will not be the same
// return mean motion
// positions and velocity computed with orbelements
// orbels computed with GM taking into account both masses
double add_pt_mass_kep(struct reb_simulation* const r, 
    int il, int ih, int ip, double m1, double r1,
    double aa,double ee, double ii,
    double longnode,double argperi,double meananom)  
{
   struct reb_particle* particles = r->particles;
   struct reb_particle pt;
   pt.m = m1;
   pt.r = r1;
   pt.ax = 0; pt.ay = 0; pt.az = 0;
   double GM,m0;
   double  x0=0.0;  double y0=0.0; double  z0=0.0;
   double vx0=0.0; double vy0=0.0; double vz0=0.0;
   if (ip<0){ // moving resolved body
     m0 = sum_mass(r,il,ih);
     subtractcom(r,il,ih);  // zero center of mass position of extended body
     subtractcov(r,il,ih);  // subtract center of mass velocity of extended body
   }
   else { // new particle has motion w.r.t to particle at ip
     m0 = particles[ip].m;
     x0 = particles[ip].x;    y0 = particles[ip].y;   z0 = particles[ip].z;
     vx0 = particles[ip].vx; vy0 = particles[ip].vy; vz0 = particles[ip].vz;
   }
   GM = r->G*(m1+m0); 
   OrbitalElements orbel;
   orbel.a = aa;
   orbel.e = ee;
   orbel.i = ii;
   orbel.argperi = argperi;
   orbel.longnode= longnode;
   orbel.meananom= meananom;
   PhaseState ps;
   cartesian(GM, orbel, &ps);
   double om_orb = sqrt(GM/aa)/aa;

   if (ip<0){  // move position of extended body
     pt.x   = 0.0; pt.y   = 0.0; pt.z   = 0.0;
     pt.vx  = 0.0; pt.vy  = 0.0; pt.vz  = 0.0; // new particle is at origin
     move_resolved(r,ps.x,ps.y,ps.z,ps.xd,ps.yd,ps.zd, il,ih);
   }
   else { // position and velocity of new particle with respected to mass at index ip
     pt.x       = ps.x + x0;
     pt.y       = ps.y + y0;
     pt.z       = ps.z + z0;
     pt.vx      = ps.xd + vx0;
     pt.vy      = ps.yd + vy0;
     pt.vz      = ps.zd + vz0;
   }
   reb_add(r,pt); // add the new particle
   return om_orb; // return mean motion
}

// return mean motion, set up a single perturber in orbit
// m1 is mass of perturber
// mball is mass of extended body  [0,N)
// use orbital elements to add m1 to particle list
// origin is extended body
double add_one_mass_kep(struct reb_simulation* const r, double m1,
    double aa,double ee, double ii,
    double longnode,double argperi,double meananom, 
    double mball, double r1)
{
   struct reb_particle pt;
   double GM = r->G*(m1+mball); // mball here is extended body
   OrbitalElements orbel;
   orbel.a = aa;
   orbel.e = ee;
   orbel.i = ii;
   orbel.argperi = argperi;
   orbel.longnode= longnode;
   orbel.meananom= meananom;
   PhaseState ps;
   cartesian(GM, orbel, &ps);
   double om_orb = sqrt(GM/aa)/aa;

   pt.m       = m1;
   pt.x       = ps.x;
   pt.y       = ps.y;
   pt.z       = ps.z;
   pt.vx      = ps.xd;
   pt.vy      = ps.yd;
   pt.vz      = ps.zd;
   pt.r       = r1;
   pt.ax = 0; pt.ay = 0; pt.az = 0;
   reb_add(r,pt);
   return om_orb;   
}


//////////////////////////////////////////////
// drift the semi-major axis of  two point masses
//     (a binary) using Beauge et al. 06's
// formulae, Beauge, Michtchenko, Ferraz-Mello 2006, MNRAS 365, 1160
// we can write evolution as dr^2/dt^2 = -v/taua - (v-v_c)/taue
// inv_taua and inv_taue are in units of 1/time 
// inv_taua = (1/a) da/dt, inv_taue = (1/e) de/dt
//  if inv_tau_a >0 drifting apart !!!!!!!!!!!!!!
//
// cc>0 and alpha=0 corresponds to no eccentricity change
// with taua^-1 = 2C(1-alpha), taue^-1 = C*alpha from Beauge et al.
// kk = taua/taue
// eqns of motion can be written
// dv/dt = - C*(1-alpha) v - E(v-vc)
//       = -invtaua*v/2 - invtaue(v-vc)
// 
// This routine drifts a binary!!!!
// The two masses have indexes: im1, im2 
// The drift is done trying to maintain momentum of center of mass of binary
// this routine will probably work for three or more body system
// and will probably work for binaries with various mass ratios
// this routine will not drift a point mass w.r.t. a resolved mass 
//     (comprised of many particles)
/////////////////////////////////////////////////////////////
void dodrift_bin(struct reb_simulation* const r, double tstep,
     double inv_taua, double inv_taue, int im1, int im2)  
{
   struct reb_particle* particles = r->particles;

   // int im1 = r->N -1; // index for primary perturber
   // int im2 = r->N -2; // index for member of binary 
   double m1 = particles[im1].m;
   double m2 = particles[im2].m;
   double GMM = r->G *(m1+m2);
   
      double x  = particles[im2].x - particles[im1].x; 
      double y  = particles[im2].y - particles[im1].y; 
      double z  = particles[im2].z - particles[im1].z;
      double vx = particles[im2].vx -particles[im1].vx; 
      double vy = particles[im2].vy -particles[im1].vy; 
      double vz = particles[im2].vz -particles[im1].vz;
      double rdotv = x*vx + y*vy + z*vz;
      double rad = sqrt(x*x+y*y+z*z);
      double r2 = rad*rad;
      double vc = sqrt(GMM/rad); // circular velocity
      //  rxL is  vector in plane of orbit, perp to r, 
         //  direction of rotation
      // r x v x r = r x -L
      // vector identity axbxc = (adotc)b-(adotb)c
      double rcrossl_x = r2*vx - rdotv*x; 
      double rcrossl_y = r2*vy - rdotv*y; 
      double rcrossl_z = r2*vz - rdotv*z;
      double vl = sqrt(rcrossl_x*rcrossl_x+ rcrossl_y*rcrossl_y+ 
                       rcrossl_z*rcrossl_z); // length of rcrossl
      rcrossl_x /= vl; rcrossl_y /= vl; rcrossl_z /= vl; // unit vector now
      double vcvec_x = vc*rcrossl_x; 
      double vcvec_y = vc*rcrossl_y; 
      double vcvec_z = vc*rcrossl_z;
      // difference between velocity and vc
      double dd_vc_x = vx - vcvec_x;
      double dd_vc_y = vy - vcvec_y;
      double dd_vc_z = vz - vcvec_z;

// compute changes in velocity 
      double dvx =  tstep*(vx*inv_taua/2.0 + dd_vc_x*inv_taue);
      double dvy =  tstep*(vy*inv_taua/2.0 + dd_vc_y*inv_taue);
      double dvz =  tstep*(vz*inv_taua/2.0 + dd_vc_z*inv_taue);

// update velocities , in such a way as to conserve
// momentum of binary, lower mass is affected much more than higher mass
      particles[im1].vx -=  m2*dvx/(m1+m2); 
      particles[im1].vy -=  m2*dvy/(m1+m2);
      particles[im1].vz -=  m2*dvz/(m1+m2);
      particles[im2].vx +=  m1*dvx/(m1+m2); 
      particles[im2].vy +=  m1*dvy/(m1+m2);
      particles[im2].vz +=  m1*dvz/(m1+m2);
    
}


//////////////////////////////////////////////
// just like the previous routine: dodrift_bin
//    but  one mass is a resolved body!  the other is a point mass
/////////////////////////////////////////////////////////////
void dodrift_res(struct reb_simulation* const r, double tstep,
     double inv_taua, double inv_taue, int im1, int il, int ih)
{
   struct reb_particle* particles = r->particles;

   double m1 = particles[im1].m;
   double m2 = sum_mass(r,il,ih);
   double xc =0.0; double yc =0.0; double zc =0.0;
   compute_com(r,il, ih, &xc, &yc, &zc);
   double vxc =0.0; double vyc =0.0; double vzc =0.0;
   compute_cov(r,il, ih, &vxc, &vyc, &vzc);


   double GMM = r->G *(m1+m2);
   
      double x  = xc               - particles[im1].x; 
      double y  = yc               - particles[im1].y; 
      double z  = zc               - particles[im1].z;
      double vx = vxc              - particles[im1].vx; 
      double vy = vyc              - particles[im1].vy; 
      double vz = vzc              - particles[im1].vz;
      double rdotv = x*vx + y*vy + z*vz;
      double rad = sqrt(x*x+y*y+z*z);
      double r2 = rad*rad;
      double vc = sqrt(GMM/rad); // circular velocity
      //  rxL is  vector in plane of orbit, perp to r, 
         //  direction of rotation
      // r x v x r = r x -L
      // vector identity axbxc = (adotc)b-(adotb)c
      double rcrossl_x = r2*vx - rdotv*x; 
      double rcrossl_y = r2*vy - rdotv*y; 
      double rcrossl_z = r2*vz - rdotv*z;
      double vl = sqrt(rcrossl_x*rcrossl_x+ rcrossl_y*rcrossl_y+ 
                       rcrossl_z*rcrossl_z); // length of rcrossl
      rcrossl_x /= vl; rcrossl_y /= vl; rcrossl_z /= vl; // unit vector now
      double vcvec_x = vc*rcrossl_x; 
      double vcvec_y = vc*rcrossl_y; 
      double vcvec_z = vc*rcrossl_z;
      // difference between velocity and vc
      double dd_vc_x = vx - vcvec_x;
      double dd_vc_y = vy - vcvec_y;
      double dd_vc_z = vz - vcvec_z;

// compute changes in velocity 
      double dvx =  tstep*(vx*inv_taua/2.0 + dd_vc_x*inv_taue);
      double dvy =  tstep*(vy*inv_taua/2.0 + dd_vc_y*inv_taue);
      double dvz =  tstep*(vz*inv_taua/2.0 + dd_vc_z*inv_taue);

// update velocities , in such a way as to conserve
// momentum of binary, lower mass is affected much more than higher mass
      particles[im1].vx -=  m2*dvx/(m1+m2); 
      particles[im1].vy -=  m2*dvy/(m1+m2);
      particles[im1].vz -=  m2*dvz/(m1+m2);
      double cvx = m1*dvx/(m1+m2);  // changes to velocity of resolved body
      double cvy = m1*dvy/(m1+m2); 
      double cvz = m1*dvz/(m1+m2); 
      move_resolved(r,0.0,0.0,0.0,cvx,cvy,cvz, il,ih); // change only velocities
    
}

// add a single extra mass  of mass m1
// use cartesian coordinates to add m1 to particle list
void add_one_mass_cartesian(struct reb_simulation* const r, double m1, double r1,
    double x,double y, double z, double vx, double vy, double vz)
{
   struct reb_particle pt;
   pt.m       = m1;
   pt.x       = x;
   pt.y       = y;
   pt.z       = z;
   pt.vx      = vx;
   pt.vy      = vy;
   pt.vz      = vz;
   pt.r       = r1;
   pt.ax = 0; pt.ay = 0; pt.az = 0;
   reb_add(r,pt);
}
