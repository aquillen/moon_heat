
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


extern int NS; // number of springs
extern int NPERT;  // number of external perturbing bodies
struct stresstensor *stressvec; // global so can be reached by all routines here

void spring_force_one();
int markfailure();

// convension should be tensile stress is negative
// update the stress tensor, create it if does not exist
// assuming total volume = 4/3pi (used to normalize) is a constant
void update_stresstensor(struct reb_simulation* const r)
{
   static int first = 0;
   if (first==0){
      first=1;
      stressvec = malloc(r->N*sizeof(struct stresstensor));
   }
   for(int i=0;i<r->N;i++){ // over nodes
      stressvec[i].sigxx=0.0;
      stressvec[i].sigyy=0.0;
      stressvec[i].sigzz=0.0;
      stressvec[i].sigxy=0.0;
      stressvec[i].sigyz=0.0;
      stressvec[i].sigxz=0.0;
      stressvec[i].eig1=0.0;
      stressvec[i].eig2=0.0;
      stressvec[i].eig3=0.0;
      stressvec[i].maxF=0.0;
      stressvec[i].s_index=-1;
      stressvec[i].fail=0;
   }

   for(int k=0;k<r->NS;k++){ // over springs
      int ii = springs[k].i; 
      int jj = springs[k].j;
      double Fx,Fy,Fz,Lx,Ly,Lz;
      spring_force_one(r,k, &Fx,&Fy,&Fz,&Lx,&Ly,&Lz); // compute spring force 
      // sum over springs onto each node
      stressvec[ii].sigxx += Fx*Lx; stressvec[jj].sigxx -= Fx*Lx;
      stressvec[ii].sigyy += Fy*Ly; stressvec[jj].sigyy -= Fy*Ly;
      stressvec[ii].sigzz += Fz*Lz; stressvec[jj].sigzz -= Fz*Lz;
      stressvec[ii].sigxy += Fx*Ly; stressvec[jj].sigxy -= Fx*Ly;
      stressvec[ii].sigyz += Fy*Lz; stressvec[jj].sigyz -= Fy*Lz;
      stressvec[ii].sigxz += Fx*Lz; stressvec[jj].sigxz -= Fx*Lz;
      // signs should give negative stress for tensile stress
      double Fmag = sqrt(Fx*Fx + Fy*Fy + Fz*Fz);
      if (Fmag > stressvec[ii].maxF) {
          stressvec[ii].maxF = Fmag;   // keep track of which spring has max force
          stressvec[ii].s_index = k;  
      }
      if (Fmag > stressvec[jj].maxF) {
          stressvec[jj].maxF = Fmag;  
          stressvec[jj].s_index = k;  
      }
   }

   double vol = 4.0*M_PI/(3.0*(double)r->N);   // volume of a node depends on number of particles
   // we have R=1 for our masses, so the volume is 4pi/3 in units of R^3
   for(int i=0;i<r->N;i++){ // normalize
      stressvec[i].sigxx /=vol;
      stressvec[i].sigyy /=vol;
      stressvec[i].sigzz /=vol;
      stressvec[i].sigxy /=vol;
      stressvec[i].sigyz /=vol;
      stressvec[i].sigxz /=vol;
   }
   for(int i=0;i<r->N;i++){ // compute and store eigenvalues!
      double eig1,eig2,eig3;
      eigenvalues(stressvec[i].sigxx,stressvec[i].sigyy,stressvec[i].sigzz,
                  stressvec[i].sigxy,stressvec[i].sigyz,stressvec[i].sigxz,
                  &eig1, &eig2, &eig3); // eig>=eig2>=eig3
      stressvec[i].eig1 =eig1;
      stressvec[i].eig2 =eig2;
      stressvec[i].eig3 =eig3;
   }
}

// check for failure at all nodes 
// sigt_I is tensile strength of interior (these are negative)
// sigt_S is tensile strength of shell
// we need to know if the particle was originally in shell or interior
// we can't use current radius as particles are going to move! 
// decide if particle is in interior by whether particle mass is below pmass_div
// return number of failed nodes
int markfailure(struct reb_simulation* const r, int npert, 
       double pmass_div, double sigt_I, double sigt_S)
{
    int il = 0;
    int ih = r->N - npert;
    int nfail = 0;
    for (int i=il;i<ih;i++){
       double m = r->particles[i].m;
       double sigt; 
       // if particle is in interior by whether particle mass is below pmass_div
       if (m > pmass_div) sigt =  sigt_I; // in that case this is failure tensile stress
       else sigt = sigt_S;

       double tau1 = stressvec[i].eig1;
       double tau3 = stressvec[i].eig3;
       int fail=0;
       double tau13 = tau1 - tau3;
       if ((tau1  <  -3.0*tau3) && (tau3< sigt)) fail = 1;
       if ((tau1  >= -3.0*tau3) && (tau13*tau13 + 8.0*sigt*tau13 >0)) fail = 1;
       if (fail==1){
          stressvec[i].fail=1;
          nfail++;
       }
    }
    if (nfail>0)
       printf("markfailure: nfail=%d\n",nfail);

    return nfail;
}


// here is the hard part
void killsprings(struct reb_simulation* const r, int npert) {
    int il = 0;
    int ih = r->N - npert;
    for (int i=il;i<ih;i++){
       if (stressvec[i].fail==1){
          int s_index = stressvec[i].s_index;  // spring with max force amplitude on it
          springs[s_index].ks=0.0;   // spring dead
          springs[s_index].gamma=0.0;  
       }
    }
}

#define L_EPS 1e-6; // softening for spring length
// return the spring force vector for spring i
// to get force on particle j multiply by mass of particle j
void spring_force_one(struct reb_simulation* const r,int i, 
   double *Fx, double *Fy, double *Fz,
   double *Lx, double *Ly, double *Lz){

     double L = spring_length(r,springs[i]) + L_EPS; // spring length
     double rs0 = springs[i].rs0; // rest length
     int ii = springs[i].i; 
     int jj = springs[i].j;
     double dx = r->particles[ii].x - r->particles[jj].x;
     double dy = r->particles[ii].y - r->particles[jj].y;
     double dz = r->particles[ii].z - r->particles[jj].z;
     double mii = r->particles[ii].m;
     double mjj = r->particles[jj].m;
     double ks = springs[i].ks;
     double fac = -ks*(L-rs0)/L; // L here to normalize direction
            // accelerations are force divided by mass
     *Fx =  fac*dx;
     *Fy =  fac*dy;
     *Fz =  fac*dz;
     *Lx = dx; *Ly = dy; *Lz = dz;

     if (springs[i].gamma>0.0) {  // damping force too!
          double dvx = r->particles[ii].vx - r->particles[jj].vx;
          double dvy = r->particles[ii].vy - r->particles[jj].vy;
          double dvz = r->particles[ii].vz - r->particles[jj].vz;
          double dLdt = (dx*dvx + dy*dvy + dz*dvz)/L;
                // divide dL/dt by L to get strain rate
          double mbar = mii*mjj/(mii+mjj); // reduced mass
          double dampfac  = springs[i].gamma*mbar*dLdt/L;
               // factor L here to normalize dx,dy,dz
          *Fx -=  dampfac*dx;
          *Fy -=  dampfac*dy;
          *Fz -=  dampfac*dz;
     }
}


// make a stress file name depending on numbers of tp (toprint)
void sfilename(struct reb_simulation* const r,char *root, double tp, char *fname){
   int xd = (int)(r->t/tp);
   char junks[20];
   sprintf(junks,"%d",xd);
   sprintf(fname,"%s_",root);
   if (xd < 100000) strcat(fname,"0");
   if (xd < 10000)  strcat(fname,"0");
   if (xd < 1000)   strcat(fname,"0");
   if (xd < 100)    strcat(fname,"0");
   if (xd < 10)     strcat(fname,"0");
   strcat(fname,junks);
   strcat(fname,"_stress.txt");
}


// print out stress file
// npert is number of point mass perturbers
void print_stress(struct reb_simulation* const r, int npert,
      char* filename)
{
   struct reb_particle* particles = r->particles;

   FILE *fpo;
   fpo = fopen(filename, "w");
   fprintf(fpo,"#%.2e\n",r->t);
   int il =0;
   int ih =r->N-npert; // NPERT?
   double xc,yc,zc;
   compute_com(r,il, ih, &xc, &yc, &zc);
   for(int i=il;i<ih;i++){
      double x = particles[i].x - xc;
      double y = particles[i].y - yc;
      double z = particles[i].z - zc;
      double m = particles[i].m;
      fprintf(fpo,"%d %.4f %.3f %.3f %.3f " ,i,m,x,y,z);
      // double sigxx = stressvec[i].sigxx;
      // double sigyy = stressvec[i].sigyy;
      // double sigzz = stressvec[i].sigzz;
      // double sigxy = stressvec[i].sigxy;
      // double sigyz = stressvec[i].sigyz;
      // double sigxz = stressvec[i].sigxz;
      double eig1  = stressvec[i].eig1;
      double eig2  = stressvec[i].eig2;
      double eig3  = stressvec[i].eig3;
      fprintf(fpo," %.3f %.3f %.3f " ,eig1,eig2,eig3);
      fprintf(fpo," %d " ,stressvec[i].fail);
      fprintf(fpo,"\n");
   }

}

