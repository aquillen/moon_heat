
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h"
#include "tools.h"
#include "output.h"
#include "spring.h"


extern int NS;

// printout out positions and velocities of only surface particles
void print_surf(struct reb_simulation* const r, int il, int ih, int *surfarr, char* filename){
   FILE *fpo;
   static int tstart =0;
   static double toff = 0.0;
   if (tstart==0){
     tstart=1;
     toff = r->t;
     printf("print_surf toff=%.6f\n",toff);
   }
   fpo = fopen(filename, "w");
   fprintf(fpo,"# t=%.8f #txyzvxvyvzaxayazrm\n",r->t); // actual time
   for(int i=0;i< r->N;i++){
      if (surfarr[i] ==1)
        fprintf(fpo,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6e\n",
          r->t - toff, 
          r->particles[i].x, r->particles[i].y, r->particles[i].z,
          r->particles[i].vx, r->particles[i].vy, r->particles[i].vz,
          r->particles[i].ax, r->particles[i].ay, r->particles[i].az,
          r->particles[i].r, r->particles[i].m); 
   }
   fclose(fpo);
   printf("\n print_surf: N=%d %s\n",r->N,filename);
}



// print out information for an extended body indices [il,ih)
void print_extended(struct reb_simulation* const r, int il, int ih, char* filename)
{
   static int first=0;
   FILE *fpo;
   if (first==0){
     first=1;
     fpo = fopen(filename, "w");
     fprintf(fpo,"# t x y z vx vy vz omx omy omz llx lly llz Ixx Iyy Izz Ixy Iyz Ixz KErot PEspr PEgrav Etot\n");
   }
   else {
     fpo = fopen(filename, "a");
   }
   fprintf(fpo,"%.3f ",r->t);
   double xc =0.0; double yc =0.0; double zc =0.0;
   compute_com(r,il, ih, &xc, &yc, &zc);
   double vxc =0.0; double vyc =0.0; double vzc =0.0;
   compute_cov(r,il, ih, &vxc, &vyc, &vzc);
   fprintf(fpo,"%.5f %.5f %.5f ",xc,yc,zc);
   fprintf(fpo,"%.5f %.5f %.5f ",vxc,vyc,vzc);

   double omx,omy,omz, Ibig, Imid, Ismall, llx,lly,llz;
   // computes omega using inverse of moment of inertia matrix and angular momentum
   body_spin(r, il, ih, &omx, &omy, &omz, &Ibig, &Imid, &Ismall);
   measure_L(r,il, ih, &llx, &lly, &llz); // computes spin angular momentum of spining body
   fprintf(fpo,"%.4e %.4e %.4e ",omx, omy, omz);
   fprintf(fpo,"%.4e %.4e %.4e ",llx,lly,llz);
   // fprintf(fpo,"%.4e %.4e %.4e ",Ibig, Imid, Ismall);
   double Ixx,Iyy,Izz,Ixy,Iyz,Ixz;
   mom_inertia(r,il,ih, &Ixx, &Iyy, &Izz,&Ixy, &Iyz, &Ixz);
   fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",Ixx,Iyy,Izz,Ixy,Iyz,Ixz);

   double KErot = compute_rot_kin(r,il,ih); // rotational kinetic energy
   double pe_springs = spring_potential_energy(r); // potential energy in springs
   double pe_grav = grav_potential_energy(r,il,ih); // potential energy in gravity 
   double E_tot = KErot+pe_springs+pe_grav; 
   fprintf(fpo,"%.5e %.5e %.10e %.10e",KErot,pe_springs,pe_grav,E_tot);
   // printf("%.5e %.5e %.5e %.5e\n",KErot,pe_springs,pe_grav,E_tot); //  
   fprintf(fpo,"\n");
   fclose(fpo);
}

// print out information for an extended body indices [il,ih)
void print_extended_simp(struct reb_simulation* const r, int il, int ih, char* filename)
{
   static int first=0;
   FILE *fpo;
   if (first==0){
     first=1;
     fpo = fopen(filename, "w");
     fprintf(fpo,"# t x y z vx vy vz omx omy omz llx lly llz Ixx Iyy Izz Ixy Iyz Ixz \n");
   }
   else {
     fpo = fopen(filename, "a");
   }
   fprintf(fpo,"%.3f ",r->t);
   double xc =0.0; double yc =0.0; double zc =0.0;
   compute_com(r,il, ih, &xc, &yc, &zc);
   double vxc =0.0; double vyc =0.0; double vzc =0.0;
   compute_cov(r,il, ih, &vxc, &vyc, &vzc);
   fprintf(fpo,"%.5f %.5f %.5f ",xc,yc,zc);
   fprintf(fpo,"%.5f %.5f %.5f ",vxc,vyc,vzc);

   double omx,omy,omz, Ibig, Imid, Ismall, llx,lly,llz;
   // computes omega using inverse of moment of inertia matrix and angular momentum
   body_spin(r, il, ih, &omx, &omy, &omz, &Ibig, &Imid, &Ismall);
   measure_L(r,il, ih, &llx, &lly, &llz); // computes spin angular momentum of spining body
   fprintf(fpo,"%.4e %.4e %.4e ",omx, omy, omz);
   fprintf(fpo,"%.4e %.4e %.4e ",llx,lly,llz);
   // fprintf(fpo,"%.4e %.4e %.4e ",Ibig, Imid, Ismall);
   double Ixx,Iyy,Izz,Ixy,Iyz,Ixz;
   mom_inertia(r,il,ih, &Ixx, &Iyy, &Izz,&Ixy, &Iyz, &Ixz);
   fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",Ixx,Iyy,Izz,Ixy,Iyz,Ixz);

   fprintf(fpo,"\n");
   fclose(fpo);
}

// print out information for two specific nodes in an extended body indices [il,ih)
// the first node has largest x value at initial time
// the second node has largest x value at initial time
// same nodes are printed each time after they are chosen
void print_extended_2nodes(struct reb_simulation* const r, int il, int ih, 
    // struct node* nodevec, 
    char* filename)
{
   static int first=0;
   static int iz=0; // index of particle with largest z value
   static int ix=0; // index of particle with largest x value
   FILE *fpo;
   if (first==0){
     first=1;
     fpo = fopen(filename, "w");
     fprintf(fpo,"#t xyzvxvyvz(xnode) xyzvxvyvz(znode) \n");
     iz = nearest_to_shape(r,il,ih, 0.0,0.0,10.0);
     ix = nearest_to_shape(r,il,ih,10.0,0.0, 0.0);
   }
   else {
     fpo = fopen(filename, "a");
   }
   fprintf(fpo,"%.3f ",r->t);
   fprintf(fpo,"%.5f %.5f %.5f ",
        r->particles[ix].x,r->particles[ix].y,r->particles[ix].z);
   fprintf(fpo,"%.5f %.5f %.5f ",
        r->particles[ix].vx,r->particles[ix].vy,r->particles[ix].vz);
   fprintf(fpo,"%.5f %.5f %.5f ",
        r->particles[iz].x,r->particles[iz].y,r->particles[iz].z);
   fprintf(fpo,"%.5f %.5f %.5f ",
        r->particles[iz].vx,r->particles[iz].vy,r->particles[iz].vz);
   fprintf(fpo,"\n");
   fclose(fpo);
}

#define NPMAXE 100
// print out information for a point mass with index ip
// which point mass is ipert  (0 being the central one)
// but ip is its index in the particle list
void print_pm(struct reb_simulation* const r, int ip, int ipert, char* filename)
{
   struct reb_particle* particles = r->particles;
   static int firstarr[NPMAXE];
   static int first=0;
   if (first==0){
      first=1;
      for (int i=0;i<NPMAXE;i++) firstarr[i]=0;
   }
   FILE *fpo;
   if (firstarr[ipert]==0){  
     firstarr[ipert]=1;
     fpo = fopen(filename, "w");                     
     fprintf(fpo,"# t x y z vx vy vz m\n");
     printf("%s w openned\n",filename);
   }
   else {
     fpo = fopen(filename, "a");
   }
   double x = particles[ip].x;
   double y = particles[ip].y;
   double z = particles[ip].z;
   double vx = particles[ip].vx;
   double vy = particles[ip].vy;
   double vz = particles[ip].vz;
   fprintf(fpo,"%.3f ",r->t);
   fprintf(fpo,"%.5f %.5f %.5f ",x,y,z);
   fprintf(fpo,"%.5f %.5f %.5f ",vx,vy,vz);
   fprintf(fpo,"%.3e ",particles[ip].m);
   fprintf(fpo,"\n");
   fclose(fpo);
}

void print_tab(struct reb_simulation* const r, int npert, char* filename)
{
   // struct reb_particle* particles = r->particles;
   static int first=0;
   FILE *fpo;

   if (first==0){
     first=1;
     fpo = fopen(filename, "w");
     fprintf(fpo,"# t a n e i omx omy omz A B C E lx ly lz ang lox loy loz Ixx Iyy Izz Ixy Iyz Ixz\n");
   }
   else {
     fpo = fopen(filename, "a");
   }
   int il = 0;
   int ih = r->N - npert;
   double a,meanmo,ecc,incl,Lorb;
   if (npert ==1) {
      compute_semi(r, il, ih, r->N-1, &a, &meanmo, &ecc, &incl, &Lorb);
   }
   else
      compute_semi_bin(r, il, ih, npert, &a, &meanmo, &ecc, &incl, &Lorb);
   fprintf(fpo,"%.3f %.6e %.4e %.4f %.4e ",r->t, a, meanmo, ecc, incl);
   double omx,omy,omz, Ibig, Imid, Ismall;
   // computes sping using inverse of moment of inertia matrix and angular momentum vec
   body_spin(r, il, ih, &omx, &omy, &omz, &Ibig, &Imid, &Ismall);
   fprintf(fpo,"%.4e %.4e %.4e %.3f %.3f %.3f ",omx, omy, omz, Ibig, Imid, Ismall);
   double E = Young_mush(r, il, ih, 0.0, 0.5);
   fprintf(fpo,"%.3f ",E);
   double llx,lly,llz;
   measure_L(r,il, ih, &llx, &lly, &llz); // spin angular momentum of spining body
   fprintf(fpo,"%.5f %.5f %.5f ",llx,lly,llz);
   double llx_o,lly_o,llz_o;
   compute_Lorb(r, il,ih, npert, &llx_o, &lly_o, &llz_o);  // total orbital angular momentum
   double L_spin = sqrt(llx*llx + lly*lly + llz*llz);
   double L_orb = sqrt(llx_o*llx_o + lly_o*lly_o + llz_o*llz_o);
   double LsdotLo = llx*llx_o + lly*lly_o + llz*llz_o; 
   double angle = acos(LsdotLo/(L_spin*L_orb));  // angular difference  between 
            // spin angular momentu and orbital angular momentum
   fprintf(fpo,"%.5f ",angle);
   fprintf(fpo,"%.5f %.5f %.5f ",llx_o,lly_o,llz_o); // orbital angular momentum
   double Ixx,Iyy,Izz,Ixy,Iyz,Ixz;
   mom_inertia(r,il,ih, &Ixx, &Iyy, &Izz,&Ixy, &Iyz, &Ixz);
   fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",Ixx,Iyy,Izz,Ixy,Iyz,Ixz);


   fprintf(fpo,"\n");
   fclose(fpo);
}

void print_bin(struct reb_simulation* const r, int npert, char* filename)
{
   struct reb_particle* particles = r->particles;
   static int mark =-1;
   static int first=0;
   FILE *fpo;
   if (first==0){
     first=1;
     fpo = fopen(filename, "w");
     fprintf(fpo,"# t xyz,vxyz(res) xyz,vxyz(PT) xyz,vxyz(PC)  ");
     fprintf(fpo," lx ly lz Ixx Iyy Izz Ixy Iyz Ixz\n  ");
     double r2max = 0.0;
     for(int i=0;i<r->N-npert;i++){
       double rc2 = particles[i].x*particles[i].x + particles[i].y*particles[i].y
                  + particles[i].z*particles[i].z;
       if (rc2 > r2max){   // mark particle with largest initial r2 value 
         mark = i;
         r2max = rc2;
       }
     }
   }
   else {
     fpo = fopen(filename, "a");
   }

   fprintf(fpo,"%.3f ",r->t);

   int il = 0; // index range for resolved body
   int ih = r->N - npert;
   double xc =0.0; double yc =0.0; double zc =0.0;
   double vxc =0.0; double vyc =0.0; double vzc =0.0;
   compute_com(r,il, ih, &xc, &yc, &zc); // center of mass of resolved body
   compute_cov(r,il, ih, &vxc, &vyc, &vzc); // center of velocity of resolved body
   int iml = r->N -npert; // index range for perturbing masses, binary
   int imh = r->N; 
   double xb=0.0; double yb=0.0; double zb=0.0; 
   double vxb=0.0; double vyb=0.0; double vzb = 0.0;
   compute_com(r,iml, imh, &xb, &yb, &zb); // center of mass of binary 
   compute_cov(r,iml, imh, &vxb, &vyb, &vzb); // center of velocity of binary  

   double xr = xc-xb; double vxr = vxc-vxb;  // of resolved w.r.t binary center
   double yr = yc-yb; double vyr = vyc-vyb;
   double zr = zc-zb; double vzr = vzc-vzb;
   // store resolved center vs bin center
   fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",xr,yr,zr,vxr,vyr,vzr);

   double xs = particles[mark].x - xc;   // marked w.r.t to resolved center
   double ys = particles[mark].y - yc; 
   double zs = particles[mark].z - zc; 
   double vxs = particles[mark].vx - vxc; 
   double vys = particles[mark].vy - vyc; 
   double vzs = particles[mark].vz - vzc; 
   // store marked w.r.t. to resolved center 
   fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",xs,ys,zs,vxs,vys,vzs);

   if (npert==2){
     double x_CP  = particles[iml].x  - particles[iml+1].x;  // charon w.r.t pluto
     double y_CP  = particles[iml].y  - particles[iml+1].y; 
     double z_CP  = particles[iml].z  - particles[iml+1].z; 
     double vx_CP = particles[iml].vx - particles[iml+1].vx; 
     double vy_CP = particles[iml].vy - particles[iml+1].vy; 
     double vz_CP = particles[iml].vz - particles[iml+1].vz; 
     // store Charon w.r.t Pluto
     fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",x_CP,y_CP,z_CP,vx_CP,vy_CP,vz_CP);
   }
   else {
     fprintf(fpo,"0 0 0 0 0 0 ");
   }

   double llx,lly,llz;
   measure_L(r,il, ih, &llx, &lly, &llz); // spin angular momentum of spining body
   fprintf(fpo,"%.5f %.5f %.5f ",llx,lly,llz);

   //compute moment of inertia matrix
   double Ixx,Iyy,Izz,Ixy,Iyz,Ixz;
   mom_inertia(r,il,ih, &Ixx, &Iyy, &Izz,&Ixy, &Iyz, &Ixz);
   fprintf(fpo,"%.5f %.5f %.5f %.5f %.5f %.5f ",Ixx,Iyy,Izz,Ixy,Iyz,Ixz);

   fprintf(fpo,"\n");
   fclose(fpo);
   //  xxxxxx
}


// write out springs to a file
void write_springs(struct reb_simulation* const r,char *fileroot, int index){
   FILE *fpo;
   char filename[100];
   char istring[100]; 
   toistring(istring, index);
   strcpy(filename,fileroot);
   strcat(filename,"_");
   strcat(filename,istring);
   strcat(filename,"_springs.txt");
   fpo = fopen(filename,"w");
   for(int i=0;i<NS;i++){
      double dr = spring_length(r,springs[i]);
      fprintf(fpo,"%d %d %.6f %.6f %.6f %.6f ",
        springs[i].i, springs[i].j, 
        springs[i].ks, springs[i].rs0, 
        springs[i].gamma, springs[i].k_heat);
      fprintf(fpo,"%.6f %.6f\n",
         strain(r,springs[i]), dr);
   }
   fclose(fpo);
   printf("\n write_springs: NS=%d %s\n",NS,filename);
}


void read_springs(struct reb_simulation* const r,char *fileroot, int index){
   char filename[100];
   char istring[100]; 
   toistring(istring, index);
   strcpy(filename,fileroot);
   strcat(filename,"_");
   strcat(filename,istring);
   strcat(filename,"_springs.txt");
   printf("\n reading in springs %s\n",filename);
   FILE *fpi;
   fpi = fopen(filename,"r");
   char string[300];
   struct spring spr;
   int  i,j; 
   double ks,rs0,gamma,k_heat;
   while(fgets(string,300,fpi) != NULL){
      sscanf(string,"%d %d %lf %lf %lf %lf",
        &i,&j,&ks,&rs0,&gamma,&k_heat);
      spr.i=i; spr.j=j;
      spr.ks=ks; spr.rs0=rs0; spr.gamma=gamma; spr.k_heat=k_heat; 
      springs_add(r,spr);
   }
   fclose(fpi);
   printf("read_springs: NS=%d\n",NS);
}


void read_particles(struct reb_simulation* const r,char *fileroot, int index){
   char filename[100];
   char istring[100]; 
   toistring(istring, index);
   strcpy(filename,fileroot);
   strcat(filename,"_");
   strcat(filename,istring);
   strcat(filename,"_particles.txt");
   printf("\n reading in particles %s\n",filename);
   FILE *fpi;
   fpi = fopen(filename,"r");
   char string[300];
   struct reb_particle pt;
   pt.ax = 0.0; pt.ay = 0.0; pt.az = 0.0;
   double x,y,z,vx,vy,vz,m,rad;
   while(fgets(string,300,fpi) != NULL){
      sscanf(string,"%lf %lf %lf %lf %lf %lf %lf %lf\n",
        &x, &y, &z, &vx, &vy, &vz, &rad, &m);
      pt.x = x;   pt.y = y;   pt.z = z;
      pt.vx = vx; pt.vy = vy; pt.vz = vz;
      pt.r = rad; pt.m = m; 
      reb_add(r,pt);
   }
   fclose(fpi);
   printf("read_particles: N=%d\n",r->N);
}

void write_particles(struct reb_simulation* const r,char *fileroot, int index){
   FILE *fpo;
   char filename[100];
   char istring[100]; 
   toistring(istring, index);
   strcpy(filename,fileroot);
   strcat(filename,"_");
   strcat(filename,istring);
   strcat(filename,"_particles.txt");
   fpo = fopen(filename,"w");
   for(int i=0;i< r->N;i++){
      fprintf(fpo,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6e\n",
        r->particles[i].x, r->particles[i].y, r->particles[i].z,
        r->particles[i].vx, r->particles[i].vy, r->particles[i].vz,
        r->particles[i].r, r->particles[i].m); 
   }
   fclose(fpo);
   printf("\n write_particles: N=%d %s\n",r->N,filename);
}

// give me an integer string
void toistring(char *istring, int i){
   char junks[20];
   sprintf(junks,"%d",i);
   if (i < 100000) strcpy(istring,"0");
   if (i < 10000)  strcat(istring,"0");
   if (i < 1000)   strcat(istring,"0");
   if (i < 100)    strcat(istring,"0");
   if (i < 10)     strcat(istring,"0");
   strcat(istring,junks);
}

