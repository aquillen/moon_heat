/**
 * resolved mass spring model
 * using the leap frog integrator. 
 */
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


int NS; 
struct spring* springs;
void reb_springs();// to pass springs to display
#define NPMAX 10  // maximum number of point masses
double itaua[NPMAX],itaue[NPMAX];
int *surfp; // is allocated with marksurface routine!
struct node *nodevec;

double gamma_fac; // for adjustment of gamma of all springs
double t_damp;    // end faster damping, relaxation
double t_print;   // for table printout 
double t_heat;    // for heat  printout 
char froot[30];   // output files
int npert; // number of point mass perturbers
double powerfac; // fraction of heat to center of spring rather than nodes
double dote_radg;  // radiogenic heating rate  
                        // energy per unit mass per unit time
void print_run_double();

int icentral=-1; // central mass location

void heartbeat(struct reb_simulation* const r);

void additional_forces(struct reb_simulation* r){
   spring_forces(r); // spring forces
}

struct spring spring_mush_hot; // spring parameters for mush
struct spring spring_mush_cold; // spring parameters for mush
        double Ttrans;  // Transition temperature 

int main(int argc, char* argv[]){
	struct reb_simulation* const r = reb_create_simulation();
	// Setup constants
	r->integrator	= REB_INTEGRATOR_LEAPFROG;
	r->gravity	= REB_GRAVITY_BASIC;
	r->boundary	= REB_BOUNDARY_NONE;
	r->G 		= 1;		
        r->additional_forces = additional_forces;  // setup callback function for additional forces
        double mball = 1.0;          // total mass of ball
        double rball = 1.0;          // radius of a ball
        double tmax = 0.0;  // if 0 integrate forever

// things to set! ////////////////////// could be read in with parameter file
        double dt, b_distance,omegaz,mush_fac,surfdist;
        double ratio1,ratio2,obliq_deg;
        double Tsurf=0.0;   // surface temperature
        double Tinit =0.0;  // initial interior node temperature
        dote_radg = 0.0;   // radiogenic heating rate
        int lattice_type;
        double rad[NPMAX],mp[NPMAX]; // orbit stuff
        double aa[NPMAX],ee[NPMAX],ii[NPMAX];
        double longnode[NPMAX],argperi[NPMAX],meananom[NPMAX];
        int npointmass=0;
        npert=0;
        double k_heat_hot=1.0;     // thermal transport coeffs
        double k_heat_cold = 0.0;  
        double gamma_hot= 0.0;     // damping parm
        double gamma_cold= 0.0;
        double ks_hot = 0.0;       // spring constant
        double ks_cold = 0.0;
        double cp_hot=1.0;         // heat capacity
        double cp_cold = 1.0;
        Ttrans=0.0;                // transition temperature two state model

    if (argc ==1){
        strcpy(froot,"t1");   // to make output files
	dt	   = 1e-3;    // Timestep
	lattice_type       =0;        // 0=rand 1=hcp
        b_distance = 0.15;    // for creating random sphere, min separation between particles
        mush_fac    = 2.3;    // ratio of smallest spring distance to minimum interparticle dist
        omegaz     = 0.2;     // initial spin
        // spring damping
        gamma_fac   = 1.0;    // initial factor for initial damping value for springs
        gamma_hot   = 1.0;    // final damping coeff
        t_damp      = 1.0;    // gamma to final values for all springs at this time
        ks_hot      = 8e-2;   // spring constant

        ratio1 =1.0; // shape of resolved body  y/x b/a
        ratio2 =0.95; // z/x c/a
        t_print =  1.0;  // printouts for table
        t_heat =  10000.0;  // heat printouts 
        powerfac = 1.0; // fraction of heat to center of springs
        obliq_deg=0.0; // obliquity in degrees
        surfdist=0.1; // for identifying surface particles
     }
     else{
        FILE *fpi;
        fpi = fopen(argv[1],"r");
        char line[300];
        fgets(line,300,fpi);  sscanf(line,"%s",froot);
        fgets(line,300,fpi);  sscanf(line,"%lf",&dt);
        fgets(line,300,fpi);  sscanf(line,"%lf",&tmax);
        fgets(line,300,fpi);  sscanf(line,"%lf",&t_print);
        fgets(line,300,fpi);  sscanf(line,"%lf",&t_heat);
        fgets(line,300,fpi);  sscanf(line,"%lf",&powerfac);
        fgets(line,300,fpi);  sscanf(line,"%d" ,&lattice_type);
        fgets(line,300,fpi);  sscanf(line,"%lf",&b_distance);
        fgets(line,300,fpi);  sscanf(line,"%lf",&mush_fac);
        fgets(line,300,fpi);  sscanf(line,"%lf",&ratio1);
        fgets(line,300,fpi);  sscanf(line,"%lf",&ratio2);
        fgets(line,300,fpi);  sscanf(line,"%lf",&omegaz);
        fgets(line,300,fpi);  sscanf(line,"%lf",&obliq_deg);
        fgets(line,300,fpi);  sscanf(line,"%lf",&gamma_fac);
        fgets(line,300,fpi);  sscanf(line,"%lf",&t_damp);
        fgets(line,300,fpi);  sscanf(line,"%lf",&surfdist);
        

        fgets(line,300,fpi);  sscanf(line,"%lf",&ks_hot);
        fgets(line,300,fpi);  sscanf(line,"%lf",&ks_cold);
        fgets(line,300,fpi);  sscanf(line,"%lf",&gamma_hot);
        fgets(line,300,fpi);  sscanf(line,"%lf",&gamma_cold);
        fgets(line,300,fpi);  sscanf(line,"%lf",&cp_hot);
        fgets(line,300,fpi);  sscanf(line,"%lf",&cp_cold);
        fgets(line,300,fpi);  sscanf(line,"%lf",&k_heat_hot);
        fgets(line,300,fpi);  sscanf(line,"%lf",&k_heat_cold);
        fgets(line,300,fpi);  sscanf(line,"%lf",&Ttrans);
        fgets(line,300,fpi);  sscanf(line,"%lf",&Tinit);
        fgets(line,300,fpi);  sscanf(line,"%lf",&dote_radg);

        fgets(line,300,fpi);  sscanf(line,"%d",&npointmass);
        for (int ip=0;ip<npointmass;ip++){
           fgets(line,300,fpi);  sscanf(line,"%lf %lf %lf %lf",
             mp+ip,rad+ip,itaua+ip,itaue+ip);
           fgets(line,300,fpi);  sscanf(line,"%lf %lf %lf %lf %lf %lf",
             aa+ip,ee+ip,ii+ip,longnode+ip,argperi+ip,meananom+ip);
        }
     printf("parm file read in\n");

     }
     double obliquity = obliq_deg*M_PI/180.0;
     if (powerfac >1.0) powerfac = 1.0;
     if (powerfac <0.0) powerfac = 1.0;
     

/// end of things to set /////////////////////////

        r->dt=dt; // set integration timestep
	const double boxsize = 3.2*rball;    // display
	reb_configure_box(r,boxsize,1,1,1);
	r->softening      = b_distance/100.0;	// Gravitational softening length


   // properties of springs
   spring_mush_hot.gamma          = gamma_fac*gamma_hot; // initial damping coefficient
   spring_mush_cold.gamma         = gamma_cold; // damping coefficient
   spring_mush_hot.ks             = ks_hot; // spring constant
   spring_mush_cold.ks            = ks_cold; 
   spring_mush_hot.k_heat         = k_heat_hot; // heat transport coefficient
   spring_mush_cold.k_heat        = k_heat_cold;
   double mush_distance=b_distance*mush_fac; 
       // distance for connecting and reconnecting springs

   FILE *fpr;
   char fname[200];
   sprintf(fname,"%s_run.txt",froot);
   fpr = fopen(fname,"w");

   NS=0; // start with no springs 

// do you want volume to be the same? yes, adjusting here!!!!
   double volume_ratio = pow(rball,3.0)*ratio1*ratio2;  // neglecting 4pi/3 factor
   double vol_radius = pow(volume_ratio,1.0/3.0);

   rball /= vol_radius; // volume radius used to compute semi-major axis
// assuming that body semi-major axis is rball
   fprintf(fpr,"a %.3f\n",rball); 
   fprintf(fpr,"b %.3f\n",rball*ratio1); 
   fprintf(fpr,"c %.3f\n",rball*ratio2); 
   volume_ratio = pow(rball,3.0)*ratio1*ratio2;  // neglecting 4pi/3 factor
   fprintf(fpr,"vol_ratio %.6f\n",volume_ratio); // with respect to 4pi/3 
   // printf("vol_ratio %.6f\n",volume_ratio); // with respect to 4pi/3 
   // so I can check that it is set to 1
   // working in units of volume equivalent sphere

   // create particle distribution
   if (lattice_type==0){
      rand_football(r,b_distance,rball,rball*ratio1, rball*ratio2,mball );
      surfp = marksurface_football(r, surfdist, rball,rball*ratio1,rball*ratio2);
      nodevec = mknodevec(r,cp_hot,Tinit); // temperatures on nodes
      surface_nodes(surfp, nodevec, r->N,Tsurf); // mark surface on nodevec
                 // also set temperature on surface nodes to Tsurf
   }
   else{
      exit(0);
   }

   int il=0;
   int ih=r->N;
   centerbody(r,il,ih);  // move reference frame to resolved body 
   subtractcov(r,il,ih);  // subtract center of velocity
   // spin it
   spin(r,il, ih, 0.0, 0.0, omegaz);  // you can change one of these to tilt!
   subtractcov(r,il,ih);  // subtract center of velocity
   double speriod  = fabs(2.0*M_PI/omegaz);
   print_run_double(speriod, "spin period", fpr);
   print_run_double(omegaz, "omegaz", fpr);
   if (obliquity != 0.0)
     rotate_body(r, il, ih, 0.0, obliquity, 0.0); // tilt by obliquity in radians

   // make springs, all pairs connected within interparticle distance mush_distance
   connect_springs_dist(r,mush_distance, 0, r->N, spring_mush_hot);
   addto_heatvec(); // allocate heat vector, there should be no heat in it

   // assume minor semi is rball*ratio2
   double ddr = rball*ratio2 - 0.5*mush_distance;
   ddr = 0.5;
   double Emush = Young_mush(r,il,ih, 0.0, ddr);
   double Emush_big = Young_mush_big(r,il,ih);
   print_run_double(mush_distance, "max spring length", fpr);
   print_run_double(ddr, "ddr", fpr);
   print_run_double(Emush, "Young's modulus hot", fpr);
   print_run_double(Emush_big, "Young's modulus big hot", fpr);
   print_run_double(Emush*ks_cold/ks_hot, "Young's modulus cold", fpr);
   print_run_double(Emush_big*ks_cold/ks_hot, "Young's modulus big cold", fpr);

   double K_T = Kthermal_mush(r,il,ih,0.0,ddr); // depends on k_heat

   print_run_double(K_T, "Thermal conductivity hot", fpr);
   print_run_double(K_T*k_heat_cold/k_heat_hot,"Thermal conductivity cold", fpr);

   print_run_double(mean_L(r),"Mean spring length", fpr);

   double om = 0.0; // set up the perturbing central mass
   if (npointmass >0){
      // set up central star
      int ip=0;
      om = add_pt_mass_kep(r, il, ih, -1, mp[ip], rad[ip],
           aa[ip],ee[ip], ii[ip], longnode[ip],argperi[ip],meananom[ip]);
      fprintf(fpr,"resbody mm=%.3f period=%.2f\n",om,2.0*M_PI/om);
      printf("resbody mm=%.3f period=%.2f\n",om,2.0*M_PI/om);
      icentral = ih;
      double na = om*aa[ip];
      double adot = 3.0*mp[ip]*na/pow(aa[ip],5.0); // should approximately be adot
      fprintf(fpr,"adot %.3e\n",adot);
      // set up rest of point masses
      for(int ipp = 1;ipp<npointmass;ipp++){
          double omp; // central mass assumed to be first one ipp=ih = icentral
          omp = add_pt_mass_kep(r, il, ih, icentral, mp[ipp], rad[ipp],
             aa[ipp],ee[ipp], ii[ipp], longnode[ipp],argperi[ipp],meananom[ipp]);
          fprintf(fpr,"pointm %d mm=%.3f period=%.2f\n",ipp,omp,2.0*M_PI/omp);
          printf("pointm %d mm=%.3f period=%.2f\n",ipp,omp,2.0*M_PI/omp);
      }
      npert = npointmass;
   }

   // factor of 0.5 is due to reduced mass being used in calculation
   double tau_relax_hot = 1.0*gamma_hot*0.5*(mball/(r->N -npert))
                   /spring_mush_hot.ks; // Kelvin Voigt relaxation time
   double tau_relax_cold = 1.0*gamma_cold*0.5*(mball/(r->N -npert))
                   /spring_mush_cold.ks; // Kelvin Voigt relaxation time
   print_run_double(tau_relax_hot , "relaxation time hot", fpr);
   print_run_double(tau_relax_cold, "relaxation time cold", fpr);

   double barchi_hot = 2.0*fabs(om - omegaz)*tau_relax_hot;  // initial value of barchi
   double barchi_cold = 2.0*fabs(om - omegaz)*tau_relax_cold;  // initial value of barchi
   print_run_double(barchi_hot , "barchi hot", fpr);
   print_run_double(barchi_cold, "barchi cold", fpr);

   double Nratio = (double)NS/(double)r->N;
   printf("N=%d  NS=%d NS/N=%.1f\n", r->N, NS, Nratio);
   fprintf(fpr,"N=%d  NS=%d NS/N=%.1f\n", r->N, NS,Nratio);
   fclose(fpr);

   reb_springs(r); // pass spring index list to display

   r->heartbeat = heartbeat;
#ifdef LIBPNG
// system("mkdir png");
#endif // LIBPNG

   if (tmax ==0.0)
      reb_integrate(r, INFINITY);
   else
      reb_integrate(r, tmax);
}


#define NSPACE 50
void heartbeat(struct reb_simulation* const r){
        char hfile[100];
        char nfile[100];
        static int first=0;
        static char extendedfile[50];
        static char pointmassfile[NPMAX*NSPACE];
        static int j0 = 0;
        if (first==0){
           first=1;
           sprintf(extendedfile,"%s_ext.txt",froot);
           for(int i=0;i<npert;i++){
             sprintf(pointmassfile+i*NSPACE,"%s_pm%d.txt",froot,i);
           }
           centerbody(r,0,r->N-npert);  // move reference frame, position only
           j0 = nearest_to_shape(r,0,r->N-npert,0.0,0.0,0.0); // find particle nearest origin
           // printf("j0=%d\n",j0);
        }

	if (reb_output_check(r,10.0*r->dt)){
		reb_output_timing(r,0);
	}
        if (fabs(r->t - t_damp) < 0.9*r->dt) set_gamma_fac(gamma_fac); 
            // damp initial bounce only 
            // reset gamma only at t near t_damp
	
         // stuff to do every timestep
         centerbody(r,0,r->N-npert);  // move reference frame, position only
         // don't store heat until after dampdown!
         if (r->t > t_damp){ 
            addto_heatvec(r); // store heat dE/dt accumulated each timestep
                              // in heatvec
            // the heat vector is zeroed and normalized only when printed out

            heat_nodes_tidal(r, npert, nodevec, r->dt);
              // stores tidal heat every timestep
              // applies 0.5*dEdt from each spring to the node (does not use heatvec)
              // tidal heating only right now
              // nodevec[j0].temp = 1.0; // fix temperature of node near origin

            transport_heat(r, npert, nodevec, r->dt); // adjust temps over network
            
         }
	 if (reb_output_check(r,t_heat)) { // heat files
               int ndt = (int)(t_heat/r->dt);  // how long tidal heat is stored up
               hfilename(r,froot, t_heat, hfile);
               print_heat(r,npert,hfile,ndt,powerfac); // heat info printed out!

               nfilename(r,froot, t_heat, nfile);
               print_node(r,npert,nodevec,nfile); // temperature info printed out!
               if (r->t > t_damp) 
                   adjust_spring_temp_ts(r, nodevec, 
                      spring_mush_hot, spring_mush_cold, Ttrans);
         }

         if (reb_output_check(r,t_print)) {
            print_extended(r,0,r->N-npert,extendedfile); // orbital info and stuff
            if (npert>0)
               for(int i=0;i<npert;i++){
                  int ip = icentral+i;
                  print_pm(r,ip,i,pointmassfile+i*NSPACE);
               }
         }



}

// make a spring index list
void reb_springs(struct reb_simulation* const r){
   r->NS = NS;
   r->springs_ii = malloc(NS*sizeof(int));
   r->springs_jj = malloc(NS*sizeof(int));
   for(int i=0;i<NS;i++){
     r->springs_ii[i] = springs[i].i;
     r->springs_jj[i] = springs[i].j;
   }
}


void print_run_double(double quantity, char* lstring, FILE *fpr)
{
   double lq = log10(quantity);
   if ((lq>4) || (lq>-4)){
      printf("%s %.4f\n",lstring,quantity);
      fprintf(fpr,"%s %.4f\n",lstring,quantity);
   }
   else {
      printf("%s %.3e\n",lstring,quantity);
      fprintf(fpr,"%s %.3e\n",lstring,quantity);
   }
}

