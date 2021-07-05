/* VFSA Depth velocity gradient inversion - First iteration of stereoniptomo program

This program uses Very Fast Simulated Annealing (VFSA) to obtain the depth velocity gradient for a linear velocity model that velocity varies with depth. This is the background model used as the first iteration for sfstereoniptomo program.

*/

#include <math.h>
#include <rsf.h>
#include "tomography.h"
#include "zgradvfsa_lib.h"
#define MAX_ITERATIONS 30
#define temp0 5
#define c0 0.1

int main(int argc, char* argv[])
{
	bool verb; // Verbose parameter
	int n[2]; // Velocity grid dimensions n[0]=n1, n[1]=n2
	float d[2]; // Velocity grid sampling d[0]=d1, d[1]=d2
	float o[2]; // Velocity grid origin o[0]=o1, o[1]=o2
	float** s; // NIP source position (z,x)
	float gznew=0.0; // Disturbed grad z of the current iteration
	float gzots; // Optimal Grad z (result)
	float tmis0=100; // Time misfit of each VFSA iteration (convergence criteria)
	float otmis=0; // Best time misfit of all iterations
	float deltaE; // Difference between the time misfit of this iteration and the previous
	float Em0=0; // Energy function of VFSA
	float PM; // Probability function of the Metrópolis criteria of VFSA
	float temp=1; // This iteration temperature VFSA
	float u=0; // random number between 0 and 1
	int ndim; // n1 dimension in shotsfile, should be equal 2
	int nshot; // n2 dimensions in shotsfile
	int nm; // Number of samples in velocity grid
	float* a; // Normal Ray initial angle
	float* slow; // slowness
	int im; // loop counter
	float v; // Velocity
	float v0; // Near surface velocity
	int ns; // Number of NIP sources
	int q; // Loop counter for VFSA iteration
	float tmis; // data time misfit value
	float *m0; // Central CMPs
	float *t0; // Normal ray traveltime for each m0
	float *RNIP; // RNIP parameter for each m0
	float *BETA; // BETA parameter for each m0
	float gz; // Velocity gradient in depth (z)
	int i, j; // loop counters
	sf_file shots; // NIP sources for tomography
	sf_file vel; // initial velocity model
	sf_file velinv; // inverted model
	sf_file angles; // NIP sources normal rays angles (Degrees)
	sf_file m0s; // Central CMPs
	sf_file t0s; // Normal ray traveltimes
	sf_file rnips; // RNIP parameter for each m0
	sf_file betas; // BETA parameter for each m0
	sf_file gradz; // Velocity gradient in depth

	sf_init(argc,argv);

	shots = sf_input("shotsfile");
	vel = sf_input("in");
	velinv = sf_output("out");
	angles = sf_input("anglefile");
	m0s = sf_input("m0s");
	t0s = sf_input("t0s");
	rnips = sf_input("rnips");
	betas = sf_input("betas");
	gradz = sf_output("gz");

	/* Velocity model: get 2D grid parameters */
	if(!sf_histint(vel,"n1",n)) sf_error("No n1= in input");
	if(!sf_histint(vel,"n2",n+1)) sf_error("No n2= in input");
	if(!sf_histfloat(vel,"d1",d)) sf_error("No d1= in input");
	if(!sf_histfloat(vel,"d2",d+1)) sf_error("No d2= in input");
	if(!sf_histfloat(vel,"o1",o)) o[0]=0.;
	if(!sf_histfloat(vel,"o2",o+1)) o[1]=0.;
	
	if(!sf_getbool("verb",&verb)) verb=true;
	/* verbose */

	if(!sf_getfloat("v0",&v0)) v0=1.5;
	/* Near surface velocity (Km/s) */

	/* Shotsfile: get m0 shot points */
	if(!sf_histint(shots,"n1",&ndim) || 2 != ndim)
		sf_error("Must have n1=2 in shotsfile");
	if(!sf_histint(shots,"n2",&nshot)) sf_error("No n2= in shotfile");
	s = sf_floatalloc2(ndim,nshot);
	sf_floatread(s[0],ndim*nshot,shots);
	gz=0.;
	sf_fileclose(shots);

	/* Anglefile: get initial emergence angle */
	if(!sf_histint(angles,"n1",&ns)) sf_error("No n1= in anglefile");
	a = sf_floatalloc(ns);
	sf_floatread(a,ns,angles);
	if(ns!=nshot) sf_error("n1 in anglefile should be equal to n2 in shotsfile!");

	/* allocate parameters vectors */
	m0 = sf_floatalloc(ns);
	sf_floatread(m0,ns,m0s);
	t0 = sf_floatalloc(ns);
	sf_floatread(t0,ns,t0s);
	RNIP = sf_floatalloc(ns);
	sf_floatread(RNIP,ns,rnips);
	BETA = sf_floatalloc(ns);
	sf_floatread(BETA,ns,betas);

	/* get slowness squared */
	nm = n[0]*n[1];
	slow =  sf_floatalloc(nm);
	sf_floatread(slow,nm,vel);

	for(im=0;im<nm;im++){
		v = slow[im];
		slow[im] = 1./(v*v);
	}

	if(verb){
		sf_warning("Input file (Velocity model)");
		sf_warning("n1=%d d1=%f o1=%f",*n,*d,*o);
		sf_warning("n2=%d d2=%f o2=%f",*(n+1),*(d+1),*(o+1));
		sf_warning("Input file (shotsfile)");
		sf_warning("n1=%d",ndim);
		sf_warning("n2=%d",nshot);
		sf_warning("Input file (anglefile)");
		sf_warning("n1=%d",ns);
	}

	/* Velocity model from inversion */
	sf_putint(velinv,"n1",n[0]);
	sf_putint(velinv,"n2",n[1]);
	sf_putint(velinv,"n3",1);
	sf_putfloat(velinv,"d1",d[0]);
	sf_putfloat(velinv,"d2",d[1]);
	sf_putfloat(velinv,"o1",o[0]);
	sf_putfloat(velinv,"o2",o[1]);
	sf_putfloat(velinv,"d3",1);
	sf_putfloat(velinv,"o3",0);

	for (q=0; q <MAX_ITERATIONS; q++){
	
		/* calculate VFSA temperature for this iteration */
		temp=getVfsaIterationTemperature(q,c0,temp0);
						
		/* parameter disturbance */
		disturbGradZ(temp,&gznew,gz,0.001);

		/* Function to update velocity gradient */
		updateGzVelModel(slow, n, o, d, v0, gznew);

		tmis=0;
	
		/* Calculate time missfit through forward modeling */		
		tmis=calculateTimeMissfit(s,v0,t0,m0,RNIP,BETA,n,o,d,slow,a,nshot);

		if(fabs(tmis) < fabs(tmis0) ){
			otmis = fabs(tmis);
			/* optimized parameters matrix */
			gzots=gznew;
			tmis0 = fabs(tmis);			
		}

		/* VFSA parameters update condition */
		deltaE = -fabs(tmis) - Em0;
		
		/* Metrópolis criteria */
		PM = expf(-deltaE/temp);
		
		if (deltaE<=0){
			gz=gznew;
			Em0 = -fabs(tmis);
		} else {
			u=getRandomNumberBetween0and1();
			if (PM > u){
				gz=gznew;
				Em0 = -fabs(tmis);
			}	
		}	
			
		sf_warning("%f => %d/%d (%f)",gz,q,MAX_ITERATIONS,otmis);	

	} /* loop over VFSA iterations */

	/* Print optimal velocity gradient */
	sf_warning("(%f)=> vgrad=%f v0=%f",tmis0,gzots,v0);	
	
	/* Generate optimal velocity model */
	for(i=0;i<n[0];i++){

		v = gzots*(i*d[0]+o[0])+v0;

		for(j=0;j<n[1];j++){
			slow[j*n[0]+i]=v;
		}
	}
	sf_floatwrite(slow,nm,velinv);

	/* Output gradient gz */
	sf_putint(gradz,"n1",1);
	sf_putint(gradz,"n2",1);
	sf_putint(gradz,"n3",1);
	sf_floatwrite(&gzots,1,gradz);
}
