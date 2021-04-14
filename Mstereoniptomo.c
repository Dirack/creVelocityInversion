/* Ray tracer based on stereotomography strategy

Trace rays from NIP sources to acquisition surface in order to get traveltime curves used in NIP tomography.

*/

#include <math.h>
#include <rsf.h>
#include "tomography.h"
#include "vfsacrsnh_lib.h"
#define MAX_ITERATIONS 100
#define temp0 5
#define c0 0.1

int main(int argc, char* argv[])
{
	bool verb; // Verbose parameter
	int n[2]; // Velocity grid dimensions n[0]=n1, n[1]=n2
	float d[2]; // Velocity grid sampling d[0]=d1, d[1]=d2
	float o[2]; // Velocity grid origin o[0]=o1, o[1]=o2
	float** s; // NIP source position (z,x)
	float* cnew;
	float* ots;
	float tmis0=100, otmis=0, deltaE, Em0=0, PM, temp=1, u=0;
	int ndim; // n1 dimension in shotsfile, should be equal 2
	int nshot; // n2 dimensions in shotsfile
	int nm; // Number of samples in velocity grid
	float* a; // Normal Ray initial angle
	float* slow; // slowness
	int im; // loop counter
	float v; // Velocity
	float v0; // Near surface velocity
	int ns; // Number of NIP sources
	int is; // Loop counter for NIP sources
	int q; // Loop counter for VFSA iteration
	float tmis; // data time misfit value
	float *m0, *t0, *RNIP, *BETA;
	float x[2];
	int iv, nv=10;
	float* gv;
	int i, j;
	sf_file out, shots, vel, velinv, angles, m0s, t0s, rnips, betas;

	sf_init(argc,argv);

	shots = sf_input("shotsfile");
	vel = sf_input("in");
	out = sf_output("out");
	velinv = sf_output("velinv");
	angles = sf_input("anglefile");
	m0s = sf_input("m0s");
	t0s = sf_input("t0s");
	rnips = sf_input("rnips");
	betas = sf_input("betas");

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
	cnew = sf_floatalloc(1);
	ots = sf_floatalloc(1);
	gv = sf_floatalloc(1);
	gv[0]=0.01;
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

	/* Optimized parameters */
	sf_putint(out,"n1",ndim);
	sf_putint(out,"n2",ns);
	sf_putint(out,"n3",1);
	sf_putfloat(out,"d1",1);
	sf_putfloat(out,"o1",0);
	sf_putfloat(out,"d2",1);
	sf_putfloat(out,"o2",0);
	sf_putfloat(out,"d3",1);
	sf_putfloat(out,"o3",0);

	/* Velocity models from inversion */
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
		disturbParameters(temp,cnew,gv,1,0.001);

		/* Function to update velocity gradient */
		updatevelmodel(slow, n, o, d, v0, gv[0]);

		tmis=0;
	
		/* Calculate time missfit through forward modeling */		
		tmis=calculateTimeMissfit(s,v0,t0,m0,RNIP,BETA,n,o,d,slow,a,nshot);

		if(fabs(tmis) < fabs(tmis0) ){
			otmis = fabs(tmis);
			/* optimized parameters matrix */
			ots[0]=cnew[0];
			tmis0 = fabs(tmis);			
		}

		/* VFSA parameters update condition */
		deltaE = -fabs(tmis) - Em0;
		
		/* Metrópolis criteria */
		PM = expf(-deltaE/temp);
		
		if (deltaE<=0){
			gv[0]=cnew[0];
			Em0 = -fabs(tmis);
		} else {
			u=getRandomNumberBetween0and1();
			if (PM > u){
				gv[0]=cnew[0];
				Em0 = -fabs(tmis);
			}	
		}	
			
		sf_warning("%f => %d/%d (%f)",gv[0],q,MAX_ITERATIONS,otmis);	

	} /* loop over iterations */

	/* Print optimal velocity gradient */
	sf_warning("=> vgrad=%f v0=%f",ots[0],v0);	
	
	/* Generate optimal velocity model */
	for(i=0;i<n[0];i++){

		v = ots[0]*(i*d[0]+o[0])+v0;

		for(j=0;j<n[1];j++){
			slow[j*n[0]+i]=v;
		}
	}
	sf_floatwrite(slow,nm,velinv);

	/* NIP sources position */	
	sf_floatwrite(s,ndim*nshot,out);
}
