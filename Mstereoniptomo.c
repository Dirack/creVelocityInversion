/* VFSA velocity inversion based on stereotomography and NIP tomography strategies

The initial velocity model and NIP sources position used in this program is set up using sfnipmodsetup. This program does the forward modelling by ray tracying from NIP sources to surface and gets reflection traveltime.

The time misfit is calculated by the difference between the reflection traveltime obtained in the forward modelling and the traveltime calculated by CRE traveltime approximation formula for RNIP and BETA parameters given. This time misfit is used as a convergence criteria for VFSA global optimization algorithm to obtain optimized velocity model.

*/

#include <math.h>
#include <rsf.h>
#include "tomography.h"
#include "vfsacrsnh_lib.h"
#define N_STRIPES 50

void (*updateVelocityModel)(int*, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float*, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float*, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float*, /* Velocity model disturbance */
			   float*, /* Depth coordinates of sv vector */
			   float*, /* Velocity model */
			   int, /* sv n1 dimwnsion */
			   int, /* sv n2 dimension */
			   float, /* Near surface velocity */
			   float/* Depth velocity gradient */);

int main(int argc, char* argv[])
{
	bool verb; // Verbose parameter
	int n[2]; // Velocity grid dimensions n[0]=n1, n[1]=n2
	float d[2]; // Velocity grid sampling d[0]=d1, d[1]=d2
	float o[2]; // Velocity grid origin o[0]=o1, o[1]=o2
	float** s; // NIP source position (z,x)
	float* cnew; // Temporary parameters vector used in VFSA
	float* ots; // Optimized parameters vector
	float tmis0=100; // Best time misfit
	float otmis=0; // Best time misfit
	float deltaE; // Delta (Metrópolis criteria in VFSA)
	float Em0=0; // Energy (VFSA algorithm)
	float PM; // Metrópolis criteria
	float temp=1; // Temperature for VFSA algorithm
	float u=0; // Random number between 0 and 1
	int nit; // Number of VFSA iterations
	float temp0; // Initial temperature for VFSA
	float c0; // Damping factor for VFSA
	int ndim; // n1 dimension in shotsfile, should be equal 2
	int nshot; // n2 dimensions in shotsfile, number of shots
	int nm; // Number of samples in velocity grid n1*n2
	float* a; // Normal Ray initial angle for each NIP source
	float* slow; // slowness model
	int im, k, i; // loop counter
	float v; // Velocity temporary variable
	float v0; // Near surface velocity
	int ns; // Number of NIP sources
	int q; // Loop counter for VFSA iteration
	float tmis; // data time misfit value
	float *m0; // CMP's for normal rays
	float *t0; // t0's for normal rays
	float *RNIP; // Rnip parameters vector
	float *BETA; // Beta parameters vector
	float* sz; // Depth coordinates of the spline velocity function
	int nsz; // Dimension of sz vector
	float* sv; // Velocity coordinates of the spline velocity function
	float* gz; // Depth velocity gradient for backgorund velocity model
	sf_file shots; // NIP sources (z,x)
	sf_file vel; // background velocity model
	sf_file velinv; // Inverted velocity model
	sf_file angles; // Normal ray angles (degrees)
	sf_file m0s; // Central CMPs m0
	sf_file t0s; // Normal ray traveltimes
	sf_file rnips; // RNIP parameter for each m0
	sf_file betas; // BETA parameter for each m0
	sf_file sz_file; // z coordinates of the cubic spline functions
	sf_file gradz; // Depth velocity gradient for background model
	sf_file vspline; // Cubic spline velocity model

	updateVelocityModel = interpolateSlowModel;

	sf_init(argc,argv);

	shots = sf_input("shotsfile");
	vel = sf_input("in");
	velinv = sf_output("out");
	vspline = sf_output("vspline");
	angles = sf_input("anglefile");
	m0s = sf_input("m0s");
	t0s = sf_input("t0s");
	rnips = sf_input("rnips");
	betas = sf_input("betas");
	sz_file = sf_input("sz");
	gradz = sf_input("gz");

	/* Velocity model: get 2D grid parameters */
	if(!sf_histint(vel,"n1",n)) sf_error("No n1= in input");
	if(!sf_histint(vel,"n2",n+1)) sf_error("No n2= in input");
	if(!sf_histfloat(vel,"d1",d)) sf_error("No d1= in input");
	if(!sf_histfloat(vel,"d2",d+1)) sf_error("No d2= in input");
	if(!sf_histfloat(vel,"o1",o)) o[0]=0.;
	if(!sf_histfloat(vel,"o2",o+1)) o[1]=0.;
	
	if(!sf_getbool("verb",&verb)) verb=true;
	/* verbose parameter (y/n) */

	if(!sf_getfloat("v0",&v0)) v0=1.5;
	/* Near surface velocity (Km/s) */

	if(!sf_getint("nit",&nit)) nit=1;
	/* Number of VFSA iterations */

	if(!sf_getfloat("temp0",&temp0)) temp0=5;
	/* Initial temperature for VFSA algorithm */

	if(!sf_getfloat("c0",&c0)) c0=0.1;
	/* Damping factor for VFSA algorithm */

	/* Shotsfile: get shot points */
	if(!sf_histint(shots,"n1",&ndim) || 2 != ndim)
		sf_error("Must have n1=2 in shotsfile");
	if(!sf_histint(shots,"n2",&nshot)) sf_error("No n2= in shotfile");
	s = sf_floatalloc2(ndim,nshot);
	sf_floatread(s[0],ndim*nshot,shots);
	sf_fileclose(shots);

	/* Cubic spline vectors */
	if(!sf_histint(sz_file,"n1",&nsz)) sf_error("No n1= in sz file");
	
	/* Build cubic spline velocity matrix */
	gz = sf_floatalloc(1);
	sf_floatread(gz,1,gradz);
	sv = sf_floatalloc(N_STRIPES*nsz);

	sz = sf_floatalloc(nsz);
	sf_floatread(sz,nsz,sz_file);

	/* Initialize sv, the disturbance matrix in background model */
	for(k=0;k<N_STRIPES;k++){
		for(i=0;i<nsz;i++){
			sv[(k*nsz)+i]=0.0;
		}
	}

	/* VFSA parameters vectors */
	cnew = sf_floatalloc(N_STRIPES*nsz);
	ots = sf_floatalloc(N_STRIPES*nsz);

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

	/* get slowness squared (Background model) */
	nm = n[0]*n[1];
	slow =  sf_floatalloc(nm);
	sf_floatread(slow,nm,vel);

	for(im=0;im<nm;im++){
		v = slow[im];
		slow[im] = 1./(v*v);
	}

	if(verb){
		sf_warning("Command line Parameters");
		sf_warning("v0=%f nit=%d temp0=%f c0=%f",v0,nit,temp0,c0);
		sf_warning("Input file (Velocity model)");
		sf_warning("n1=%d d1=%f o1=%f",*n,*d,*o);
		sf_warning("n2=%d d2=%f o2=%f",*(n+1),*(d+1),*(o+1));
		sf_warning("Input file (shotsfile)");
		sf_warning("n1=%d",ndim);
		sf_warning("n2=%d",nshot);
		sf_warning("Input file (anglefile)");
		sf_warning("n1=%d",ns);
		sf_warning("Input file (sz)");
		sf_warning("nz=%d",nsz);
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

	/* cubic spline velocity matrix */
	sf_putint(vspline,"n1",nsz);
	sf_putint(vspline,"n2",N_STRIPES);

	/* Very Fast Simulated Annealing (VFSA) algorithm */
	for (q=0; q<nit; q++){
	
		/* calculate VFSA temperature for this iteration */
		temp=getVfsaIterationTemperature(q,c0,temp0);
						
		/* parameter disturbance */
		disturbParameters(temp,cnew,sv,nsz*N_STRIPES,0.001);

		/* Function to update velocity model */
		updateVelocityModel(n,o,d,sv,sz,slow,nsz,N_STRIPES,v0,gz[0]);
		//interpolateSlowModel(n, o, d,sv,sz,slow,nsz,N_STRIPES,v0,gz[0]);

		tmis=0;
	
		/* Calculate time missfit through forward modeling */		
		tmis=calculateTimeMissfit(s,v0,t0,m0,RNIP,BETA,n,o,d,slow,a,nshot);

		if(fabs(tmis) < fabs(tmis0) ){
			otmis = fabs(tmis);
			/* optimized parameters */
			for(im=0;im<nsz*N_STRIPES;im++)
				ots[im]=cnew[im];
			tmis0 = fabs(tmis);
		}

		/* VFSA parameters update condition */
		deltaE = -fabs(tmis) - Em0;
		
		/* Metrópolis criteria */
		PM = expf(-deltaE/temp);
		
		if (deltaE<=0){
			for(im=0;im<nsz*N_STRIPES;im++)
				sv[im]=cnew[im];
			Em0 = -fabs(tmis);
		} else {
			u=getRandomNumberBetween0and1();
			if (PM > u){
				for(im=0;im<nsz*N_STRIPES;im++)
					sv[im]=cnew[im];
				Em0 = -fabs(tmis);
			}	
		}	
			
		sf_warning("%d/%d => (%f)",q+1,nit,otmis);

	} /* loop over VFSA iterations */

	/* Generate optimal velocity model */
	interpolateVelModel(n, o, d,sv,sz,slow,nsz,N_STRIPES,v0,gz[0]);

	/* Write velocity model file */
	sf_floatwrite(slow,nm,velinv);

	/* Write velocity cubic spline function */
	sf_floatwrite(ots,nsz*N_STRIPES,vspline);
}
