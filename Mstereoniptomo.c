/* Ray tracer based on stereotomography strategy

Trace rays from NIP sources to acquisition surface in order to get traveltime curves used in NIP tomography.

*/

#include <math.h>
#include <rsf.h>
#include "tomography.h"

int main(int argc, char* argv[])
{
	bool verb; // Verbose parameter
	int n[2]; // Velocity grid dimensions n[0]=n1, n[1]=n2
	float d[2]; // Velocity grid sampling d[0]=d1, d[1]=d2
	float o[2]; // Velocity grid origin o[0]=o1, o[1]=o2
	float** s; // NIP source position (z,x)
	int ndim; // n1 dimension in shotsfile, should be equal 2
	int nshot; // n2 dimensions in shotsfile
	int nm; // Number of samples in velocity grid
	float* a; // Normal Ray initial angle
	int nr=2; // Number of rays to trace
	float* slow; // slowness
	int im; // loop counter
	float v; // Velocity
	float v0; // Near surface velocity
	int ns; // Number of NIP sources
	int is; // Loop counter for NIP sources
	int ir; // Loop counter of reflection rays
	int i; // Loop over iterations
	float x[2]; // Ray initial position
	float tmis[1]; // data misfit t vector
	float dmis;
	float *ts, *tr, *xs, *xr, offset, cmp;
	float *m0, *t0, *RNIP, *BETA;
	sf_file shots, vel, angles, timeCurve, xCurve, m0s, t0s, rnips, betas;

	sf_init(argc,argv);

	shots = sf_input("shotsfile");
	vel = sf_input("in");
	timeCurve = sf_output("out");
	xCurve = sf_output("x");
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

	sf_putint(xCurve,"n1",2);
	sf_putint(xCurve,"n2",ns);
	sf_putint(timeCurve,"n1",2);
	sf_putint(timeCurve,"n2",ns);

	ts = sf_floatalloc(nr);
	xs = sf_floatalloc(nr);
	tr = sf_floatalloc(nr);
	xr = sf_floatalloc(nr);

	for(i=0; i<MAX_ITERATIONS;i++){

		for(is=0; is<ns; is++){

			/* initialize position */
			x[0] = s[is][0];
			x[1] = s[is][1];

			raystraveltimes(ts,tr,xs,xr,x,a[is],n,o,d,slow,nr);

			for(ir=0;ir<nr;ir++){

				/* Calculte data misfit (t,x) */
				offset=(xr[ir]-xs[ir])/2.;
				cmp=(xr[ir]+xs[ir])/2.;
				tmis[is] = creTimeApproximation(offset,cmp,v0,t0[is],m0[is],RNIP[is],BETA[is],0) - (ts[ir] + tr[ir]);
				sf_warning("t=%f tcre=%f tmis=%f xs=%f xr=%f",(ts[ir]+tr[ir]),tmis[is]+(ts[ir]+tr[ir]),
				tmis[is],xs[ir],xr[ir]);
				dmis += tmis[is];
			} /*Loop over reflection rays */


		} /* Loop over NIP sources */

		/* TODO */
		updatevelmodel(x,slow,nm,dmis,i);
		dmis = 0;

	} /* Loop over iterations */
}