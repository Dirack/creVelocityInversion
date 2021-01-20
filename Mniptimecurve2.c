/* New version of NIP Point source modeling ray tracer 

Trace rays from NIP sources to acquisition surface in order to get traveltime curves used in NIP tomography.

*/

#include <math.h>

#include <rsf.h>

#include "raytrace.h"

#define DANGLE 15.0
#define INITIAL_ANGLE 45.0

int main(int argc, char* argv[])
{
	bool verb; // Verbose parameter
	int n[2]; // Velocity grid dimensions n[0]=n1, n[1]=n2
	float d[2]; // Velocity grid sampling d[0]=d1, d[1]=d2
	float o[2]; // Velocity grid origin o[0]=o1, o[1]=o2
	float t; // Ray traveltime
	float** s; // NIP source position (z,x)
	int ndim; // n1 dimension in shotsfile, should be equal 2
	int nshot; // n2 dimensions in shotsfile
	int nm; // Number of samples in velocity grid
	float* a; // Normal Ray initial angle
	int nr=7; // Number of rays to trace
	float* slow; // slowness
	int im; // loop counter
	float v0; // Velocity
	float dt=0.001; // ray time sampling
	int nt=1000; // Numer of time samples for each ray
	int ns; // Number of NIP sources
	int it,ir,is; // loop counters
	float** traj; // Ray trajectory (z,x)
	float x[2]; // Ray initial position
	float p[2]; // Ray initial slowness vector
	float currentRayAngle; // Ray initial angle
	sf_file shots, vel, angles, timeCurve,xCurve;
	raytrace rt;

	sf_init(argc,argv);

	shots = sf_input("shotsfile");
	vel = sf_input("in");
	timeCurve = sf_output("out");
	xCurve = sf_output("x");
	angles = sf_input("anglefile");

	/* Velocity model: get 2D grid parameters */
	if(!sf_histint(vel,"n1",n)) sf_error("No n1= in input");
	if(!sf_histint(vel,"n2",n+1)) sf_error("No n2= in input");
	if(!sf_histfloat(vel,"d1",d)) sf_error("No d1= in input");
	if(!sf_histfloat(vel,"d2",d+1)) sf_error("No d2= in input");
	if(!sf_histfloat(vel,"o1",o)) o[0]=0.;
	if(!sf_histfloat(vel,"o2",o+1)) o[1]=0.;
	
	if(!sf_getbool("verb",&verb)) verb=true;
	/* verbose */

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


	/* get slowness squared */
	nm = n[0]*n[1];
	slow =  sf_floatalloc(nm);
	sf_floatread(slow,nm,vel);

	for(im=0;im<nm;im++){
		v0 = slow[im];
		slow[im] = 1./(v0*v0);
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

	sf_putint(xCurve,"n1",nr);
	sf_putint(xCurve,"n2",ns);
	sf_putint(timeCurve,"n1",nr);
	sf_putint(timeCurve,"n2",ns);

	for(is=0; is<ns; is++){
		for(ir=0; ir<nr; ir++){

			/* initialize ray tracing object */
			rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);

			/* Ray tracing */
			traj = sf_floatalloc2(ndim,nt+1);
			
			/* initialize position */
			x[0] = s[is][0]; 
			x[1] = s[is][1];

			/* initialize direction */
			currentRayAngle=(-INITIAL_ANGLE+ir*DANGLE+a[is])*DEG2RAD;
			p[0] = -cosf(currentRayAngle);
			p[1] = sinf(currentRayAngle);

			it = trace_ray (rt, x, p, traj);

			if(it>0){
				t = it*dt;
				sf_floatwrite(&traj[it][1],1,xCurve);
			}else if(it == 0){
				t = abs(nt)*dt;
				sf_floatwrite(&traj[nt][1],1,xCurve);
			}else{
				t = -it * dt;
				sf_floatwrite(&traj[-it][1],1,xCurve);
			}

			sf_floatwrite(&t,1,timeCurve);

			/* Raytrace close */
			raytrace_close(rt);
			free(traj);
		}
	}

}
