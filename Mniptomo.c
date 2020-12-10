/* NIP tomography ray tracer */

#include <math.h>

#include <rsf.h>

#include "raytrace.h"

int main(int argc, char* argv[])
{
	bool verb;
	int n[2];
	float d[2];
	float o[2];
	float t;
	float** s;
	int ndim;
	int nshot;
	int nm;
	float* a;
	int nr;
	float* slow;
	int im;
	float v0;
	int nt0;
	float dt0;
	float* t0;
	int nt;
	int it;
	float** traj;
	float x[2];
	float p[2];
	int i;
	sf_file shots, vel, rays, angles, t0s;
	raytrace rt;

	sf_init(argc,argv);

	shots = sf_input("shotsfile");
	vel = sf_input("in");
	rays = sf_output("out");
	t0s = sf_input("t0s");
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
	if(sf_histfloat(shots,"o2",&t)) sf_putfloat(rays,"o3",t);
	if(sf_histfloat(shots,"d2",&t)) sf_putfloat(rays,"d3",t);
	s = sf_floatalloc2(ndim,nshot);
	sf_floatread(s[0],ndim*nshot,shots);
	sf_fileclose(shots);

	/* Anglefile: get initial emergence angle */
	if(!sf_histint(angles,"n1",&nr)) sf_error("No n1= in anglefile");
	a = sf_floatalloc(nr);
	sf_floatread(a,nr,angles);
	a[0]=a[0]*DEG2RAD;

	/* T0s file: File with ray's traveltime */
	if(!sf_histint(t0s,"n1",&nt0)) sf_error("No n1= in t0s file");
	t0 = sf_floatalloc(nt0);
	sf_floatread(t0,nt0,t0s);
	dt0=0.001;
	sf_putfloat(rays,"d1",dt0);
	sf_putfloat(rays,"o1",o[0]);
	sf_putfloat(rays,"o2",180.);
	sf_putfloat(rays,"d2",1.);

	/* get slowness squared */
	nm = n[0]*n[1];
	slow =  sf_floatalloc(nm);
	sf_floatread(slow,nm,vel);

	for(im=0;im<nm;im++){
		v0 = slow[im];
		slow[im] = 1./(v0*v0);
	}

	sf_putint(rays,"n1",nr);
	sf_putint(rays,"n2",1);
	sf_settype(rays,SF_COMPLEX);
	sf_putstring(rays,"label1","time");
	sf_putstring(rays,"unit1","s");
	sf_putstring(rays,"label2","Degrees");
	sf_putstring(rays,"unit2","Angle");
	sf_fileflush(rays,NULL);
	sf_settype(rays,SF_FLOAT);

	if(verb){
		sf_warning("Input file (Velocity model)");
		sf_warning("n1=%d d1=%f o1=%f",*n,*d,*o);
		sf_warning("n2=%d d2=%f o2=%f",*(n+1),*(d+1),*(o+1));
		sf_warning("Command line parameters");
		sf_warning("");
	}

	for(i=0; i<nr; i++){

		/* initialize ray tracing object */
		nt = (int) t0[0]/dt0;
		rt = raytrace_init(2,true,nt,dt0,n,o,d,slow,ORDER);

		/* Ray tracing */
		traj = sf_floatalloc2(ndim,nt+1);
		
		/* initialize position */
		x[0] = s[0][0]; 
		x[1] = s[0][1];

		/* initialize direction */
		p[0] = -cosf(a[0]);
		p[1] = sinf(a[0]);

		it = trace_ray (rt, x, p, traj);

		sf_floatwrite (traj[nt-1],ndim,rays);

		if(it>0){
			i = it >= 2 ? it - 2 : it - 1;
			/* Escape vector */
			x[0]=traj[it][0];
			x[1]=traj[it][1];
			x[0]-=traj[i][0];
			x[1]-=traj[i][1];
			/* Dot product with unit vector pointing upward */
			t = sqrt(x[0]*x[0]+x[1]*x[1]); /* Length */
			t = acos(x[0]/t);
			if(x[1]<0) t = -t;
		}else{
			t = a[0];
			t /= DEG2RAD;
		}

		/* Raytrace close */
		raytrace_close(rt);
	}
}
