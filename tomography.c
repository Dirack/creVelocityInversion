#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <rsf.h>
#include "raytrace.h"
#include "tomography.h"

#ifndef GDB_DEBUG
	#define DSLOW 0.04
	#define DANGLE 5.0
	#define MAX_ITERATIONS 2
#else
	#define DSSLOW 0.04
	#define DANGLE 5.0
	#define MAX_ITERATIONS 1
#endif
/*^*/

void updatevelmodel(float* x, float* slow, int nm, float dmis,int i){
/*< TODO update velocity model >*/
	float v;
	int im;
	for(im=0;im<nm;im++){
		slow[im]-=0.01;
	}
}

float creTimeApproximation(float h, 
			 float m,
			 float v0,
			 float t0,
			 float m0,
			 float RNIP,
			 float BETA,
			 bool cds){ 
/*< CRE traveltime approximation >*/
	float alpha;
	float d = m-m0;
	float c1;
	float c2;
	float a1, a2, b2, b1, Fd, Fd1, Fd2;
	float t;

	if(!cds){
		c1 = (d+h)/RNIP;
		c2 = (d-h)/RNIP;
		alpha = sin(BETA)/RNIP;
		t = (t0-2*RNIP/v0)+(RNIP/v0)*sqrt(1-2*alpha*(d+h)+c1*c1)+(RNIP/v0)*sqrt(1-2*alpha*(d-h)+c2*c2);
	}else{

		a1=(2*sin(BETA))/(v0);
		a2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RNIP);
		b2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RNIP);
		b1=2*b2+a1*a1-a2;
		Fd=(t0+a1*d)*(t0+a1*d)+a2*d*d;
		Fd2=(t0+a1*(d-h))*(t0+a1*(d-h))+a2*(d-h)*(d-h);
		Fd1=(t0+a1*(d+h))*(t0+a1*(d+h))+a2*(d+h)*(d+h);
		t=sqrt((Fd+b1*h*h+sqrt(Fd2*Fd1))*0.5);
	}
	return t;
}

void raystraveltimes(
		float *ts, /* sources-NIP rays traveltime */
		float *tr, /* NIP-receivers rays traveltimes */
		float *xs, /* sources-NIP rays position */
		float *xr, /* NIP-receivers rays position */
		float x[2], /* NIP position */
		float nrdeg, /* Normal ray angle in degrees */
		int n[2],
		float o[2],
		float d[2],
		float* slow,
		int nr /* number of reflection rays */
		    )
/*< Return traveltimes of source-NIP-receiver rays and endpoints >*/
{

	float currentRayAngle;
	int i, ir, it;
	float p[2], s[2], t;
	float nt=5000;
	float dt=0.001;
	raytrace rt;
	float** traj; // Ray trajectory (z,x)

	s[0] = x[0];
	s[1] = x[1];

	#ifdef GDB_DEBUG
	nr=2;
	#endif
	for(ir=0;ir<nr;ir++){

		for(i=0; i<2; i++){

			/* initialize ray tracing object */
			rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);

			/* Ray tracing */
			traj = sf_floatalloc2(2,nt+1);
			
			/* initialize ray direction */
			currentRayAngle=(i==0)?(nrdeg-(ir+1)*DANGLE)*DEG2RAD:(nrdeg+(ir+1)*DANGLE)*DEG2RAD;

			p[0] = -cosf(currentRayAngle);
			p[1] = sinf(currentRayAngle);

			it = trace_ray (rt, x, p, traj);

			if(it>0){
				t = it*dt;
				if(i==0){
					ts[ir]=t;
					xs[ir]=x[1];
					sf_warning("ts=%f",ts[ir]);
				}else{ 
					tr[ir]=t;
					xr[ir]=x[1];
					sf_warning("tr=%f",tr[ir]);
				}
			}else if(it == 0){
				t = abs(nt)*dt;
				nt += 1000;
			}else{
				sf_error("Bad angle, ray get to the model side/bottom");
			}

			/* Raytrace close */
			raytrace_close(rt);
			free(traj);

			x[0] = s[0];
			x[1] = s[1];
		} /* Loop over source-NIP-receiver rays */
	} /* Loop over reflection rays */
}
