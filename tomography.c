#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <rsf.h>
#include "raytrace.h"
#include "tomography.h"

#ifndef GDB_DEBUG
	#define DSLOW 0.04
	#define DANGLE 1.0
#else
	#define DSSLOW 0.04
	#define DANGLE 1.0
#endif
/*^*/

void updatevelmodel(float* slow, /* Slowness vector */
		    int* n, /* n[0]=n1 n2=n[1] */
		    float* o, /* o[0]=o1 o[1]=o2 */
		    float* d, /* d[0]=d1 d[1]=d2 */
		    float v0, /* First velocity is the near surface velocity */
		    float grad /* Velocity gradient */)
/*< Funcion to update constant velocity gradient:
Note:
This is a scratch of the function to update the velocity model,
it uses a constant gradient velocity model, and update the 
gradient in each iteration of the process.

The purpose is to show that NIP sources will converge to the
reflector interface with the "right" gradient used.

TODO: Modify this function to VFSA optimization of the velocity model.
 >*/
{
	int i, j;
	float x, v;

	for(i=0;i<n[0];i++){

		x = i*d[0]+o[0];
		v = grad*x+v0;

		for(j=0;j<n[1];j++){
			slow[j*n[0]+i]=1./(v*v);
		}
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

float calculateTimeMissfit(float** s, /* NIP sources matrix (z,x) pairs */
			   float v0, /* Near surface velocity */
			   float* t0, /* Normal ray traveltime for each NIP source */
			   float* m0, /* Central CMP for each NIP source */
			   float* RNIP, /* RNIP parameter for each NIP source */
			   float* BETA, /* BETA parameter for each NIP source */
			   int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float *slow, /* Slowness velociy model */
			   float *a, /* Normal ray angle for each NIP source */
			   int ns /* Number of NIP sources */)
/*< Return time missfit sum of source-NIP-receiver rays 

Note: This function traces nr reflection rays pairs from each NIP source
position passed though s matrix. It also calculate the difference between
the traveltime of the traced rays with calculated traveltime using CRE
traveltime approximation to calculate the time misfit returned by the function.
 >*/
{

	float currentRayAngle;
	int i, ir, it, is;
	float p[2], t, nrdeg;
	int nt=5000, nr=5; //TODO to correct nr
	float dt=0.001;
	raytrace rt;
	float** traj; // Ray trajectory (z,x)
	float m, h, tmis=0;
	float xs, xr, tr, ts, *x;

	x = sf_floatalloc(2);

	for(is=0;is<ns;is++){

		x[0]=s[is][0];
		x[1]=s[is][1];
		nrdeg = a[is]; // TODO verify is it in degree?

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
						ts=t;
						xs=x[1];
					}else{ 
						tr=t;
						xr=x[1];
					}
				}else if(it == 0){
					t = abs(nt)*dt;
					nt += 1000;
				}else{
					sf_warning("=> x=%f y=%f t=%f",s[1],s[0],t);
					sf_error("Bad angle, ray get to the model side/bottom");
				}

				/* Raytrace close */
				raytrace_close(rt);
				free(traj);

				x[0] = s[is][0];
				x[1] = s[is][1];
			} /* Loop over source-NIP-receiver rays */

			m = (xr+xs)/2.;
			h = (xr-xs)/2.;
			t = creTimeApproximation(h,m,v0,t0[is],m0[is],RNIP[is],BETA[is],true);
			tmis += fabs((ts+tr)-t);

		} /* Loop over reflection rays */

	} /* Loop over NIP sources */

	/* TODO: Evaluate the best function to calcullate the time misfit */
	tmis = (tmis*tmis)/(nr*ns);
	return tmis;
}
