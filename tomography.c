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

void updatevelmodel(float* x, float* slow, int nm, float dmis,int i){
/*< TODO update velocity model >*/
	//float v;
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

float calculateTimeMissfit(float* s, /* NIP sources matrix */
			   float v0,
			   float* t0,
			   float* m0,
			   float* RNIP,
			   float* BETA,
			   int *n, 
			   float *o,
			   float *d,
			   float *slow,
			   float *a,
			   int is)
/*< Return time missfit sum of source-NIP-receiver rays >*/
{

	float currentRayAngle;
	int i, ir, it;
	float p[2], t, nrdeg;
	int nt=5000, nr=4; //TODO to correct nr
	float dt=0.001;
	raytrace rt;
	float** traj; // Ray trajectory (z,x)
	float m, h, tmis=0;
	float xs, xr, tr, ts, *x;

	x = sf_floatalloc(2);

	//for(is=0;is<ns;is++){

		x[0]=s[0];
		x[1]=s[1];
		nrdeg = a[is]; // TODO is in degree?
		//sf_warning("=> sx=%f sy=%f sa=%f",x[1],x[0],nrdeg);

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
						//sf_warning("xs=%f ts=%f",xs,ts);
					}else{ 
						tr=t;
						xr=x[1];
						//sf_warning("xr=%f tr=%f",xr,tr);
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

				x[0] = s[0];
				x[1] = s[1];
			} /* Loop over source-NIP-receiver rays */

			m = (xr+xs)/2.;
			h = (xr-xs)/2.;
			t = creTimeApproximation(h,m,v0,t0[is],m0[is],RNIP[is],BETA[is],true);
			tmis += fabs((ts+tr)-t);
			//sf_warning("=> tc=%f t=%f tmis=%f dtmis=%f\n",t,ts+tr,ts+tr-t,tmis);

		} /* Loop over reflection rays */
	//} /* Loop over NIP sources */

	tmis = (tmis*tmis)/(2*nr);
	return tmis;
}
