#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <rsf.h>
#include "raytrace.h"
#include "tomography.h"
#define N_STRIPES 5

#ifndef GDB_DEBUG
	#define DSLOW 0.04
	#define DANGLE 1.0
#else
	#define DSSLOW 0.04
	#define DANGLE 1.0
#endif
/*^*/

void calculateSplineCoeficients(int n, /* Vectors (x,y) dimension */
				float* x, /* x coordinates */
				float** y, /* y coordinates */
				float** coef /* Spline coeficients */)
/*< Function to calculate natural cubic spline coeficients

Note: It Receives n points and two vectors x and y with n dimension.
It returns a coeficients vector with 4 coeficients for each of the
n-1 natural cubic splines, coef[n-1)*4].

IMPORTANT: The number of points must be equal or major than 3 (n>3)
and x vector must be in crescent order.

>*/
{

	float s2[n]; // Second derivatives matrix
	int i, ip1, ip2, im1, m, k; // Loop counter
	float hb, ha, deltaa, deltab, t; // temporary variables
	float e[n-2]; // hi's vector
	float dp[n-2]; // main diagonal

	/* Vectors dimension must be major than 3 */
	if(n<3){
		fprintf(stderr,"Erro, n<3\n");
		exit(-1);
	}

	/* x vector must be in crescent order */
	for(i=1;i<n;i++){
		if(x[i-1]>x[i]){
			fprintf(stderr,"Erro, vetor x deve possuir ordem crescente\n");
			exit(-2);
		}
	}

	for(k=0;k<N_STRIPES;k++){
		
		/* Simetric tridiagonal linear system build */
		ha = x[1]-x[0]; deltaa = (y[k][1]-y[k][0])/ha; m=n-2;
		for(i=0;i<m;i++){
			ip1 = i+1; ip2 = i+2;
			hb = x[ip2]-x[ip1];
			deltab = (y[k][ip2]-y[k][ip1])/hb;
			e[i] = hb; dp[i] = 2*(ha+hb);
			s2[ip1] = 6*(deltab-deltaa);
			ha=hb; deltaa=deltab;
		}

		/* Gauss elimination */
		for(i=1;i<m;i++){
			ip1=i+1; im1=i-1;
			t = e[im1]/dp[im1];
			dp[i] = dp[i]-t*e[im1];
			s2[ip1] = s2[ip1]-t*s2[i];
		}

		/* Retroactive substitutive solution */
		s2[m]=s2[m]/dp[m-1];
		for(i=m-1;i>0;i--){
			ip1=i+1; im1=i-1;
			s2[i]=(s2[i]-e[im1]*s2[ip1])/dp[im1];
		}
		s2[0]=0; s2[n-1]=0;

		/* Calculate spline coeficients */
		for(i=0;i<n-1;i++){
			ha = x[i+1]-x[i];
			coef[k][0+i*4] = (s2[i+1]-s2[i])/(6*ha);
			coef[k][1+i*4] = s2[i]/2;
			coef[k][2+i*4] = (y[k][i+1]-y[k][i])/ha-(s2[i+1]+2*s2[i])*(ha/6);
			coef[k][3+i*4] = y[k][i];
		}
	}
}

void updateSplineCubicVelModel( float* slow, /* Slowness vector */
		    		int* n, /* n[0]=n1 n2=n[1] */
		    		float* o, /* o[0]=o1 o[1]=o2 */
		    		float* d, /* d[0]=d1 d[1]=d2 */
			    	int dim, /* Dimension of (z,vz) vectors */
				float* sz, /* Spline not Depth coordinates */
				float** sv /* Spline not Velocity coordinates */)
/*< Funcion to update spline cubic velocity model:
Note:
Make a velocity varying with depth model using the spline cubic interpolation
for a set of points (z,vz) given. TODO

 >*/
{
	int i, j=0, ic, k;
	float z=0.0;
	float coef[N_STRIPES][4*(dim-1)];
	float v[N_STRIPES][n[0]];
	int app, app_len=n[1]/N_STRIPES;

	/* Calculate spline coeficients */
	//calculateSplineCoeficients(dim,sz,sv,coef);

	/* Calculate velocity function */
	for(k=0;k<N_STRIPES;k++){

		z = o[0];
		j = 0;

		for(i=1;i<dim;i++){
			
			ic = (i-1)*4;

			while(z<=sz[i]){
				z = (j*d[0]+o[0])-sz[i-1];
				if(j>=n[0]) break;
				v[k][j] = coef[k][0+ic]*z*z*z+coef[k][1+ic]*z*z+coef[k][2+ic]*z+coef[k][3+ic];
				j++;
			}
		}
	}

	/* Update slowness model */
	for(k=0;k<N_STRIPES;k++){

		app = (k*app_len);

		for(i=0;i<n[0];i++){

			for(j=app;j<((k+1)*app_len);j++){
				slow[j*n[0]+i]=1./(v[k][i]*v[k][i]);
				#ifdef GDB_DEBUG
				sf_warning("[%d][%d][%d] %f %f ",i,j,k,v[k][i],slow[j*n[0]+i]);
				#endif
			} /* Loop over distance */
		} /* Loop over depth */
	} /* Loop over cubic spline functions */
	#ifdef GDB_DEBUG
	sf_error("fim");
	#endif
}

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
	float z, v;

	for(i=0;i<n[0];i++){

		z = i*d[0]+o[0];
		v = grad*z+v0;

		for(j=0;j<n[1];j++){
			slow[j*n[0]+i]=1./(v*v);
		} /* Loop over distance */
	} /* Loop over depth */
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
