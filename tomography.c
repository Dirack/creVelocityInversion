#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <rsf.h>
#include "raytrace.h"
#include "tomography.h"

#define NR 5
#define NT 5000
#define DT 0.001
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
				float* y, /* y coordinates */
				float** coef, /* Spline coeficients */
				int n_stripes)
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

	for(k=0;k<n_stripes;k++){
		
		/* Simetric tridiagonal linear system build */
		ha = x[1]-x[0]; deltaa = (y[k*n+1]-y[k*n+0])/ha; m=n-2;
		for(i=0;i<m;i++){
			ip1 = i+1; ip2 = i+2;
			hb = x[ip2]-x[ip1];
			deltab = (y[k*n+ip2]-y[k*n+ip1])/hb;
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
			coef[k][2+i*4] = (y[k*n+i+1]-y[k*n+i])/ha-(s2[i+1]+2*s2[i])*(ha/6);
			coef[k][3+i*4] = y[k*n+i];
		}
	}
}

void updateCubicSplineVelModel( float* slow, /* Slowness vector */
		    		int* n, /* n[0]=n1 n2=n[1] */
		    		float* o, /* o[0]=o1 o[1]=o2 */
		    		float* d, /* d[0]=d1 d[1]=d2 */
			    	int dim, /* Dimension of (z,vz) vectors */
				float* sz, /* Spline not Depth coordinates */
				float* sv, /* Spline not Velocity coordinates */
				float gzbg, /* Background velocity gradient in z */
				float v0, /* Near surface velocity */
				int n_stripes)
/*< Funcion to update spline cubic velocity model:
Note:
Make a velocity varying with depth model using the spline cubic interpolation
for a set of points (z,vz) given. TODO

 >*/
{
	int i, j=0, k;
	//int ic;
	float z=0.0;
	//float** coef;
	float v[n_stripes][n[0]];
	int app, app_len=n[1]/n_stripes;

	//coef = sf_floatalloc2(4*(dim-1),n_stripes);

	/* Calculate spline coeficients */
	//calculateSplineCoeficients(dim,sz,sv,coef,n_stripes);

	/* Calculate vel(city function */
	for(k=0;k<n_stripes;k++){

		z = o[0];
		j = 0;

		for(i=1;i<dim;i++){
			
			//ic = (i-1)*4;

			while(z<=sz[i]){
				z = (j*d[0]+o[0])-sz[i-1];
				if(j>=n[0]) break;
				//v[k][j] = coef[k][0+ic]*z*z*z+coef[k][1+ic]*z*z+coef[k][2+ic]*z+coef[k][3+ic];
				v[k][j] = v0+gzbg*z+sv[(k*dim)+i-1];
				//v[k][j] = sv[(k*dim)+i-1];
				j++;
			}
		}
	}

	/* Update slowness model */
	for(k=0;k<n_stripes;k++){

		app = (k*app_len);

		for(i=0;i<n[0];i++){

			for(j=app;j<((k+1)*app_len);j++){
				/* TODO Use 2D eno interpolation to obtain velocity model*/
				slow[j*n[0]+i]=1./(v[k][i]*v[k][i]);
				//#ifdef GDB_DEBUG
				//sf_warning("[%d][%d][%d] %f %f ",i,j,k,v[k][i],slow[j*n[0]+i]);
				//#endif
			} /* Loop over distance */
		} /* Loop over depth */
	} /* Loop over cubic spline functions */
	//#ifdef GDB_DEBUG
	//sf_error("fim");
	//#endif
}


float creTimeApproximation(float h, 
			 float m,
			 float v0,
			 float t0,
			 float m0,
			 float RNIP,
			 float BETA,
			 bool cds)
/*< CRE traveltime approximation t(m,h)
Note: If cds parameter is false, it uses the CRE formula to calculate time.
If cds parameter is true, it uses the non-hyperbolic CRS formula with CDS condition (RN=RNIP)
to calculate time.
>*/
{ 
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
			   float *a, /* Normal ray angle for each NIP source (degrees) */
			   int ns /* Number of NIP sources */)
/*< Return time missfit sum of source-NIP-receiver rays 

Note: This function traces nr reflection rays pairs from each NIP source
position passed though s matrix. It also calculate the difference between
the traveltime of the traced rays with calculated traveltime using CRE
traveltime approximation to calculate the time misfit returned by the function.
 >*/
{

	float currentRayAngle; // Emergence angle from source (radians)
	int i, ir, it, is; // loop counters
	float p[2]; // slowness vector
	float t; // Ray traveltime
	float nrdeg; // Emergence angle in degrees
	int nt=NT; // number of time samples in each ray
	int nr=NR; // number of reflection ray pairs for each source
	float dt=DT; // time sampling of rays
	raytrace rt; // raytrace struct
	float** traj; // Ray trajectory (z,x)
	float m; // CMP
	float h; // half-offset
	float tmis=0; // time misfit
	float xs; // Source position
	float xr; // Receiver position
	float tr; // NIP to receiver ray traveltime
	float ts; // NIP to source ray traveltime
	float *x; // Source position (z,x)

	x = sf_floatalloc(2);

	for(is=0;is<ns;is++){

		x[0]=s[is][0];
		x[1]=s[is][1];
		nrdeg = a[is]; // angle in degree

		for(ir=0;ir<nr;ir++){

			for(i=0; i<2; i++){

				/* initialize ray tracing object */
				rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);

				traj = sf_floatalloc2(2,nt+1);
				
				/* initialize ray direction */
				if(i==0){
					// NIP to source ray
					currentRayAngle=(nrdeg-(ir+1)*DANGLE)*DEG2RAD;
				}else{
					// NIP to receiver ray
					currentRayAngle=(nrdeg+(ir+1)*DANGLE)*DEG2RAD;
				}

				p[0] = -cosf(currentRayAngle);
				p[1] = sinf(currentRayAngle);

				/* Ray tracing */
				it = trace_ray (rt, x, p, traj);

				if(it>0){
					t = it*dt;
					if(i==0){ // Keep NIP to source ray traveltime
						ts=t;
						xs=x[1];
					}else{ // Keep NIP to receiver ray traveltime
						tr=t;
						xr=x[1];
					}
				}else if(it == 0){ // Ray endpoint inside model
					t = abs(nt)*dt;
					nt += 1000;
				}else{ // Side or bottom ray
					/* TODO to correct the way you treat side rays */
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

	/* L2 norm to evaluate the time misfit */
	tmis = sqrt(tmis*tmis);
	return tmis;
}

void interpolateVelModel(  int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float *sv, /* Velocity model disturbance */
			   float *sz, /* Depth coordinates of sv vector */
			   float *slow, /* Velocity model */
			   int nsz, /* sv n1 dimwnsion */
			   int nsx, /* sv n2 dimension */
			   float v0, /* Near surface velocity */
			   float gzbg /* Depth velocity gradient */)
/*< Velocity model interpolation
Note: This function uses a sv control points grid to obtain the complete
velocity model matrix through eno 2D interpolation. The sv vector is the
velocity disturbance of a constant velocity depth gradient model, that
velocity increases linearly with depth for gzbg gradient given.
 >*/
{

	int i, j;

	/* Calculate velocity function */
        for(j=0;j<nsx;j++){

                for(i=0;i<nsz;i++){

			sv[(j*nsz)+i] = v0+gzbg*sz[i]+sv[(j*nsz)+i];
                }

	}

	/* Interpolate velocity matrix */
	enoInterpolation2d(n,o,d,sv,slow,nsz,nsx);
}

void interpolateSlowModel( int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float *sv, /* Velociy disturbance */
			   float *sz, /* Depth coordinate of disturbance */
			   float *slow, /* Slowness model */
			   int nsz, /* n1 dimension of sv */
			   int nsx, /* n2 dimension of sv */
			   float v0, /* Near surface velocity */
			   float gzbg /* Background gradient in depth */)
/*< Slowness model interpolation
Note: This function uses a sv control points grid to obtain the complete
slowness model matrix through eno 2D interpolation. The sv vector is the
velocity disturbance of a constant velocity depth gradient model, that
velocity increases linearly with depth for gzbg gradient given.
 >*/
{

	int i, nm; // Loop counters and indexes

	interpolateVelModel(n, o, d,sv,sz,slow,nsz,nsx,v0,gzbg);

	/* transform velocity to slowness */
	nm =n[0]*n[1];
	for(i=0;i<nm;i++){
			slow[i] = 1.0/(slow[i]*slow[i]);
	}
}

void enoInterpolation2d(int *n, /* Interpolated vector dimension n1=n[0] n2=n[1] */
			float *o, /* Interpolated vector  axis origin o1=o[0] o2=o[1] */
			float *d, /* Interpolated vector sampling d1=d[0] d2=d[1] */
			float *ov, /* Original vector to interpolate */
			float *iv, /* Interpolated vector */
			int nov1, /* Original vector n1 dimension */
			int nov2 /* Orignanl vector n2 dimension */)
/*< Eno interpolation 2D function
Note: This function interpolates a vector increasing the number of
samples in the interpolated vector using eno interpolation. This vector
is a 2D matrix stored in a vector ov by columns (ov[j*n1+i]) and the new
vector iv will be the interpolated vector.
 >*/
{

	sf_eno2 map;
	float f2[2];
	int i, j, i1, i2;
	float x, y;

	map = sf_eno2_init(3,nov1,nov2);

	sf_eno2_set1(map,ov);

        for(i2=0;i2<n[1];i2++){

                for(i1=0;i1<n[0];i1++){
                        x = i1*d[0]+o[0]; i=x; x -= i;
                        y = i2*d[1]+o[1]; j=y; y -= j;
                        sf_eno2_apply(map,i,j,x,y,&iv[i2*n[0]+i1],f2,FUNC);
                }
        }
        sf_eno2_close(map);
}
