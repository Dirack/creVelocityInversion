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
				int n_stripes)
/*< Funcion to update spline cubic velocity model:
Note:
Make a velocity varying with depth model using the spline cubic interpolation
for a set of points (z,vz) given. TODO

 >*/
{
	int i, j=0, ic, k;
	float z=0.0;
	float** coef;
	float v[n_stripes][n[0]];
	int app, app_len=n[1]/n_stripes;

	coef = sf_floatalloc2(4*(dim-1),n_stripes);

	/* Calculate spline coeficients */
	calculateSplineCoeficients(dim,sz,sv,coef,n_stripes);

	/* Calculate velocity function */
	for(k=0;k<n_stripes;k++){

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
	for(k=0;k<n_stripes;k++){

		app = (k*app_len);

		for(i=0;i<n[0];i++){

			for(j=app;j<((k+1)*app_len);j++){
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
	/* TODO buil another way to define nt and nr and to define its best values */
	int nt=5000; // TODO nt is the number of time samples in each ray
	int nr=5; //TODO nr is the number of ray pairs for each source
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

				traj = sf_floatalloc2(2,nt+1);
				
				/* initialize ray direction */
				/* TODO this part is confusing */
				currentRayAngle=(i==0)?(nrdeg-(ir+1)*DANGLE)*DEG2RAD:(nrdeg+(ir+1)*DANGLE)*DEG2RAD;

				p[0] = -cosf(currentRayAngle);
				p[1] = sinf(currentRayAngle);

				/* Ray tracing */
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

	/* TODO: Evaluate the best function to calcullate the time misfit */
	tmis = (tmis*tmis)/(nr*ns);
	return tmis;
}

void interpolateVelModel(  int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float *slow /* Slowness velociy model */)
/*< TODO to finish this function >*/
{

	int nt=5000;
	float dt=0.001;
	raytrace rt;
	sf_eno2 e2;
	int m, i;
	float x={1.0,1.0};

	/* initialize ray tracing object */
	//rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
	e2 = sf_eno2_init(4,n[0],n[1]);
	sf_eno2_set1(e2,slow);
	sf_eno2_close(e2);

	//sf_warning("%f",grid2_vel(rt->grd2,x));
	//m = n[0]*n[1];

	//for(i=0;i<m;i++)
	//	slow[i] = rt->grd2->pnt->ent->diff[0][i];
}
