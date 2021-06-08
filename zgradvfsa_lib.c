/*
	 zgradvfsa_lib.c (c)
	 
	 Purpose: VFSA functions library of 'Mzgradtomo.c'.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 07/06/2021

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#define MAX 0.8
#define MIN 0
#define APERTURE 0.2
#define hMAX 50
#define mMAX 50
#define ITMAX 3000
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <rsf.h>
/*^*/

#define signal(s) ((s<0)?(-1.):(1.))
/*< Signal function >*/
/*^*/

float getRandomNumberBetween0and1(){
/*< Function to get a random number between 0 and 1 >*/

	return (float)(rand()%1000)/1000;
}

float getVfsaIterationTemperature(int iteration,float dampingFactor,float inicialTemperature){
/*< Temperature function for VFSA algorithm >*/

	return inicialTemperature*expf(-dampingFactor*pow(iteration,0.25));

}

void disturbGradZ( 	float temperature, /* Temperature of this VFSA interation */
			float* gznew, /* Grad z disturbed */
			float gz, /* Grad z to disturb */
			float scale /* Scale to multiply by disturbance */)
/*< Disturb depth velocity gradient from VFSA previous iteration >*/
{

	float u;
	float disturbance;

	u=getRandomNumberBetween0and1();
			
	disturbance = signal(u - 0.5) * temperature * (pow( (1+temperature),fabs(2*u-1) )-1);

	*gznew = gz + (disturbance*scale) * (APERTURE);

	if (*gznew >= MAX) {

		*gznew = MAX - (APERTURE) * getRandomNumberBetween0and1();
		
	}

	if (*gznew <= MIN) {

		*gznew = (APERTURE) * getRandomNumberBetween0and1() + MIN;
		
	}

}

void updateGzVelModel(float* slow, /* Slowness vector */
		    int* n, /* n[0]=n1 n2=n[1] */
		    float* o, /* o[0]=o1 o[1]=o2 */
		    float* d, /* d[0]=d1 d[1]=d2 */
		    float v0, /* First velocity is the near surface velocity */
		    float grad /* Velocity gradient */)
/*< Function to update constant depth velocity gradient:
Note:
This is a function to update the velocity model,
it uses a constant gradient velocity model, and update the 
gradient in each iteration of the process.

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

