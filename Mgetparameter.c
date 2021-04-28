/* Get a parameter from parameters cube for (t0,m0) pairs given

Use beta=y if you want to get BETA parameter with automatic conversion
to degrees
 */

#include <rsf.h>

#define RAD2DEG 180./SF_PI /* Degrees to radians conversion */

int main(int argc, char* argv[])
{
	bool verb; // Verbose parameter
	bool beta;
	int nm0, nt0, np1, np2;
	float dp1, dp2, op1, op2;
	float* a, *t0, *m0;
	float** par;
	int i, it0, im0;
	sf_file in, t0s, m0s, parfile;

	sf_init(argc,argv);

	in = sf_input("in");
	m0s = sf_input("m0s");
	parfile = sf_output("out");
	t0s = sf_input("t0s");

	/* Velocity model: get 2D grid parameters */
	if(!sf_histint(m0s,"n1",&nm0)) sf_error("No n1= in m0s file");
	if(!sf_histint(t0s,"n1",&nt0)) sf_error("No n1= in t0s file");
	
	if(!sf_histint(in,"n1",&np1)) sf_error("No n1= in input");
	if(!sf_histint(in,"n2",&np2)) sf_error("No n2= in input");
	if(!sf_histfloat(in,"d1",&dp1)) sf_error("No d1= in input");
	if(!sf_histfloat(in,"d2",&dp2)) sf_error("No d2= in input");
	if(!sf_histfloat(in,"o1",&op1)) sf_error("No o1= in input");
	if(!sf_histfloat(in,"o2",&op2)) sf_error("No o2= in input");

	if(!sf_getbool("verb",&verb)) verb=true;
	/* verbose */

	if(!sf_getbool("beta",&beta)) beta=false;
	/* Get beta parameter */

	t0 = sf_floatalloc(nt0);
	sf_floatread(t0,nt0,t0s);
	m0 = sf_floatalloc(nm0);
	sf_floatread(m0,nm0,m0s);
	a = sf_floatalloc(nt0);
	par = sf_floatalloc2(np1,np2);
	sf_floatread(par[0],np1*np2,in);

	sf_putint(parfile,"n1",nt0);
	sf_putint(parfile,"n2",1);

	for(i=0;i<nt0;i++){
		it0 = (int)((t0[i]-op1)/dp1);
		im0 = (int)((m0[i]-op2)/dp2);
		
		if(beta){
			a[i] = (par[im0][it0])*RAD2DEG;
			a[i]=180-a[i];
		}else{
			a[i] = par[im0][it0];
		}
	}

	sf_floatwrite(a,nt0,parfile);
}
