/*
 *
 *	a set of functions to print the P-wave radiation pattern
 *		on a line printer
 *
 *	After a code from Jeoren Ritsema and Hong Kie Thio
 *
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
 
#define sqrt2 1.414213562
#define x_to_y_ratio 1.85
#define deg_to_rad 0.017453292
#define rad_to_deg 57.29577951
#define kMAXchar 72
#define kBoundary 0.95
#define kNrows 29

#define N 0
#define E 1
#define D 2

void GetlprPradiation(FILE *outStream,char thePattern[kMAXchar][kMAXchar],float Mij[3][3],int *,int *,float plnt,float stkt,float plnp,float stkp);
void GetAmplitudes(float Mij[3][3],float ,float ,float *,float *,float *);
void PrintPradiation(FILE *outStream, float Mij[3][3], float plnt,float stkt,float plnp,float stkp);
void PrintMTensor(FILE *outStream, float Mij[3][3]);
void CompactPattern(char thePattern[kMAXchar][kMAXchar],int *,int *);

void PrintPradiation(FILE *outStream, float Mij[3][3],float plnt,float stkt,float plnp,float stkp)
{
	char thePattern[kMAXchar][kMAXchar];
	int i,j,nlines,ncols;

	nlines = kNrows;

	GetlprPradiation(outStream,thePattern,Mij,&nlines,&ncols, plnt, stkt, plnp, stkp);
	
	for(i=0;i<nlines;i++)
	{
		for(j=0;j<ncols;j++)
		{
			fprintf(outStream,"%c",thePattern[i][j]);
		}
		fprintf(outStream,"\n");
	}



}

/*
 *	Function to return a character array containing the P-wave
 *		radiation pattern for the input moment tensor.
 *
 *	nlines is the number of lines the radiation pattern will fill.
 *		the width of the pattern is returned in ncols
 *
 * 	the character array is filled with two symbols:
 *  	# -> positive amplitudes
 *		- -> negative ampltiudes
 *
 *
 */

#define ABS(a)   ( (a) > (0) ? (a): -(a) )
void
GetlprPradiation(FILE *outStream, char thePattern[kMAXchar][kMAXchar],float Mij[3][3], int *nlines, int *ncols,float plnt,float stkt,float plnp,float stkp)
{
 
	int nr,nc,i,j;
	float pw,svw,shw,theta,phi;
	double x0,y0,r0,xx,yy,r,temp;
	float degrad = 3.141592653/180.0 ;
	int ip, jp, it, jt;
	float misp, mist;
 	
	if(stkt > 180)stkt -= 360 ;
	if(stkp > 180)stkp -= 360 ;
	/* clear out the character image */
	for(i=0;i<kMAXchar;i++)
	{
		for(j=0;j<kMAXchar;j++)
		{
			thePattern[i][j] = ' ';
		}
	}

	misp = 1.0e+38;
	mist = 1.0e+38;
	
	/* 
	 *	set up the scale for the radiation pattern image 	
	 * 		used the fixed aspect ratio 
	 */
	nr = *nlines;
	nc = (int) ((float) nr * x_to_y_ratio);
	/* make sure the number of columns is odd */
	nc = (2 * (int)(nc/2) + 1);
	*ncols = nc;
	
	/* set up the origin of the focal sphere */
	x0 = (((double) (nc) * 0.5 + 1) / x_to_y_ratio);
	y0 =  ((double) (nr) * 0.5);
	/* set the radius of the focal sphere */
	r0 = sqrt(x0*x0 + y0*y0) / 1.85;
	
	/* now fill in the negative polarity amplitudes */
	for(i=0;i<nr;i++)
	{
		 yy = y0 - (double) (i);
		 
		 for(j=0;j<nc;j++)
		 {
		 	xx = ((double) (j) / x_to_y_ratio) - x0;
		 	
		 	r = sqrt(xx*xx + yy*yy) / r0;
		 	
		 	if(r < 1)
		 	{
		 		/*
		 		 * 	phi is the azimuth, theta is the takeoff angle
		 		 */
		 		phi = atan2(xx,yy) * rad_to_deg;
		 		
		 		/* for a stereographic projection */
		 		/* theta = 2*atan(r) * rad_to_deg; */
		 		/* for an equal area projection */
		 		theta = 2*asin(r/sqrt2) * rad_to_deg;
				if(ABS(phi - stkt)+ABS(90 - theta - plnt) < mist){
					mist = ABS(phi - stkt)+ABS(90-theta - plnt);
					it = i;
					jt = j;
				}
				if(ABS(phi - stkp)+ABS(90 - theta - plnp) < misp){
					misp = ABS(phi - stkp)+ABS(90-theta - plnp);
					ip = i;
					jp = j;
				}
		 				 		
		 		/* get the body wave radiation amplitudes */
		 		GetAmplitudes(Mij,phi,theta,&pw,&svw,&shw);
		 		
		 		/* assign a character based on the P-wave polarity */
		 		if(pw > 0)
		 			thePattern[i][j] = '#';
		 		else
		 			thePattern[i][j] = '-';
		 	}
		 }
	}
	for(i = ip -1 ; i < ip +2; i++)
		for(j=jp -1 ; j < jp +2 ; j++)
			thePattern[i][j] = ' ';
	thePattern[ip][jp] = 'P';
	for(i = it -1 ; i < it +2; i++)
		for(j=jt -1 ; j < jt +2 ; j++)
			thePattern[i][j] = ' ';
	thePattern[it][jt] = 'T';
	/* now handle the P and T axes */
	
	
}
/*
 *	Given a moment tensor stored in Mij
 *
 *	& a tekeoff angle theta and an azimuth phi (in degrees)
 *
 *	Return the P,SV,SH amplitude.
 *
 */

void
GetAmplitudes(float Mij[3][3],float phi,float theta,float *pw,float *svw,float *shw)
{
	float sth, cth, cphi, sphi, p_dir[3],sv_dir[3],sh_dir[3];
	int i,j;
	
	
	/* first compute some direction cosines */
	
	sth = sin(theta*deg_to_rad);
	cth = cos(theta*deg_to_rad);
	cphi = cos(phi*deg_to_rad);
	sphi = sin(phi*deg_to_rad);
	
	/* now construct the ray vectors */
	p_dir[N] = sth*cphi;
	p_dir[E] = sth*sphi;
	p_dir[D] = cth;

	sv_dir[N] = cth*cphi;
	sv_dir[E] = cth*sphi;
	sv_dir[D] = -sth;	

	sh_dir[N] = -sphi;
	sh_dir[E] = cphi;
	sh_dir[D] = 0;
	
	/* now compute the radiated amplitudes */
	
	*pw = 0; *svw = 0; *shw = 0;
	
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			*pw  += Mij[i][j]*p_dir[i]* p_dir[j];
			*svw += Mij[i][j]*p_dir[i]*sv_dir[j];
			*shw += Mij[i][j]*p_dir[i]*sh_dir[j];
		}
	}

}


/*
 *	Function to print out the moment tensor
 *
 */
 
void PrintMTensor(FILE *outStream, float Mij[3][3])
{

	int i,j;
	
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			fprintf(outStream,"%10.2e ",Mij[i][j]);
		}
		fprintf(outStream,"\n");
	}
		
	fprintf(outStream,"\n");
	

}
#undef sqrt2
#undef x_to_y_ratio
#undef rad_to_deg
#undef deg_to_rad
#undef kMAXchar
#undef kBoundary
#undef kNrows
#undef N
#undef E
#undef D

