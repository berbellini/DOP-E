

/* originally developed by C J Ammon Pennsylvania State University
   modified by RBHermann of Saint Louis University to make more compatible with
   Berkeley annouincement */

#include <stdio.h>
#include <math.h>
void PrintPradiation(FILE *outStream, float Mij[3][3], float plnt,float stkt,float plnp,float stkp);
void phase(float x,float y,float z,float *del,float *bet);
void trans(float dip,float stk,float slip,float *stkt,float *plnt,float *stkp,float *plnp, float *plnn, float *stkn);
void tpdss(float stkt,float plnt,float stkp,float plnp,
	float *stk0, float *dip0, float *rak0,
	float *stk1, float *dip1, float *rak1) ;

main(ac,av)
int ac;
char *av[];
{
	float theStrike,theDip,theRake,m[3][3];
	float theMw, theDepth;
	float Mo;
	int i,j;
	
   if(ac != 7)
	 {
	    fprintf(stdout,"usage: lprmech TITLE strike dip rake MW depth\n\n");
			exit(1);
	 }
	 
	 sscanf(av[2],"%f",&theStrike);
	 sscanf(av[3],"%f",&theDip);
	 sscanf(av[4],"%f",&theRake);
	 sscanf(av[5],"%f",&theMw);
	 sscanf(av[6],"%f",&theDepth);
	 Mo = 1.5*theMw + 16.05;
	 Mo = pow (10.0, Mo);

	 fprintf(stdout," SLU Moment Tensor Solution\n");
         fprintf(stdout," %s\n",av[1]);
	 fprintf(stdout," \n");
	 fprintf(stdout," Best Fitting Double Couple\n");
	 fprintf(stdout,"    Mo = %8.2e dyne-cm\n",Mo);
	 fprintf(stdout,"    Mw = %4.2f \n",theMw);
	 fprintf(stdout,"    Z  = %d km\n",(int)theDepth);

	 float stkt, stkp, stkn, plnt, plnp, plnn;
	 float stk0, dip0, rak0;
	 float stk1, dip1, rak1;
trans( theDip, theStrike, theRake, &stkt, &plnt, &stkp, &plnp,  &plnn,  &stkn);
tpdss( stkt, plnt, stkp, plnp,
	 &stk0,  &dip0,  &rak0,
	 &stk1,  &dip1,  &rak1) ;
	fprintf(stdout,"     Plane   Strike  Dip  Rake\n");
	fprintf(stdout,"      NP1      %3.0f   %3.0f   %3.0f\n",stk0,dip0,rak0);
	fprintf(stdout,"      NP2      %3.0f   %3.0f   %3.0f\n",stk1,dip1,rak1);
	fprintf(stdout," Principal Axes:\n");
	fprintf(stdout,"   Axis    Value   Plunge  Azimuth\n");
	fprintf(stdout,"     T  %9.2e    %3.0f     %3.0f\n",Mo,plnt,stkt);
	fprintf(stdout,"     N  %9.2e    %3.0f     %3.0f\n",0.0,plnn,stkn);
	fprintf(stdout,"     P  %9.2e    %3.0f     %3.0f\n",- Mo,plnp,stkp);

	 
   sdr_to_mij(theStrike,theDip,theRake,m);
	 for(i=0 ; i < 3 ; i++)
		 for(j=0;j<3;j++)
			 m[i][j] *= Mo;

	 fprintf(stdout,"\n");
	 fprintf(stdout,"\n");
	 fprintf(stdout,"\n");
	 fprintf(stdout," Moment Tensor: (dyne-cm)\n");
	 fprintf(stdout,"    Component  Value\n");
	 fprintf(stdout,"       Mxx    %9.2e\n",m[0][0]);
	 fprintf(stdout,"       Mxy    %9.2e\n",m[0][1]);
	 fprintf(stdout,"       Mxz    %9.2e\n",m[0][2]);
	 fprintf(stdout,"       Myy    %9.2e\n",m[1][1]);
	 fprintf(stdout,"       Myz    %9.2e\n",m[1][2]);
	 fprintf(stdout,"       Mzz    %9.2e\n",m[2][2]);

	 
	 PrintPradiation(stdout, m,plnt,stkt,plnp,stkp);

   ar_to_hrv_mij(m);
	 
	 fprintf(stdout,"\n Harvard Convention\n Moment Tensor:\n");
	 fprintf(stdout,"      R          T          F\n");
	 PrintMTensor(stdout, m);


}
#define ABS(a)   ( (a) > (0) ? (a): -(a) )
#define SGN(a,b) ( (b) > (0) ? (ABS(a)): -(ABS(a)))

void trans(float dip,float stk,float slip,float *stkt,float *plnt,float *stkp,float *plnp, float *plnn, float *stkn)
{
/*
c-----
c     This takes strike, dip, and slip angles and defines a coordinate
c     transformation in which the slip vector, fault normal vector
c     and the normal vector of the auxiliary plane in terms of the
c     local NS, EW, and vertical cartesian coordinate system.
c-----
c     dip   angle of dip measured from horizontal plane.
c           0 < dip < 90
c     stk   direction of nodal plane strike. When looking
c           along stirke, fault plane dips downward to the
c           right
c           0 < stk < 360
c     slip  direction of movement along nodal plane.
c           If -180 < slip < 180, then the P-wave
c           vertically downward will a first motion polarity
c           equal to the sign of slip
c-----
c
c     REFERENCE:
c
c     Herrmann, R. B. (1975). A student's guide to the use of
c     P and S wave data, Earthquake Notes 46, 29-39
c
c-----
c     The X, Y and Z axes form a right handed coordinate
c     system, such that the compressional quadrant is
c     occurs whenever the xy > 0. The Z axis is the
c     null axis. The Y-axis is normal to the fault
c     plane and the X-axis is normal to the auxiliary plane
c     Note that for the input convention on dip, the
c     Y-axis will always point in the negative z-direction
c
c-----
*/
	float degrad, sins, coss, sind, cosd, sinf, cosf;
	float a11, a12, a13, a21, a22, a23, a31, a32, a33;
	float t1, t2, t3, p1, p2, p3;
	float xy, z3;
	float tmp1, tmp2, tmp3;
	float delx, dely, delz, betx, bety, betz;
	float delt, delp, bett, betp;
	float pnorm, tnorm, xnorm, ynorm, znorm;
	degrad=0.017452329;
	sins=sin(stk*degrad);
	coss=cos(stk*degrad);
	sind=sin(dip*degrad);
	cosd=cos(dip*degrad);
	sinf=sin(slip*degrad);
	cosf=cos(slip*degrad);
	/*
		X-axis
	*/
	a11=cosf*coss+sinf*cosd*sins;
	a12=cosf*sins-sinf*cosd*coss;
	a13=-sinf*sind;
	/*
		Y-axis
	*/
	a21=-sins*sind;
	a22=coss*sind;
	a23=-cosd;
	/*
		Z-axis
	*/
	a31=coss*sinf-cosd*cosf*sins;
	a32=cosd*cosf*coss+sinf*sins;
	a33=cosf*sind;
	t1=a11+a21;
	t2=a12+a22;
	t3=a13+a23;
	tnorm=sqrt(2.);
	if(t3 <= 0.0) tnorm=-tnorm;
	t1=t1/tnorm;
	t2=t2/tnorm;
	t3=t3/tnorm;
	p1=a11-a21;
	p2=a12-a22;
	p3=a13-a23;
	pnorm=sqrt(2.);
	if(p3 <= 0.0) pnorm=-pnorm;
	p1=p1/pnorm;
	p2=p2/pnorm;
	p3=p3/pnorm;
	xnorm=1.0;
	if(a13 < 0.0) xnorm=-1.0;
	ynorm=1.0;
	if(a23 < 0.0) ynorm=-1.0;
	znorm=1.0;
	if(a33 < 0.0) znorm=-1.0;
	a11=a11/xnorm;
	a12=a12/xnorm;
	a13=a13/xnorm;
	a21=a21/ynorm;
	a22=a22/ynorm;
	a23=a23/ynorm;
	a31=a31/znorm;
	a32=a32/znorm;
	a33=a33/znorm;
/*
c------
c     now all vectors point into the lower hemisphere.
c     To get the correct orientations, we require
c     that if the center of the focal sphere is a compression
c     that the X,Y, Z axes form a right handed coordinate system in the
c     lower hemisphere, otherwise it will form a left handed
c     coordinate system
c-----
c     obtain P-wave polarity at the center of the focal
c     sphere
c-----
*/
	xy = a13*a23;
	/*
		determine if right handed or left handed coordinate system
	*/
	z3 = a11*a22 - a12*a21;
	/*
		make right handed coordinate system
	*/
	if(z3 < 0.0){
		tmp1=a11 ;
		tmp2=a12 ;
		tmp3=a13 ;
		a11=a21 ;
		a12=a22 ;
		a13=a23 ;
		a21=tmp1 ;
		a22=tmp2 ;
		a23=tmp3 ;
	}
	if(SGN(1.0,xy) != SGN(1.0,slip)){
		tmp1=a11 ;
		tmp2=a12 ;
		tmp3=a13 ;
		a11=a21 ;
		a12=a22 ;
		a13=a23 ;
		a21=tmp1 ;
		a22=tmp2 ;
		a23=tmp3 ;
	}
	phase(a11,a12,a13,&delx,&betx);
	phase(a21,a22,a23,&dely,&bety);
	phase(a31,a32,a33,&delz,&betz);
	phase(t1,t2,t3,&delt,&bett);
	phase(p1,p2,p3,&delp,&betp);
/*	plx = 90-betx;
	ply = 90-bety;
	plz = 90-betz;
	plt=90-bett;
	plp=90-betp;
*/
	*stkt=delt;
	*plnt=90.0-bett;
	*stkp=delp;
	*plnp=90.0-betp;
	*stkn=delz;
	*plnn=90.0-betz;
}

void phase(float x,float y,float z,float *del,float *bet)
{
/*
c-----
c     This expresses a three component vector in terms
c     of two angles in a spherical coordinate system
c-----
c     x     x-coordinate of vector
c     y     y-coordinate of vector
c     z     z-coordinate of vector
c-----
c     del   direction of projection of the vector on a
c           horizontal plane, measured positively, clockwise
c           from north - the trend of the vector.
c     bet   angle that the 3-D vector makes with the
c           downward vertical. This is just 90 - plunge
c-----
*/
	float degrad;
	degrad=0.017452329;
	*bet=acos(z);
	*bet/=degrad;
	if(x < 0.0){
		*del=atan((y/x))/degrad+180. ;
	} else if(x == 0.0){
		if(y < 0.0)
			*del = 270;
		else if(y == 0.0)
			*del = 0;
		else
			*del = 90.;
	} else {
		if(y < 0.0)
			*del = atan((y/x))/degrad + 360.0;
		else if(y >= 0.0)
			*del=atan((y/x))/degrad;
		
	}
}

void tpdss(float stkt,float plnt,float stkp,float plnp,
	float *stk0, float *dip0, float *rak0,
	float *stk1, float *dip1, float *rak1) 
{
/*
	This takes the orientations of the tension and pressure 
	axes, and determines the two nodal planes 
	convert (P,T) to (dip,slip,stk) and (X,Y,Z). 
*/
	float t1,t2,t3,p1,p2,p3,f1,f2,f3,xn1,xn2,xn3 ;
	float deg;
	float stk,slip,dip,stkk,slipp,dipp,xx,yy ;
	deg=3.141592653/180.0 ;
/*
	strike angle is measured from the north clockwise. 
	plunge angle is measured from the horizon downward. 
*/
	stkp=stkp*deg ;
	plnp=plnp*deg ;
	stkt=stkt*deg ;
	plnt=plnt*deg ;
	t1=cos(plnt)*cos(stkt) ;
	t2=cos(plnt)*sin(stkt) ;
	t3=sin(plnt) ;
	p1=cos(plnp)*cos(stkp) ;
	p2=cos(plnp)*sin(stkp) ;
	p3=sin(plnp) ;
	f1=t1+p1 ;
	f2=t2+p2 ;
	f3=t3+p3 ;
	yy=sqrt(f1*f1+f2*f2+f3*f3) ;
	f1=f1/yy ;
	f2=f2/yy ;
	f3=f3/yy ;
	xn1=t1-p1 ;
	xn2=t2-p2 ;
	xn3=t3-p3 ;
	yy=sqrt(xn1*xn1+xn2*xn2+xn3*xn3) ;
	xn1=xn1/yy ;
	xn2=xn2/yy ;
	xn3=xn3/yy ;
	if(xn3 > 0.0){ 
		xn1=-xn1 ;
		xn2=-xn2 ;
		xn3=-xn3 ;
		f1=-f1 ;
		f2=-f2 ;
		f3=-f3 ;
	}
	dip=acos(-xn3) ;
	stk=atan2(-xn1,xn2) ;
	xx=f1*xn2-f2*xn1 ;
	slip=atan2(-f3,xx) ;
	*dip0=dip/deg ;
	*rak0=slip/deg ;
	if(*rak0 < -180.)*rak0+=180. ;
	if(*rak0 > 180.)*rak0-=360. ;
	*stk0=stk/deg ;
	if(*stk0 < 0.0) *stk0+=360. ;

	if(f3 > 0.0){
		xn1=-xn1 ;
		xn2=-xn2 ;
		xn3=-xn3 ;
		f1=-f1 ;
		f2=-f2 ;
		f3=-f3 ;
	}
	dip=acos(-f3) ;
	stk=atan2(-f1,f2) ;
	xx=f2*xn1-f1*xn2 ;
	slip=atan2(-xn3,xx) ;
	*dip1=dip/deg ;
	*rak1=slip/deg ;
	if(*rak1 < -180.)*rak1+=180. ;
	if(*rak1 > 180.)*rak1-=360. ;
	*stk1=stk/deg ;
	if(*stk1 < 0.0) *stk1+=360. ;
}
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
#define kNrows 15

#define N 0
#define E 1
#define D 2

void GetlprPradiation(FILE *outStream,char thePattern[kMAXchar][kMAXchar],float Mij[3][3],int *,int *);
void GetAmplitudes(float Mij[3][3],float ,float ,float *,float *,float *);
void PrintPradiation(FILE *outStream, float Mij[3][3]);
void PrintMTensor(FILE *outStream, float Mij[3][3]);
void CompactPattern(char thePattern[kMAXchar][kMAXchar],int *,int *);

void PrintPradiation(FILE *outStream, float Mij[3][3])
{
	char thePattern[kMAXchar][kMAXchar];
	int i,j,nlines,ncols;

	nlines = kNrows;

	GetlprPradiation(outStream,thePattern,Mij,&nlines,&ncols);
	
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

void
GetlprPradiation(FILE *outStream, char thePattern[kMAXchar][kMAXchar],float Mij[3][3], int *nlines, int *ncols)
{
 
	int nr,nc,i,j;
	float pw,svw,shw,theta,phi;
	double x0,y0,r0,xx,yy,r,temp;
 	
	/* clear out the character image */
	for(i=0;i<kMAXchar;i++)
	{
		for(j=0;j<kMAXchar;j++)
		{
			thePattern[i][j] = ' ';
		}
	}
	
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
		 				 		
		 		/* get the body wave radiation amplitudes */
		 		GetAmplitudes(Mij,phi,theta,&pw,&svw,&shw);
		 		
		 		/*printf("phi= %.0f		theta = %.0f		pw = %.2f\n",phi,theta,pw);*/
		 		/* assign a character based on the P-wave polarity */
		 		if(pw > 0)
		 			thePattern[i][j] = '#';
		 		else
		 			thePattern[i][j] = '-';
		 	}
		 }
	}
	
	
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
			fprintf(outStream,"%5.2f ",Mij[i][j]);
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

#include <stdio.h>
#include <math.h>


sdr_to_mij(s,d,r,m)
float s,d,r,m[3][3];
{
      float d2r,sd,cd,s2d,c2d,ss,cs,c2s,s2s,sr,cr;

      d2r = acos(-1) / 180;
					
      sd = sin(d2r * d);
      cd = cos(d2r * d);
      s2d = sin(d2r * 2 * d);
      c2d = cos(d2r * 2 * d);
      
      ss = sin(d2r * s);
      cs = cos(d2r * s);
      s2s = sin(d2r * 2 * s);
      c2s = cos(d2r * 2 * s);
      
      sr = sin(d2r * r);
      cr = cos(d2r * r);
      
      m[0][0] = -sd*cr*s2s - s2d*sr*ss*ss;
      m[1][1] =  sd*cr*s2s - s2d*sr*cs*cs;
			m[2][2] =  s2d*sr;
      m[0][1] =  sd*cr*c2s + 0.5*s2d*sr*s2s;
      m[0][2] = -cd*cr*cs  - c2d*sr*ss;
      m[1][2] = -cd*cr*ss  + c2d*sr*cs;
      
			m[1][0] = m[0][1];
			m[2][1] = m[1][2];
			m[2][0] = m[0][2];

      return;
}

ar_to_hrv_mij(m)
float m[3][3];
{
      float tmp[3][3];
      int x=0,y=1,z=2,r=0,t=1,f=2, i,j;
      
      tmp[r][r] =  m[z][z];
			tmp[r][t] =  m[z][x];
			tmp[r][f] = -m[z][y];
			
			tmp[t][r] =  tmp[r][t];
			tmp[t][t] =  m[x][x];
			tmp[t][f] = -m[x][y];
			
			tmp[f][r] =  tmp[r][f];
			tmp[f][t] =  tmp[t][f];
			tmp[f][f] =  m[y][y];
			
			for(i=0;i<3;i++)
			 for(j=0;j<3;j++)
				 m[i][j] = tmp[i][j];
			
      return;
}

hrv_to_ar_mij(m)
float m[3][3];
{
      float tmp[3][3];
      int x=0,y=1,z=2,r=0,t=1,f=2, i,j;
			
      tmp[x][x] =  m[2][2];
			tmp[x][y] = -m[t][f];
			tmp[x][z] =  m[t][r];
			
			tmp[y][x] =  tmp[x][y];
			tmp[y][y] =  m[f][f];
			tmp[y][z] = -m[f][r];
			
			tmp[z][x] =  tmp[x][z];
			tmp[z][y] =  tmp[y][z];
			tmp[z][z] =  m[r][r];
			
			for(i=0;i<3;i++)
			 for(j=0;j<3;j++)
				 m[i][j] = tmp[i][j];
			
      return;
}
