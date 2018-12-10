/* system includes */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* local includes for Computer Programs in Seismology */
#include "sacsubc.h"

/* define macross */
#define ABS(a  ) ( (a) >  0  ? (a):-(a))
#define YES 1
#define NO 0

/* procedure prototypes */
void saciterd_gfilter(float *x, float gwidth , int n, int nn, float dt,float *z1);
void saciterd_getres(float *x, float *y, int n, float *r, float *sumsq);
void saciterd_usage (void);
void saciterd_zero(float *x, int n );
void saciterd_getmax(float *x, int n, float * maxvalue, int *maxindex);
void saciterd_getabsmax(float *x,int n, float *thevalue,int *maxindex);
void gcmdln(int argc, char **argv,  char **fnum,char **fden,int *niter, 
	float *delay, float *error, int *dopos, float *gwidth,int *verbose, 
	int *dotwice, float *rayp, int *quiet) ;
void saciterd_build_decon(float *amps, int *shifts, int nshifts, float *p, int n, int nn, float gwidth,
	float dt, float *z1);
void saciterd_convolve(float *x, float *y, int n, int nn, float dt,float *z1,float *z2);
void saciterd_fcorrelate(float *x, float *y, int n, int nn, float dt, float *z1, float *z2);
void saciterd_phs_shift(float *x,float theshift,int  n,int  nn,float dt, float *z1);
void saciterd_tdomain(char *fnum,char *fden,float gwidth,float theshift,
          int verbose,int maxbumps,int lpositive,float tol,int dotwice,float rayp, int quiet);



/* the following are taken from gsac code */
int npow2(int n);
void four(float data[], int nn, int isign, float *dt, float *df);

/*
c---------------------------------------------------------------------c
c                                                                    c
c     COMPUTER PROGRAMS IN SEISMOLOGY                                c
c     VOLUME V                                                       c
c                                                                    c
c     PROGRAM: SACITERD                                              c
c                                                                    c
c     COPYRIGHT 1996 1998 CHARELS J AMMON                            c
c                                                                    c
c     Department of Earth and Atmospheric Sciences                   c
c     Saint Louis University                                         c
c     221 North Grand Boulevard                                      c
c     St. Louis, Missouri 63103                                      c
c     U. S. A.                                                       c
c     VERSION 1.04                                                   c
c                                                                    c
c---------------------------------------------------------------------c
c      Changes:
c
c      22 JAN 2002 - add -RAYP rayp to the command line so that
c          the ray parameter is placed in the USER4 SAC variable
c      21 NOV 2002 - modified KEVNM header as per 
c          Mejian An Dept Geophysics, IAG, Sao Paulo, Br   
c      09 JAN 2005 - nbumps undefined in  subroutine saciterd_tdomain, changed
c          to nshift. tshift undefined 
c          in subroutine saciterd_tdomain  changed the theshift, baker@usgs.gov
c      23 APR 2007 - repaired the doubling construct at line  line 490
c          which was incorrectly implemented on early 2007 - the code now
c          recreates the iamges of 2006-11-10
c      24 APR 2009 - permit ray parameter in E format
c      23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c      12 FEB 2011  - converted FORTRAN to C
c      15 AUG 2011 - introduced -Q (quiet flag) to suppress interation output
c-----
c
c    Based on the Kikuchi and Kanamori (1982) iterative
c     deconvolution algorithm. The output of this code
c     is the deconvolution of the "denominator" from
c     the "numerator"
c    Kikuchi, M., and H. Kanamori (1982). Inversion of complex body waves, 
c    Bull. Seism. Soc. Am. 72, 491-506.
c
c    Ligorria, J. P., and C. A. Ammon (1999). Iterative deconvolution and 
c    receiver-function estimation,  Bull. Seism. Soc. Am. 89, 1395-1400.
c
c    The final deconvolution is called "decon.out"
c      the observed time series (original convolved with Gaussian)
c       is called "observed"
c      the predicted is call "predicted"
c
c    Header values from the "numerator" are copied into the
c      decon.out file. The gwidth is stored in user0.
c
c-----
c  INTERNAL PARAMETERS
c
c  rhdr    = SAC waveform real header fields (70).
c  ihdr    = SAC waveform integer header fields (40).
c  chdr    = SAC waveform character header fields (48).
c  x(i)    = array containing the waveform data.
c  delta   = waveform sampling interval in seconds.
c  btime   = waveform begin time.
c-----
*/
main(int argc, char **argv)
{
	char *fnum, *fden;
	int maxbumps, lpositve, verbose, dotwice, lpositive, quiet;
	float theshift, tol, gwidth, rayp;
	/* parse the command line */
	gcmdln(argc,argv,&fnum,&fden,&maxbumps,&theshift,&tol,&lpositive,&gwidth,
		&verbose,&dotwice,&rayp,&quiet);
	printf("Numerator file   : %s\n",fnum);
	printf("Denominator file : %s\n",fden);
	printf("Maxbumps         : %d\n",maxbumps);
	printf("Gwidth           : %f\n",gwidth);
	printf("Time shift       : %f\n",theshift);
	printf("Tolerance        : %f\n",tol     );
	printf("Positivity       : %d\n",lpositive     );
	printf("Double length    : %d\n",dotwice     );
	printf("Ray parameter    : %f\n",rayp     );
	printf("Verbos           : %d\n",verbose);
	/* do the deconvolution */
	saciterd_tdomain(fnum,fden,gwidth,theshift,
          verbose,maxbumps,lpositive,tol,dotwice,rayp,quiet);

}

void saciterd_build_decon(float *amps, int *shifts, int nshifts, float *p, int n, int nn, float gwidth,
	float dt, float *z1)
{
/*
       compute the predicted time series from a set of
       amplitudes and shifts

       amps    R*4 - array of amplitudes
       shifts  R*4 - array of time shifts
       p   R*4 - array generated
       n   I*4 - number of points in array
       nn  I*4 - number of points rounded up to nearest power of 2
       gwidth  R*4 - filter parameter
       dt  R*4 - sampling interval
*/
        int i;

        saciterd_zero(p,nn);
        for(i=0;i< nshifts;i++){
            p[shifts[i]] = p[shifts[i]] + amps[i] ;
	}

}

char *saciterd_usage_help [] = {
"Usage: saciterd  -FN file_num -FD file_den [-E error] [-N niter]\n",
"    [-POS] [-D delay] [-ALP alpha] [-2] [-RAYP rayp] [-V] [-?] [-h]\n",
" -FN file_num (required ) numerator SAC binary\n",
" -FD file_den (required ) denominator SAC binary\n",
" -E error     (default 0.001) convergence criteria\n",
" -ALP alpha   (default 1.0) Gaussian Filter Width\n",
"                  H(f) = exp( - (pi freq/alpha)**2) \n",
"                  Filter corner ~ alpha/pi \n",
" -N niter     (default 100, 1000 max) Number iterations/bumps\n",
" -D delay     (default 5 sec) Begin output delay sec before t=0\n",
" -POS         (default false) Only permit positive amplitudes\n",
" -2           (default false) use double length FFT to \n",
"                  avoid FFT wrap around in convolution\n",
" -RAYP rayp   (default -12345.0) ray parameter in \n",
"                  sec/km used by rftn96/joint96\n",
" -V           Verbose output \n",
" -?           This help message\n",
" -h           This help message\n",
"Output files:\n",
"   observed   : original numerator convolved with Gaussian\n",
"   numerator  : original numerator convolved with Gaussian\n",
"   denominator: original numerator convolved with Gaussian\n",
"   decon.out  : Receiver function for Gaussian\n",
"   predicted  : Receiver function for Gaussian\n",
" \n",
"SAC header values set\n",
" B     :  delay                          \n",
" USERO :  gwidth        KUSER0:  Rftn    \n",
" USER4 :  rayp (sec/km) \n",
" USER5 :  fit in %      \n",
" KEVNM :  Rftn          KUSER1:  IT_DECON\n",
""
} ;

void saciterd_usage(void){
int i;
	for(i = 0 ; strlen(saciterd_usage_help[i]) > 0 ; i++ )
		printf("%s",saciterd_usage_help[i]);
	exit(0) ;
}

void saciterd_zero(float *x, int n ){
/*
	zero a real array 
	x   float - array to be zeroed
	n   int   - number of points in array
*/
	int i;
        for(i=0; i < n ; i++)
		x[i] = 0.0;
        
}

void saciterd_getmax(float *x, int n, float * maxvalue, int *maxindex)
{

/*
       get maximum value in array and return index to array position
       x   float * - array
       n   int    - number of points in array
       maxvalue float *    - maximum value
       maxindex int  *     - position of maximum value in array
*/
	int i;
        *maxvalue = x[0];
        *maxindex = 1 ;
	for(i=1;i< n; i++){
		if(x[i] > *maxvalue) {
			*maxvalue = x[i];
			*maxindex = i;
		}
	}
}

void saciterd_getabsmax(float *x,int n, float *thevalue,int *maxindex)
{
/*
       get maximum value in array and return index to array position
       x *   float - array
       n     int   - number of points in array
       thevalue float  *   - maximum value
       maxindex int  *     - position of maximum value in array

*/
	float maxvalue ;
	int i;
        *thevalue = x[0];
        maxvalue = ABS(x[0]);
        *maxindex = 1 ;
	for(i=1;i< n; i++){
		if(ABS(x[i]) > maxvalue) {
			*thevalue = x[i];
			maxvalue = ABS(x[i]);
			*maxindex = i;
		}
	}
}



void saciterd_getres(float *x, float *y, int n, float *r, float *sumsq)
{
/*
      get sum square residuals between x and y arrays

      x   float - array
      y   float - array
      n   int   - number of points
      r   float - array of residuals
      sumsq   float - sum of residuals squared
*/
	int i;
        
        *sumsq = 0 ;
        for(i=0 ; i < n ; i++){
            r[i] = x[i] - y[i] ;
            *sumsq +=  r[i]*r[i] ;
        }
}

void gcmdln(int argc, char **argv,  char **fnum,char **fden,int *niter, 
	float *delay, float *error, int *dopos, float *gwidth,int *verbose, 
	int *dotwice, float *rayp, int *quiet) {
/*
      parse command line arguments

      requires subroutine targ() and funtion mnmarg()

      fnum    Ch* - file name for numerator
      fden    Ch* - file name for denominator
      niter   I*4 - numebr of iterations
      delay   R*4 - time shift (seconds before P)
      error   R*4 - convergence criteria
      dopos   L   - .true. permit only positive bumps
      gwidth   R   - Gaussian filter parameter
      sacbin  L   - .true. if SAC binary output, else ascii
      verbose L   - .true. output many intermediate files
      dotwice L   - .true. use double length FFT to avoid wrap around
      rayp    R   - ray parameter in sec/km, default of -12345.
*/
	int nfile, i;

	/* set default values */
        *niter = 100 ;
        *delay = 5.0 ;
        *error = 0.001 ;
        *dopos = NO ;
        *gwidth = 1.0 ;
        *verbose = NO ;
        *quiet = NO ;
        *dotwice = NO ;
        *rayp = -12345. ;

	/* initialize internal parameters */
        nfile = 0 ;

	char *cp;
 	while(argc-- > 1) {
		cp = argv[1];
		if(strncmp(cp,"-FN",3) == 0 || strncmp(cp,"-fn",3) == 0){
			nfile++;
			argv++;
			argc--;
			*fnum = argv[1];
		} else if(strncmp(cp,"-FD",3) == 0 || strncmp(cp,"-fd",3) == 0){
			nfile++;
			argv++;
			argc--;
			*fden = argv[1];
		} else if(strncmp(cp,"-ALP",4) == 0 ){
			argv++;
			argc--;
			*gwidth = atof(argv[1]);
		} else if(strncmp(cp,"-E",2) == 0 ){
			argv++;
			argc--;
			*error = atof(argv[1]);
		} else if(strncmp(cp,"-D",2) == 0 ){
			argv++;
			argc--;
			*delay = atof(argv[1]);
		} else if(strncmp(cp,"-N",2) == 0 ){
			argv++;
			argc--;
			*niter = atoi(argv[1]);
		} else if(strncmp(cp,"-RAYP",5) == 0 ){
			argv++;
			argc--;
			*rayp = atof(argv[1]);
		} else if(strncmp(cp,"-POS",4) == 0 ){
			*dopos = YES;
		} else if(strncmp(cp,"-V",2) == 0 ){
			*verbose = YES;
		} else if(strncmp(cp,"-Q",2) == 0 ){
			*quiet = YES;
		} else if(strncmp(cp,"-2",2) == 0 ){
			*dotwice = YES;
		} else if(strcmp("h",cp) == 0 || strcmp("?",cp)==0){
			saciterd_usage();
		}
		argv++;
	}
	/* safety checks */
	if(nfile != 2)saciterd_usage();
}



void saciterd_convolve(float *x, float *y, int n, int nn, float dt,float *z1,float *z2)
{

/*
       convolve x and y, replacing the x array

       x   float - array
       y   float - array
       n   int   - number of points in time series
       nn  int   - number of points rounded up to next power of 2
       dt  float - sampling interval
       z1  float - complex work array
       z2  float - complex work array
*/

	/* internal variables */
	int i,j,n2, jr, ji, kr, ki;
	float df;
	float tmpr, tmpi;


	/* convolve  by frequency domain multiplication */
	for(i=0,j=0; i< nn; i++,j+=2){
		if(i < n){
			z1[j  ] = x[i] ;	/* real */
			z1[j+1] = 0.0 ;		/* imaginary */
			z2[j  ] = y[i] ;	/* real */
			z2[j+1] = 0.0;		/* imaginary */
		} else {
			z1[j  ] = 0.0;	/* real */
			z1[j+1] = 0.0 ;		/* imaginary */
			z2[j  ] = 0.0;	/* real */
			z2[j+1] = 0.0;		/* imaginary */
		}
	}
	/* get Fourier transforms */
        four(z1,nn,-1,&dt,&df) ;
        four(z2,nn,-1,&dt,&df) ;

	/* convolution = F  G */
        n2 = nn / 2 ;
	for(j=0,i=0; i <= n2 ; i++, j+=2){
		/* this can be simplified since the z1 and z2 arrays are both real 
			but full complex arithmetic is used here
		*/
		jr = j   ;
		ji = j+1 ;
		tmpr = z2[jr]*z1[jr] - z2[ji]*z1[ji];
		tmpi = z2[jr]*z1[ji] + z2[ji]*z1[jr];
                z1[jr] = tmpr;
		z1[ji] = tmpi;
		if(i > 0){
			/* apply symmetry to make the convolution real */
			kr = 2*nn +2 - jr -2;
			ki = 2*nn +2 - jr -1   ;
			z1[kr] =   z1[jr] ;
			z1[ki] = - z1[ji] ;
		}
	}
	/* ensure Nyquist frequency element is real */
	/*  preserve the real part, change sign of imaginary part
	z1[n2-1] =   z1[n2-1] ;
		*/
	z1[n2  ] = - z1[n2  ] ;

	/*      compute inverse Fourier transform */
        four(z1,nn,+1,&dt,&df) ;

	/* save the correlated function */
	for(i=0,j=0; i< nn; i++, j+=2){
		x[i] = z1[j];
	}

	/* free temporary arrays */
}

void saciterd_gfilter(float *x, float gwidth , int n, int nn, float dt, float *z1){
/*

       convolve a function with a unit-area Gaussian filter
       and overwrite original array

       x        float * - array to be filtered - this is overwritter
       Gwidth   float   - Filter is exp( - ( pi freq/gwidth)**2)
       n   int          - number of points in time series
       nn  int          - number of points at next largest power of 2
       dt  float        - sampling interval
       z1  float        - work array for complex analysis
*/

	/* internal variables */
	int i,j,n2, jr, ji, kr, ki;
	float df;
	float tmpr, tmpi;
	float fac, freq;
	
	/* define temporrary arrays */


	for(i=0,j=0; i< nn; i++,j+=2){
		if(i < n){
			z1[j  ] = x[i] ;	/* real */
			z1[j+1] = 0.0 ;		/* imaginary */
		} else {
			z1[j  ] = 0.0;	/* real */
			z1[j+1] = 0.0 ;		/* imaginary */
		}
	}
	/* get Fourier transforms */
        four(z1,nn,-1,&dt,&df) ;

        /* Gaussian filter */
	n2 = nn / 2 ;
	for(i=0,j=0; i <= n2 ; i++, j+=2){
		freq =  i  * df ;
		fac = 3.1415927*freq/gwidth ;
		fac = fac * fac ;
		/* get exp(-fac) but ensure no underflow */
		if(fac > 50){
			fac = 0.0 ;
		} else {
			fac = exp(-fac);
		}
		jr = j;
		ji = j + 1;
		z1[jr] *= fac;
		z1[ji] *= fac;
		if(i > 0){
			/* apply symmetry to make the convolution real */
			kr = 2*nn +2 - jr -2;
			ki = 2*nn +2 - jr -1   ;
			z1[kr] =   z1[jr] ;
			z1[ki] = - z1[ji] ;
		}
	}
	/* ensure Nyquist frequency element is real */
	/*  preserve the real part, change sign of imaginary part
	z1[n2-1] =   z1[n2-1] ;
		*/
	z1[n2  ] = - z1[n2  ] ;

	/*      compute inverse Fourier transform */
        four(z1,nn,+1,&dt,&df) ;

	/* reconstitute the real series from the compelx array */
	for(i=0,j=0; i< nn; i++, j+=2){
		x[i] = z1[j];
	}

}

void saciterd_fcorrelate(float *x, float *y, int n, int nn, float dt, float *z1, float *z2)
{

/*
	correlate x and y, replacing the y array with the correlation
        which is norrmalized by the autocorrelation of g


       x   float - array
       y   float - array
       n   int   - number of points in time series
       nn  int   - number of points rounded up to next power of 2
       dt  float - sampling interval
       z1  float - complex work array
       z2  float - complex work array
*/

	/* internal variables */
	int i,j,n2, jr, ji, kr, ki;
	float df;
	float tmpr, tmpi;
	float sum0;



	/* compute the zero lag auto-correlation of y */
	sum0 = 0.0 ;
	for(i=0 ; i < n ; i++)
		sum0 += y[i] * y[i] ;
	sum0 *= dt;

	/* convolve  by frequency domain multiplication */
	for(i=0,j=0; i< nn; i++,j+=2){
		if(i < n){
			z1[j  ] = x[i] ;	/* real */
			z1[j+1] = 0.0 ;		/* imaginary */
			z2[j  ] = y[i] ;	/* real */
			z2[j+1] = 0.0;		/* imaginary */
		} else {
			z1[j  ] = 0.0;	/* real */
			z1[j+1] = 0.0 ;		/* imaginary */
			z2[j  ] = 0.0;	/* real */
			z2[j+1] = 0.0;		/* imaginary */
		}
	}
	/* get Fourier transforms */
        four(z1,nn,-1,&dt,&df) ;
        four(z2,nn,-1,&dt,&df) ;

	/* cross correlation = F  G */
        n2 = nn / 2 ;
	for(j=0,i=0; i <= n2 ; i++, j+=2){
		/* this can be simplified since the z1 and z2 arrays are both real 
			but full complex arithmetic is used here
		*/
		jr = j   ;
		ji = j+1 ;
		/* we need z1 conjg(z2) */
		tmpr = z2[jr]*z1[jr] + z2[ji]*z1[ji];
		tmpi = z2[jr]*z1[ji] - z2[ji]*z1[jr];
                z1[jr] = tmpr;
		z1[ji] = tmpi;
		if(i > 0){
			/* apply symmetry to make the convolution real */
			kr = 2*nn +2 - jr -2;
			ki = 2*nn +2 - jr -1   ;
			z1[kr] =   z1[jr] ;
			z1[ki] = - z1[ji] ;
		}
	}
	/* ensure Nyquist frequency element is real */
	/*  preserve the real part, change sign of imaginary part
	z1[n2-1] =   z1[n2-1] ;
		*/
	z1[n2  ] = - z1[n2  ] ;

	/*      compute inverse Fourier transform */
        four(z1,nn,+1,&dt,&df) ;

	/* update the y array and apply the normalization */
	for(i=0,j=0; i< nn; i++, j+=2){
		y[i] = z1[j] / sum0;
	}

	/* free temporary arrays */
}

void saciterd_phs_shift(float *x,float theshift,int  n,int  nn,float dt, float *z1)
{
/*
       time shift a signal

       X   float - signal to be shifted
       theshift R*4    - time shift in seconds
       n   int   - length of signal
       nn  int   - length of signal to nearest power of 2
       dt  float - sampling interval
       z1  float - complex work array
*/
        float  pi, two_pi, d_omega, df;
	float c, s, tmpr, tmpi, ang;
	int i, j, jr, ji, kr, ki , n2;


	for(i=0,j=0; i< nn; i++,j+=2){
		if(i < n){
			z1[j  ] = x[i] ;	/* real */
			z1[j+1] = 0.0 ;		/* imaginary */
		} else {
			z1[j  ] = 0.0;	/* real */
			z1[j+1] = 0.0 ;		/* imaginary */
		}
	}
	/* get Fourier transforms */
        four(z1,nn,-1,&dt,&df) ;

        pi = 3.1415927 ;
        two_pi = 2 * pi ;
        d_omega = two_pi * df ;

        /* time shift in frequency domain  */
	n2 = nn / 2 ;
	for(i=0,j=0; i <= n2 ; i++, j+=2){
		ang =  i  * d_omega * theshift ;
		c = cos(ang) ;
		s = sin(ang) ;
		jr = j;
		ji = j + 1;
                tmpr =   z1[jr]*c + z1[ji]*s ;
		tmpi = - z1[jr]*s + z1[ji]*c ;
		z1[jr] = tmpr ;
		z1[ji] = tmpi ;
		if(i > 0){
			/* apply symmetry to make the convolution real */
			kr = 2*nn +2 - jr -2;
			ki = 2*nn +2 - jr -1   ;
			z1[kr] =   z1[jr] ;
			z1[ki] = - z1[ji] ;
		}
	}

	/* ensure Nyquist frequency element is real */
	/*  preserve the real part, change sign of imaginary part
	z1[n2-1] =   z1[n2-1] ;
		*/
	z1[n2  ] = - z1[n2  ] ;

	/*      compute inverse Fourier transform */
        four(z1,nn,+1,&dt,&df) ;

	/* reconstitute the real series from the compelx array */
	for(i=0,j=0; i< nn; i++, j+=2){
		x[i] = z1[j];
	}
}

void saciterd_tdomain(char *fnum,char *fden,float gwidth,float theshift,
          int verbose,int maxbumps,int lpositive,float tol,int dotwice,float rayp,int quiet){

/*
	time domain deconvolution

	fnum    char* - file name for numerator
	fden    char* - file name for denominator
	f       float - numerator array
	g       float - denominator array
	p       float - array
	r       float - array
	nn      int   - power of 2 for fft
	npts    int   - number of points in initial time series
	dt      float - sampling interval
	gwidth  float - filter width
	theshift float    - time delay for output of result
	verbose  int    - 1 for a lot of intermediate files
	maxbumps int       - number of iterators
	lpositive int    1 only ouput positive pulses
	tol	float - convergence tolerance
	dotwice	int    - 1 use double length FFT to avoid wrap around
	rayp	float - ray parameter in sec/km
*/
/*      purpose: given N and D find the filter F such that N = F * D
		in the time domain
	method:  Kikuchi and Kanomore iterative deconvolution
	routine:
		1. read in N and D as xn and xd
		2. Gaussian filter each and store as fn and fd
		3. Initialize first bump
		   a.Cross correlate fn and fd and determine the
			maximum and lag.  The model is that
			fn = A fd. Then the estimate of A is
			A = ( fn x fd)    / ( fd x fd )
			              max              0
                        where 'x' is the correlation operator and
			'max' is the maximum value of the quantity
			and the '0' is the zero lag value
		   b.Since the forward model will consist of a
			sequence of impulses, the impulse will have
			a height of A/dt to insure that the proper
			amplitude in the spectral domain
		
		   c.Compute the predicted trace from the convolution of the
			impulse series with the filtered denominator, e.g.,
			p = A exp ( - i omega delay) x fd
                   d.Compute the residual
			r = fn - p
		5. Begin loop
		   Repeat steps 3a-3d starting with the residual file r 
			instead of the fn

		6. Output results
*/

	float *xn, *xd;		/* trace arrays calloc'd by the brsac */
	float *f, *g, *p, *r;
	float *fn, *fd;
	float *z1, *z2;			/* temporary complex arrays */
	float btime_num, btime_den;
	int nn, npts, nout, nwrout;
	int nnum, nden;			/* number of points in num/den */
	float dt, cmpdt, power, sumsq_ip1, sumsq_i, d_error ;
	float depmax, depmin, depmen ;
	float b, beg, delta, theend;
	float fit;
	int i, maxlag, nshifts;
	int maxpts;  	/* a dummy variable for compatibility with FORTRAN
		must be quite large */
	int ndummy;
	int indmax, indmin;


#define MAXG 1000	/* maximum number of bumps */
        float amps[MAXG] ;
        int shifts[MAXG] ;
#define NSAMP 131072   /* this is a FORTRAN ism that is built into brsac */
	

	char resfile[13], filename[13];
	int nerr, ierr;

	maxpts = NSAMP;

	/* parameter check */
        if(maxbumps > MAXG ){
            maxbumps = MAXG ;
	}
        if(maxbumps < 0){
            maxbumps = MAXG ;
        }

	/* READ THE NUMERATOR AND DENOMINATOR FILES */
	/* open the sac file for the numerator */
	maxpts = NSAMP;
        brsac(maxpts,fnum,&xn,&nerr) ;  /* brsac performs a calloc for f */
        if(nerr < 0){
            fprintf(stderr,"Problem reading the numerator   file: nerr=%d\n",nerr);
            saciterd_usage() ;
        }
	/* get header values */
        getfhv("DELTA   ",&cmpdt,&nerr) ;
        getfhv("B       ",&btime_num,&nerr) ;
        getnhv("NPTS    ",&npts, &nerr) ;
	nnum = npts; 
        dt = cmpdt ;
        nn = npts ;
	/* work on the denominator */
	maxpts = NSAMP;
	brsac(maxpts,fden,&xd,&nerr) ;
	if(nerr < 0){
	    fprintf(stderr,"Problem reading the denominator file: nerr=%d\n",nerr);
	    saciterd_usage() ;
	}
	/* get header values for denominator */
        getnhv("NPTS    ",&nden, &nerr) ;
        getfhv("B       ",&btime_den,&nerr) ;
	
	/* END OF WAVEFORM INPUT */

	/* INITIALIZE PROCESSING ARRAYS */
	/* find power of 2 that >= number of observed points */
        nn = npow2(npts) ;
	/*
		if use double length FFT , nn -> 2*nn, but comparison is
		for residuals and output is based on nn before doubling
	*/
        if(dotwice == YES){
            nwrout = nn ;
            nn = 2*nn ;
            nout = nn ;
	} else {
            nout = nn ;
            nwrout = nn ;
        }
	/* allocate space for real arrays f, g,  r and p */
        f = (float *)calloc(nn,sizeof(float));
        g = (float *)calloc(nn,sizeof(float));
        p = (float *)calloc(nn,sizeof(float));
        r = (float *)calloc(nn,sizeof(float));
        z1 = (float *)calloc(2*nn,sizeof(float));
        z2 = (float *)calloc(2*nn,sizeof(float));
	fn = (float *)calloc(nn,sizeof(float));
        fd = (float *)calloc(nn,sizeof(float));


	/* CREATE THE FILTERED WAVEFORMS */

	/* copy the xn array into f work array */
	saciterd_zero(f,nn);
	for(i=0 ; i < nnum;i++)
		f[i] = xn[i];
	
        scmxmn(f,npts,&depmax,&depmin,&depmen,&indmax,&indmin) ;
        printf("Max Min amplitudes of Numerator  : %g %g\n",depmin,depmax );
	/* Filter and save the numerator */
        saciterd_gfilter(f,gwidth,npts,nn,dt,z1) ;

        setnhv("NPTS    ",nout, &nerr) ;
        scmxmn(f,nout,&depmax,&depmin,&depmen,&indmax,&indmin) ;
        setfhv("TIMMAX  ",btime_num+indmax*dt, &nerr) ;
        setfhv("TIMMIN  ",btime_num+indmin*dt, &nerr) ;
        setfhv("DEPMAX  ",depmax, &nerr) ;
        setfhv("DEPMIN  ",depmin, &nerr) ;
        setfhv("DEPMEN  ",depmen, &nerr) ;
        bwsac(nout,"numerator",f) ;
        bwsac(nout,"observed" ,f) ;
	/* save the filtered numerator */
	for(i=0;i<nn;i++)fn[i] = f[i];

	/* copy the xn array into g work array */
	saciterd_zero(g,nn);
	/* copy the xd array into g */
	for(i=0 ; i < nden;i++)
		g[i] = xd[i];

        scmxmn(g,npts,&depmax,&depmin,&depmen,&indmax,&indmin) ;
        printf("Max Min amplitudes of Denominator: %g %g\n",depmin,depmax );
	/* filter and save the denominator */
        saciterd_gfilter(g,gwidth,npts,nn,dt,z1) ;
        setnhv("NPTS    ",nout, &nerr) ;
        scmxmn(g,nout,&depmax,&depmin,&depmen,&indmax,&indmin) ;
        setfhv("TIMMAX  ",btime_den+indmax*dt, &nerr) ;
        setfhv("TIMMIN  ",btime_den+indmin*dt, &nerr) ;
        setfhv("DEPMAX  ",depmax, &nerr) ;
        setfhv("DEPMIN  ",depmin, &nerr) ;
        setfhv("DEPMEN  ",depmen, &nerr) ;
        bwsac(nout,"denominator",g) ;
	/* save the filtered denominator */
	for(i=0;i<nn;i++)fd[i] = g[i];
	/* compute the power in the numerator for error scaling */
        power = 0.0 ;
	for(i=0 ; i < nn ; i++){
            power +=  fn[i]*fn[i] ;
	}
	
	/* free the original trace arrays */
	free(xn);
	free(xd);

	/* ESTIMATE AMPLITUDE AND DELAY OF THE FIRST PULSE */
	/* correlate the signals */
	saciterd_fcorrelate(f,g,nn,nn,dt, z1, z2) ;
	/* fine the peak in the correlation */
        maxlag = nout/2 ;
        printf("\nThe maximum spike delay is  %f\n",maxlag * dt );
	if(lpositive == YES) {
		saciterd_getmax(g,maxlag,&amps[0],&shifts[0]) ;
	} else {
		saciterd_getabsmax(g,maxlag,&amps[0],&shifts[0]) ;
	}
        printf(" maxlag,amps,shifts: %d %f %d\n",maxlag,amps[0],shifts[0] );
        amps[0] /=  dt ;

        nshifts = 0 ;
	/* read in the signals again */
	saciterd_zero(f,nn);
	for(i=0 ; i < nnum;i++)
		f[i] = fn[i];
	saciterd_zero(g,nn);
	for(i=0 ; i < nden;i++)
		g[i] = fd[i];

        npts = nn ;

	/* compute the predicted deconvolution result
		Note that this will be a sequence if
		impulses */

        saciterd_zero(p,nn) ;
        saciterd_build_decon(amps,shifts,nshifts+1,p,npts,nn,gwidth,dt,z1) ;
        setnhv("NPTS    ",nn, &nerr) ;
        if(verbose) {
		saciterd_phs_shift(p,theshift,nn,nn,dt,z1) ;     
		wsac1("d001",p,npts,-theshift,dt,&nerr) ;
		saciterd_phs_shift(p,-theshift,nn,nn,dt,z1) ;     
        }

	/* convolve the prediction with the denominator signal */
        saciterd_convolve(p,g,npts,nn,dt,z1,z2) ;

        if(verbose) {
           wsac1("p001",p,npts,beg,dt,&nerr) ;
        }


        if(verbose){
	  /* the initial residual is just the numerator */
	  i = 0;
          sprintf(resfile,"r%3.3d",i);
          wsac1(resfile,f,nout,beg,dt,&nerr) ;
        }
	/* compute the residual (initial error is 1.0 ) */
        saciterd_getres(f,p,nout,r,&sumsq_ip1) ;

        sumsq_i = 1.0 ;
        sumsq_ip1 /=  power ;
        d_error = 100*(sumsq_i - sumsq_ip1)  ;

	
	i = 1;
	sprintf(resfile,"r%3.3d",i);
        if(verbose){
          wsac1(resfile,r,nout,beg,dt,&nerr) ;
        }
    
        if( quiet != YES ){
		printf("\n File         Spike amplitude   Spike delay   Misfit   Improvement\n") ;
		printf("%10s  %16.9e  %10.3f   %7.2f%%   %9.3f%%\n",
			resfile, dt*amps[0],shifts[0]*dt,100*sumsq_ip1,d_error);
	}

        while(d_error > tol && nshifts < maxbumps -1) {
		nshifts++  ;

		sumsq_i = sumsq_ip1 ;

		for(i=0 ; i < nn;i++)
			g[i] = fd[i];
		saciterd_fcorrelate(r,g,nn,nn,dt, z1, z2) ;
		if(lpositive == YES){
			saciterd_getmax(g,maxlag,&amps[nshifts],&shifts[nshifts]) ;
		} else {
			saciterd_getabsmax(g,maxlag,&amps[nshifts],&shifts[nshifts]) ;
		}
		amps[nshifts] /=  dt ;

		saciterd_zero(p,nn) ;
		saciterd_build_decon(amps,shifts,nshifts+1,p,npts,nn,gwidth,dt,z1) ;
		if(verbose){
			sprintf(filename,"d%3.3d",nshifts);
				saciterd_phs_shift(p,theshift,nn,nn,dt,z1)   ;
			wsac1(filename,p,nout,-theshift,dt,&nerr) ;
			saciterd_phs_shift(p,-theshift,nn,nn,dt,z1) ;
		}

		for(i=0 ; i < nn;i++)
			g[i] = fd[i];
		saciterd_convolve(p,g,npts,nn,dt,z1,z2) ;
		if(verbose){
			sprintf(filename,"p%3.3d",nshifts);
			wsac1(filename,p,nout,beg,dt,&nerr) ;
		}
               
		for(i=0 ; i < nnum;i++)
			f[i] = fn[i];
		saciterd_getres(f,p,nout,r,&sumsq_ip1) ;
          
		sumsq_ip1 /=  power ;
	
		sprintf(resfile,"r%3.3d",nshifts+1);
		if(verbose){
			wsac1(resfile,r,nout,beg,dt,&nerr);
		}
		d_error = 100*(sumsq_i - sumsq_ip1) ;
        
                if( quiet != YES ) {
			printf("%10s  %16.9e  %10.3f   %7.2f%%   %9.3f%%\n",
				resfile, dt*amps[nshifts],shifts[nshifts]*dt,100*sumsq_ip1,d_error);
		}
	}

	/* ITERATIVE DECONVOLUTION IS FINISHED OUTPUT THE FINAL RESULTS */
	printf("\nLast Error Change = %9.4f%%\n",d_error) ;
	/* if the last change made no difference, drop it */
        fit = 100 - 100*sumsq_ip1 ;

        if(d_error <= tol){
           nshifts = nshifts - 1 ;
           fit = 100 - 100*sumsq_i ;
           printf("Hit the min improvement tolerance - halting.\n");
        }

        if(nshifts > maxbumps){
           printf("Hit the max number of bumps - halting\n");
        }
        printf("Number of bumps in final result: %d\n", nshifts );
	printf("The final deconvolution reproduces %6.1f%% for the signal\n",fit);
	/* compute the final prediction */
        saciterd_zero(p,nn) ;
        saciterd_build_decon(amps,shifts,nshifts+1,p,npts,nn,gwidth,dt,z1) ;
        saciterd_zero(g,nn) ;
       	brsach(fden,&nerr) ;  /* brsac performs a calloc for f */
	getnhv("NPTS    ",&ndummy,&ierr) ;
	getfhv("DELTA   ",&dt  ,&ierr) ;
	getfhv("B       ",&b  ,&ierr) ;
		saciterd_zero(g,nn);
		for(i=0 ; i < ndummy;i++)
			g[i] = fd[i];
        saciterd_convolve(p,g,npts,nn,dt,z1,z2) ;
        wsac1("predicted",p,nout,beg,dt,&nerr) ;
        saciterd_zero(g,nn) ;

	/* write out the answer */
        saciterd_zero(p,nn) ;
        saciterd_build_decon(amps,shifts,nshifts+1,p,npts,nn,gwidth,dt,z1) ;
        saciterd_gfilter(p,gwidth,npts,nn,dt,z1) ;
        saciterd_phs_shift(p,theshift,nn,nn,dt,z1) ;
	/*
	we use the original numerator file to get header values,
	related to event location and station, and then modify the fields
	accordingly
	*/
        newhdr() ;
       	brsach(fnum,&nerr) ;  /* brsac performs a calloc for f */
        scmxmn(p,ndummy,&depmax,&depmin,&depmen,&indmax,&indmin) ;
	getnhv("NPTS    ",&ndummy,&ierr) ;
	getfhv("DELTA   ",&dt  ,&ierr) ;
	getfhv("B       ",&b  ,&ierr) ;
        setnhv("NPTS    ",nwrout, &nerr) ;
        setfhv("B       ",-theshift, &nerr) ;
        setfhv("TIMMAX  ",-theshift+indmax*dt, &nerr) ;
        setfhv("TIMMIN  ",-theshift+indmin*dt, &nerr) ;
        setfhv("DEPMAX  ",depmax, &nerr) ;
        setfhv("DEPMIN  ",depmin, &nerr) ;
        setfhv("DEPMEN  ",depmen, &nerr) ;
        theend = -theshift + (nn-1)*dt ;
        setfhv("E       ",theend, &nerr) ;
        setnhv("NZSEC   ",-12345, &nerr) ;
        setfhv("USER0   ",gwidth, &nerr) ;
        setkhv("KUSER0  ","Rftn    ", &nerr) ;
        setkhv("KUSER1  ","IT_DECON", &nerr) ;
        setkhv("KEVNM   ","Rftn    ", &nerr) ;
        setkhv("KEVNMC  ","        ", &nerr) ;
        setfhv("USER5   ",fit,&nerr) ;
        if(rayp > 0.0){
            setfhv("USER4   ", rayp, &nerr) ;
        }
	bwsac(ndummy,"decon.out",p) ;

	/* write out the gaussian filter */
	if(verbose == YES){
          newhdr() ;
          saciterd_zero(p,nn) ;
          p[1] = 1 / dt ;
          saciterd_phs_shift(p,theshift,nn,nn,dt,z1) ;
          saciterd_gfilter(p,gwidth,npts,nn,dt,z1) ;
          wsac1("thefilter",p,nn,beg,dt,&nerr) ;
	}
	/* clean up temporary arrays */
	free(f);
	free(g);
	free(p);
	free(r);
	free(z1);
	free(z2);
}

/* 
	the following code is from PROGRAMS.330/VOLVIII/gsac.src/gsac_sub.c
*/

int npow2(int n)
{
	/* return integer that is
		>= n, and is a power of 2
	*/
        int k, nn;
        k = 1; nn = 2;
        while(nn < n){
                nn *=2;
        }
        return(nn);

}

void four(float data[], int nn, int isign, float *dt, float *df){
/*
	subroutine four(data,nn,isign,dt,df)
c-----
c     THE input is a complex array
c     which has numbers stored in memory as
c     R1, I1, R2, I2, ..., Rnn, Inn
c     where nn must be a power of 2 R and I are the real and imaginary
c     parts of the complex number
c
c     For isign -1 this is a complex time series
c     For isign +1 this is a complex frequency series with
c        index 1 (in fortran corresponding to f=0
c              2                              f=df
c            nn/2 + 1                         f = 1/2dt = Nyquist
c            nn - 1                           f = -2df
c             nn                              f = -df
c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(data(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  data is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     data are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array data, replacing the input data.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c-----
      real*4 data(*)
*/
	int n;
	int i, j, m, mmax, iiii, istep;
	float tempr, tempi;
	float wr, wi;
	float wstpr, wstpi;
	float sinth, theta;
	float dtt, dff;

	dtt = (*dt);
	dff = (*df);

	n = 2 * nn;
	if((dtt) == 0.0) dtt = 1./(nn*(dff)) ;
	if((dff) == 0.0) dff = 1./(nn*(dtt)) ;
	if((dtt) != (nn*(dff))) dff = 1./(nn*(dtt)) ;
	*dt = dtt;
	*df = dff;
	j = 1;
	for (i=1;i<= n ; i+=2){
		if(i < j){
			tempr = data[j-1];
			tempi = data[j  ];
			data[j-1] = data[i-1];
			data[j  ]=data[i  ];
			data[i-1] = tempr;
			data[i  ] = tempi;
		}
		m = n/2;
statement3:		if(j <= m) goto statement4 ;
			j = j-m;
			m = m/2;
			if(m >= 2)goto statement3 ;
statement4:		
		j=j+m;
    	}
	mmax = 2 ;
	while(mmax < n ){
		istep= 2 *mmax;
		theta = 6.283185307/(float)(isign*mmax);
		sinth=sin(theta/2.);
		wstpr=-2.*sinth*sinth;
		wstpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1; m <= mmax ; m +=2){
				for(i = m ; i <= n ; i+=istep){
					j=i+mmax;
					tempr=wr*data[j-1]-wi*data[j  ];
					tempi=wr*data[j  ]+wi*data[j-1];
					data[j-1]=data[i-1]-tempr;
					data[j  ]=data[i  ]-tempi;
					data[i-1]=data[i-1]+tempr;
					data[i  ]=data[i  ]+tempi;
				}
				tempr = wr;
				wr = wr*wstpr-wi*wstpi + wr;
				wi = wi*wstpr+tempr*wstpi + wi;
		}
		mmax = istep;
	}
 /*
  * 	get the correct dimensions for a Fourier Transform 
  * 	from the Discrete Fourier Transform
*/ 
	if(isign > 0){
		/*
		frequency to time domain
		*/
		for (iiii= 0 ; iiii < n ; iiii++){
			data[iiii] = data[iiii] * (*df);
		}
	} else {
		/*
		time to frequency domain
		*/
		for (iiii= 0 ; iiii < n ; iiii++){
			data[iiii] = data[iiii] * (*dt);
		}
 }
}
