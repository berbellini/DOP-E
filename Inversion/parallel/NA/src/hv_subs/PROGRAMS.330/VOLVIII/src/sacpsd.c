/* sacpsd
	Changes:
	28 JAN 2010 - created
		Certain routines are taken from gsac
	28 FEB 2010 - fixed error in demean_trend() by
		replacing an instance of npts by n
        20 SEP 2010  - error in input at lines 999, 1005
			while(fscanf(pfileptr,"%s %d %s\n",&pfile)!=EOF)
                changed to
			while(fscanf(pfileptr,"%s\n",&pfile)!=EOF)
                cleanup up code after using -Wall compile flag

References:

Ifeachor, E. C., and B. W. Jervis (1993). Digital Signal
	Processing: A Practical Approach, Addison Wesley,
	Wokingham, England, 760pp

Windowing corrections:
	given s(n), the corrected time series that
accounts for the decrease in amplitude due to the windowing function w(n)
and the mean level introducted by windowing is S(n) where

	S(n) = k  w(n)[s(n) - k  ]       (10.16)
                2              1

where k   is the mean windowed amplitude
       1

     k  = SUM [ w(n) s(n) ] / SUM w(n)   (10.17)
      1

and   2            2
     k  = N / SUM w (n)                  (10.20)
      2

and N is the number of samples. 
To save computations (10.16) can be written as

        S(n) = k  w(n) s(n)   - k  k  w(n)
                2                2  1
The windowing functions return w(n) s(n), k  and k
                                           1      2
The windowing routine will return the windowed time series w(n) s(n), 
the windowing function w(n) and the k2 and k1 constants.  
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sacsubc.h"
#include "grphsubc.h"
#include "calplot.h"
#include "noisemodel.h"
#include "csstime.h"
#include <ctype.h>

/* define function prototypes */
void gcmdln(int argc, char **argv);
void usage(void);
void four(float data[], int nn, int isign, float *dt, float *df);
int npow2(int n);
int sacdatchk(char *filename);
/* windowing routines */
void hann(float *w, float *y, int n);
void sin10(float *w, float *y, int n);
/* demean and linear trend routine */
void dmean_trend(float *x, int n, int *nerr);
/* smoothing routine */
void smooth_5(double *x, double *y, int np);
/* main routine for performing the PSD */
void dopsd (int *n21, float *df);
void deblank(char *s, int n);

/* routines from gsac */
void arr_locate(float *x, int n, float xval, int *jind);
void gsac_evaamp(float freq,float *famp);
void doplot(int n21, float df);
int gsac_set_eval(char *ampname);
int  sfgetline(FILE *fp, char s[], int lim);
void evarray(float xval, double *yval, float *x, double *y, int n);

/* graphics routines */
void doplot_axes(float x0,float y0,float xlen,float ylen,float xmin,
	float xmax,float ymin,float ymax);
void doplot_obspsd(float x0,float y0,float xlen,float ylen,float xmin,
	float xmax,float ymin,float ymax,  int n, float df, int do_smooth);
void do_plotnm(float x0,float y0,float xlen,float ylen,float xmin,
	float xmax,float ymin,float ymax,struct _PXY *nm);
void do_plotfile(float x0,float y0,float xlen,float ylen,float xmin,
	float xmax,float ymin,float ymax, int kolor, float width, 
	char *title, char *fname);

/* define useful macros and defines */
#define MIN(a,b) ( (b) > (a) ? (a):(b) )
#define MAX(a,b) ( (b) < (a) ? (a):(b) )
#define ABS(a  ) ( (a) >  0  ? (a):-(a))
#define SIGN(a ) ( (a) >  0  ? (1):(-1))
#define YES 1
#define NO  0

/* CALPLOT COLORS */
#define BLACK    1 /* these two are reversed if use -I on plotnps, plotxvig */
#define WHITE    0
#define RED      2
#define GREEN    3
#define BLUE     4
#define ORANGE   5
#define BLUGREEN 6
#define YELLOW   7
/* 1000 = red, 1100 = blue for progression*/

/* global parameters */
FILE  *fampfile;
int have_sacfile, have_ampfile, have_otherpsd, have_anotherpsd;
char *sacfile,  *ampfile;
char *listfile;
char *alistfile;
int npts;
float dt;
char kstnm[9], kcmpnm[9], khole[9], knetwk[9];
int nzyear, nzjday, nzmon, nzday, nzhour, nzmin, nzsec, nzmsec;
float window;
struct date_time t;
char smoothedpsdname[32];
int nft, nseg;		/* length of FFT use number of segments averaged */


float *x, *y, *w; /* storage for trace, temporary array for trace segment 
	and window  function */
float *z; /* storage for FFT */

/* EVAL  */
#define NEVAL 50000
int evnamp;            /* number of amplitude entries */
float evfamp[NEVAL];   /* each files has its own frequency column */
float evaamp[NEVAL];   /* each files has its own frequency column */

#define NCHAR 200
static char input_lineptr[NCHAR] ;

/* 4096 8192 16384 32768 */
#define NARRAY 16384
double psdave[NARRAY], period[NARRAY];	/* storage arrays for stacked PSD */
double psdsmooth[NARRAY];
float frqarr[NARRAY];
float pmin, pmax;			/* min and max periods from dt */

/* this is used to plot user input spectra */
FILE *ofileptr;
char pfile[NCHAR];
char ptitle[NCHAR];
FILE *pfileptr;
int do_title;
float position_y;

/* set a flag to select the type of windowing */
#define HANNING 0
#define SIN10   1
int window_type = SIN10;

/* windowing constants */
double k22, k1 ;  /* k2 squared and k1 */
int do_smooth_5 = NO;



main(int argc, char **argv)
{

	float df;
	int n21;
	/* parse command line arguments */
	gcmdln(argc, argv);


	n21 = 0;
	if(have_sacfile && have_ampfile)
		dopsd(&n21, &df);
	/* now do the graphics */
	doplot(n21,df);
	return 0;
}
 

void dopsd(int *nn21, float *ddf)
{
	int nerr;
	int i,j,offset,n21;
	float dt, df, freq, fnyq;
	double psd;
	double facpsd;
	double dr, di;
	float amp;
	FILE *fpsdout;
	int imin, imax;
	double period_mult, period_current;
	FILE* smoothedpsdfd;

	float depmax, depmin, depmen;

        printf("Processing %s %s\n",sacfile,ampfile);

	/* get the acceleration sensitivity */
	if(gsac_set_eval(ampfile) == NO)
		return;

	/* get the trace and get trace parameters */

	if(sacdatchk(sacfile) == NO)
		return;
	/* allocate memory 
		note that the x array is allocated in the brsac routine
		of sacsubc.c - we can free it later though */
	y = (float *)calloc(NARRAY, sizeof(float));
	w = (float *)calloc(NARRAY, sizeof(float));
	brsac(npts,sacfile,&x,&nerr);
        getfhv("DELTA",&dt,&nerr);
	window = (npts -1 ) * dt;
	fnyq = 0.5/dt;



	/* initialize the PSD average */
	for(i=0 ; i < NARRAY; i++)
		psdave[i] = 0.0;

	if(npts > NARRAY){
		nft = NARRAY;
	} else {
		/* get power of 2 >= npts */
		nft = npow2(npts);
		nft /= 2;
	}
	n21 = nft/2;
	z = (float *)calloc(2.*nft, sizeof(float));


	/* march along the time series to get
		overlying segments - shifting by NARRAY/2 each time */
	for(offset = 0, nseg=0; offset< npts - nft;offset += nft/2, nseg++){
		for(i=0;i<nft;i++){
			y[i] = x[i+offset];
		}
		/* demean and taper */
        	dmean_trend(y,nft,&nerr);

		if(nerr >= 0){
			if(window_type == HANNING)
				hann(w,y,nft);
			else if(window_type == SIN10)
	        		sin10(w,y,nft);
			
			/* initialize FFT array accounting for
				the mean correction (10.16) */
			for(i=0,j=0;i<nft;i++){
				z[j++] = y[i] - k1 * w[i];
				z[j++] = 0.0;
			}
			/* must define df to some initial value */
			df = 0.0;
			four(z,nft,-1,&dt,&df);
			/* now output the PSD */
			facpsd = 2.0 * k22 /(nft*dt);
			/* safety define the first element */
			frqarr[0] = 0.0;
			psdave[0] = 0.0;
			for(i=1;i<n21;i++){
				freq = i*df;
				gsac_evaamp(freq,&amp);
				j = 2*i;
				dr = z[j];
				di = z[j+1];
				psd = facpsd * (dr*dr + di*di )
					/(amp*amp);
				psdave[i] += psd;
				frqarr[i] = freq;
		
				
			}
		}
	}
	/* determine the avarage PSD */
	for(i=1;i<n21;i++){
			psdave[i] /= nseg;
	}

		/* apply a five point smoothing */	
	if(do_smooth_5 == YES)
		smooth_5(&psdave[1],&psdsmooth[1],n21-1);

	/* determine the average PSD  and output in the file
		sacpsd.out */
	/* we can now output the selected periods */
	pmax = 0.25*nft*dt;
	pmin = 1.90/fnyq;
	for(i=1,imin=n21,imax=0;i<n21;i++){
		freq = i*df;
		if(freq <= 0.53*fnyq && freq >= 4.0/(nft*dt) ){
			if(i < imin) imin = i;
			if(i > imax) imax = i;
			if(psdave[i] > 0.0){
				psdave[i] = 10.0*log10(psdave[i]);
				period[i] = 1./freq;
			}
		} else { 
				period[i] = -1.0 ; 
				psdave[i] = 0.0 ;
		}
	}
	fpsdout = fopen("sacpsd.out", "w+");
	for(i=imin;i < imax; i++){
		fprintf(fpsdout,"%9.4f %6.0f\n",period[i],psdave[i]);
	}
	fclose(fpsdout);

	/* output periods at 1/8 octave intervals */
	smoothedpsdfd = fopen(smoothedpsdname,"w+");
	period_mult = pow(2.0, 0.125);
	period_current = 1.0 * 1024 ;
	/* start with some outrageously large period */
	for(;period_current > pmin/period_mult;period_current /= period_mult){
		if(period_current >= pmin && period_current <= pmax){
			freq = 1.0/period_current;
			 evarray(freq, &psd, &frqarr[0],  &psdsmooth[0], n21);
			if(psd > 0.0){ 
			psd = 10.0*log10(psd);
fprintf(smoothedpsdfd,"%12.6f %10.3f\n",period_current,psd);
			}
		}
	}
	fclose(smoothedpsdfd);
	

	/* clean up on termination */
	free(w);
	free(x);
	free(y);
	free(z);
	*nn21 = n21;
	*ddf = df;
}
		

int sacdatchk(char *filename)
{
	int nerr;
	brsach(filename, &nerr);
        if(nerr < 0)
		return NO;
	getnhv("NPTS    ",&npts,&nerr);
	getnhv("NZYEAR    ",&nzyear,&nerr);
	getnhv("NZJDAY    ",&nzjday,&nerr);
	getnhv("NZHOUR    ",&nzhour,&nerr);
	getnhv("NZMIN    ",&nzmin,&nerr);
	getnhv("NZSEC    ",&nzsec,&nerr);
	getnhv("NZMSEC    ",&nzmsec,&nerr);
	getfhv("DELTA   ",&dt,&nerr);
 	getkhv("KSTNM",kstnm,&nerr);
 	getkhv("KCMPNM",kcmpnm,&nerr);
 	getkhv("KHOLE",khole,&nerr);
 	getkhv("KNETWK",knetwk,&nerr);
	if(strncmp(khole,"-12345",6) == 0  || strncmp(khole,"  ",2) == 0 )
		strcpy(khole,"__");
	if(strncmp(knetwk,"-12345",6) == 0 || strncmp(knetwk,"  ",2) == 0 )
		strcpy(knetwk,"__");
	deblank(kstnm,5);
	deblank(knetwk,2);
	deblank(khole,2);
	deblank(kcmpnm,3);
	t.year = nzyear;
	t.doy = (long)nzjday;
	t.hour = nzhour;
	t.minute = nzmin;
	t.second = (float)nzsec + 0.001*(float)nzmsec;
	/* convert from human to epoch */
	month_day(&t);
	mdtodate(&t);

	/* define the name of the smoothed PSD output file of the form
		NNSSSSSCCCLL.YYYY.DDD.HH.MM.psd where
	NN    is network code
	SSSSS is station code
	CCC   is component code
	LL    is location code
	YYYY  is year
	DDD   is day of year
	HH    is hour 
	MM  is minute
	*/
	sprintf(smoothedpsdname,"%2.2s%5.5s%3.3s%2.2s.%4.4d.%3.3d.%2.2d.%2.2d.psd",
		knetwk,kstnm,kcmpnm,khole,t.year,(int)t.doy,t.hour,t.minute);
printf("%s\n",smoothedpsdname);

	
	if(nerr < 0)
		return NO;
	else
		return YES;
}

int npow2(int n)
{
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

void gcmdln(int argc, char **argv)
{
/* parse command line arguments, e.g.,
	sacpsd -f sacfile -r response
   where
	sacfile is a sac binary trace
	response is the acceleration sensitivity given as a table 
	of frequency-response pairs created by invoking evalresp with 
	the -u "acc" flag to give the amplitude response in
	units of COUNTS/M/S**2  as a function of frequency (Hz)
*/
	char *cp;


	have_sacfile = NO;
	have_ampfile = NO;
	have_otherpsd = NO;
	have_anotherpsd = NO;
	do_smooth_5 = YES;
	do_title = YES;
	listfile = (char *)NULL;
	alistfile = (char *)NULL;
	

	while(argc-- > 1 ) {
		if(*argv[1] == '-'){
			cp = argv[1];
			cp++;
			if(strcmp("f",cp) == 0 || strcmp("F",cp)==0){
				argv++;
				argc--;
				sacfile = argv[1];
				have_sacfile++;
			} else if(strcmp("r",cp) == 0 || strcmp("R",cp)==0){
				argv++;
				argc--;
				ampfile = argv[1];
				have_ampfile++;
			} else if(strcmp("L",cp) == 0 || strcmp("l",cp)==0){
				argv++;
				argc--;
				listfile = argv[1];
				have_otherpsd++;
			} else if(strcmp("A",cp) == 0 || strcmp("l",cp)==0){
				argv++;
				argc--;
				alistfile = argv[1];
				have_anotherpsd++;
			} else if(strcmp("H",cp)==0){
				window_type = HANNING;
			} else if(strcmp("S",cp)==0){
				window_type = SIN10;
			} else if(strcmp("5",cp)==0){
				do_smooth_5 = YES;
			} else if(strcmp("NT",cp)==0){
				do_title = NO;
			} else if(strcmp("h",cp) == 0 || strcmp("?",cp)==0){
				usage();
			}
			argv++;
		}
	}
	/* put in some syntax checks */
	/* if no other plots, a sacfile requires a response amplitude file */
	/* a sacfile requires a response amplitude file */
		if(have_sacfile == 0 && have_ampfile > 0) usage();
	/* a response amplitude file requires a sacfile */
		if(have_ampfile == 0 && have_sacfile > 0) usage();
}

void usage(void)
{
	fprintf(stderr,
	"Usage: sacpsd -f sacfile -r response -H -S -L listfile -A afilelist -NT\n");
	fprintf(stderr," -f sacfile   (default none) sac binary trace\n");
	fprintf(stderr," -r response  (default none) file giving table of \n");
	fprintf(stderr,
	"           frequency-response pairs created by invoking\n");
	fprintf(stderr,
	"           evalresp with the -u \"acc\" flag to give the \n");
	fprintf(stderr,
	"           amplitude response in units of COUNTS/M/S**2  vs f(Hz)\n");
	fprintf(stderr," -H           (default false ) use Hanning window\n");
	fprintf(stderr," -S           (default true ) use 10%% sin window\n");
	fprintf(stderr,
	" -L listfile plots according to list which has entries  \n");
	fprintf(stderr,"     file1 kolor1 title1\n");
	fprintf(stderr,"     file2 kolor2 title2\n");
	fprintf(stderr,"  where file1 is listing of period-psd(db) pairs\n");
	fprintf(stderr,
	"  kolor is a CALPLOT color, e.g., 1000 to 1100 for red ->blue\n");
	fprintf(stderr,"  title is a string without spaces\n");
	fprintf(stderr,
	" -A alistfile plots according to list which has entries  \n");
	fprintf(stderr,"     file1 \n");
	fprintf(stderr,"     file2 \n");
	fprintf(stderr," The color is automatically generated\n");
	fprintf(stderr," -NT    (default  false ) do not annotate with file name\n");
	fprintf(stderr,
	" -5  (number 5 not ess)  (default true ) apply 5 point PSD smoothing\n");
	fprintf(stderr," -h           (default false) online help\n");
	exit(EXIT_SUCCESS);
}


void hann(float *w, float *x, int n)
{
	int m, i;
	double weight;
	double sum_wx, sum_ww;
	double   fac;
	fac = 6.2831853/(float)n;

	sum_wx = 0.0;
	sum_ww = 0.0;

        for(m=0;m<n;m++){
		i = (-n/2) + m;
		weight = 0.5*(1.0+cos(i*fac));
		w[m] = weight;
		x[m] *= weight;
		sum_wx += x[m];
		sum_ww += weight*weight;
	}
	/* (10.17) and (10.20) */
	k1 = sum_wx/(double)n;
	k22 = (double) n/sum_ww;
}

#define TAPER_PCT 10
void sin10(float *w, float *x, int n)
{
	int j, ffl;
	double weight;
	double sum_wx, sum_ww;
	double   fac;

	sum_wx = 0.0;
	sum_ww = 0.0;

	ffl = (int) ((float) (n*(TAPER_PCT/100.)));
	/* VERY IMPORTANT THE CONSTANT HERE IS 5 PI which is 
		only used for 10 PERCENT TAPER so that 5 PI/n * ffl = PI/2
	*/
	fac = TAPER_PCT*0.5*3.1415927/n ;
	/* fac = 15.70796/n; */
	for(j=0;j<n;j++){
		if(j < ffl)
			weight = sin(j*fac);
		else if( j > n-ffl-1)
			weight = sin((n-j-1)*fac);
		else
			weight = 1.0;
		w[j] = weight;
		x[j] *= weight; 
		sum_wx += (double)x[j];
		sum_ww += weight*weight;
	}
	/* (10.17) and (10.20) */
	k1 = sum_wx/(double)n;
	k22 = (double) n/sum_ww;
}

void dmean_trend(float *y, int n, int *nerr){
/*
	remove the mean and trend from the evently spaced
	time series

	Given regression model  y = A + B x, the L2-norm estimates
	of A is a, and B is b where

        | a |             1            | SUM x^2   -SUM x | | SUM y  |
        |   | = ---------------------- |                  | |        |
        | b |   N SUM x^2 - (SUM x )^2 | -SUM x       N   | | SUM xy |

	which can be simplified if x = [ 1,2,3,...,N]^T  since

        SUM x   = N(N+1)/2
        SUM x^2 = N(N+1)(2N+1)/6

        so

        N SUM x^2 - (SUM x)^2 = N^2 (N+1)(N-1)/12

        However since N can be very very large such an analytic 
	expression have problems. Instead we normalize the x to
	vary between 0 and 1 , which is also donw in gsac

       Note this will usually remove the DC offset of the
       trace, but for a one sided pulse this may seriously
       affect the trace and its spectrum
       y   - time series array
       n   - number of samples	
       nerr - <0 failure due to a glitch or a non-valid point 
*/
        int i, k, ntrc;
        double sumn, sumx, sumxx, sumy, sumxy;
        float depmax, depmin, depmen;
        double yval;
        double det, a, b;
	double x;

	*nerr = 0;
	/* remove the linear trend */
        if(n > 1){
		sumn  = 0.0;
		sumx  = 0.0;
		sumy  = 0.0;
		sumxx = 0.0;
		sumxy = 0.0;
		for(i=0; i < n ; i++){
		x = (double)i/(double)n;
		yval = y[i];
		sumy  += yval;
		sumxy += x*yval;
		sumn  += 1.0;
		sumx  += x;
		sumxx += x*x;
		}
		det = sumn*sumxx - sumx *sumx;
		if(det == 0.0) {
			a = 0.0;
			b = 0.0;
		} else {
			a = ( sumxx * sumy -  sumx*sumxy)/det;
			b = (-sumx  * sumy +  sumn*sumxy)/det;
		}
		for(i=0; i < n ; i++){
			x = (double)i/(double)n;
			y[i] -= (a + b*x);
		}
	}
}


void arr_locate(float *x, int n, float xval, int *jind)
{
	int increase;
	int jlow, jup, jmid;
/*
 	Written by RB Herrmann Saint Louis University 2003
 	Note: the array must be ordered, however it can be ordered
 	increasing or decreasing
 
 	The technique essentially uses an interval halving technique 
 	which means that the number of comparisons is on the order of log_2(n)
 	instead of n for a linear search
 
 	Given an ordered array x[], find jind such that
 	xval belongs to ( x[jind] x[jind+1] )
 	legitimate values are 0 <= jind <= n-2
 	jind = -1 and jind = n-1 indicate failure
*/
/* do the arrays increase or decrease ? */
	
	if(x[n-1] > x[0])
		increase = YES;
	else
		increase = NO;
	/* do end member  test */
	if(increase){
		if(xval < x[0]){
			*jind = -1 ;
			return ;
		}
		if(xval > x[n-1]){
			*jind = n ;
			return ;
		}
	} else {
		if(xval > x[0]){
			*jind = 0 ;
			return ;
		}
		if(xval < x[n-1]){
			*jind = n ;
			return ;
		}
	}
	/* safety */
	if(xval == x[0]){
		*jind = 0 ;
		return ;
	}
	if(xval == x[n-1]){
		*jind = n-1 ;
		return ;
	}
	/* 
	jlow and jup are the current extremal bounds
	jmid is the current test
	*/
	jlow = 0 ;
	jup = n -1 ;
	while( ABS(jup-jlow) != 1){
		jmid = ( jlow + jup)/2 ;
		if(increase){
			if(xval > x[jmid]){
				jlow = jmid ;
			} else {
				jup  = jmid ;
			}
		} else {
			if(xval > x[jmid]){
				jup = jmid ;
			} else { 
				jlow  = jmid ;
			}
		}
	}
	*jind = jlow ;
	return ;
}

void gsac_evaamp(float freq,float *famp)
{
	int jind;
	float p;
	float amp;
	
	*famp = 1.0;
	arr_locate(evfamp, evnamp, freq, &jind);
	if(jind < 0) {
		amp = evaamp[0];
	} else if(jind >= 0 && jind < evnamp-1){
		p = (freq - evfamp[jind])/(evfamp[jind+1] - evfamp[jind]);
		amp = (1.0 -p)*evaamp[jind] + p*evaamp[jind+1];
	} else {
		amp = evaamp[evnamp-1];
	}
	*famp = amp;
}


int gsac_set_eval(char *ampname)
{
	FILE *fin;
	char *v;
	char *p;
	char s1[100], s2[100];
	double v1, v2;
	if((fin = fopen(ampname, "r")) ==NULL)
		return NO;	
	evnamp = 0;
	while(sfgetline(fin, input_lineptr, NCHAR) != -1 ){
		p = &input_lineptr[0];
		/* get rid of initial blanks */
	/* HACK 17 JAN 2007 since old Sparcs did not have isblank
			while (*p && isblank(*p) )
	*/
		while (*p && isspace(*p) )
			p++;
		sscanf(p,"%s %s",s1,s2);
		v1 = strtod(s1,&v);
		if(*v != '\0'){
			printf("Error in string %s of input %s %s\n",s1,s1,s2);
			return NO;
			}
		v2 = strtod(s2,&v);
		if(*v != '\0'){
			printf("Error in string %s of input %s %s\n",s2,s1,s2);
			return NO;
			}
		evfamp[evnamp] = v1;
		evaamp[evnamp] = v2;
		evnamp++;
		if(evnamp == NEVAL){
			printf("Number of lines in %s exceeds limit of %d\n",
					ampname,NEVAL);
			return NO;
		}
	}
	fclose(fin);

	return YES;
}

/* I need to read a line and parse tokens in the gsac_trans.c However I
	do not require the Prompts 
*/
/* http://www.math.ncu.edu.tw/~shann/Teach/comp/c-samp/getline.htm */
int sfgetline(FILE *fp, char s[], int len) {
    char *t;
    int c, lim;

    lim=len;

    t = s;
    while (--lim>1 && (c=getc(fp)) != EOF && c != '\n'){
        *s++ = c;
    }
    if (c == '\n')
        *s++ = '\0';
    else if (lim == 1) {
	*s++ = '\0';
	fprintf(stderr, "WARNING. fgetline: Line too long, split.\n");
    } else if(c == EOF) {
	    return EOF;
    }
    *s = '\0';
    return s - t;
}

void smooth_5(double *x, double *y, int np){
/* apply a five point smoother */
int nx, mx, mxx, m;
	/* special case first 5 points */
	for(m=0;m<2;m++){
		y[m] = (x[m] + x[m+1] + x[m+2] )/3.0;
	}
	/* smoothing interior points */
	for(m=2; m<np -2 ;m++){
		y[m] = (x[m-2] + x[m-1] + x[m] + x[m+1] + x[m+2])/5.0;
	}
	/* special case last 5 points */
	for(m=np-2; m< np;m++){
		y[m] = (x[m-2] +x[m-1] + x[m])/3.0;
	}
}

void doplot(int n, float df)
{
/* plot the observed spectra, the new low and high noise models, and
   noise spectra from the list file */
	int i,kolor,nkolor;
	float freq;
	float x0, y0, xx, yy;
	float x,y;
	float xlen, ylen;
	float xmin, ymin, xmax, ymax;
	int ipen;
	float per;

	float width = 0.025;

	/* plot position coordinates and limits in user space */
	x0 = 1.0;
	y0 = 1.0;
	xlen = 8.0;
	ylen = 6.0;
	xmin = 0.01 ;
	xmax = 200;
	ymin = -200;
	ymax = -50;

	/* initiate graphics */
	pinitf("SACPSD.PLT");

	/* plot the axes */
	doplot_axes(x0,y0,xlen,ylen,xmin,xmax,ymin,ymax);
	/* define a clip region */
	gclip("ON" ,x0,y0,x0+xlen,y0+ylen);

	/* plot the PSD if computed */
	doplot_obspsd(x0,y0,xlen,ylen,xmin,xmax,ymin,ymax,n,df,do_smooth_5);
	/* plot least squares smoothed spectrum */
        position_y = y0 + ylen -0.3;
	if(n > 0){
		position_y -= 0.40;
		do_plotfile( x0, y0, xlen, ylen, xmin, xmax, ymin, 
				ymax,  BLACK,  0.025, " ", smoothedpsdname);
		position_y -= 0.10;
		}


	/* draw the lower an upper limits of the noise models */
	do_plotnm(x0,y0,xlen,ylen,xmin,xmax,ymin,ymax,nlnm);
	do_plotnm(x0,y0,xlen,ylen,xmin,xmax,ymin,ymax,nhnm);

	/* plot other spectra */
	/* next do those contained in the listfile */
	if(have_otherpsd == YES){
		if((pfileptr = fopen(listfile,"r") )!=NULL){
			while(fscanf(pfileptr,"%s %d %s\n",pfile, &kolor, ptitle)!=EOF)
			{
				do_plotfile( x0, y0, xlen, ylen, xmin, xmax, ymin, 
					ymax,  kolor,  width, ptitle, pfile);
				position_y -= 0.10;
			}
			fclose(pfileptr);
		}
	}
	if(have_anotherpsd == YES){
		if((pfileptr = fopen(alistfile,"r") )!=NULL){
			nkolor = 0;
			/* read the list to get the number of entries for the kolor */
			while(fscanf(pfileptr,"%s\n",pfile)!=EOF)
			{
				nkolor++;
			}
			rewind(pfileptr);
			i=0;
			while(fscanf(pfileptr,"%s\n",pfile)!=EOF)
			{
				ptitle[0] = '\0';
				if(nkolor == 1){
					kolor = 1000;
				} else {
				kolor = 1000 + (int)(100.0*i/(float)(nkolor-1) );
				}
				if(kolor > 1100)kolor = 1100;
				do_plotfile( x0, y0, xlen, ylen, xmin, xmax, ymin, 
					ymax,  kolor,  width, ptitle, pfile);
				position_y -= 0.10;
				i++;
			}
			fclose(pfileptr);
		}
	}
	
	/* reset the clip region */
	gclip("OFF",x0,y0,x0+xlen,y0+ylen);
	/* reset color since graphics can be concatenated */
	newpen(BLACK);
	/* terminate graphics */
	pend();
}

void doplot_axes(float x0,float y0,float xlen,float ylen,float xmin,
	float xmax,float ymin,float ymax)
{
/* draw axes */
	newpen(BLACK);
	gbox(x0,y0,x0+xlen,y0+ylen);
	dologx(x0,y0+ylen,xlen,xmax,xmin,0.1,NO ,NO ,NO , 0,"          ");
	dologx(x0,y0     ,xlen,xmax,xmin,0.1,YES,NO ,YES,10,"Period (s)");
	doliny(x0     ,y0,ylen,ymax,ymin,0.1,NO ,YES,YES,34,
		" Power[10 log10(m**2/s**4/Hz)](db)");
	doliny(x0+xlen,y0,ylen,ymax,ymin,0.1,YES,YES,NO , 0," ");
}

void doplot_obspsd(float x0,float y0,float xlen,float ylen,float xmin,
	float xmax,float ymin,float ymax, int n, float df, int do_smooth)
{
/* plot the raw spectra in RED - these may or may not have been
	smoothed */
	int ipen, i;
	float per, xx, yy;
	char ostr[50];


	/* draw the noise spectrum */
	if(n > 0){
		newpen(RED);
		ipen = 3;
		for (i=1;i<n;i++){
			per = 1.0/(i*df);
			if(per >= pmin &&per <= pmax){
				if(per >= xmin && per <= xmax && period[i] > 0.0){
					xx = x0 + xlen*log10(period[i]/xmin)
						/log10(xmax/xmin);
					yy = y0 + ylen*(psdave[i] - ymin)/
						(ymax-ymin);
					plot(xx,yy,ipen);
					ipen = 2;
				}
			}
		}
		/* lift the pen */
		plot(xx,yy,3);
		gwidth(0.0);
		newpen(BLACK);

	
		/* annotate the plot with information about the station etc */
		sprintf(ostr,"%4.4d%2.2d%2.2d%2.2d%2.2d%2.2d.%3.3d\n",
			t.year,t.month,t.day,t.hour,t.minute,nzsec,nzmsec);
		gleft(x0+0.2 ,y0+ylen -0.30,0.10,"Station",0.0);
		gleft(x0+1.0 ,y0+ylen -0.30,0.10,"Channel",0.0);
		gleft(x0+1.9 ,y0+ylen -0.30,0.10,"Net",0.0);
		gleft(x0+2.3 ,y0+ylen -0.30,0.10,"Loc",0.0);
		gleft(x0+3.0 ,y0+ylen -0.30,0.10,"Date_Time",0.0);
		gleft(x0+4.7 ,y0+ylen -0.30,0.10,"Doy",0.0);
		gleft(x0+5.2 ,y0+ylen -0.30,0.10,"dt",0.0);
		gleft(x0+5.9 ,y0+ylen -0.30,0.10,"Window",0.0);
  	
		gleft(x0+0.3 ,y0+ylen -0.45,0.10,kstnm,0.0);
		gleft(x0+1.2 ,y0+ylen -0.45,0.10,kcmpnm,0.0);
		gleft(x0+1.9 ,y0+ylen -0.45,0.10,knetwk,0.0);
		gleft(x0+2.3 ,y0+ylen -0.45,0.10,khole,0.0);
		gleft(x0+2.7 ,y0+ylen -0.45,0.10,ostr,0.0);
        	number(x0+4.7,y0+ylen -0.45,0.10, nzjday,0.0,-1);
        	number(x0+5.1,y0+ylen -0.45,0.10, dt    ,0.0,3);
        	number(x0+5.8,y0+ylen -0.45,0.10, window,0.0,3);

		gleft(x0+0.2 ,y0      +0.75,0.10,"nft",0.0);
		gleft(x0+1.2 ,y0      +0.75,0.10,"nseg",0.0);

        	number(x0+0.2,y0      +0.60,0.10, nft,0.0,-1);
        	number(x0+1.2,y0      +0.60,0.10, nseg,0.0,-1);

	}
}

void do_plotnm(float x0,float y0,float xlen,float ylen,float xmin,
	float xmax,float ymin,float ymax,struct _PXY *nm){
/* plot the noise models which are contained in a data structure */
	int i, ipen;
	float xx,yy;
	newpen(BLACK);
	for(i=0 ,ipen=3 ; nm[i].per > 0 ; i++){
			xx = x0 + xlen*log10(nm[i].per/xmin)/log10(xmax/xmin);
			yy = y0 + ylen*(nm[i].psd - ymin)/(ymax-ymin);
			plot(xx,yy,ipen);
			ipen = 2;
	}
	/* lift the pen */
	plot(xx,yy,3);
}

void do_plotfile(float x0,float y0,float xlen,float ylen,float xmin, float xmax,	float ymin,float ymax, int kolor, float width, char *title, char *pfile)
{
/* open the file pfile which has two entires per line
	period PSD(db)
   plot using the color kolor with the line width gwidth
*/
	float x, y, xx, yy;
	int ipen;
	if((ofileptr = fopen(pfile,"r") )!=NULL){
		printf("plotting %s\n",pfile);
		ipen = 3;
		newpen(kolor);
		gwidth(width);
		while(fscanf(ofileptr,"%f %f\n",&x,&y)!=EOF){
			if(x >= xmin && x <= xmax){
				xx = x0 + xlen*log10(x/xmin)/log10(xmax/xmin);
				yy = y0 + ylen*(y - ymin)/(ymax-ymin);
				plot(xx,yy,ipen);
				ipen = 2;
			}
		}
		/* lift the pen */
		plot(xx,yy,3);
		if(do_title){
			plot(x0+xlen-2.8,position_y+0.03,3);
			plot(x0+xlen-2.5,position_y+0.03,2);
			plot(x0+xlen-2.5,position_y+0.03,3);
		}
		gwidth(0.0);
		newpen(BLACK);
		if(do_title)
			gleft(x0+xlen-2.4,position_y,0.07,pfile,0.0);
		fclose(ofileptr);
	} else {
		printf("cannot open %s\n",pfile);
	}
}

void deblank(char *s, int n)
{
	/* replace blanks with underscore
		also ensure that n+1 element is null */
	s[n] = '\0';
	while (*s != '\0'){
		if(*s == ' ')*s = '_';
		s++;
	}
}


void evarray(float xval, double *yval, float *x, double *y, int n)
{
	int jind;
	float p;
	double amp;
	
	*yval = 1.0;
	arr_locate(x, n, xval, &jind);
	if(jind < 0) {
		amp = y[0];
	} else if(jind >= 0 && jind < n-1){
		p = (xval - x[jind])/(x[jind+1] - x[jind]);
		amp = (1.0 -p)*y[jind] + p*y[jind+1];
	} else {
		amp = y[n-1];
	}
	*yval = amp;
}

