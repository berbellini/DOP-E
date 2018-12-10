/* Changes:
	20 JUL 2007 - place the time of the maximum of the cross-correlation
		into T9
*/

#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include	"csstim.h"


extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	CORR_DFLT	0
#define	CORR_MAST	1
#define	CORR_LENG	2
#define	CORR_SUF	3
#define	CORR_TWO	4
#define	CORR_NUMB	5
#define	CORR_REV	6

static int   corr_master = 0;
static float corr_length;
static int   corr_dolength = NO;
static int   corr_number;
static int   corr_donumber = NO;
static int   corr_dotwo = NO;
static char  corr_suffix[80];
static int   corr_dorev = NO;
static int   corr_dosuffix = NO;
static float *tx = (float *)NULL;
static float *ty = (float *)NULL;
static float *datax = (float *)NULL;
static float *datay = (float *)NULL;

void gsac_corr(float *x, float dtx, int nptsx, float *y, float dty, int nptsy, 
		float **yout, int *nptsout);


struct arghdr corrarg[] = {
	{CORR_DFLT, "DEFAULT", IHDR, 0, 0, NO, "",  1},
	{CORR_SUF , "SUFFIX" , CHDR, 0, 1, NO, "SUFFIX suffix",  1},
	{CORR_MAST, "MASTER" , IHDR, 0, 1, NO, "MASTER n",  1},
	{CORR_LENG, "LENGTH" , RHDR, 0, 1, NO, "LENGTH window ", 1},
	{CORR_NUMB, "NUMBER" , IHDR, 0, 1, NO, "NUMBER n ", 1},
	{CORR_TWO , "2"      , IHDR, 0, 0, NO, "2  ", 1},
	{CORR_REV , "REVERSE", IHDR, 0, 0, NO, "Reverse  ", 1},
	{0,	""		, IHDR, 0, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float corr_real[10];
int   corr_int [10];
int   corr_yn;
int   corr_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_corr(int ncmd, char **cmdstr)
{
	int i;
	char instr[80];
	/* set up default */
	corr_dosuffix = YES; 
	corr_suffix[0]='\0';
		strcpy(corr_suffix,".cor");

	/* parsing code here */


	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, corrarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; corrarg[i].key[0] != '\0' ; i++){
		if(corrarg[i].used > 0){
			if(corrarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, corrarg[i].key, 
					corrarg[i].mfit,corrarg[i].narg, corr_real);
			} else if(corrarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, corrarg[i].key, 
					corrarg[i].mfit,corrarg[i].narg, corr_int );
			} else if(corrarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, corrarg[i].key, 
					corrarg[i].mfit,corrarg[i].narg, &corr_yn );
			} else if(corrarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, corrarg[i].key, 
					corrarg[i].mfit,corrarg[i].narg, &corr_num );
			} else if(corrarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, corrarg[i].key, 
					corrarg[i].mfit, corrarg[i].narg, instr );
			}
			switch(corrarg[i].id){
				case CORR_SUF:
					if(strlen(instr) < 80){
						strcpy(corr_suffix, instr);
						corr_dosuffix = YES;
					}
					break;
				case CORR_DFLT:
					corr_dosuffix = NO;
					corr_master = 0;
					corr_dolength = NO;
					corr_dotwo = NO;
					corr_donumber = NO;
					corr_dorev = NO;
					break;
				case CORR_MAST:
					corr_master = corr_int[0];
					break;
				case CORR_LENG:
					corr_length = corr_real[0];
					corr_dolength = YES;
					break;
				case CORR_NUMB:
					corr_number = corr_int[0];
					if(corr_number > 0)
						corr_donumber = YES;
					else
						corr_donumber = NO ;
					break;
				case CORR_TWO:
					corr_dotwo = YES;
					break;
				case CORR_REV:
					corr_dorev = YES;
					break;


			}
		}
	}
			
}


void gsac_exec_corr(void)
{
	int k, ntrc, nptsx, nptsy, npts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float dtx, dty;
	int i;
	float stlax, stlox;
	float stlay, stloy;
	float stla, stlo, evla, evlo;
	float dist, az, baz, gcarc;
	int month, day;
	double tzbegx, tzbeg;
	float tdiff;



	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
printf("CORRELATE MASTER %d DOTWO %s DOREVERSE %s LENGTH %s NUMBER %s \n",
		corr_master,
		corr_dotwo == YES ? "YES" : "NO",
		corr_dorev == YES ? "YES" : "NO",
		corr_dolength == YES ? "YES" : "NO",
		corr_donumber == YES ? "YES" : "NO");
	if(corr_master < 0 || corr_master > ntrc){
		printf("corr_master = %d not in range of [0,%d]\n",
				corr_master,ntrc-1);
		return;
	}


	/* CHECK TO SEE THAT ALL DT ARE THE SAME ELSE RETURN */
	/* define the master trace */
	nptsx = sacdata[corr_master].sachdr.ihdr[H_NPTS];
	dtx = sacdata[corr_master].sachdr.rhdr[H_DELTA];
	for ( k=0 ; k < ntrc ; k ++){
		dty = sacdata[k].sachdr.rhdr[H_DELTA];
		if(ABS (dtx - dty) > 0.01 * ABS(dtx)){
			printf("DT not equal - correlation not done\n");
			printf("Files %d %d with %f %f, respectively\n",
				corr_master,k,dtx,dty);
			return;
		}
	}
	/* set timing parameters for plot */
	gsac_control.begmin= 1.0e+30;
	gsac_control.endmax=-1.0e+30;
	gsac_control.fft = NO;

	/* SAVE THE MASTER TRACE INFORMATION SINCE THIS WILL BE 
	 * REQUIRED MORE THAN ONCE */
	if(tx == (float *)NULL)
		tx = (float *)calloc(nptsx,sizeof(float));
	else
		tx = (float *)realloc(tx,nptsx*sizeof(float));
	for(i=0;i< nptsx; i++)
		tx[i] = sacdata[corr_master].sac_data[i];
	stlax = sacdata[corr_master].sachdr.rhdr[H_STLA];
	stlox = sacdata[corr_master].sachdr.rhdr[H_STLO];
	tzbegx = sacdata[corr_master].tzbeg;
	/* process the traces */
	for ( k=0 ; k < ntrc ; k ++){
		nptsy = sacdata[k].sachdr.ihdr[H_NPTS];
		dty = sacdata[k].sachdr.rhdr[H_DELTA];
		stlay = sacdata[k].sachdr.rhdr[H_STLA];
		stloy = sacdata[k].sachdr.rhdr[H_STLO];
		tzbeg = sacdata[k].tzbeg;
		if(ty == (float *)NULL)
			ty = (float *)calloc(nptsy,sizeof(float));
		else
			ty = (float *)realloc(ty,nptsy*sizeof(float));
		for(i=0;i< nptsy; i++)
			ty[i] = sacdata[k].sac_data[i];
		/* filter the trace */
		gsac_corr(tx, dtx, nptsx, ty, dty, nptsy, 
				&sacdata[k].sac_data, &npts);
		/*
		*/
		/* update the header values */
		sacdata[k].sachdr.ihdr[H_NPTS] = npts ;
		/*
		*/
		sacdata[k].sachdr.rhdr[H_EVLA] = stlax ;
		sacdata[k].sachdr.rhdr[H_EVLO] = stlox ;
		stla = sacdata[k].sachdr.rhdr[H_STLA];
		stlo = sacdata[k].sachdr.rhdr[H_STLO];
		evla = sacdata[k].sachdr.rhdr[H_EVLA];
		evlo = sacdata[k].sachdr.rhdr[H_EVLO];
		if(stla != -12345. && stlo != -12345. 
			&& evla != -12345. && evlo != -12345.){
			delaz( evla,  evlo,  stla,  stlo,  
				&gcarc,  &az,  &baz,  &dist);
			sacdata[k].sachdr.rhdr[H_DIST] = dist;
			sacdata[k].sachdr.rhdr[H_AZ] = az;
			sacdata[k].sachdr.rhdr[H_BAZ] = baz;
			sacdata[k].sachdr.rhdr[H_GCARC] = gcarc;
		}
		getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
		sacdata[k].sachdr.ihdr[H_NPTS] = npts;
		if(corr_dorev == YES)
			tdiff = tzbegx - tzbeg;
		else
			tdiff = tzbeg - tzbegx;
		sacdata[k].sachdr.rhdr[H_B] = (- npts/2 )*dty
			+ tdiff ;
		sacdata[k].sachdr.rhdr[H_E] = sacdata[k].sachdr.rhdr[H_B] 
			+ (npts -1 )*dty;
		sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_O]  =     0.0 ;
		sacdata[k].sachdr.rhdr[H_A]  = -12345.0;
		sacdata[k].sachdr.rhdr[H_T0] = -12345.0;
		sacdata[k].sachdr.rhdr[H_T1] = -12345.0;
		sacdata[k].sachdr.rhdr[H_T2] = -12345.0;
		sacdata[k].sachdr.rhdr[H_T3] = -12345.0;
		sacdata[k].sachdr.rhdr[H_T4] = -12345.0;
		sacdata[k].sachdr.rhdr[H_T5] = -12345.0;
		sacdata[k].sachdr.rhdr[H_T6] = -12345.0;
		sacdata[k].sachdr.rhdr[H_T7] = -12345.0;
		sacdata[k].sachdr.rhdr[H_T8] = -12345.0;
		sacdata[k].sachdr.rhdr[H_T9] = -12345.0;
		sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
		sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
		sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		/* find point in array correspoding to depmax 
			and place into H_T9 */
		for(i=0;i<npts;i++){
			if(sacdata[k].sac_data[i] == depmax){
				sacdata[k].sachdr.rhdr[H_T9] = 
					sacdata[k].sachdr.rhdr[H_B] 
					+ (i  )*dty;
			}
		}
		/*
		 * NOTE HERE GET THE EPOCH TIME AND THEN SUBTRACT 
		*/
		sacdata[k].tzref = 0.0;
		sacdata[k].tzbeg = sacdata[k].tzref + 
			sacdata[k].sachdr.rhdr[H_B] ;
		sacdata[k].tzend = sacdata[k].tzref + 
			sacdata[k].sachdr.rhdr[H_E] ;
		sacdata[k].tzbegx = sacdata[k].tzbeg;
		sacdata[k].tzendx = sacdata[k].tzend;
		if(sacdata[k].tzbeg < gsac_control.begmin)
			gsac_control.begmin = sacdata[k].tzbeg;
		if(sacdata[k].tzend > gsac_control.endmax)
			gsac_control.endmax = sacdata[k].tzend;


		/* for cross correlation - reference time is always 0 lag */
		/* fro cross correlation origin time is always zero lag */
		etoh(sacdata[k].tzref, &sacdata[k].sachdr.ihdr[H_NZYEAR], 
			&sacdata[k].sachdr.ihdr[H_NZJDAY], &month, &day,
			&sacdata[k].sachdr.ihdr[H_NZHOUR], 
			&sacdata[k].sachdr.ihdr[H_NZMIN],
			&sacdata[k].sachdr.ihdr[H_NZSEC], 
			&sacdata[k].sachdr.ihdr[H_NZMSEC]);
		/* create the name of the output file */
		chofname(sacdata[corr_master].schdr[H_KSTNM],
			sacdata[corr_master].schdr[H_KCMPNM],
			sacdata[k].schdr[H_KSTNM],
			sacdata[k].schdr[H_KCMPNM],
			sacdata[k].sac_ofile_name,
			corr_dosuffix, corr_suffix);


	}
}

void gsac_corr(float *x, float dtx, int nptsx, float *y, float dty, int nptsy, 
		float **yout,  int *nptsout)
{
	int n2, n22, n2out;
	int nseg, ns;
	int i, j, ii;
	int jr, ji, kr, ki;
	float df;
	float xr, xi, yr, yi;
	float fac;


	/* first determine the size of the FFT to be used */
	if(corr_dolength == YES || corr_donumber == YES){
		if(corr_dolength == YES){
			n2 = npow2(corr_length/dtx);
			*nptsout = n2 ; 
			n2out = n2;
			if(corr_dotwo == YES)n2 *=2;
			nseg = nptsx/n2 + 1;
		} else if(corr_donumber == YES){
			n2 = npow2(MAX(nptsx/corr_number, nptsy/corr_number));
			*nptsout = n2 ; 
			n2out = n2;
			if(corr_dotwo == YES)n2 *=2;
			nseg = nptsx/n2 + 1;
		}
	} else {
		n2 = npow2(MAX(nptsx, nptsy));
		n2out = n2;
		*nptsout = n2 ; 
		if(corr_dotwo == YES)n2 *=2;
		nseg = 1;
	}
	/* initialise the output stream - note that yout already exists */
	*yout = (float *)realloc( *yout,n2out*sizeof(float));
	for(i=0;i < n2out ; i++)
		(*yout)[i] = 0.0;

	n22 = n2 / 2;
	/* define the arrays for the FFT */
	if(datax == (float *)NULL)
		datax = (float *)calloc(2*n2,sizeof(float));
	else
		datax = (float *)realloc(datax,2*n2*sizeof(float));
	if(datay == (float *)NULL)
		datay = (float *)calloc(2*n2,sizeof(float));
	else
		datay = (float *)realloc(datay,2*n2*sizeof(float));
	/* now fill the arrays and get the FFT */
/*
printf("nseg %d n2 %d nptxx %d nptsy %d nseg*n2 %d\n",nseg,n2,nptsx,nptsy,nseg*n2);
*/
	ii = 0;
	for(ns=0 ; ns < nseg ; ns++){
		printf("SEGMENT ns=%d of %d, data(%d) to data(%d) \n",ns,nseg,ii,MIN(ii+n2-1, nptsy));
		for(i=0 , j = 0 ; i < n2 ; i++){
			if((i+ii) < nptsx )
				datax[j] = x[i + ii];
			else
				datax[j] = 0.0;
			if((i+ii) < nptsy )
				datay[j] = y[i + ii];
			else
				datay[j] = 0.0;
			/* imaginary part of the time series */
			j++;
			datax[j] = 0.0;
			datay[j] = 0.0;
			j++;
		}
		ii += n2;
		/* get the DFT  with the proper dimensions */
		df = 0.0;
		four(datax,n2,-1,&dtx,&df);
		four(datay,n2,-1,&dtx,&df);
		/* form the cross correlation  in the frequency domain 
		 * conjg X  Y
		 * 
		 * also apply the centering operator 
		 * 	This is done as follows: If x(k) <=> X(n) for 
		 * 	0 <=k < N  and 0 <=n < N, then
		 * 	x(k - N/2) <=> X(n) exp (- j pi n) This is done here
		 * 	by the artifice of multiplying by fac and then
		 * 	changing the sign of fac 
		 * 	*/
		fac = 1.0;
		for(i=0, j = 0; i <= n22 ; i++){
			jr = j;
			ji = j +1;
			xr = datax[jr];
			xi = datax[ji];
			yr = datay[jr];
			yi = datay[ji];
			datax[jr] = (xr*yr + xi*yi)*fac;
			if(corr_dorev)
				datax[ji] = - (xr*yi - yr*xi)*fac;
			else
				datax[ji] =   (xr*yi - yr*xi)*fac;
			fac = -1.0 * fac;
			/* ensure a real time series which says
				real part even, imaginary part odd */
			if(i>0){
				if(i == n22){
					datax[n2+1] = 0.0;
				} else {
					kr = 2*(n2 - i );
					ki = kr + 1;
					datax[kr] =   datax[jr] ;
					datax[ki] = - datax[ji] ;
				}
			}
			j+=2;
		}
		/* apply the inverse FFT */
		four(datax,n2,+1,&dtx,&df);
		/* stack the results */
		if(corr_dotwo == YES){
			/* we do not want first and last fourths , just the
				center */
			for(i=0 , j=n2out ; i <  n2out ; i++){
				(*yout)[i] += datax[j];
				j+=2 ;
			}
		} else {
			for(i=0 , j=0 ; i < n2out ; i++){
				(*yout)[i] += datax[j];
				j+=2 ;
			}
		}

	}
}

