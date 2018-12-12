/* implement TRANSFER and FILTER commands 
 * The transfer is actually formed by parsing the command line
 * and then invoking filter twice to remove and then apply a 
 * response. THUS FILTER is a new GSAC command */

/* CHANGES
	15 JUN 2006 - error at line 933
		Change from	for(i=0 ; i < npoles ; i++){
		Change to	for(i=0 ; i < nzeros ; i++){
			error at 966-967
		Change from	zreal[ip] = v1;
				zimag[ip] = v2;
		Change to	zreal[iz] = v1;
				zimag[iz] = v2;
	16 SEP 2006 Build in a skip blank space in all response files
		using isblank
	28 NOV 2006 put in warning that transfer requires both a FROM and a TO
*/


#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int  sfgetline(FILE *fp, char s[], int lim);

#define TR_FROM 0
#define TR_TO	1
#define TR_FREQ	2
#define TR_PZ	3
#define TR_FAP	4
#define TR_EVAL	5
#define TR_APPLY 6
#define TR_REMOVE 7
#define TR_NONE 8
#define TR_DISP 9
#define TR_VEL 10
#define TR_ACC 11
#define TR_NOFREQ	12

#define FROM 0
#define TO   1

#define NONE 0
#define DISP 0
#define VEL  1
#define ACC 2
#define POLEZERO 3
#define EVAL 4
#define FAP	5
#define SUBTYPE 6


struct arghdr trarg[] = {
	{TR_FROM, "FROM"	, XHDR, 0, 0, YES, "FROM",-1},
	{TR_TO  , "TO"	        , XHDR, 0, 0, YES, "TO",-1},
	{TR_FREQ, "FREQLIMITS"	, RHDR, 0, 4, YES, "FREQLIMITS F1 F2 F3 F4",4},
	{TR_PZ  , "POLEZERO"	, CHDR, 0, 1, YES, "POLEZERO pzfile ", 1},
	{TR_FAP , "FAPFILE"	, CHDR, 0, 1, YES, "FAPFILE fapfile ", 2},
	{TR_EVAL, "EVAL"	, CHDR, 0, 2, YES, "EVAL afile pfile ", 1},
	{TR_REMOVE, "REMOVE"	, XHDR, 0, 0, YES, "REMOVE", 1},
	{TR_APPLY  , "APPLY"	        , XHDR, 0, 0, YES, "APPLY", 1},
	{TR_NONE  , "NONE"	        , XHDR, 0, 0, YES, "NONE",-1},
	{TR_DISP  , "DISP"	        , XHDR, 0, 0, YES, "DISP",-1},
	{TR_VEL  , "VEL"	        , XHDR, 0, 0, YES, "VEL",-1},
	{TR_ACC  , "ACC"	        , XHDR, 0, 0, YES, "ACC",-1},
	{TR_NOFREQ, "NOFREQLIMITS"	, RHDR, 0, 0, YES, "NOFREQLIMITS",3},
	{-10, ""        , CHDR, 0, 0, YES, "" ,-1}
};

/* frequencies for freqlimit control important for deconvolution */
static float tr_f1 = -10.;
static float tr_f2 =  -5.;
static float tr_f3 =  1.0e6;
static float tr_f4 =  1.0e7;

static float *y = (float *)NULL;

static int gsac_instype;	/* array to define FROM and TO parts of
				   instrument response */
/* data structures for each instrument response type 
 * the first subscript refers to the FROM and the second to the TO */

/* EVAL  */
#define NEVAL 400000
char aname[1000];
char pname[1000];
int evnamp;	       /* number of amplitude entries */
int evnpha;	       /* number of amplitude entries */
float evaamp[NEVAL];   /* each files has its own frequency column */
float evfamp[NEVAL];   /* each files has its own frequency column */
float evapha[NEVAL];   /* each files has its own frequency column */
float evfpha[NEVAL];   /* each files has its own frequency column */

/* FAP  */
/* FAP not implemented yet */
#define NFAP 100000
char fapname[1000];
int fapnum;	       /* number of amplitude-phase entries */
float fapamp[NEVAL];   /* file amplitude response column */
float fapfrq[NEVAL];   /* file frequency  column */
float fappha[NEVAL];   /* file phase response column */

/* POLEZERO
 * The ridiculously long numbers of real and imaginary parts 
 * here 500 are used to avoid having to worry about dynamic memory
 * the overhead head is not that large*/
#define NPAZ 500
static char pzname[1000];	/* name of the pole zero file */
static int  npoles;		/* number of zeros */
static int  nzeros;		/* number of poles */
static double zreal[NPAZ];	/* complex zero specification */
static double zimag[NPAZ];
static double preal[NPAZ];	/* complex pole specificaiton */
static double pimag[NPAZ];
static double pzconst;

/* special names for TRANSFER which will be implemented using
 * two FILTER calls
 * */
int gstr(int ncmd, char **cmdstr, char *pat, int n);
static int stage_instype[2] ;
static char fpz[2][1000];
static char fap[2][1000];
static char fevala[2][1000];
static char fevalp[2][1000];

static char tmpstr[1000];
static char *Top[2] = { "FROM", "TO" };
struct Second_ {
	int id ; char *sec ; int narg ; };
static struct Second_ Second[] = {
	{ ACC, "ACC", 0},
	{ VEL, "VEL", 0},
	{ DISP, "DISP", 0},
	{ NONE, "NONE", 0},
	{ POLEZERO, "POLEZERO", 1},
	{ EVAL, "EVAL", 1},
	{ FAP, "FAP", 1},
	{ FAP, "FAPFILE", 1},
	{ -1, "", 0}
};



/* these are temporary variables only used here */
float tr_real[10];
int do_tr;
static int from_to;
int do_taper;

#define NCHAR 200
static char input_lineptr[NCHAR] ;


void gsac_dofilt(float **x, float dt, int npts);
void gsac_dotrans(float **x, float dt, int npts);
void gsac_getresp(float freq, float *hr, float *hi, int from_to);
float gsac_taper3(float xl, float xh, float x);
void gsac_evalpz(float freq,float *fr,float *fi);
void gsac_evaleval(float freq,float *fr,float *fi);
void gsac_evalfap(float freq,float *fr,float *fi);
int gsac_set_pz(char *pzname);
int gsac_set_eval(char *ampname, char *phaname);
int gsac_set_fap(char *fapname);
void arr_locate(float *x, int n, float xval, int *jind);


void gsac_set_param_filter(int ncmd, char **cmdstr)
{
	int i;

	/* note when the testrg routine is used, if the argument is
		NO then you must use internal variables to define the 
		state of the operation - if you use YES, then things are
		not changed until the input is proven correct. An exmple of
		this concept with YES is the following:
		Assume we wish aa LP filter with fc 1 np 2 p 1 
		If we enter  fc 2 np2   there is a syntax error and we
		should not chnge the fc since the np2 is wrong. One way to
		do this in the code would be to do two calls

			if(testarc,ncmd, cmdstr, cmdargs, YES) is OK
			then
				testarc,ncmd, cmdstr, cmdargs, NO)
		*/

	if(testarg(ncmd, cmdstr, trarg, NO, YES))
		return;
	do_tr = NO;
	do_taper = NO;
	for(i=0 ; trarg[i].key[0] != '\0' ; i++){
		if(trarg[i].used > 0){
			if(trarg[i].ricell == RHDR){
			getargr(ncmd, cmdstr, trarg[i].key, 
				trarg[i].mfit, trarg[i].narg, tr_real );
			}
			switch(trarg[i].id){
				case TR_FREQ:
					tr_f1 = tr_real[0];
					tr_f2 = tr_real[1];
					tr_f3 = tr_real[2];
					tr_f4 = tr_real[3];
					do_taper = YES;
					break;
				case TR_NOFREQ:
					do_taper = NO;
					break;
				case TR_PZ:
					getargs(ncmd, cmdstr, trarg[i].key, 
						trarg[i].mfit,trarg[i].narg, pzname );
					gsac_instype = POLEZERO;
					break;
				case TR_EVAL:
					getargs2(ncmd, cmdstr, trarg[i].key, 
						trarg[i].mfit,trarg[i].narg, aname, pname );
					gsac_instype = EVAL;
					break;
				case TR_FAP:
					getargs(ncmd, cmdstr, trarg[i].key, 
						trarg[i].mfit,trarg[i].narg, fapname );
					gsac_instype = FAP;
					break;
				case TR_NONE:
					gsac_instype = NONE;
					break;
				case TR_DISP:
					gsac_instype = NONE;
					break;
				case TR_VEL:
					gsac_instype = VEL;
					break;
				case TR_ACC:
					gsac_instype = ACC;
					break;
				case TR_FROM:
				case TR_REMOVE:
					from_to = FROM;
					break;
				case TR_TO:
				case TR_APPLY:
					from_to = TO;
					break;

			}

		}
	}
	/* set up the filter parameters */
	switch(gsac_instype){
		case POLEZERO:
			do_tr = gsac_set_pz(pzname);
			break;
		case EVAL:
			do_tr = gsac_set_eval(aname,pname);
			break;
		case FAP:
			do_tr = gsac_set_fap(fapname);
			break;
	}
			
}

void gsac_set_param_trans(int ncmd, char **cmdstr)
{
        int i, ind;

        for(i=0; i < ncmd; i++)
                printf("%s ",cmdstr[i]);
        printf("\n");
	
	int from_to;
	/* lets parse */
	do_tr = NO ;
	do_taper = NO;

	for(from_to = 0 ; from_to < 2 ; from_to++ ){
		ind = gstr( ncmd, cmdstr, Top[from_to],strlen(Top[from_to]));
		if(ind < 0 || (ind+1) >= ncmd){
			do_tr = NO ;
printf("Did not find  a %s: Syntax is  TRANSFER FROM ... TO ... \n",Top[from_to]);
			return;
		}
		for(i=0; Second[i].sec[0] != '\0'; i++){
			strcpy(tmpstr,cmdstr[ind+1]);
			gsac_strupr(tmpstr);
			if(strcmp(tmpstr,Second[i].sec)==0){
				switch(Second[i].id){
					case ACC:
						stage_instype[from_to] = ACC;
						do_tr = YES;
						break;
					case VEL:
						stage_instype[from_to] = VEL;
						do_tr = YES;
						break;
					case DISP:
						stage_instype[from_to] = NONE;
						do_tr = YES;
						break;
					/* for the next we must have a
					 * subtype keyword, followed by
					 * one or more file names */
					case POLEZERO:
						strcpy(tmpstr,cmdstr[ind+2]);
						gsac_strupr(tmpstr);
						if(strcmp(tmpstr,"SUBTYPE")==0){
						strcpy(fpz[from_to],cmdstr[ind+3]);
							do_tr = YES;
						} else {
							do_tr = NO;
printf("Error in polezero syntax: polezero subtype polezero_file\n");
						}
						stage_instype[from_to] = POLEZERO;
						break;
					case EVAL:
						strcpy(tmpstr,cmdstr[ind+2]);
						gsac_strupr(tmpstr);
						if(strcmp(tmpstr,"SUBTYPE")==0){
						strcpy(fevala[from_to],cmdstr[ind+3]);
						strcpy(fevalp[from_to],cmdstr[ind+4]);
							do_tr = YES;
						} else {
							do_tr = NO;
printf("Error in evalresp syntax: eval subtype amp_file phase_file\n");
						}
						stage_instype[from_to] = EVAL;
						break;
					case FAP:
						strcpy(tmpstr,cmdstr[ind+2]);
						gsac_strupr(tmpstr);
						if(strcmp(tmpstr,"SUBTYPE")==0){
						strcpy(fap[from_to],cmdstr[ind+3]);
							do_tr = YES;
						} else {
							do_tr = NO;
printf("Error in fap syntax: fap subtype fap_file \n");
						}
						stage_instype[from_to] = FAP;
						break;
				}
			}
		}
	}
	/* now look for FREQLIMITS */
	ind = gstr( ncmd, cmdstr, "FREQL",5);
	if( (ind+4 ) >= ncmd) {
		do_tr = NO ;
		printf("Wrong syntax for FREQLIMITS f1 f2 f3 f4 \n");
		return;
	}
	if( ind < 0 )
		do_taper = NO;
	else {
		if(getargr(ncmd, cmdstr, "FREQLIMITS", 4,  4, tr_real ) > 0){
			tr_f1 = tr_real[0];
			tr_f2 = tr_real[1];
			tr_f3 = tr_real[2];
			tr_f4 = tr_real[3];

			do_taper = YES;
		} else {
			do_tr = NO ;
			printf("Wrong syntax for FREQLIMITS f1 f2 f3 f4 \n");
			return;
		}
	}
}

void gsac_exec_filter(void)
{
	int k, ntrc, npts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float dt;
	float permin, permax;

	/* only apply the filtering if permitted */
	if(do_tr == NO)
		return;


	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;

	/* process the traces */
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		dt = sacdata[k].sachdr.rhdr[H_DELTA];
		/* filter the trace */
	printf("filter processing %s npts %d dt %f\n",sacdata[k].sac_ifile_name,npts,dt);
		gsac_dofilt(&sacdata[k].sac_data, dt, npts);
		/* update the header values */
		getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
		sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
		sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
		sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		/* set USER1 = permin USER2 = permax */
		if(do_taper){
			/* if filter bound is set at default, finally use */
			if(tr_f3 > 0.0){
				permin = 1.0/tr_f3;
				/* if filter bound is set at default, finally use */
				if(sacdata[k].permin == -12345.)
					sacdata[k].permin = permin;
				/* never expand range of periods */
				if(sacdata[k].permin < permin)
					sacdata[k].permin = permin;
			}
			if(tr_f2 > 0.0){
				permax = 1.0/tr_f2;
				/* if filter bound is set at default, finally use */
				if(sacdata[k].permax == -12345.)
					sacdata[k].permax = permax;
				/* never expand range of periods */
				if(sacdata[k].permax > permax)
					sacdata[k].permax = permax;
			}
		}

	}
}

void gsac_exec_trans(void)
{
	int k, ntrc, npts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float dt;
	float permin, permax;

	/* only apply the filtering if permitted */
	if(do_tr == NO)
		return;


	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;

	/* process the traces */
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		dt = sacdata[k].sachdr.rhdr[H_DELTA];
		/* filter the trace */
/* DEBUG 
printf("gsac_dotrans k=%d of %d\n",k,ntrc);
printf("FROM %s %d\n",fpz[0],stage_instype[0]);
printf("TO   %s %d\n",fpz[1],stage_instype[1]);
 DEBUG */
		gsac_dotrans(&sacdata[k].sac_data, dt, npts);
		/* update the header values */
		getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
		sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
		sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
		sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		/* set USER1 = permin USER2 = permax */
		if(do_taper){
			/* if filter bound is set at default, finally use */
			if(tr_f3 > 0.0){
				permin = 1.0/tr_f3;
				/* if filter bound is set at default, finally use */
				if(sacdata[k].permin == -12345.)
					sacdata[k].permin = permin;
				/* never expand range of periods */
				if(sacdata[k].permin < permin)
					sacdata[k].permin = permin;
			}
			if(tr_f2 > 0.0){
				permax = 1.0/tr_f2;
				/* if filter bound is set at default, finally use */
				if(sacdata[k].permax == -12345.)
					sacdata[k].permax = permax;
				/* never expand range of periods */
				if(sacdata[k].permax > permax)
					sacdata[k].permax = permax;
			}
		}

	}
}

void gsac_dofilt(float **x, float dt, int npts)
{
	int i, j;
	int n2, n22;
	int jr, ji, kr, ki;
	float df, freq, f;
	float xr, xi, hr, hi, yr, yi;
	float taper;
	n2 = npow2(npts);
	/* define Nyquist index */
	n22 = n2 / 2;
	/* ensure that the temporary array is of the proper size 
	 * The y array will initially represent the complex
	 * time series, and then the complex spectra */

	if(y == (float *)NULL)
		y = (float *)calloc(2*n2,sizeof(float));
	else
		y = (float *)realloc(y,2*n2*sizeof(float));
	for(i=0, j=0; i < n2 ; i++){
		/* real part of time series */
		if(i < npts )
			y[j] = (*x)[i];
		else
			y[j] = 0.0;
		/* imaginary part of time series */
		j++;
		y[j++] = 0.0;
	}
	/* get the DFT  with the proper dimensions */
	df = 0.0;
	four(y,n2,-1,&dt,&df);

	/* now apply the instrument response 
	 * For a real times series, the spectra of the real part is
	 * even in frequency and the imaginary part is odd in frequency
	 * This is how the spectra of the negative frequencies is
	 * constructed 
	 * i : frequency index
	 * j : index into linear array. j = even is real, j= odd is imag */

	for(i=0, j = 0; i <= n22 ; i++){
		freq = i*df;
		jr = j   ;
		ji = j+1 ;
		xr = y[jr];
		xi = y[ji];
		/* safety for inverse */
		if(freq == 0.0)
			f = 0.01 * df;
		else
			f = freq ;
		gsac_getresp(f, &hr, &hi, from_to);
		/* taper */
		if(do_taper){
			if(freq < tr_f2)
				taper = gsac_taper3(tr_f1, tr_f2, freq);
			else if(freq > tr_f3)
				taper = gsac_taper3(tr_f4, tr_f3, freq);
			else
				taper = 1.0;
		} else {
			taper = 1.0;
		}
		hr *= taper;
		hi *= taper;
		yr = xr*hr - xi*hi ;
		yi = xr*hi + xi*hr ;
		y[jr] = yr;
		y[ji] = yi;
		/* ensure a real time series */
		if(i > 0){
			kr = 2*(n2 - i );
			ki = kr + 1;
			y[kr] =   y[jr] ;
			y[ki] = - y[ji] ;
		}
		j+=2 ;
	}
	/* take care of Nyquist frequency n2 is the real part at the Nyquist
	 * n2+1 is the imaginary part */

	y[n2+1] = 0.0;

	/* now obtain the inverse DFT */
	four(y,n2,+1,&dt,&df);
	/* now update the time series */
	for(i=0 , j=0 ; i < npts ; i++){
		(*x)[i] = y[j];
		j+=2 ;
	}
}

void gsac_dotrans(float **x, float dt, int npts)
{
	int i, j;
	int n2, n22;
	int jr, ji, kr, ki;
	float df, freq, f;
	float xr, xi, hr, hi, yr, yi;
	float taper;
	n2 = npow2(npts);
	n22 = n2 / 2;
	/* ensure that the temporary array is of the proper size 
	 * The y array will initially represent the complex
	 * time series, and then the complex spectra */

	if(y == (float *)NULL)
		y = (float *)calloc(2*n2,sizeof(float));
	else
		y = (float *)realloc(y,2*n2*sizeof(float));
	for(i=0, j=0; i < n2 ; i++){
		/* real part of time series */
		if(i < npts )
			y[j] = (*x)[i];
		else
			y[j] = 0.0;
		/* imaginary part of time series */
		j++;
		y[j++] = 0.0;
	}
	/* get the DFT  with the proper dimensions */
	df = 0.0;
	four(y,n2,-1,&dt,&df);

	/* now apply the instrument response 
	 * For a real times series, the spectra of the real part is
	 * even in frequency and the imaginary part is odd in frequency
	 * This is how the spectra of the negative frequencies is
	 * constructed 
	 * i : frequency index
	 * j : index into linear array. j = even is real, j= odd is imag */

	/* we apply this in stages 
	 * FIRST we do the FROM, the optional taper and then the TO */

	/* FROM 
	 * define the response */
	gsac_instype = stage_instype[FROM];
	switch(gsac_instype){
		case POLEZERO:
			do_tr = gsac_set_pz(fpz[FROM]);
			break;
		case FAP:
			do_tr = gsac_set_fap(fap[FROM]);
			break;
		case EVAL:
			do_tr = gsac_set_eval(fevala[FROM],fevalp[FROM]);
			break;
	}
	/* now implement FROM */
	for(i=0, j = 0; i <= n22 ; i++){
		freq = i*df;
		jr = j   ;
		ji = j+1 ;
		xr = y[jr];
		xi = y[ji];
		/* safety for inverse */
		if(freq == 0.0)
			f = 0.01 * df;
		else
			f = freq ;
		gsac_getresp(f, &hr, &hi, FROM);
		if(do_taper){
			if(freq < tr_f2)
				taper = gsac_taper3(tr_f1, tr_f2, freq);
			else if(freq > tr_f3)
				taper = gsac_taper3(tr_f4, tr_f3, freq);
			else
				taper = 1.0;
		} else {
			taper = 1.0;
		}
		hr *= taper;
		hi *= taper;
		yr = xr*hr - xi*hi ;
		yi = xr*hi + xi*hr ;
		y[jr] = yr;
		y[ji] = yi;
		if(i > 0){
			kr = 2*(n2 - i );
			ki = kr + 1;
			y[kr] =   y[jr] ;
			y[ki] = - y[ji] ;
		}
		j+=2 ;
	}
	/* TO 
	 * define the response */
	gsac_instype = stage_instype[TO];
	switch(gsac_instype){
		case POLEZERO:
			do_tr = gsac_set_pz(fpz[TO]);
			break;
		case FAP:
			do_tr = gsac_set_fap(fap[TO]);
			break;
		case EVAL:
			do_tr = gsac_set_eval(fevala[TO],fevalp[TO]);
			break;
	}
	/* now implement TO */
	for(i=0, j = 0; i <= n22 ; i++){
		freq = i*df;
		jr = j   ;
		ji = j+1 ;
		xr = y[jr];
		xi = y[ji];
		f = freq ;
		gsac_getresp(f, &hr, &hi, TO);
		yr = xr*hr - xi*hi ;
		yi = xr*hi + xi*hr ;
		y[jr] = yr;
		y[ji] = yi;
		if(i > 0){
			kr = 2*(n2 - i );
			ki = kr + 1;
			y[kr] =   y[jr] ;
			y[ki] = - y[ji] ;
		}
		j+=2 ;
	}
	/* take care of Nyquist frequency */
	y[n2+1] = 0.0;

	/* now obtain the inverse DFT */
	four(y,n2,+1,&dt,&df);
	/* now update the time series */
	for(i=0 , j=0 ; i < npts ; i++){
		(*x)[i] = y[j];
		j+=2 ;
	}
}

void gsac_getresp(float freq, float *hr, float *hi, int from_to)
{
	float fr, fi;
	float fact;
	switch (gsac_instype){
		case DISP: 
			fr = 1.0 ; 
			fi = 0.0;
			break;
		case VEL: 
			fr = 0.0 ; 
			fi = 6.2831853*freq ;
			break;
		case ACC: 
			fr = -(6.2831853*freq)*(6.2831853*freq) ;
			fi =  0.0;
			break;
		case POLEZERO: 
			gsac_evalpz(freq,&fr,&fi);
			break;
		case EVAL: 
			gsac_evaleval(freq,&fr,&fi);
			break;
		case FAP: 
			gsac_evalfap(freq,&fr,&fi);
			break;

	}
	if(from_to == FROM){
		fact = fr*fr + fi*fi;
		*hr =  fr/fact;
		*hi = -fi/fact;

	} else {
		*hr = fr;
		*hi = fi;
	}
}

float gsac_taper3(float xl, float xh, float x)
{
	/* a cubic taper operator
	 * gsac_taper3    x
	 *     0          xl    
	 *     0.5	0.5(xl+xh)
	 *     1.0	  xh
	 *     */
	float p;  /* mapping p=0 -> xl, p=1 -> xh */
	/* safety */
	if(xl == xh)
		return 1.0;

	p = (x - xl)/(xh - xl);
	p = (2*p - 1 );
	if(p <= -1.0)
		return (0.0) ;
	else if(p > -1.0 && p < 1.0)
		return (0.5 + 0.75*p*(1.0 - p*p/3.0));
	else
		return (1.0);

}

/* these routines evaluate the response for the different user filters */
void gsac_evalpz(float freq,float *fr,float *fi)
{
	/* lifted from sacfilt.f in VOLV/src */
	int np, i;
	double zr, zi;
	double sr, si;
	double tr, ti;
	double tmpr, tmpi;
	double fac;

	sr = 0.0 ;
	si = 6.2831853*freq ;
	if(npoles < nzeros)
		np = nzeros;
	else
		np = npoles;
	zr = 1.0;
	zi = 0.0;
	zr *= pzconst;
	for(i=0 ; i < np ; i++){
		if(i < nzeros){
			tr = sr - zreal[i];
			ti = si - zimag[i];
			tmpr = zr * tr - zi * ti;
			tmpi = zr * ti + zi * tr;
			zr = tmpr ;
			zi = tmpi;
		}
		if(i < npoles){
			tr = sr - preal[i];
			ti = si - pimag[i];
			fac = tr * tr + ti * ti ;
			tmpr =  zr * tr + zi * ti ;
			tmpi = -zr * ti + zi * tr ;
			zr = tmpr/fac ;
			zi = tmpi/fac;
		}
	}
	*fr = zr ; 
	*fi = zi ;
}

void gsac_evaleval(float freq,float *fr,float *fi)
{
	float degrad ;
	int jind;
	float p;
	float amp, pha;
	degrad = 3.1415927/180.0;
	/* for safety we process the amplitude and phase separately */
	*fr = 1.0;
	*fi = 0.0;
	arr_locate(evfamp, evnamp, freq, &jind);
	if(jind < 0) {
		amp = evaamp[0];
	} else if(jind >= 0 && jind < evnamp-1){
		p = (freq - evfamp[jind])/(evfamp[jind+1] - evfamp[jind]);
		amp = (1.0 -p)*evaamp[jind] + p*evaamp[jind+1];
	} else {
		amp = evaamp[evnamp-1];
	}
	arr_locate(evfpha, evnpha, freq, &jind);
	if(jind < 0) {
		pha = evapha[0];
	} else if(jind >= 0 && jind < evnpha-1){
		p = (freq - evfpha[jind])/(evfpha[jind+1] - evfpha[jind]);
		pha = (1.0 -p)*evapha[jind] + p*evapha[jind+1];
	} else {
		pha = evapha[evnpha-1];
	}


	*fr = amp*cos(pha*degrad);
	*fi = amp*sin(pha*degrad);
}

void gsac_evalfap(float freq,float *fr,float *fi)
{

	float degrad ;
	int jind;
	float p;
	float amp, pha;
	degrad = 3.1415927/180.0;
	/* for safety we process the amplitude and phase separately */
	*fr = 1.0;
	*fi = 0.0;
	arr_locate(fapfrq, fapnum, freq, &jind);
	if(jind < 0) {
		amp = fapamp[0];
		pha = fappha[0];
	} else if(jind >= 0 && jind < fapnum-1){
		p = (freq - fapfrq[jind])/(fapfrq[jind+1] - fapfrq[jind]);
		amp = (1.0 -p)*fapamp[jind] + p*fapamp[jind+1];
		pha = (1.0 -p)*fappha[jind] + p*fappha[jind+1];
	} else {
		amp = fapamp[fapnum-1];
		pha = fappha[fapnum-1];
	}


	*fr = amp*cos(pha*degrad);
	*fi = amp*sin(pha*degrad);
}

/* these routines read in the filter parameters */
int gsac_set_pz(char *pzname)
{
	FILE *fin;
	int i;
	int ip, iz;  /* index for pole or zero */
	char *v;
	char *p;
	char s1[100], s2[100];
	double v1, v2;
	int nread;
	int porz; /* -1 = pole, 0 = neither, 1 = zero */
	if((fin = fopen(pzname, "r")) ==NULL)
		return NO;	
	porz = 0;
	ip = 0;
	iz = 0;
	while((nread=sfgetline(fin, input_lineptr, NCHAR)) != -1 ){
		p = &input_lineptr[0];
		/* get rid of initial blanks */
/* HACK 17 JAN 2007 since old Sparcs did not havw isblank
		while (*p && isblank(*p) )
*/
		while (*p && isspace(*p) )
			p++;
		if(strncmp(p,"CON",3) == 0 ||
			strncmp(p,"con",3) == 0 ){
				if(sscanf(p,"%s %s",s1,s2) ==2 ){
				pzconst = strtod(s2,&v);
				/*
					printf("CON s1 %s s2 %s\n",s1,s2);
					*/
				if(*v != '\0'){
					printf("Error in string %s of input %s %s\n",s2,s1,s2);
					return NO;
					}
				}
			}
		else if(strncmp(p,"P",1) == 0 ||
			strncmp(p,"p",1) == 0 ){
				porz = -1 ;
				if(sscanf(p,"%s %s",s1,s2)==2){
					printf("P s1 %s s2 %s\n",s1,s2);
				/*
					*/
				npoles = atoi(s2);
				/* initialize */
				if(npoles >= NPAZ){
					printf("Number of poles in %s exceeds internal limit of %d\n",pzname,NPAZ);
					return NO ;
				}
				for(i=0 ; i < npoles ; i++){
					preal[i] = 0.0;
					pimag[i] = 0.0;
				}
				}
			}
		else if(strncmp(p,"Z",1) == 0 ||
			strncmp(p,"z",1) == 0 ){
				porz = 1;
				if(sscanf(p,"%s %s",s1,s2) == 2){
					printf("Z s1 %s s2 %s\n",s1,s2);
					/*
					*/
				nzeros = atoi(s2);
				if(nzeros >= NPAZ){
					printf("Number of zeros in %s exceeds internal limit of %d\n",pzname,NPAZ);
					return NO ;
				}
				/* initialize */
				for(i=0 ; i < nzeros ; i++){
					zreal[i] = 0.0;
					zimag[i] = 0.0;
				}
				}
			}
		else if(strncmp(p,"*",1) == 0 ){
			}
		else if(strncmp(p,"#",1) == 0 ){
			}
		else {
				if(sscanf(p,"%s %s",s1,s2) == 2){
					/*
					printf("D s1 %s s2 %s\n",s1,s2);
					*/
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
				if(porz == -1){
					preal[ip] = v1;
					pimag[ip] = v2;
					ip++;
				} else if(porz == 1 ){
					zreal[iz] = v1;
					zimag[iz] = v2;
					iz++;
				}
				}
		}
	}
	fclose (fin);
	printf("POLES %d\n",npoles);
		for(i=0;i< npoles;i++)
			printf("   %f %f\n",preal[i],pimag[i]);
	printf("ZEROS %d\n",nzeros);
		for(i=0;i< nzeros;i++)
			printf("   %f %f\n",zreal[i],zimag[i]);
	printf("CONSTANT %g\n",pzconst);
/* DEBUG 
 end of debug */
	return YES;
}

int gsac_set_eval(char *ampname, char *phaname)
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
/* HACK 17 JAN 2007 since old Sparcs did not havw isblank
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

	if((fin = fopen(phaname, "r")) ==NULL)
		return NO;	
	evnpha = 0;
	while(sfgetline(fin, p, NCHAR) != -1 ){
		p = &input_lineptr[0];
		/* get rid of initial blanks */
/* HACK 17 JAN 2007 since old Sparcs did not havw isblank
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
		evfpha[evnpha] = v1;
		evapha[evnpha] = v2;
		evnpha++;
		if(evnamp == NEVAL){
			printf("Number of lines in %s exceeds limit of %d\n",
					phaname,NEVAL);
			return NO;
		}
	}
	fclose(fin);
	return YES;
}

int gsac_set_fap(char *fapname)
{
/*
#
 # Displacement response for Array
 #
 # Example:  ST01 z
 #
 # Geotech 23900 seismometer
 #
 # Phase unwrapped
 #
  theoretical  0   instrument    fap Organization
 40
 0.100000   1.576582e-01   -52.923801  0.000000  0.000000
 0.125990   3.511520e-01   -61.669102  0.000000  0.000000
 0.200000   1.634426e+00   -79.966599  0.000000  0.000000
 0.368400   1.171214e+01  -107.522003  0.000000  0.000000
 0.500000   3.135000e+01  -126.447998  0.000000  0.000000
 0.683990   8.322500e+01  -155.035004  0.000000  0.000000
 0.800000   1.273452e+02  -174.207001  0.000000  0.000000
*/
/* steps
	1. skip leading blanks
	2. if then see # this is a comment
	3. If see theoretical then this is a key word
	4. If see one entry then that is the number
		skip this
	Eventually key to see one entry (e.g., parse into words)
		and it is a number
	5. read to end of file
*/



	FILE *fin;
	char *v;
	char *p;
	char s1[100], s2[100], s3[100];
	double v1, v2,v3;
	if((fin = fopen(fapname, "r")) ==NULL)
		return NO;	
	fapnum = 0;
	while(sfgetline(fin, input_lineptr, NCHAR) != -1 ){ 
		p = &input_lineptr[0];
		/* get rid of initial blanks */
/* HACK 17 JAN 2007 since old Sparcs did not havw isblank
		while (*p && isblank(*p) )
*/
		while (*p && isspace(*p) )
			p++;
		/* look for a # */
		if(strchr(p,'#')==NULL) {
		if(strchr(p, 't') == NULL &&
			strchr(p, 'T') == NULL){
			/* skip the next line */
			sfgetline(fin, p, NCHAR);
		sscanf(p,"%s %s %s",s1,s2,s3);
		v1 = strtod(s1,&v);
		if(*v != '\0'){
			printf("Error in string %s of input %s %s %s\n",s1,s1,s2,s3);
			return NO;
			}
		v2 = strtod(s2,&v);
		if(*v != '\0'){
			printf("Error in string %s of input %s %s %s\n",s2,s1,s2,s3);
			return NO;
			}
		v3 = strtod(s3,&v);
		if(*v != '\0'){
			printf("Error in string %s of input %s %s %s\n",s3,s1,s2,s3);
			return NO;
			}
		fapfrq[fapnum] = v1;
		fapamp[fapnum] = v2;
		fappha[fapnum] = v3;
		fapnum++;
		if(fapnum == NFAP){
			printf("Number of lines in %s exceeds limit of %d\n",
					fapname,NFAP);
			return NO;
		}
	}
	}
	}
	fclose(fin);

	return YES;
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

int gstr(int ncmd, char **cmdstr, char *pat, int n)
{
	int i;
	for(i=1;i<ncmd;i++){
		strcpy(tmpstr,cmdstr[i]);
		gsac_strupr(tmpstr);
		if(strncmp(tmpstr,pat,n) == 0)
			return i;
	}
	return -1;
}
