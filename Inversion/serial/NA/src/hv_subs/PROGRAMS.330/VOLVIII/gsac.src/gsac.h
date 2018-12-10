/* Changes:
	13 JUN 2010: Added defines for IFTYPE -> ITIME == 1
						 IXY   == 4
                                       IZTYPE -> IB    == 9
						 IO    == 11
                                                 IA    == 12
						 ITn (n=0,9) 13+n
				       LEVEN  -> TRUE == 1
	11 JUL 2010 - add support for the ~ screendump
*/

#ifndef _GSAC_H
#define _GSAC_H
#include        <math.h>
#include        <stdio.h>

#ifdef MAIN
#define EXTERN
#else
#define EXTERN extern
#endif  /* MAIN */


/* headers and master control statements */



#define	PLT 1
#define WIN 0
#define YES 1
#define NO  0

/* define shorthand for accessing SAC header variables so that
 * code is slightly easier to read
 * */
/* real header */
#define H_DELTA   0
#define H_DEPMIN  1
#define H_DEPMAX  2
#define H_B       5
#define H_E       6
#define H_O       7
#define H_A       8
#define H_T0     10
#define H_T1     11
#define H_T2     12
#define H_T3     13
#define H_T4     14
#define H_T5     15
#define H_T6     16
#define H_T7     17
#define H_T8     18
#define H_T9     19
#define H_F      20
#define H_STLA   31
#define H_STLO   32
#define H_STEL   33
#define H_EVLA   35
#define H_EVLO   36
#define H_EVDP   38
#define H_MAG     39
#define H_USER0   40
#define H_USER1   41
#define H_USER2   42
#define H_USER3   43
#define H_USER4   44
#define H_USER5   45
#define H_USER6   46
#define H_USER7   47
#define H_USER8   48
#define H_USER9   49
#define H_DIST   50
#define H_AZ     51
#define H_BAZ    52
#define H_GCARC  53
#define H_DEPMEN 56
#define H_CMPAZ 57
#define H_CMPINC 58
#define H_TIMMAX 64
#define H_TIMMIN 65

/* integer header */
#define H_NZYEAR 0
#define H_NZJDAY 1
#define H_NZHOUR 2
#define H_NZMIN  3
#define H_NZSEC  4
#define H_NZMSEC 5
#define H_NVHDR  6
#define H_NPTS 9
#define H_IFTYPE  15
#define H_IZTYPE  17
#define H_IHDR11  25
#define H_IHDR20  34
#define H_LPSPOL  36
#define H_LOVROK  37
#define H_LCALDA  38
#define H_LHDR5  39

/* character header */
#define H_KSTNM   0
#define H_KO       4
#define H_KA      5
#define H_KT0     6
#define H_KT1     7
#define H_KT2     8
#define H_KT3     9
#define H_KT4     10
#define H_KT5     11
#define H_KT6     12
#define H_KT7     13
#define H_KT8     14
#define H_KT9     15
#define H_KF      16
#define H_KCMPNM  20

/* enmerated header values */
#define ENUM_ITIME  1
#define ENUM_IXY    4
#define ENUM_IB     9
#define ENUM_IO    11
#define ENUM_IA    12
#define ENUM_IT0   13
#define ENUM_IT1   14
#define ENUM_IT2   15
#define ENUM_IT3   16
#define ENUM_IT4   17
#define ENUM_IT5   18
#define ENUM_IT6   19
#define ENUM_IT7   20
#define ENUM_IT8   21
#define ENUM_IT9   22

#define MIN(a,b) ( (b) > (a) ? (a):(b) )
#define MAX(a,b) ( (b) < (a) ? (a):(b) )
#define ABS(a  ) ( (a) >  0  ? (a):-(a))
#define SIGN(a ) ( (a) >  0  ? (1):(-1))

struct sac_control_ {
	int number_itraces;	/* current number of traces in memory */
	int number_iheaders;	/* current number of headers in memory */
	int number_otraces;	/* current number of traces for write */
	int max_number_traces;	/* maximum number of traces every read */
	int more;
	int plotdevice;
	int plotinit;
	int plotdevicechange;
	int plotchange;
	int plotcount;
	int plotcount_prs;
	int plotcount_refr;
	int plotcount_dump;
	int everinteractive;
	double begmin;		/* minimum epoch time for all data */
	double endmax;		/* maximum epoch time for all data */
	float fmax;		/* maximum Nyquist frequency */
	int docut;		/* cut YES or NO */
	int doxlim;		/* xlim YES or NO */
	int cutint[2];		/* pointer to header structure for
				   reference to cut value */
	float cutoff[2];	/* offset with respect to cut value */
	double cutepoch[2];	/* for use with CAL or GMT cut */
	char cutkey[2][20];
	int xlimint[2];		/* pointer to header structure for
				   reference to xlim value */
	float xlimoff[2];	/* offset with respect to xlim value */
	char xlimkey[2][3];
	int qdp;		/* QDP decimation factor only for WINDOW plot 
					-1 QDP off
					 0 automatic determination QDP on
					>0 decimate value is specified by user
				*/
	int local;		/* 0 PRESERVE original byte order.
				   1 Convert to LOCAL byte order */
	float x0;		/* (x0,y0 ) is the lower left corner of the */
	float y0;		/* psp, prs and p1 plots */
	float xlen;		/* xlen and ylen are the axis lengths */
	float ylen;
	int grid;		/* positioning grid switch */
	float fft;		/* FFT done ? - reset on each read */
	int hold;		/* if >0 continue plot on current frame */
	int inpltmode;		/* used for hold synchronization for PLT */
	double begminx;		/* minimum epoch time for all data for XLIM */
	double endmaxx;		/* maximum epoch time for all data for XLIM */
	int kolor;		/* device color capability */
	int black;		/* device definition of black */
	int xvigenv;	/* has PLOTXVIG environment been set */
	float uxcen;
	float uxmul;
	float ylim_low;
	float ylim_high;
	float ylim_ctrl;
	float ylim_rsc;
	int plotlinx;	/* for plot and plotpk YES=linear, NO=LOG */
	int plotliny;
	float XmaxDev;	/* for window resolution */
	float YmaxDev;	/* for window resolution */
	int prs;	/* 0 = prs, 1 = prs-refraction, 2= prs-reflection */
	FILE *prshist;	/* file pointer for prs history - file starts
		afresh with each read */
	FILE *refrpick;	/* file pointer for prs refraction picks */
	int xgrid;	/* 0 off 1 on */
	int xgrid_type;	/* 1 solid, 2 dotted */
	int xgrid_color;	/* CALPLOT color */
	int xgrid_minor; /* YES or NO */
	int ygrid;	/* 0 off 1 on */
	int ygrid_type;	/* 1 solid, 2 dotted */
	int ygrid_color;	/* CALPLOT color */
	int ygrid_minor; /* YES or NO */
	int background ; /* YES or NO */
	int background_color ; /* back ground color if >= 0 */
};
EXTERN struct sac_control_ gsac_control;


/* function prototypes */
int gsac_parse_command(void);
void gsac_init(void);
char *gsac_strupr(char *s);
int parsecommand(char *str);
void gsac_filt(float fl, float fh, int np, int p, int lhp, int type, float cheb_eps);
void gsac_pulconv(float t1, float t2, float t3);
void getmxmn(float *x, int npts, float *depmax, float *depmin, float *depmen, int *indmax, int *indmin);
void delaz(float elat, float elon, float slat, float slon, float *deldeg, float *az, float *baz, float *delkm);
void gsac_setcolor(int onoff, int k, int ntrc);
void gsac_alloc_trace(int oldmax);
void fillstr1(char *instr, char *str1);
void fillstr12(char *instr, char *str1, char *str2);
void chcmpnm(char *kstnm, char *kcmpnm, char *cmpstr, char *kfilename, int dosuffix, char *suffix);
int findblank(char *str);
void chofname(char *kstnm1, char *kcmpnm1, char *kstnm2, char *kcmpnm2, char *kfilename, int dosuffix, char *suffix);
int npow2(int n);
void four(float data[], int nn, int isign, float *dt, float *df);
int gsac_countgmt(int ncmd, char **cmdstr, char *pat);
int gsac_findgmt(int ncmd, char **cmdstr, char *pat);
void gsac_mergegmt(int *ncmd, char **cmdstr, char *pat);
void cleanstring(char *str);






#endif
