/* Bug Fixes
 * 11 JAN 2010 - added Chebyshev Type I filter
*/
#include        <stdio.h>
#include        "gsac.h"
#include        "gsac_sac.h"
#include "gsac_docommand.h"
#include "gsac_arg.h"
#include "gsac_pz.h"


#define LP_BUTTER	0
#define LP_CORNER	1
#define LP_NPOLES	2
#define LP_PASSES	3
#define LP_BESSEL       5
#define LP_C1           6
#define LP_TRANBW      7
#define LP_ATTEN        8

#define LOWPASS 3
#define BUTTERWORTH 0
#define BESSEL 1
#define CHEBYSHEVI 2
static char *lpstr[] = {"Butterworth", "Bessel", "Chebyshev-I"};


struct arghdr lparg[] = {
	{LP_BESSEL, "BESSEL"	, IHDR, 0, 0, NO, "BESSEL", 2},
	{LP_BUTTER, "BUTTER"	, IHDR, 0, 0, NO, "BUTTER", 2},
	{LP_BUTTER, "B"		, IHDR, 0, 0, NO, "B ",-1},
	{LP_C1    , "C1"	, IHDR, 0, 0, NO, "C1",2},
	{LP_CORNER, "CORNER"	, RHDR, 0, 1, NO, "CORNER  fh", 2},
	{LP_CORNER, "C"	        , RHDR, 0, 1, NO, "CORNER  fh", -1},
	{LP_NPOLES, "NPOLES"	, IHDR, 0, 1, NO, "NPOLES np", 1},
	{LP_PASSES, "PASSES"	, IHDR, 0, 1, NO, "PASSES num - 1 or 2", 1},
	{LP_TRANBW, "TRANBW"	, RHDR, 0, 1, NO, "TRANBW tranbw", 1},
	{LP_ATTEN, "ATTEN"	, RHDR, 0, 1, NO, "ATTEN atten", 1},
	{0,	""		, IHDR, 0, 0, NO, "",-1}
};


/* these are temporary variables only used here */
float lp_real[10];
int   lp_int [10];

static float lpfilt_fl = 0.0;
static float lpfilt_fh = 1.0e+10;
static int lpfilt_np = 1;
static int lpfilt_p = 1;
static int lpfilt_lhp = LOWPASS;
static int lp_filt_type = BUTTERWORTH;
static float lp_cheb_eps;
static float lp_cheb_tranbw = 0.3;
static float lp_cheb_atten = 30.;


void gsac_set_param_lp(int ncmd, char **cmdstr)
{
	int i;

	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, lparg, NO, YES))
	       	return	;
	/* systematically go through the arguments to see if they have
	 * been correctly invoked. If the have, then get the values */
	for(i=0 ; lparg[i].key[0] != '\0' ; i++){
		if(lparg[i].used > 0){
			if(lparg[i].ricell == RHDR)
				getargr(ncmd, cmdstr, lparg[i].key, 
					lparg[i].mfit, lparg[i].narg, lp_real);
			else if(lparg[i].ricell == IHDR)
				getargi(ncmd, cmdstr, lparg[i].key, 
					lparg[i].mfit, lparg[i].narg, lp_int );
			switch(lparg[i].id){
				case LP_BESSEL:
					lp_filt_type = BESSEL;
					break;
				case LP_BUTTER:
					lp_filt_type = BUTTERWORTH;
					break;
				case LP_C1:
					lp_filt_type = CHEBYSHEVI;
					break;
				case LP_NPOLES:
					if(lp_int[0] < 1)
						lpfilt_np = 1;
					else if(lp_int[0] > 10)
						lpfilt_np = 10;
					else
						lpfilt_np = lp_int[0];
					break;
				case LP_PASSES:
					if(lp_int[0] < 1)
						lpfilt_p = 1;
					else if(lp_int[0] > 2)
						lpfilt_p = 2;
					else
						lpfilt_p = lp_int[0];
					break;
				case LP_CORNER:
					if(lp_real[0] < 0.0)
						lpfilt_fl = 10000.0;
					lpfilt_fh = lp_real[0];
					break;
				case LP_ATTEN:
					if(lp_real[0] < 1.0)
						lp_cheb_atten = 1.0;
					lp_cheb_atten = lp_real[0];
					break;
				case LP_TRANBW:
					if(lp_real[0] < 0.1)
						lp_cheb_tranbw = 0.1;
					lp_cheb_tranbw = lp_real[0];
					break;

			}
		}
	}
	lpfilt_lhp = LOWPASS;
}

void gsac_exec_lp(void)
{
float rat, fac, x, ripple;

if(lp_filt_type == CHEBYSHEVI){
	rat = 1. + lp_cheb_tranbw;
	fac = (rat + sqrt(rat*rat -1 ));
	fac = pow(fac,lpfilt_np);
	x = (fac*fac +1.0)/(2.0*fac);
	lp_cheb_eps = sqrt(lp_cheb_atten*lp_cheb_atten -1.0)/x;
	ripple = 1./sqrt(1.0+lp_cheb_eps*lp_cheb_eps);
printf("LP: corner fc %f  npoles %d pass %d tranbw %f atten %f eps %f %s\n",lpfilt_fh ,lpfilt_np, lpfilt_p, lp_cheb_tranbw, lp_cheb_atten, lp_cheb_eps, lpstr[lp_filt_type]);
} else {
printf("LP: corner fc %f  npoles %d pass %d %s\n",lpfilt_fh ,lpfilt_np, lpfilt_p, lpstr[lp_filt_type]);
}

gsac_filt(0.0, lpfilt_fh, lpfilt_np, lpfilt_p, lpfilt_lhp, lp_filt_type, lp_cheb_eps);

	if(gsac_control.prs > 0){
		if(gsac_control.prshist == NULL)
			gsac_control.prshist = fopen("prshist.tmp","w+");
		fprintf(gsac_control.prshist,"lp c fc %f n %d p %d %s\n",
			lpfilt_fh ,lpfilt_np, lpfilt_p,lpstr[lp_filt_type]);
		fflush(gsac_control.prshist);
	}

}
