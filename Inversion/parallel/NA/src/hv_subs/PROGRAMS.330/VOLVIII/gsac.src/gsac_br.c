/* Bug Fixes
 * 11 JAN 2010 - added Chebyshev Type I filter
*/
#include "gsac_docommand.h"
#include	<stdio.h>
#include "gsac.h"
#include "gsac_arg.h"


#define BR_BUTTER	0
#define BR_CORNER	1
#define BR_NPOLES	2
#define BR_PASSES	3
#define BR_BESSEL	5
#define BR_C1           6
#define BR_TRANBW      7
#define BR_ATTEN        8

#define BANDREJECT 4 
#define BUTTERWORTH 0
#define BESSEL 1
#define CHEBYSHEVI 2
static char *brstr[] = {"Butterworth", "Bessel", "Chebyshev-I"};

struct arghdr brarg[] = {
	{BR_BESSEL, "BESSEL"	, IHDR, 0, 0, NO, "BESSEL", 2},
	{BR_BUTTER, "BUTTER"	, IHDR, 0, 0, NO, "BUTTER", 2},
	{BR_BUTTER, "B"		, IHDR, 0, 0, NO, "B ",-1},
	{BR_C1    , "C1"	, IHDR, 0, 0, NO, "C1 ",2},
	{BR_CORNER, "CORNER"	, RHDR, 0, 2, NO, "CORNER fl fh", 2},
	{BR_CORNER, "C"	        , RHDR, 0, 2, NO, "CORNER fl fh", -1},
	{BR_NPOLES, "NPOLES"	, IHDR, 0, 1, NO, "NPOLES np", 1},
	{BR_PASSES, "PASSES"	, IHDR, 0, 1, NO, "PASSES num - 1 or 2", 1},
        {BR_TRANBW, "TRANBW"    , RHDR, 0, 1, NO, "TRANBW tranbw", 1},
        {BR_ATTEN, "ATTEN"      , RHDR, 0, 1, NO, "ATTEN atten", 1},
	{0,	""		, IHDR, 0, 0, NO, "",-1}
};

/* these are temporary variables only used here */
float br_real[10];
int   br_int [10];

static float brfilt_fl = 0.0;
static float brfilt_fh = 1.0e+10;
static int brfilt_np = 1;
static int brfilt_p = 1;
static int brfilt_lhp = BANDREJECT;
static int br_filt_type = BUTTERWORTH;
static float br_cheb_eps;
static float br_cheb_tranbw = 0.3;
static float br_cheb_atten = 30.;


void gsac_set_param_br(int ncmd, char **cmdstr)
{
	int i;
	float tmp;
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, brarg, NO, YES))
	       	return	;
	for(i=0 ; brarg[i].key[0] != '\0' ; i++){
		if(brarg[i].used > 0){
			if(brarg[i].ricell == RHDR)
				getargr(ncmd, cmdstr, brarg[i].key, brarg[i].mfit, brarg[i].narg, br_real);
			else if(brarg[i].ricell == IHDR)
				getargi(ncmd, cmdstr, brarg[i].key, brarg[i].mfit, brarg[i].narg, br_int );
			switch(brarg[i].id){
				case BR_BESSEL:
					br_filt_type = BESSEL;
					break;
				case BR_BUTTER:
					br_filt_type = BUTTERWORTH;
					break;
				case BR_C1:
					br_filt_type = CHEBYSHEVI;
					break;
				case BR_NPOLES:
					if(br_int[0] < 1)
						brfilt_np = 1;
					else if(br_int[0] > 10)
						brfilt_np = 10;
					else
						brfilt_np = br_int[0];
					break;
				case BR_PASSES:
					if(br_int[0] < 1)
						brfilt_p = 1;
					else if(br_int[0] > 2)
						brfilt_p = 2;
					else
						brfilt_p = br_int[0];
					break;
				case BR_CORNER:
					if(br_real[0] > br_real[1]){
						tmp = br_real[0];
						br_real[0] = br_real[1];
						br_real[1] = tmp;
					}
					if(br_real[0] < 0.0)
						brfilt_fl = 0.0;
					brfilt_fl = br_real[0];
					brfilt_fh = br_real[1];
					break;
                                case BR_ATTEN:
                                        if(br_real[0] < 1.0)
                                                br_cheb_atten = 1.0;
                                        br_cheb_atten = br_real[0];
                                        break;
                                case BR_TRANBW:
                                        if(br_real[0] < 0.1)
                                                br_cheb_tranbw = 0.1;
                                        br_cheb_tranbw = br_real[0];
                                        break;


			}
		}
	}
	brfilt_lhp = BANDREJECT;
}

void gsac_exec_br(void)
{
float rat, fac, x, ripple;
	if(brfilt_fl  == brfilt_fh){
		printf("Filter corners must be different: %f %f\n",brfilt_fl, brfilt_fh);
		return;
	}

if(br_filt_type == CHEBYSHEVI){
	rat = 1. + br_cheb_tranbw;
	fac = (rat + sqrt(rat*rat -1 ));
	fac = pow(fac,brfilt_np);
	x = (fac*fac +1.0)/(2.0*fac);
	br_cheb_eps = sqrt(br_cheb_atten*br_cheb_atten -1.0)/x;
	ripple = 1./sqrt(1.0+br_cheb_eps*br_cheb_eps);
printf("BR: corners fl %f fh %f npoles %d pass %d tranbw %f atten %f eps %f %s\n",brfilt_fl, brfilt_fh, brfilt_np, brfilt_p, br_cheb_tranbw, br_cheb_atten, br_cheb_eps,brstr[br_filt_type]);
} else {
printf("BR: corners fl %f fh %f npoles %d pass %d %s\n",brfilt_fl, brfilt_fh ,brfilt_np, brfilt_p, brstr[br_filt_type]);
}

gsac_filt(brfilt_fl, brfilt_fh, brfilt_np, brfilt_p, brfilt_lhp, br_filt_type, br_cheb_eps);

	if(gsac_control.prs > 0){
		if(gsac_control.prshist == NULL)
			gsac_control.prshist = fopen("prshist.tmp","w+");
		fprintf(gsac_control.prshist,"br c fl %f fh %f n %d p %d %s\n",
		brfilt_fl, brfilt_fh, brfilt_np, brfilt_p,brstr[br_filt_type]);
		fflush(gsac_control.prshist);
	}
}
