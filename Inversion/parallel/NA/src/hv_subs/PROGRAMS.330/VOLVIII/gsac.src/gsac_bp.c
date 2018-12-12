/* Bug Fixes
 * 11 JAN 2010 - added Chebyshev Type I filter
*/
#include "gsac_docommand.h"
#include	<stdio.h>
#include "gsac.h"
#include "gsac_arg.h"


#define BP_BUTTER	0
#define BP_CORNER	1
#define BP_NPOLES	2
#define BP_PASSES	3
#define BP_BESSEL	5
#define BP_C1           6
#define BP_TRANBW      7
#define BP_ATTEN        8

#define BANDPASS 1
#define BUTTERWORTH 0
#define BESSEL 1
#define CHEBYSHEVI 2
static char *bpstr[] = {"Butterworth", "Bessel", "Chebyshev-I"};


struct arghdr bparg[] = {
	{BP_BESSEL, "BESSEL"	, IHDR, 0, 0, NO, "BESSEL", 2},
	{BP_BUTTER, "BUTTER"	, IHDR, 0, 0, NO, "BUTTER", 2},
	{BP_BUTTER, "B"		, IHDR, 0, 0, NO, "B ",-1},
	{BP_C1    , "C1"	, IHDR, 0, 0, NO, "C1 ",2},
	{BP_CORNER, "CORNER"	, RHDR, 0, 2, NO, "CORNER fl fh",2},
	{BP_CORNER, "C"	        , RHDR, 0, 2, NO, "CORNER fl fh",-1},
	{BP_NPOLES, "NPOLES"	, IHDR, 0, 1, NO, "NPOLES np", 1},
	{BP_PASSES, "PASSES"	, IHDR, 0, 1, NO, "PASSES num - 1 or 2",1},
        {BP_TRANBW, "TRANBW"    , RHDR, 0, 1, NO, "TRANBW tranbw", 1},
        {BP_ATTEN, "ATTEN"      , RHDR, 0, 1, NO, "ATTEN atten", 1},

	{0,	""		, IHDR, 0, 0, NO, "",-1}
};

/* these are temporary variables only used here */
float bp_real[10];
int   bp_int [10];

static float bpfilt_fl = 0.0;
static float bpfilt_fh = 1.0e+10;
static int bpfilt_np = 1;
static int bpfilt_p = 1;
static int bpfilt_lhp = BANDPASS;
static int bp_filt_type = BUTTERWORTH;
static float bp_cheb_eps;
static float bp_cheb_tranbw = 0.3;
static float bp_cheb_atten = 30.;


void gsac_set_param_bp(int ncmd, char **cmdstr)
{
	int i;
	float tmp;

	if(ncmd == 1)
		return;
	/* check for command syntax */
	if(testarg(ncmd, cmdstr, bparg, NO, YES))
	       	return	;
	for(i=0 ; bparg[i].key[0] != '\0' ; i++){
		if(bparg[i].used > 0){
			if(bparg[i].ricell == RHDR)
				getargr(ncmd, cmdstr, bparg[i].key, bparg[i].mfit, bparg[i].narg, bp_real);
			else if(bparg[i].ricell == IHDR)
				getargi(ncmd, cmdstr, bparg[i].key, bparg[i].mfit, bparg[i].narg, bp_int );
			switch(bparg[i].id){
				case BP_BESSEL:
					bp_filt_type = BESSEL;
					break;
				case BP_BUTTER:
					bp_filt_type = BUTTERWORTH;
					break;
				case BP_C1:
					bp_filt_type = CHEBYSHEVI;
					break;
				case BP_NPOLES:
					if(bp_int[0] < 1)
						bpfilt_np = 1;
					else if(bp_int[0] > 10)
						bpfilt_np = 10;
					else
						bpfilt_np = bp_int[0];
					break;
				case BP_PASSES:
					if(bp_int[0] < 1)
						bpfilt_p = 1;
					else if(bp_int[0] > 2)
						bpfilt_p = 2;
					else
						bpfilt_p = bp_int[0];
					break;
				case BP_CORNER:
					if(bp_real[0] > bp_real[1]){
						tmp = bp_real[0];
						bp_real[0] = bp_real[1];
						bp_real[1] = tmp;
					}
					if(bp_real[0] < 0.0)
						bpfilt_fl = 0.0;
					bpfilt_fl = bp_real[0];
					bpfilt_fh = bp_real[1];
					break;
                                case BP_ATTEN:
                                        if(bp_real[0] < 1.0)
                                                bp_cheb_atten = 1.0;
                                        bp_cheb_atten = bp_real[0];
                                        break;
                                case BP_TRANBW:
                                        if(bp_real[0] < 0.1)
                                                bp_cheb_tranbw = 0.1;
                                        bp_cheb_tranbw = bp_real[0];
                                        break;


			}
		}
	}
	bpfilt_lhp = BANDPASS;
}

void gsac_exec_bp(void)
{
float rat, fac, x, ripple;
	if(bpfilt_fl  == bpfilt_fh){
		printf("Filter corners must be different: %f %f\n",bpfilt_fl, bpfilt_fh);
		return;
	}
if(bp_filt_type == CHEBYSHEVI){
	rat = 1. + bp_cheb_tranbw;
	fac = (rat + sqrt(rat*rat -1 ));
	fac = pow(fac,bpfilt_np);
	x = (fac*fac +1.0)/(2.0*fac);
	bp_cheb_eps = sqrt(bp_cheb_atten*bp_cheb_atten -1.0)/x;
	ripple = 1./sqrt(1.0+bp_cheb_eps*bp_cheb_eps);
printf("BP: corners fl %f fh %f npoles %d pass %d tranbw %f atten %f eps %f %s\n",bpfilt_fl, bpfilt_fh, bpfilt_np,  bpfilt_p, bp_cheb_tranbw, bp_cheb_atten, bp_cheb_eps, bpstr[bp_filt_type]);
} else {
printf("BP: corners fl %f fh %f npoles %d pass %d %s\n",bpfilt_fl, bpfilt_fh, bpfilt_np, bpfilt_p,bpstr[bp_filt_type]);
}
gsac_filt(bpfilt_fl, bpfilt_fh, bpfilt_np, bpfilt_p, bpfilt_lhp, bp_filt_type, bp_cheb_eps);

	if(gsac_control.prs > 0){
		if(gsac_control.prshist == NULL)
			gsac_control.prshist = fopen("prshist.tmp","w+");
		fprintf(gsac_control.prshist,"bp c fl %f fh %f n %d p %d %s\n",
		bpfilt_fl, bpfilt_fh, bpfilt_np, bpfilt_p,bpstr[bp_filt_type]); 
		fflush(gsac_control.prshist);
	}
}
