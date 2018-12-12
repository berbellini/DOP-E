/* Bug Fixes
 * 11 JAN 2010 - added Chebyshev Type I filter
*/
#include        <stdio.h>
#include        "gsac.h"
#include        "gsac_sac.h"
#include	"gsac_docommand.h"
#include	"gsac_arg.h"

extern struct sacfile_ *sacdata;

#define HP_BUTTER	0
#define HP_CORNER	1
#define HP_NPOLES	2
#define HP_PASSES	3
#define HP_BESSEL	5
#define HP_C1           6
#define HP_TRANBW      7
#define HP_ATTEN        8

#define HIGHPASS 2
#define BUTTERWORTH 0
#define BESSEL 1
#define CHEBYSHEVI 2
static char *hpstr[] = {"Butterworth", "Bessel", "Chebyshev-I"};

struct arghdr hparg[] = {
	{HP_BESSEL, "BESSEL"	, IHDR, 0, 0, NO, "BESSEL", 2},
	{HP_BUTTER, "BUTTER"	, IHDR, 0, 0, NO, "BUTTER", 2},
	{HP_BUTTER, "B"		, IHDR, 0, 0, NO, "B ",-1},
	{HP_C1    , "C1"	, IHDR, 0, 0, NO, "C1",2},
	{HP_CORNER, "CORNER"	, RHDR, 0, 1, NO, "CORNER  fl", 2},
	{HP_CORNER, "C"	        , RHDR, 0, 1, NO, "CORNER  fl", -1},
	{HP_NPOLES, "NPOLES"	, IHDR, 0, 1, NO, "NPOLES np", 1},
	{HP_PASSES, "PASSES"	, IHDR, 0, 1, NO, "PASSES num - 1 or 2", 1},
        {HP_TRANBW, "TRANBW"    , RHDR, 0, 1, NO, "TRANBW tranbw", 1},
        {HP_ATTEN, "ATTEN"      , RHDR, 0, 1, NO, "ATTEN atten", 1},

	{0,	""		, IHDR, 0, 0, NO, "",-1}
};


/* these are temporary variables only used here */
float hp_real[10];
int   hp_int [10];

static float hpfilt_fl = 0.0;
static int hpfilt_np = 1;
static int hpfilt_p = 1;
static int hpfilt_lhp = HIGHPASS;
static int hp_filt_type = BUTTERWORTH;
static float hp_cheb_eps;
static float hp_cheb_tranbw = 0.3;
static float hp_cheb_atten = 30.;


void gsac_set_param_hp(int ncmd, char **cmdstr)
{
	int i;
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, hparg, NO, YES))
	       return	;
	for(i=0 ; hparg[i].key[0] != '\0' ; i++){
		if(hparg[i].used > 0){
			if(hparg[i].ricell == RHDR)
				getargr(ncmd, cmdstr, hparg[i].key, 
					hparg[i].mfit, hparg[i].narg, hp_real);
			else if(hparg[i].ricell == IHDR)
				getargi(ncmd, cmdstr, hparg[i].key, 
					hparg[i].mfit, hparg[i].narg, hp_int );
			switch(hparg[i].id){
				case HP_BESSEL:
					hp_filt_type = BESSEL;
					break;
				case HP_BUTTER:
					hp_filt_type = BUTTERWORTH;
					break;
				case HP_C1:
					hp_filt_type = CHEBYSHEVI;
					break;
				case HP_NPOLES:
					if(hp_int[0] < 1)
						hpfilt_np = 1;
					else if(hp_int[0] > 10)
						hpfilt_np = 10;
					else
						hpfilt_np = hp_int[0];
					break;
				case HP_PASSES:
					if(hp_int[0] < 1)
						hpfilt_p = 1;
					else if(hp_int[0] > 2)
						hpfilt_p = 2;
					else
						hpfilt_p = hp_int[0];
					break;
				case HP_CORNER:
					if(hp_real[0] < 0.0)
						hpfilt_fl = 0.0;
					hpfilt_fl = hp_real[0];
					break;
                                case HP_ATTEN:
                                        if(hp_real[0] < 1.0)
                                                hp_cheb_atten = 1.0;
                                        hp_cheb_atten = hp_real[0];
                                        break;
                                case HP_TRANBW:
                                        if(hp_real[0] < 0.1)
                                                hp_cheb_tranbw = 0.1;
                                        hp_cheb_tranbw = hp_real[0];
                                        break;

			}
		}
	}
	hpfilt_lhp = HIGHPASS;
}

void gsac_exec_hp(void)
{
float rat, fac, x, ripple;

if(hp_filt_type == CHEBYSHEVI){
	rat = 1. + hp_cheb_tranbw;
	fac = (rat + sqrt(rat*rat -1 ));
	fac = pow(fac,hpfilt_np);
	x = (fac*fac +1.0)/(2.0*fac);
	hp_cheb_eps = sqrt(hp_cheb_atten*hp_cheb_atten -1.0)/x;
	ripple = 1./sqrt(1.0+hp_cheb_eps*hp_cheb_eps);
printf("HP: corner fc %f  npoles %d pass %d tranbw %f atten %f eps %f %s\n",hpfilt_fl ,hpfilt_np,  hpfilt_p,hp_cheb_tranbw, hp_cheb_atten, hp_cheb_eps,hpstr[hp_filt_type]);
} else {
printf("HP: corner fc %f  npoles %d pass %d %s\n",hpfilt_fl ,hpfilt_np, hpfilt_p,hpstr[hp_filt_type]);
}

gsac_filt(hpfilt_fl, 1.0e+10, hpfilt_np, hpfilt_p, hpfilt_lhp, hp_filt_type, hp_cheb_eps);

	if(gsac_control.prs > 0){
		if(gsac_control.prshist == NULL)
			gsac_control.prshist = fopen("prshist.tmp","w+");
		fprintf(gsac_control.prshist,"hp c fc %f  n %d p %d %s\n",
			hpfilt_fl, hpfilt_np, hpfilt_p,hpstr[hp_filt_type]);
		fflush(gsac_control.prshist);
	}
}
