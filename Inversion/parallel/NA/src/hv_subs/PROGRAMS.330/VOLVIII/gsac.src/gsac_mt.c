/* CHANGES
	22 MAY 2009 - created
	19 MAY 2010 - added N E and ZNE options
        24 MAY 2010 - found parsing error - gives wrong result

	 	to R R rake

		will give rake = 0 since the parsing actually goes by the pattern
		and not by position.
            
        17 JUN 2010 - if the TSS TDS do not exist this bombs
                      get around this since we may wish to generate
                      only Love or Rayleigh Green and still apply mech
        17 JUN 2010 - add explosion to the source
        08 JUN 2011 - fixed code so tha the t0 and t1 arrival times of SV and SH are 
                      correctly read and then saved in the headers
	13 OCT 2011 - now supports general moment tensor source
*/
#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	MT_DFLT	0
#define	MT_AZ	1
#define	MT_BAZ	2
#define	MT_STK	3
#define	MT_DIP	4
#define	MT_RAKE	5
#define	MT_MW	6
#define	MT_FILE	7
#define MT_TO   8
#define MT_ISO  9
#define MT_FN   10
#define MT_FE   11
#define MT_FD   12
#define MT_MXX   13
#define MT_MXY   14
#define MT_MXZ   15
#define MT_MYY   16
#define MT_MYZ   17
#define MT_MZZ   18
static float mt_az;
static float mt_baz;
static float mt_stk;
static float mt_dip;
static float mt_rake;
static float mt_mw;
static float mt_fx, mt_fy, mt_fz;
static float mt_mxx, mt_mxy, mt_mxz, mt_myy, mt_myz, mt_mzz;
static char *mt_proto;
static int mt_do_z;
static int mt_do_r;
static int mt_do_t;
static int mt_do_n;
static int mt_do_e;
static int mt_source;  /* 0 = DC, 1 = ISO ; 2 = Force */


struct arghdr mtarg[] = {
	{MT_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{MT_AZ  , "AZ"  , RHDR, NO, 1, NO, "Az az", 1},
	{MT_BAZ , "BAZ" , RHDR, NO, 1, NO, "Baz baz", 1},
	{MT_STK , "STK" , RHDR, NO, 1, NO, "STK stk",-1},
	{MT_DIP , "DIP" , RHDR, NO, 1, NO, "DIP dip",-1},
	{MT_RAKE, "RAKE", RHDR, NO, 1, NO, "RAKE rake",-1},
	{MT_MW  , "MW"  , RHDR, NO, 1, NO, "MW mw", -1},
	{MT_MXX , "MXX" , RHDR, NO, 1, NO, "MXX mxx", -1},
	{MT_MXY , "MXY" , RHDR, NO, 1, NO, "MXY mxy", -1},
	{MT_MXZ , "MXZ" , RHDR, NO, 1, NO, "MXZ mxz", -1},
	{MT_MYY , "MYY" , RHDR, NO, 1, NO, "MYY myy", -1},
	{MT_MYZ , "MYZ" , RHDR, NO, 1, NO, "MYZ myz", -1},
	{MT_MZZ , "MZZ" , RHDR, NO, 1, NO, "MZZ mzz", -1},
	{MT_ISO , "ISO" , RHDR, NO, 0, NO, "ISO ", -1},
	{MT_FN  , "FN"  , RHDR, NO, 1, NO, "FN fn ", -1},
	{MT_FE  , "FE"  , RHDR, NO, 1, NO, "FE fe ", -1},
	{MT_FD  , "FD"  , RHDR, NO, 1, NO, "FD fd ", -1},
	{MT_FILE, "FILE", CHDR, NO, 1, NO, "File fileproto ", -1},
	{MT_TO  , "TO"  , CHDR, NO, 1, YES, "TO [UZ|UR|UT|ZRT|UN|UE|ZNE]", 1},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float mt_real[10];
int   mt_int [10];
int   mt_yn;
int   mt_num;

static char instr[1000];
static char *tstr;

/* these are prototypes for global variables to be used by the routine */
static void gsac_mt_get(char *grn, int *tds,int *sum_t);

void gsac_set_param_mt(int ncmd, char **cmdstr)
{
	int i;
	int ls;
	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, mtarg, NO, YES))
		return;
	/* initialize */
	mt_do_z = NO;
	mt_do_r = NO;
	mt_do_t = NO;
	mt_do_n = NO;
	mt_do_e = NO;
	mt_mw = 2.6;
	mt_az = 0;
	mt_baz = -12345;
	mt_stk = 0.0;
	mt_rake = 0.0;
	mt_dip = 0.0;
	mt_fx = 0.0;
	mt_fy = 0.0;
	mt_fz = 0.0;
	mt_mxx = 0.0;
	mt_mxy = 0.0;
	mt_mxz = 0.0;
	mt_myy = 0.0;
	mt_myz = 0.0;
	mt_mzz = 0.0;
        mt_source = 0;
	/* parse commands */
	for(i=0 ; mtarg[i].key[0] != '\0' ; i++){
		if(mtarg[i].used > 0){
			if(mtarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, mtarg[i].key, 
					mtarg[i].mfit,mtarg[i].narg, mt_real);
			} else if(mtarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, mtarg[i].key, 
					mtarg[i].mfit,mtarg[i].narg, mt_int );
			} else if(mtarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, mtarg[i].key, 
					mtarg[i].mfit,mtarg[i].narg, &mt_yn );
			} else if(mtarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, mtarg[i].key, 
					mtarg[i].mfit,mtarg[i].narg, &mt_num );
			} else if(mtarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, mtarg[i].key, 
					mtarg[i].mfit,mtarg[i].narg, instr );
			}
			switch(mtarg[i].id){
				case MT_AZ:
					mt_az = mt_real[0];
					break;
				case MT_BAZ:
					mt_baz = mt_real[0];
					break;
				case MT_STK:
					mt_stk = mt_real[0];
					mt_source = 0;
					break;
				case MT_DIP:
					mt_dip = mt_real[0];
					mt_source = 0;
					break;
				case MT_RAKE:
					mt_rake = mt_real[0];
					mt_source = 0;
					break;
				case MT_FN:
					mt_fx = mt_real[0];
                                        mt_source = 2;
					break;
				case MT_FE:
					mt_fy = mt_real[0];
                                        mt_source = 2;
					break;
				case MT_FD:
					mt_fz = mt_real[0];
                                        mt_source = 2;
					break;
				case MT_MXX:
					mt_mxx = mt_real[0];
                                        mt_source = 3;
					break;
				case MT_MXY:
					mt_mxy = mt_real[0];
                                        mt_source = 3;
					break;
				case MT_MXZ:
					mt_mxz = mt_real[0];
                                        mt_source = 3;
					break;
				case MT_MYY:
					mt_myy = mt_real[0];
                                        mt_source = 3;
					break;
				case MT_MYZ:
					mt_myz = mt_real[0];
                                        mt_source = 3;
					break;
				case MT_MZZ:
					mt_mzz = mt_real[0];
                                        mt_source = 3;
					break;
				case MT_ISO:
                                        mt_source = 1;
					break;
				case MT_MW:
					mt_mw = mt_real[0];
					break;
				case MT_TO:
					if(strcmp(instr,"ZRT") == 0 || strcmp(instr,"zrt") == 0 ){
						mt_do_z = YES;
						mt_do_r = YES;
						mt_do_t = YES;
					} else if(strcmp(instr,"UZ") == 0 || strcmp(instr,"Z") == 0 ){
                                                mt_do_z = YES;
                                        } else if(strcmp(instr,"UR") == 0 || strcmp(instr,"R") == 0 ){
                                                mt_do_r = YES;
                                        } else if(strcmp(instr,"UT") == 0 || strcmp(instr,"T") == 0 ){
                                                mt_do_t = YES;
                                        } else if(strcmp(instr,"UN") == 0 || strcmp(instr,"N") == 0 ){
						mt_do_n = YES;
                                        } else if(strcmp(instr,"UE") == 0 || strcmp(instr,"E") == 0 ){
						mt_do_e = YES;
					} else if(strcmp(instr,"ZNE") == 0 || strcmp(instr,"zne") == 0 ){
						mt_do_z = YES;
						mt_do_n = YES;
						mt_do_e = YES;
                                        } ;
					break;
				case MT_FILE:
					ls = strlen(instr);
					mt_proto = calloc(ls+5,sizeof(char));
					strcpy(mt_proto,instr);
					break;

			}
		}
	}
	/* safety check on BAZ */
	if(mt_baz == -12345.)
		mt_baz = fmod(mt_az+180.0,360.0);
}

void gsac_exec_mt(void)
{
	float m[6];  /* really the 6 components of the moment tensor */
	float degrad;
	int ixx, iyy, ixy, ixz, iyz, izz;
	float sind, cosd, sin2d, cos2d;
	float sinr, cosr;
	float sins, coss, sin2s, cos2s;
	int ls;
	int oldmax;

	int i, j, k;
	int iret;
	int npts;
	float delta, stla, stlo, evla , evlo, az, baz, dist, gcarc;
	float  depmax, depmin, depmen;
	int indmax, indmin;
	float sinaz, cosaz, sin2az, cos2az;
	float fmax;
	float mom_fac; /* adjustment factor to achieve synthetic for the
		desired Mw. This accounts for the fact that the Green's functions
		are to 10^20 dyne-cm and cm/sec and that we wish to convert to m/s */
	float force_fac; /* adjustment factor to achieve synthetic for the so that	
		for a force in dynes we get output in m/s */
	float z, r, t, n, e ; /* these are individual sums that make the time series
		since we will overwrite in place */

	int sum_r, sum_z, sum_t;
	int zdd, zds, zss, zex, rdd, rds, rss, rex, tds, tss; /* these are pointers
		to the individual sac files containing the Green s functions */
	int zvf, rvf, zhf, rhf, thf;
	int kout; /* counter to inticate which input is over written */
	int kz, kr, kt, kn, ke;

	float fz1,fz2,fz3,fz4, fr1, fr2, fr3, fr4, ft1, ft2, ft3, ft4;
	float fzh, fth;
	float t0save, t1save;
	float sbaz, cbaz;
	

	degrad = 3.1415927/180.0 ;
	sbaz = sin(degrad*mt_baz);
	cbaz = cos(degrad*mt_baz);
	mom_fac = pow(10.0, 1.5*(mt_mw-2.60));	/* get moment as a multiple of 10^20 dyne-cm */
	mom_fac /= 100.0 ; 	/* add conversion factor to correct to m/s from cm/sec */
	force_fac = 1./100.0;  /* add conversion factor to correct to m/s from cm/sec */

printf("Creating synthetic seismogram:\n");
if(mt_source == 0 ) {
	printf("     stk %f dip %f rake %f mw %f\n",mt_stk,mt_dip,mt_rake,mt_mw);
} else if(mt_source == 1){
	printf("     explosion mw %f\n",mt_mw);
} else if(mt_source == 3){
	printf("Mxx=%10.3e\n",mt_mxx);
	printf("Mxy=%10.3e\n",mt_mxy);
	printf("Mxz=%10.3e\n",mt_mxz);
	printf("Myy=%10.3e\n",mt_myy);
	printf("Myz=%10.3e\n",mt_myz);
	printf("Mzz=%10.3e\n",mt_mzz);
}
	printf("     Z[%d] R[%d] T[%d] N[%d] E[%d]\n",mt_do_z, mt_do_r, mt_do_t, mt_do_n, mt_do_e);
printf("     Az %f\n",mt_az);
printf("     File pointer: %s\n",mt_proto);
printf("     mom_fac : %f\n",mom_fac);
	/* define the transformation matrix. This is taken from wvfmch96 */
	sind = sin(mt_dip * degrad);
	cosd = cos(mt_dip * degrad);
	sin2d = sin(2.0 *mt_dip * degrad);
	cos2d = cos(2.0 *mt_dip * degrad);
	sinr = sin(mt_rake * degrad);
	cosr = cos(mt_rake * degrad);
	sins = sin(mt_stk * degrad);
	coss = cos(mt_stk * degrad);
	sin2s = sin(2.0 *mt_stk * degrad);
	cos2s = cos(2.0 *mt_stk * degrad);

	ixx = 0 ;
	iyy = 1 ;
	ixy = 2 ;
	ixz = 3 ;
	iyz = 4 ;
	izz = 5 ;

	if(mt_source == 0){
		mom_fac = pow(10.0, 1.5*(mt_mw-2.60));	/* get moment as a multiple of 10^20 dyne-cm */
		mom_fac /= 100.0 ; 	/* add conversion factor to correct to m/s from cm/sec */
		sind = sin(mt_dip * degrad);
		cosd = cos(mt_dip * degrad);
		sin2d = sin(2.0 *mt_dip * degrad);
		cos2d = cos(2.0 *mt_dip * degrad);
		sinr = sin(mt_rake * degrad);
		cosr = cos(mt_rake * degrad);
		sins = sin(mt_stk * degrad);
		coss = cos(mt_stk * degrad);
		sin2s = sin(2.0 *mt_stk * degrad);
		cos2s = cos(2.0 *mt_stk * degrad);

		m[ixx] = -sind*cosr*sin2s - sin2d*sinr*sins*sins ;
		m[iyy] =  sind*cosr*sin2s - sin2d*sinr*coss*coss ;
		m[ixy] =  sind*cosr*cos2s + 0.5*sin2d*sinr*sin2s ;
		m[ixz] = -cosd*cosr*coss  - cos2d*sinr*sins ;
		m[iyz] = -cosd*cosr*sins  + cos2d*sinr*coss ;
		m[izz] =  sin2d*sinr ;

	} else if (mt_source == 1){
		mom_fac = pow(10.0, 1.5*(mt_mw-2.60));	/* get moment as a multiple of 10^20 dyne-cm */
		mom_fac /= 100.0 ; 	/* add conversion factor to correct to m/s from cm/sec */
		m[ixx] = 1.0 ;
		m[iyy] = 1.0 ;
		m[ixy] = 0.0 ;
		m[ixz] = 0.0 ;
		m[iyz] = 0.0 ;
		m[izz] = 1.0 ;
	} else if (mt_source == 3){
		mom_fac = 1.0e-20;
		mom_fac /= 100.0 ; 	/* add conversion factor to correct to m/s from cm/sec */
		m[ixx] = mt_mxx ;
		m[iyy] = mt_myy ;
		m[ixy] = mt_mxy ;
		m[ixz] = mt_mxz ;
		m[iyz] = mt_myz ;
		m[izz] = mt_mzz ;
	} else if (mt_source == 2){
		force_fac = 1./100.0;  /* add conversion factor to correct to m/s from cm/sec */
	}

        

	sinaz = sin(mt_az * degrad);
	cosaz = cos(mt_az * degrad);
	sin2az = sin(2.0 *mt_az * degrad);
	cos2az = cos(2.0 *mt_az * degrad);

	/* now create a string that is one of what we need  */
	ls = strlen(mt_proto);
	tstr = calloc(ls+5,sizeof(char));
	/* now go through the Green's */
	gsac_control.begmin= 1.0e+30;
	gsac_control.endmax=-1.0e+30;
	gsac_control.fft = NO;
	gsac_control.number_itraces = 0;
	gsac_control.number_iheaders= 0;
	gsac_control.number_otraces = 0;
	oldmax = gsac_control.max_number_traces ;
	gsac_control.max_number_traces = 0;

	zdd = -1 ; zds = -1 ; zss = -1 ; zex = -1;
	rdd = -1 ; rds = -1 ; rss = -1 ; rex = -1;
	tds = -1 ; tss = -1 ;

	sum_z = 0;
	sum_r = 0;
	sum_t = 0;

	if(mt_do_z == YES){
		if(mt_source == 0 || mt_source == 1 || mt_source == 3){
			gsac_mt_get(".ZDD",&zdd,&sum_z);
			gsac_mt_get(".ZDS",&zds,&sum_z);
			gsac_mt_get(".ZSS",&zss,&sum_z);
			gsac_mt_get(".ZEX",&zex,&sum_z);
		} else {
			gsac_mt_get(".ZVF",&zvf,&sum_z);
			gsac_mt_get(".ZHF",&zhf,&sum_z);
		}
	};
	if(mt_do_r == YES){
		if(mt_source == 0 || mt_source == 1 || mt_source == 3){
			gsac_mt_get(".RDD",&rdd,&sum_r);
			gsac_mt_get(".RDS",&rds,&sum_r);
			gsac_mt_get(".RSS",&rss,&sum_r);
			gsac_mt_get(".REX",&rex,&sum_r);
		} else {
			gsac_mt_get(".RVF",&rvf,&sum_z);
			gsac_mt_get(".RHF",&rhf,&sum_z);
		}
	};
	if(mt_do_t == YES){
		if(mt_source == 0 || mt_source == 1 || mt_source == 3){
			gsac_mt_get(".TDS",&tds,&sum_t);
			gsac_mt_get(".TSS",&tss,&sum_t);
		} else {
			gsac_mt_get(".THF",&thf,&sum_z);
		}
	};
	if(mt_do_n == YES || mt_do_e == YES){
		if(mt_source == 0 || mt_source == 1 || mt_source == 3){
			gsac_mt_get(".RDD",&rdd,&sum_r);
			gsac_mt_get(".RDS",&rds,&sum_r);
			gsac_mt_get(".RSS",&rss,&sum_r);
			gsac_mt_get(".REX",&rex,&sum_r);
			gsac_mt_get(".TDS",&tds,&sum_t);
			gsac_mt_get(".TSS",&tss,&sum_t);
		} else {
			gsac_mt_get(".RVF",&rvf,&sum_z);
			gsac_mt_get(".RHF",&rhf,&sum_z);
			gsac_mt_get(".THF",&thf,&sum_z);
		}
	}
	gsac_exec_read();
	/* now we have the Green's functions */
	/*
	put in error checking
	1 = zdd
	2 = zds,&t0save,0
	3 = zss
	4 = zex
            	fz(1) = -(xmt(1,1)+xmt(2,2))/6.0 + xmt(3,3)/3.0
            	fr(1) = fz(1)
            	fz(2) = xmt(1,3)*cosa + xmt(2,3)*sina
            	fz(3) = 0.5*(xmt(1,1)-xmt(2,2))*cos2a + xmt(1,2)*sin2a
            	fz(4) = (xmt(1,1)+xmt(2,2)+xmt(3,3))/3.0
            	fr(2) = fz(2)
            	fr(3) = fz(3)
            	fr(4) = fz(4)
            	ft(1) = 0.0
            	ft(2) = -xmt(2,3)*cosa + xmt(1,3)*sina
            	ft(3) = 0.5*(xmt(1,1)-xmt(2,2))*sin2a - xmt(1,2)*cos2a
            	ft(4) = 0.0
	
	*/
	fz1 = (- (m[ixx]+ m[iyy]) +2.*m[izz])/6.0;
	fz2 = (m[ixz]*cosaz + m[iyz]*sinaz) ;
	fz3 = 0.5*(m[ixx]-m[iyy])*cos2az + m[ixy]*sin2az;
	fz4 = ( m[ixx] + m[iyy] + m[izz] )/3.0;
	fr1 = fz1 ;
	fr2 = fz2 ;
	fr3 = fz3 ;
	fr4 = fz4 ;
	ft2 = m[ixz]*sinaz-m[iyz]*cosaz ;
	ft3 = 0.5*(m[ixx]- m[iyy])*sin2az - m[ixy]*cos2az ;

	fzh = mt_fx*cosaz + mt_fy*sinaz;
	fth = mt_fx*sinaz - mt_fy*cosaz;
	npts = sacdata[0].sachdr.ihdr[H_NPTS];
        for(i=0 ; i < npts; i++){
		kout = 0;
		if(mt_do_z == YES){
			if(mt_source == 0 || mt_source == 1 || mt_source == 3){
	        		z = (
	              		fz1*sacdata[zdd].sac_data[i]
	            		+ fz2*sacdata[zds].sac_data[i] 
		    		+ fz3*sacdata[zss].sac_data[i]
	            		+ fz4*sacdata[zex].sac_data[i] 
				) ;
				sacdata[kout].sac_data[i] = z*mom_fac;
				kz = kout;
				kout++;
				/* only have to look at one Green function to get S time */
                                t0save =  sacdata[zdd].sachdr.rhdr[H_T0];
			} else {
	        		z = (
	              		mt_fz*sacdata[zvf].sac_data[i]
				+ fzh*sacdata[zhf].sac_data[i]
				);
				sacdata[kout].sac_data[i] = z*force_fac;
				kz = kout;
				kout++;
                                t0save =  sacdata[zvf].sachdr.rhdr[H_T0];
			}
		}
		if(mt_do_r == YES){
			if(mt_source == 0 || mt_source == 1 || mt_source == 3){
	        		r = (
	              		fr1*sacdata[rdd].sac_data[i]
	            		+ fr2*sacdata[rds].sac_data[i] 
		    		+ fr3*sacdata[rss].sac_data[i]
	            		+ fr4*sacdata[rex].sac_data[i] 
				) ;
				sacdata[kout].sac_data[i] = r*mom_fac;
				kr = kout;
				kout++;
                                t0save =  sacdata[rdd].sachdr.rhdr[H_T0];
			} else {
	        		r = (
	              		mt_fz*sacdata[rvf].sac_data[i]
				+ fzh*sacdata[rhf].sac_data[i]
				);
				sacdata[kout].sac_data[i] = r*force_fac;
				kr = kout;
				kout++;
                                t0save =  sacdata[rvf].sachdr.rhdr[H_T0];
			}
		}
		if(mt_do_t == YES){
			if(mt_source == 0 || mt_source == 1 || mt_source == 3){
				t = (
		    		ft2*sacdata[tds].sac_data[i]
		  		+ ft3*sacdata[tss].sac_data[i]
				) ;
				sacdata[kout].sac_data[i] = t*mom_fac;
				kt = kout;
                                t1save =  sacdata[tss].sachdr.rhdr[H_T1];
				kout++;
			} else {
	        		t = (
	              		fth*sacdata[thf].sac_data[i]
				);
				sacdata[kout].sac_data[i] = t*force_fac;
				kt = kout;
				kout++;
                                t1save =  sacdata[thf].sachdr.rhdr[H_T1];
			}
		}
		if(mt_do_n == YES || mt_do_e ==YES){
			if(mt_source == 0 || mt_source == 1 || mt_source == 3){
	        		r = (
	              		fr1*sacdata[rdd].sac_data[i]
	            		+ fr2*sacdata[rds].sac_data[i] 
		    		+ fr3*sacdata[rss].sac_data[i]
	            		+ fr4*sacdata[rex].sac_data[i] 
				) ;
				t = (
		    		ft2*sacdata[tds].sac_data[i]
		  		+ ft3*sacdata[tss].sac_data[i]
				) ;
				/* get n e components */
				n = -cbaz*r + sbaz*t ;
				e = -sbaz*r - cbaz*t ;
				if(mt_do_n == YES){
					sacdata[kout].sac_data[i] = n*mom_fac;
					kn = kout;
                                t0save =  sacdata[rss].sachdr.rhdr[H_T0];
                                t1save =  sacdata[tss].sachdr.rhdr[H_T1];
					kout++;
				}
				if(mt_do_e == YES){
					sacdata[kout].sac_data[i] = e*mom_fac;
					ke = kout;
                                t0save =  sacdata[rss].sachdr.rhdr[H_T0];
                                t1save =  sacdata[tss].sachdr.rhdr[H_T1];
					kout++;
				}
			} else {
	        		r = (
	              		mt_fz*sacdata[rvf].sac_data[i]
				+ fzh*sacdata[rhf].sac_data[i]
				);
	        		t = (
	              		fth*sacdata[thf].sac_data[i]
				);
				/* get n e components */
				n = -cbaz*r + sbaz*t ;
				e = -sbaz*r - cbaz*t ;
				if(mt_do_n == YES){
					sacdata[kout].sac_data[i] = n*force_fac;
					kn = kout;
                                	t0save =  sacdata[rhf].sachdr.rhdr[H_T0];
                                	t1save =  sacdata[thf].sachdr.rhdr[H_T1];
					kout++;
				}
				if(mt_do_e == YES){
					sacdata[kout].sac_data[i] = e*force_fac;
					ke = kout;
                                	t0save =  sacdata[rhf].sachdr.rhdr[H_T0];
                                	t1save =  sacdata[thf].sachdr.rhdr[H_T1];
					kout++;
				}
			}
		}
	}
	/* reset the header variables */
	for(k=0 ; k < kout ; k ++){
	getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
	sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
	sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
	sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
	sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
	sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
	}
	gsac_control.number_otraces = kout ;
	if(mt_do_z){
		strcpy(sacdata[kz].sac_ofile_name,"T.Z");
		strncpy(sacdata[kz].sachdr.chdr[H_KSTNM], "GRN     ",8);
		strncpy(sacdata[kz].sachdr.chdr[H_KCMPNM], "BHZ     ",8);
		strncpy(sacdata[kz].schdr[H_KCMPNM], "BHZ     ",8);
		strncpy(sacdata[kz].sachdr.chdr[H_KT0], "SV      ",8);
		strncpy(sacdata[kz].sachdr.chdr[H_KT1], "-12345  ",8);
		sacdata[kz].sachdr.rhdr[H_CMPINC] = 0;
		sacdata[kz].sachdr.rhdr[H_CMPAZ] = 0;
		sacdata[kz].sachdr.rhdr[H_AZ] = mt_az;
		sacdata[kz].sachdr.rhdr[H_BAZ] = fmod(mt_baz,360.0);
		sacdata[kz].sachdr.rhdr[H_T0] = t0save;
		sacdata[kz].sachdr.rhdr[H_T1] = -12345.;
	};
	if(mt_do_r){
		strcpy(sacdata[kr].sac_ofile_name,"T.R");
		strncpy(sacdata[kr].sachdr.chdr[H_KSTNM], "GRN     ",8);
		strncpy(sacdata[kr].sachdr.chdr[H_KCMPNM], "BHR     ",8);
		strncpy(sacdata[kr].schdr[H_KCMPNM], "BHR     ",8);
		strncpy(sacdata[kr].sachdr.chdr[H_KT0], "SV      ",8);
		strncpy(sacdata[kr].sachdr.chdr[H_KT1], "-12345  ",8);
		sacdata[kr].sachdr.rhdr[H_CMPINC] = 90;
		sacdata[kr].sachdr.rhdr[H_CMPAZ] = fmod(mt_baz+180.0,360.);
		sacdata[kr].sachdr.rhdr[H_AZ] = mt_az;
		sacdata[kr].sachdr.rhdr[H_BAZ] = fmod(mt_baz,360.0);
		sacdata[kr].sachdr.rhdr[H_T0] = t0save;
		sacdata[kr].sachdr.rhdr[H_T1] = -12345.;
	};
	if(mt_do_t){
		strcpy(sacdata[kt].sac_ofile_name,"T.T");
		strncpy(sacdata[kt].sachdr.chdr[H_KSTNM], "GRN     ",8);
		strncpy(sacdata[kt].sachdr.chdr[H_KCMPNM], "BHT     ",8);
		strncpy(sacdata[kt].sachdr.chdr[H_KT1], "SH      ",8);
		strncpy(sacdata[kt].sachdr.chdr[H_KT0], "-12345  ",8);
		strncpy(sacdata[kt].schdr[H_KCMPNM], "BHT     ",8);
		sacdata[kt].sachdr.rhdr[H_CMPINC] = 90;
		sacdata[kt].sachdr.rhdr[H_CMPAZ] = fmod(mt_baz-90, 360.);
		sacdata[kt].sachdr.rhdr[H_AZ] = mt_az;
		sacdata[kt].sachdr.rhdr[H_BAZ] = fmod(mt_baz,360.0);
		sacdata[kt].sachdr.rhdr[H_T1] = t1save;
		sacdata[kt].sachdr.rhdr[H_T0] = -12345.;
	};
	if(mt_do_n){
		strcpy(sacdata[kn].sac_ofile_name,"T.N");
		strncpy(sacdata[kn].sachdr.chdr[H_KSTNM], "GRN     ",8);
		strncpy(sacdata[kn].sachdr.chdr[H_KCMPNM], "BHN     ",8);
		strncpy(sacdata[kn].sachdr.chdr[H_KT1], "SV      ",8);
		strncpy(sacdata[kn].sachdr.chdr[H_KT0], "-12345  ",8);
		strncpy(sacdata[kn].schdr[H_KCMPNM], "BHN     ",8);
		sacdata[kn].sachdr.rhdr[H_CMPINC] = 90;
		sacdata[kn].sachdr.rhdr[H_CMPAZ] = 0;
		sacdata[kn].sachdr.rhdr[H_AZ] = mt_az;
		sacdata[kn].sachdr.rhdr[H_BAZ] = fmod(mt_baz,360.0);
		sacdata[kn].sachdr.rhdr[H_T0] = t0save;
		sacdata[kn].sachdr.rhdr[H_T1] = t1save;
	};
	if(mt_do_e){
		strcpy(sacdata[ke].sac_ofile_name,"T.E");
		strncpy(sacdata[ke].sachdr.chdr[H_KSTNM], "GRN     ",8);
		strncpy(sacdata[ke].sachdr.chdr[H_KCMPNM], "BHE     ",8);
		strncpy(sacdata[ke].sachdr.chdr[H_KT1], "SV      ",8);
		strncpy(sacdata[ke].sachdr.chdr[H_KT0], "-12345  ",8);
		strncpy(sacdata[ke].schdr[H_KCMPNM], "BHE     ",8);
		sacdata[ke].sachdr.rhdr[H_CMPINC] = 90;
		sacdata[ke].sachdr.rhdr[H_CMPAZ] = 90;
		sacdata[ke].sachdr.rhdr[H_AZ] = mt_az;
		sacdata[ke].sachdr.rhdr[H_BAZ] = fmod(mt_baz,360.0);
		sacdata[ke].sachdr.rhdr[H_T0] = t0save;
		sacdata[ke].sachdr.rhdr[H_T1] = t1save;
	};



	free (tstr);

}
static void gsac_mt_get(char *grn, int *tds,int *sum_t)
{
	int k;
	int oldmax;
	oldmax = gsac_control.max_number_traces ;
		strcpy(tstr,mt_proto);strcat(tstr,grn);
		if (gsac_valid_sacfile(tstr) > 0){
			/* allocate a data structure */
			gsac_alloc_trace(oldmax);
			/* now get the data */
			k = gsac_control.max_number_traces -1;
			*tds = k;
			/* beware of size limitations */
			strcpy(sacdata[k].sac_ifile_name , tstr);
			(*sum_t)++;
		}
}
