#include	<stdio.h>
#include	<string.h>
#include "gsac.h"
#include "gsac_docommand.h"
#include "gsac_sac.h"
#include "gsac_sachdr.h"
#include "gsac_arg.h"
#include "csstim.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


/* variables */

static float  fhdr_default = -12345.;
static int    ihdr_default = -12345 ;
static char  *chdr_default = "-12345  ";

char  outstr[80];

static int lh_listcolumns = YES;


#define DEFAULT -1
#define COLUMNS -2
#define KZDATE  100
#define KZTIME  101
#define KEVNM   102

/* this the is order in which they are listed 
 * int id
 * char *key
 * int ricell
 * int used
 * int narg
 * int show
 * char *errormessage;
 *
 * */
struct arghdr lharg[] = {
	{ DEFAULT, "DEFAULT" , IHDR,  0, 0,  NO, "" ,-1},
	{ COLUMNS, "COLUMNS" , IHDR,  0, 1,  NO, "COLUMNS 1 or 2" , 2},
	{  9, "NPTS"    , IHDR, YES, 0, YES, "" ,-1},
	{  5, "B"       , RHDR, YES, 0, YES, "" ,-1},
	{  6, "E"       , RHDR, YES, 0, YES, "" ,-1},
	{  0, "DELTA"   , RHDR, YES, 0, YES, "" ,-1},
	{  2, "DEPMAX"  , RHDR, YES, 0, YES, "" ,-1},
	{  1, "DEPMIN"  , RHDR, YES, 0, YES, "" ,-1},
	{ 56, "DEPMEN"  , RHDR, YES, 0, YES, "" ,-1},
	{  0, "NZYEAR"  , IHDR, YES, 0, YES, "" ,-1},
	{  1, "NZJDAY"  , IHDR, YES, 0, YES, "" ,-1},
	{  2, "NZHOUR"  , IHDR, YES, 0, YES, "" ,-1},
	{  3, "NZMIN"   , IHDR, YES, 0, YES, "" ,-1},
	{  4, "NZSEC"   , IHDR, YES, 0, YES, "" ,-1},
	{  5, "NZMSEC"  , IHDR, YES, 0, YES, "" ,-1},
	{ KZDATE, "KZDATE"  , XHDR, YES, 0, YES, "" ,-1},
	{ KZTIME, "KZTIME"  , XHDR, YES, 0, YES, "" ,-1},
	{  7, "O"       , RHDR, YES, 0, YES, "" ,-1},
	{  8, "A"       , RHDR, YES, 0, YES, "" ,-1},
	{ 10, "T0"      , RHDR, YES, 0, YES, "" ,-1},
	{  0, "KSTNM"   , CHDR, YES, 0, YES, "" ,-1},
	{ 20, "KCMPNM"  , CHDR, YES, 0, YES, "" ,-1},
	{ 31, "STLA"    , RHDR, YES, 0, YES, "" ,-1},
	{ 32, "STLO"    , RHDR, YES, 0, YES, "" ,-1},
	{ 33, "STEL"    , RHDR, YES, 0, YES, "" ,-1},
	{ 35, "EVLA"    , RHDR, YES, 0, YES, "" ,-1},
	{ 36, "EVLO"    , RHDR, YES, 0, YES, "" ,-1},
	{ 37, "EVEL"    , RHDR, YES, 0, YES, "" ,-1},
	{ 50, "DIST"    , RHDR, YES, 0, YES, "" ,-1},
	{ 51, "AZ"      , RHDR, YES, 0, YES, "" ,-1},
	{ 52, "BAZ"     , RHDR, YES, 0, YES, "" ,-1},
	{ 53, "GCARC"   , RHDR, YES, 0, YES, "" ,-1},
	{  3, "SCALE"   , RHDR, YES, 0, YES, "" ,-1},
	{  4, "ODELTA"  , RHDR, YES, 0, YES, "" ,-1},
	{  9, "FMT"     , RHDR, YES, 0, YES, "" ,-1},
	{ 11, "T1"      , RHDR, YES, 0, YES, "" ,-1},
	{ 12, "T2"      , RHDR, YES, 0, YES, "" ,-1},
	{ 13, "T3"      , RHDR, YES, 0, YES, "" ,-1},
	{ 14, "T4"      , RHDR, YES, 0, YES, "" ,-1},
	{ 15, "T5"      , RHDR, YES, 0, YES, "" ,-1},
	{ 16, "T6"      , RHDR, YES, 0, YES, "" ,-1},
	{ 17, "T7"      , RHDR, YES, 0, YES, "" ,-1},
	{ 18, "T8"      , RHDR, YES, 0, YES, "" ,-1},
	{ 19, "T9"      , RHDR, YES, 0, YES, "" ,-1},
	{ 20, "F"       , RHDR, YES, 0, YES, "" ,-1},
	{ 21, "RESP0"   , RHDR, YES, 0, YES, "" ,-1},
	{ 22, "RESP1"   , RHDR, YES, 0, YES, "" ,-1},
	{ 23, "RESP2"   , RHDR, YES, 0, YES, "" ,-1},
	{ 24, "RESP3"   , RHDR, YES, 0, YES, "" ,-1},
	{ 25, "RESP4"   , RHDR, YES, 0, YES, "" ,-1},
	{ 26, "RESP5"   , RHDR, YES, 0, YES, "" ,-1},
	{ 27, "RESP6"   , RHDR, YES, 0, YES, "" ,-1},
	{ 28, "RESP7"   , RHDR, YES, 0, YES, "" ,-1},
	{ 29, "RESP8"   , RHDR, YES, 0, YES, "" ,-1},
	{ 30, "RESP9"   , RHDR, YES, 0, YES, "" ,-1},
	{ 34, "STDP"    , RHDR, YES, 0, YES, "" ,-1},
	{ 38, "EVDP"    , RHDR, YES, 0, YES, "" ,-1},
	{ 39, "MAG"     , RHDR, YES, 0, YES, "" ,-1},
	{ 40, "USER0"   , RHDR, YES, 0, YES, "" ,-1},
	{ 41, "USER1"   , RHDR, YES, 0, YES, "" ,-1},
	{ 42, "USER2"   , RHDR, YES, 0, YES, "" ,-1},
	{ 43, "USER3"   , RHDR, YES, 0, YES, "" ,-1},
	{ 44, "USER4"   , RHDR, YES, 0, YES, "" ,-1},
	{ 45, "USER5"   , RHDR, YES, 0, YES, "" ,-1},
	{ 46, "USER6"   , RHDR, YES, 0, YES, "" ,-1},
	{ 47, "USER7"   , RHDR, YES, 0, YES, "" ,-1},
	{ 48, "USER8"   , RHDR, YES, 0, YES, "" ,-1},
	{ 49, "USER9"   , RHDR, YES, 0, YES, "" ,-1},
	{ 54, "SB"      , RHDR, YES, 0, YES, "" ,-1},
	{ 55, "SDELTA"  , RHDR, YES, 0, YES, "" ,-1},
	{ 57, "CMPAZ"   , RHDR, YES, 0, YES, "" ,-1},
	{ 58, "CMPINC"  , RHDR, YES, 0, YES, "" ,-1},
	{ 59, "XMINIMUM", RHDR, YES, 0, YES, "" ,-1},
	{ 60, "XMAXIMUM", RHDR, YES, 0, YES, "" ,-1},
	{ 61, "YMINIMUM", RHDR, YES, 0, YES, "" ,-1},
	{ 62, "YMAXIMUM", RHDR, YES, 0, YES, "" ,-1},
	{ 63, "ADJTM"   , RHDR, YES, 0, YES, "" ,-1},
	{ 64, "TIMMAX"  , RHDR, YES, 0, YES, "" ,-1},
	{ 65, "TIMMIN"  , RHDR, YES, 0, YES, "" ,-1},
	{ 66, "FHDR67"  , RHDR, YES, 0, YES, "" ,-1},
	{ 67, "FHDR68"  , RHDR, YES, 0, YES, "" ,-1},
	{ 68, "FHDR69"  , RHDR, YES, 0, YES, "" ,-1},
	{ 69, "FHDR70"  , RHDR, YES, 0, YES, "" ,-1},
	{  6, "NVHDR"   , IHDR, YES, 0, YES, "" ,-1},
	{  7, "NINF"    , IHDR, YES, 0, YES, "" ,-1},
	{  8, "NHST"    , IHDR, YES, 0, YES, "" ,-1},
	{ 10, "NSNPTS"  , IHDR, YES, 0, YES, "" ,-1},
	{ 11, "NSN"     , IHDR, YES, 0, YES, "" ,-1},
	{ 12, "NXSIZE"  , IHDR, YES, 0, YES, "" ,-1},
	{ 13, "NYSIZE"  , IHDR, YES, 0, YES, "" ,-1},
	{ 14, "NHDR15"  , IHDR, YES, 0, YES, "" ,-1},
	{ 15, "IFTYPE"  , EHDR, YES, 0, YES, "" ,-1},
	{ 16, "IDEP"    , EHDR, YES, 0, YES, "" ,-1},
	{ 17, "IZTYPE"  , EHDR, YES, 0, YES, "" ,-1},
	{ 18, "IHDR4"   , EHDR, YES, 0, YES, "" ,-1},
	{ 19, "IINST"   , EHDR, YES, 0, YES, "" ,-1},
	{ 20, "ISTREG"  , IHDR, YES, 0, YES, "" ,-1},
	{ 21, "IEVREG"  , IHDR, YES, 0, YES, "" ,-1},
	{ 22, "IEVTYP"  , EHDR, YES, 0, YES, "" ,-1},
	{ 23, "IQUAL"   , EHDR, YES, 0, YES, "" ,-1},
	{ 24, "ISYNTH"  , EHDR, YES, 0, YES, "" ,-1},
	{ 25, "IHDR11"  , IHDR, YES, 0, YES, "" ,-1},
	{ 26, "IHDR12"  , IHDR, YES, 0, YES, "" ,-1},
	{ 27, "IHDR13"  , IHDR, YES, 0, YES, "" ,-1},
	{ 28, "IHDR14"  , IHDR, YES, 0, YES, "" ,-1},
	{ 29, "IHDR15"  , IHDR, YES, 0, YES, "" ,-1},
	{ 30, "IHDR16"  , IHDR, YES, 0, YES, "" ,-1},
	{ 31, "IHDR17"  , IHDR, YES, 0, YES, "" ,-1},
	{ 32, "IHDR18"  , IHDR, YES, 0, YES, "" ,-1},
	{ 33, "IHDR19"  , IHDR, YES, 0, YES, "" ,-1},
	{ 34, "IHDR20"  , IHDR, YES, 0, YES, "" ,-1},
	{ 36, "LPSPOL"  , LHDR, YES, 0, YES, "" ,-1},
	{ 37, "LOVROK"  , LHDR, YES, 0, YES, "" ,-1},
	{ 38, "LCALDA"  , LHDR, YES, 0, YES, "" ,-1},
	{ 39, "LHDR5"   , LHDR, YES, 0, YES, "" ,-1},
	{  1,"KEVNM" , CHDL, YES, 0, YES, "" ,-1},
	{  3, "KHOLE"   , CHDR, YES, 0, YES, "" ,-1},
	{  4, "KO"      , CHDR, YES, 0, YES, "" ,-1},
	{  5, "KA"      , CHDR, YES, 0, YES, "" ,-1},
	{  6, "KT0"     , CHDR, YES, 0, YES, "" ,-1},
	{  7, "KT1"     , CHDR, YES, 0, YES, "" ,-1},
	{  8, "KT2"     , CHDR, YES, 0, YES, "" ,-1},
	{  9, "KT3"     , CHDR, YES, 0, YES, "" ,-1},
	{ 10, "KT4"     , CHDR, YES, 0, YES, "" ,-1},
	{ 11, "KT5"     , CHDR, YES, 0, YES, "" ,-1},
	{ 12, "KT6"     , CHDR, YES, 0, YES, "" ,-1},
	{ 13, "KT7"     , CHDR, YES, 0, YES, "" ,-1},
	{ 14, "KT8"     , CHDR, YES, 0, YES, "" ,-1},
	{ 15, "KT9"     , CHDR, YES, 0, YES, "" ,-1},
	{ 16, "KF"      , CHDR, YES, 0, YES, "" ,-1},
	{ 17, "KUSER0"  , CHDR, YES, 0, YES, "" ,-1},
	{ 18, "KUSER1"  , CHDR, YES, 0, YES, "" ,-1},
	{ 19, "KUSER2"  , CHDR, YES, 0, YES, "" ,-1},
	{ 21, "KNETWK"  , CHDR, YES, 0, YES, "" ,-1},
	{ 22, "KDATRD"  , CHDR, YES, 0, YES, "" ,-1},
	{ 23, "KINST"   , CHDR, YES, 0, YES, "" ,-1},
	{-10, ""        , CHDR, YES, 0, YES, "" ,-1}
};	

static int NESTR=50;
char *estr[] = {
	"ITIME   ", "IRLIM   ", "IAMPH   ", "IXY     ", "IUNKN   ", 
	"IDISP   ", "IVEL    ", "IACC    ", "IB      ", "IDAY    ", 
	"IO      ", "IA      ", "IT0     ", "IT1     ", "IT2     ", 
	"IT3     ", "IT4     ", "IT5     ", "IT6     ", "IT7     ", 
	"IT8     ", "IT9     ", "IRADNV  ", "ITANNV  ", "IRADEV  ", 
	"ITANEV  ", "INORTH  ", "IEAST   ", "IHORZA  ", "IDOWN   ", 
	"IUP     ", "ILLLBB  ", "IWWSN1  ", "IWWSN2  ", "IHGLP   ", 
	"ISRO    ", "INUCL   ", "IPREN   ", "IPOSTN  ", "IQUAKE  ", 
	"IPREQ   ", "IPOSTQ  ", "ICHEM   ", "IOTHER  ", "IGOOD   ", 
	"IGLCH   ", "IDROP   ", "ILOWSN  ", "IRLDTA  ", "IVOLTS  "
	} ;
static int NIVAL=50;
char *Istr[] = {
	"IFTYPE  ", "IFTYPE  ", "IFTYPE  ", "IFTYPE  ", "IDEP    ",
	"IDEP    ", "IDEP    ", "IDEP    ", "IZTYPE  ", "IZTYPE  ",
	"IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ",
	"IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ",
	"IZTYPE  ", "IZTYPE  ", "IRADNV  ", "ITANNV  ", "IRADEV  ",
	"ITANEV  ", "INORTH  ", "IEAST   ", "IHORZA  ", "IDOWN   ",
	"IUP     ", "ILLLBB  ", "IWWSN1  ", "IWWSN2  ", "IHGLP   ",
	"ISRO    ", "IEVTYP  ", "IEVTYP  ", "IEVTYP  ", "IEVTYP  ",
	"IEVTYP  ", "IEVTYP  ", "IEVTYP  ", "IQUAL   ", "IQUAL   ",
	"IQUAL   ", "IQUAL   ", "IQUAL   ", "ISYNTH  ", "IDEP    "    
};

int lh_int[10];
/* the memory of used */
static int *lh_used = (int *)NULL;

void gsac_set_param_lh(int ncmd, char **cmdstr)
{
	int i, setdefault;
	int anychanged ;
	/* initialization */
	if(lh_used == (int *)NULL){
		/* initialize a local array - this is not hardwired
		 * to permit expansion of the lharg structure */
		lh_used = (int *)realloc(lh_used, sizeof(int )*sizeof(lharg)/sizeof( struct arghdr ));
		/* now build up the memory */
		for(i=0; lharg[i].key[0] != '\0' ; i++){
			lh_used[i] = lharg[i].used;
		}
	}

	if(ncmd == 1)
		return;
	/* is the command syntax correct ? */
	if(testarg(ncmd, cmdstr, lharg, YES, YES))
		return;
	/* now go through the arguments to see if we want to display
	 * them. Note that the command DEFAULT means that we will
	 * reset the list of display */
	/* OK now mark all unused and reparse */
	testarg(ncmd, cmdstr, lharg, NO, YES);
	setdefault = NO; 
	anychanged = 0;
	for(i=0 ; lharg[i].key[0] != '\0' ; i++){
		/* check for special commands */
		if(lharg[i].used > 0){
			if(lharg[i].id == DEFAULT){
				setdefault = YES;
			} else if(lharg[i].id == COLUMNS){
				getargi(ncmd, cmdstr, lharg[i].key, 
					lharg[i].mfit, lharg[i].narg, lh_int );
				if(lh_int[0] == 1){
				lh_listcolumns = NO;
				} else if(lh_int[0] == 2) {
					lh_listcolumns = YES;
				}
			}
			if(lharg[i].id >= 0 )
				anychanged++;
		}
	}
	if(setdefault){
		printf("Setting default LH\n");
		for(i=0 ; lharg[i].key[0] != '\0' ; i++){
			if(lharg[i].id >= 0){
				lh_used[i] = YES;
			} else {
				lh_used[i] = NO;
			}
		}
	} else {
		/* OK if one is used reset else use old */
		if(anychanged){
		for(i=0 ; lharg[i].key[0] != '\0' ; i++){
			if(lharg[i].used){
				lh_used[i] = lharg[i].used;
			} else {
				lh_used[i] = NO;
			}
		}
		}
	}
}


void gsac_exec_lh(void)
{
int i, k, kk, ipos;
int ntrchdr, col, prn;
/* if there are no traces return */
	ntrchdr = gsac_control.number_iheaders;
	if ( ntrchdr < 1)
		return;
/* sequence through all files */
for ( kk=0 ; kk < ntrchdr ; kk ++){
	k = sortptr[kk];
	if(sacdata[k].display == YES){
	printf("%s (%d):\n",sacdata[k].sac_ofile_name,k);
	/* now that we have set everything, define the number to be listed */
	for(i=0, col = 0; lharg[i].key[0] != '\0'  ;i++){
		if(lh_used[i] == YES && lharg[i].id >=0 ){
			ipos = lharg[i].id;
			prn = NO;
			if(lharg[i].ricell == RHDR && sacdata[k].sachdr.rhdr[ipos] != fhdr_default){
				printf("     %8s %20.7g",lharg[i].key,sacdata[k].sachdr.rhdr[ipos]);
				col++;
				prn = YES;
			} else if(lharg[i].ricell == IHDR && sacdata[k].sachdr.ihdr[ipos] != ihdr_default){
				printf("     %8s %20d",lharg[i].key,sacdata[k].sachdr.ihdr[ipos]);
				col++;
				prn = YES;
			} else if(lharg[i].ricell == LHDR && sacdata[k].sachdr.ihdr[ipos] != ihdr_default){
				/* WHAT IS THIS */
				if(sacdata[k].sachdr.ihdr[ipos]==0)
				printf("     %8s                FALSE",lharg[i].key);
				else
				printf("     %8s                 TRUE",lharg[i].key);
				col++;
				prn = YES;
			} else if(lharg[i].ricell==CHDR && strncmp(sacdata[k].sachdr.chdr[ipos],chdr_default,8)!=0){
				printf("     %8s %20s",lharg[i].key,sacdata[k].schdr[ipos]);
				col++;
				prn = YES;
			} else if(lharg[i].ricell == CHDL 
				&& strncmp(sacdata[k].sachdr.chdr[1], chdr_default,8)!=0 
				&& strncmp(sacdata[k].sachdr.chdr[2], chdr_default,8)!=0){
				printf("     %8s     %8s%8s",lharg[i].key,sacdata[k].schdr[1],sacdata[k].schdr[2]);
				col++;
				prn = YES;
			} else if(lharg[i].ricell == EHDR && sacdata[k].sachdr.ihdr[ipos] != ihdr_default){
				/* we only are interested in IFTYPE (15) and IZTYPE(17) */
				/* HACK 24 APR 2006 - do not permit IEVREG - i = 97 */
				if(i != 97){
				printf("     %8s %20s",lharg[i].key,estr[sacdata[k].sachdr.ihdr[ipos]-1]);
				col++;
				prn = YES;
				} else {
				prn = NO ;
				}
			} else if(lharg[i].ricell == XHDR){
				if(lharg[i].id == KZTIME){
					printkztimestr(sacdata[k].sachdr.ihdr[H_NZHOUR],
							sacdata[k].sachdr.ihdr[H_NZMIN],
							sacdata[k].sachdr.ihdr[H_NZSEC],
							sacdata[k].sachdr.ihdr[H_NZMSEC],outstr);
					printf("       %s",outstr);
				col++;
				prn = YES;
				} else if ( lharg[i].id == KZDATE){
					printkzdatestr(sacdata[k].sachdr.ihdr[H_NZYEAR],
							sacdata[k].sachdr.ihdr[H_NZJDAY],outstr);
					printf("       %s",outstr);
				col++;
				prn = YES;

				}
			}
			/* implement single and double column */
			if(prn == YES){
				if(lh_listcolumns ==NO) {
					printf("\n") ;
				} else {
					if(col == 2){
						printf("\n");
						col = 0;
					}
				}
			}
		}
	}
	if(lh_listcolumns == YES && col == 1)
		printf("\n");
	}
	}
}


