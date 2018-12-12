/* sort according to trace header values
 * Note that only a subset of header values are subject to search
 
 09 APR 2008 - correct error that had IHDR12 as IDHR12, added IHDR11 and
    changed the array to use the ihdr instead of rhdr field
 20 AUG 2010 - permit a sort on the MAG field
 * */
#include        <stdio.h>
#include        "gsac.h"
#include        "gsac_sac.h"
#include	"gsac_docommand.h"
#include	"gsac_arg.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;
extern float *sortfloat;

#define ST_REVERSE	-1
#define ST_FORWARD	-2
#define	ST_DFLT 	-3

int sort_which;	/* define header field to search on 
			   When used in conjunction with sort_type defines
			   exactly , e.g., DIST = RHDR[50] */
int sort_type;	/* define character of sort, DEFAULT, REVERSE, 
			   RHDR IHDR */
int sort_reverse = NO;
int sort_do = NO ;
static int isort ;	/* points to sort string */
void gnomesort(int n, float ar[], int key[] ) ;
	/* http://www.cs.vu.nl/~dick/gnomesort.html
	 * */
void gnomesortint(int n, int ar[], int key[]) ;
void intrev(int *x, int n);


struct arghdr starg[] = {
	{ST_DFLT    , "OFF"	, CHDR, 0, 0, NO, "OFF ",-1},
	{ST_FORWARD , "FORWARD"	, CHDR, 0, 0, NO, "FORWARD ",-1},
	{ST_FORWARD , "UP"	, CHDR, 0, 0, NO, "FORWARD ",-1},
	{ST_FORWARD , "ASCEND"	, CHDR, 0, 0, NO, "FORWARD ",-1},
	{ST_REVERSE , "REVERSE"	, CHDR, 0, 0, NO, "REVERSE ",-1},
	{ST_REVERSE , "DOWN"	, CHDR, 0, 0, NO, "REVERSE ",-1},
	{ST_REVERSE , "DESCEND"	, CHDR, 0, 0, NO, "REVERSE ",-1},
	{ST_DFLT    , "DEFAULT"	, CHDR, 0, 0, NO, "DEFAULT",-1},
	{  9, "NPTS"    , IHDR, 0, 0, NO, "" ,-1},
	{  5, "B"       , RHDR, 0, 0, NO, "" ,-1},
	{  6, "E"       , RHDR, 0, 0, NO, "" ,-1},
	{  0, "DELTA"   , RHDR, 0, 0, NO, "" ,-1},
	{  1, "DEPMIN"  , RHDR, 0, 0, NO, "" ,-1},
	{  2, "DEPMAX"  , RHDR, 0, 0, NO, "" ,-1},
	{ 56, "DEPMEN"  , RHDR, 0, 0, NO, "" ,-1},
	{  0, "NZYEAR"  , IHDR, 0, 0, NO, "" ,-1},
	{  1, "NZJDAY"  , IHDR, 0, 0, NO, "" ,-1},
	{  2, "NZHOUR"  , IHDR, 0, 0, NO, "" ,-1},
	{  3, "NZMIN"   , IHDR, 0, 0, NO, "" ,-1},
	{  4, "NZSEC"   , IHDR, 0, 0, NO, "" ,-1},
	{  5, "NZMSEC"  , IHDR, 0, 0, NO, "" ,-1},
	{  7, "O"       , RHDR, 0, 0, NO, "" ,-1},
	{  8, "A"       , RHDR, 0, 0, NO, "" ,-1},
	{ 10, "T0"      , RHDR, 0, 0, NO, "" ,-1},
	{  0, "KSTNM"   , CHDR, 0, 0, NO, "" ,-1},
	{ 20, "KCMPNM"  , CHDR, 0, 0, NO, "" ,-1},
	{ 31, "STLA"    , RHDR, 0, 0, NO, "" ,-1},
	{ 32, "STLO"    , RHDR, 0, 0, NO, "" ,-1},
	{ 33, "STEL"    , RHDR, 0, 0, NO, "" ,-1},
	{ 35, "EVLA"    , RHDR, 0, 0, NO, "" ,-1},
	{ 36, "EVLO"    , RHDR, 0, 0, NO, "" ,-1},
	{ 37, "EVEL"    , RHDR, 0, 0, NO, "" ,-1},
	{ 50, "DIST"    , RHDR, 0, 0, NO, "" ,-1},
	{ 51, "AZ"      , RHDR, 0, 0, NO, "" ,-1},
	{ 52, "BAZ"     , RHDR, 0, 0, NO, "" ,-1},
	{ 53, "GCARC"   , RHDR, 0, 0, NO, "" ,-1},
	{  3, "SCALE"   , RHDR, 0, 0, NO, "" ,-1},
	{  4, "ODELTA"  , RHDR, 0, 0, NO, "" ,-1},
	{  9, "FMT"     , RHDR, 0, 0, NO, "" ,-1},
	{ 11, "T1"      , RHDR, 0, 0, NO, "" ,-1},
	{ 12, "T2"      , RHDR, 0, 0, NO, "" ,-1},
	{ 13, "T3"      , RHDR, 0, 0, NO, "" ,-1},
	{ 14, "T4"      , RHDR, 0, 0, NO, "" ,-1},
	{ 15, "T5"      , RHDR, 0, 0, NO, "" ,-1},
	{ 16, "T6"      , RHDR, 0, 0, NO, "" ,-1},
	{ 17, "T7"      , RHDR, 0, 0, NO, "" ,-1},
	{ 18, "T8"      , RHDR, 0, 0, NO, "" ,-1},
	{ 19, "T9"      , RHDR, 0, 0, NO, "" ,-1},
	{ 20, "F"       , RHDR, 0, 0, NO, "" ,-1},
	{ 21, "RESP0"   , RHDR, 0, 0, NO, "" ,-1},
	{ 22, "RESP1"   , RHDR, 0, 0, NO, "" ,-1},
	{ 23, "RESP2"   , RHDR, 0, 0, NO, "" ,-1},
	{ 24, "RESP3"   , RHDR, 0, 0, NO, "" ,-1},
	{ 25, "RESP4"   , RHDR, 0, 0, NO, "" ,-1},
	{ 26, "RESP5"   , RHDR, 0, 0, NO, "" ,-1},
	{ 27, "RESP6"   , RHDR, 0, 0, NO, "" ,-1},
	{ 28, "RESP7"   , RHDR, 0, 0, NO, "" ,-1},
	{ 29, "RESP8"   , RHDR, 0, 0, NO, "" ,-1},
	{ 30, "RESP9"   , RHDR, 0, 0, NO, "" ,-1},
	{ 31, "STLA"    , RHDR, 0, 0, NO, "" ,-1},
	{ 32, "STLO"    , RHDR, 0, 0, NO, "" ,-1},
	{ 33, "STEL"    , RHDR, 0, 0, NO, "" ,-1},
	{ 34, "STDP"    , RHDR, 0, 0, NO, "" ,-1},
	{ 35, "EVLA"    , RHDR, 0, 0, NO, "" ,-1},
	{ 36, "EVLO"    , RHDR, 0, 0, NO, "" ,-1},
	{ 37, "EVEL"    , RHDR, 0, 0, NO, "" ,-1},
	{ 38, "EVDP"    , RHDR, 0, 0, NO, "" ,-1},
	{ 39, "MAG"     , RHDR, 0, 0, NO, "" ,-1},
	{ 40, "USER0"   , RHDR, 0, 0, NO, "" ,-1},
	{ 41, "USER1"   , RHDR, 0, 0, NO, "" ,-1},
	{ 42, "USER2"   , RHDR, 0, 0, NO, "" ,-1},
	{ 43, "USER3"   , RHDR, 0, 0, NO, "" ,-1},
	{ 44, "USER4"   , RHDR, 0, 0, NO, "" ,-1},
	{ 45, "USER5"   , RHDR, 0, 0, NO, "" ,-1},
	{ 46, "USER6"   , RHDR, 0, 0, NO, "" ,-1},
	{ 47, "USER7"   , RHDR, 0, 0, NO, "" ,-1},
	{ 48, "USER8"   , RHDR, 0, 0, NO, "" ,-1},
	{ 49, "USER9"   , RHDR, 0, 0, NO, "" ,-1},
	{ 54, "SB"      , RHDR, 0, 0, NO, "" ,-1},
	{ 55, "SDELTA"  , RHDR, 0, 0, NO, "" ,-1},
	{ 57, "CMPAZ"   , RHDR, 0, 0, NO, "" ,-1},
	{ 58, "CMPINC"  , RHDR, 0, 0, NO, "" ,-1},
	{ 59, "XMINIMUM", RHDR, 0, 0, NO, "" ,-1},
	{ 60, "XMAXIMUM", RHDR, 0, 0, NO, "" ,-1},
	{ 61, "YMINIMUM", RHDR, 0, 0, NO, "" ,-1},
	{ 62, "YMAXIMUM", RHDR, 0, 0, NO, "" ,-1},
	{ 63, "ADJTM"   , RHDR, 0, 0, NO, "" ,-1},
	{ 64, "FHDR65"  , RHDR, 0, 0, NO, "" ,-1},
	{ 65, "FHDR66"  , RHDR, 0, 0, NO, "" ,-1},
	{ 66, "FHDR67"  , RHDR, 0, 0, NO, "" ,-1},
	{ 67, "FHDR68"  , RHDR, 0, 0, NO, "" ,-1},
	{ 68, "FHDR69"  , RHDR, 0, 0, NO, "" ,-1},
	{ 69, "FHDR70"  , RHDR, 0, 0, NO, "" ,-1},
	{  6, "NVHDR"   , IHDR, 0, 0, NO, "" ,-1},
	{  7, "NINF"    , IHDR, 0, 0, NO, "" ,-1},
	{  8, "NHST"    , IHDR, 0, 0, NO, "" ,-1},
	{ 10, "NSNPTS"  , IHDR, 0, 0, NO, "" ,-1},
	{ 11, "NSN"     , IHDR, 0, 0, NO, "" ,-1},
	{ 12, "NXSIZE"  , IHDR, 0, 0, NO, "" ,-1},
	{ 13, "NYSIZE"  , IHDR, 0, 0, NO, "" ,-1},
	{ 14, "NHDR15"  , IHDR, 0, 0, NO, "" ,-1},
	{ 16, "IDEP"    , EHDR, 0, 0, NO, "" ,-1},
	{ 17, "IZTYPE"  , EHDR, 0, 0, NO, "" ,-1},
	{ 18, "IHDR4"   , EHDR, 0, 0, NO, "" ,-1},
	{ 19, "IINST"   , EHDR, 0, 0, NO, "" ,-1},
	{ 20, "ISTREG"  , EHDR, 0, 0, NO, "" ,-1},
	{ 21, "IEVREG"  , EHDR, 0, 0, NO, "" ,-1},
	{ 22, "IEVTYP"  , EHDR, 0, 0, NO, "" ,-1},
	{ 23, "IQUAL"   , EHDR, 0, 0, NO, "" ,-1},
	{ 24, "ISYNTH"  , EHDR, 0, 0, NO, "" ,-1},
	{ 25, "IHDR11"  , IHDR, 0, 0, NO, "" ,-1},
	{ 26, "IHDR12"  , IHDR, 0, 0, NO, "" ,-1},
	{ 27, "IHDR13"  , IHDR, 0, 0, NO, "" ,-1},
	{ 28, "IHDR14"  , IHDR, 0, 0, NO, "" ,-1},
	{ 29, "IHDR15"  , IHDR, 0, 0, NO, "" ,-1},
	{ 30, "IHDR16"  , IHDR, 0, 0, NO, "" ,-1},
	{ 31, "IHDR17"  , IHDR, 0, 0, NO, "" ,-1},
	{ 32, "IHDR18"  , IHDR, 0, 0, NO, "" ,-1},
	{ 33, "IHDR19"  , IHDR, 0, 0, NO, "" ,-1},
	{ 34, "IHDR20"  , IHDR, 0, 0, NO, "" ,-1},
	{ 36, "LPSPOL"  , LHDR, 0, 0, NO, "" ,-1},
	{ 37, "LOVROK"  , LHDR, 0, 0, NO, "" ,-1},
	{ 38, "LCALDA"  , LHDR, 0, 0, NO, "" ,-1},
	{ 39, "LHDR5"   , LHDR, 0, 0, NO, "" ,-1},
	{  1, "KEVNM"   , CHDR, 0, 0, NO, "" ,-1},
	{  2, "KEVNMC"  , CHDR, 0, 0, NO, "" ,-1},
	{  3, "KHOLE"   , CHDR, 0, 0, NO, "" ,-1},
	{  4, "KO"      , CHDR, 0, 0, NO, "" ,-1},
	{  5, "KA"      , CHDR, 0, 0, NO, "" ,-1},
	{  6, "KT0"     , CHDR, 0, 0, NO, "" ,-1},
	{  7, "KT1"     , CHDR, 0, 0, NO, "" ,-1},
	{  8, "KT2"     , CHDR, 0, 0, NO, "" ,-1},
	{  9, "KT3"     , CHDR, 0, 0, NO, "" ,-1},
	{ 10, "KT4"     , CHDR, 0, 0, NO, "" ,-1},
	{ 11, "KT5"     , CHDR, 0, 0, NO, "" ,-1},
	{ 12, "KT6"     , CHDR, 0, 0, NO, "" ,-1},
	{ 13, "KT7"     , CHDR, 0, 0, NO, "" ,-1},
	{ 14, "KT8"     , CHDR, 0, 0, NO, "" ,-1},
	{ 15, "KT9"     , CHDR, 0, 0, NO, "" ,-1},
	{ 16, "KF"      , CHDR, 0, 0, NO, "" ,-1},
	{ 17, "KUSER0"  , CHDR, 0, 0, NO, "" ,-1},
	{ 18, "KUSER1"  , CHDR, 0, 0, NO, "" ,-1},
	{ 19, "KUSER2"  , CHDR, 0, 0, NO, "" ,-1},
	{ 21, "KNETWK"  , CHDR, 0, 0, NO, "" ,-1},
	{ 22, "KDATRD"  , CHDR, 0, 0, NO, "" ,-1},
	{ 23, "KINST"   , CHDR, 0, 0, NO, "" ,-1},
	{  0, ""	, IHDR,  0, 0, NO, "",-1}
};




void gsac_set_param_sort(int ncmd, char **cmdstr)
{
	int i  ;
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, starg, NO, YES))
	       return	;
	for(i=0 ; starg[i].key[0] != '\0' ; i++){
		if(starg[i].used > 0){
			if(starg[i].id >=  0 )	{
			/* get what item to sort */
				sort_which = starg[i].id;
				sort_type = starg[i].ricell;
				isort = i;
				sort_do = YES;
			} else {
				/* get sort control */
				if(starg[i].ricell == CHDR){
					switch(starg[i].id){
						case ST_REVERSE:
							sort_reverse = YES;
							break;
						case ST_FORWARD:
							sort_reverse = NO;
							break;
						case ST_DFLT:
							sort_do = NO ;
							
					}

				}
			}
		}
	}
}

void gsac_exec_sort(void)
{
	int k, ntrchdr;
	/* if there are no traces return */
	ntrchdr = gsac_control.number_iheaders ;
	if(ntrchdr < 1)
		return;
	for(k=0 ; k < ntrchdr ; k++){
		switch(sort_type){
			case ST_DFLT:
				sortfloat[k] = k;
				sortptr[k]= k;
				break;
			case RHDR:
				sortfloat[k] = sacdata[k].sachdr.rhdr[sort_which];
				sortptr[k]= k;
				break;
			case IHDR:
				sortfloat[k] = sacdata[k].sachdr.ihdr[sort_which];
				sortptr[k]= k;
				break;
			/* eventually
			case CHDR:
			 * */
		}
	}
	if(sort_do == YES){
		if(sort_reverse == YES)
			printf("Sorting on %s in descending order\n",
						starg[isort].key);
		else if(sort_reverse == NO)
			printf("Sorting on %s in ascending order \n",
						starg[isort].key);
		gnomesort(ntrchdr, sortfloat, sortptr) ;
	} else {
		printf("Sort is OFF\n");
	}
	if(sort_reverse == YES)
		intrev(sortptr, ntrchdr);
}

void intrev(int *x, int n)
{
	/* reverse a time series */
	int i;
	int tmp;
	if(n <= 1)
		return;
	for(i=0; i<= (n-1)/2 ; i++){
		tmp = x[i];
		x[i] = x[n-1-i];
		x[n-1-i]  = tmp;
	}
}

/* Gnome Sort - The Simplest Sort Algorithm
 * http://www.cs.vu.nl/~dick/gnomesort.html
 * */
void gnomesort(int n, float ar[], int key[]) {
	int i;
	float tmp;
	int   kmp;
	i = 1;
	while (i < n) {
		if( i == 0){
		       i++;
		} else if( ar[i-1] <= ar[i]) {
			i++;
		}
		else {
			tmp = ar[i]; ar[i] = ar[i-1]; ar[i-1] = tmp;
			kmp = key[i]; key[i] = key[i-1]; key[i-1] = kmp;
			i--;
		}
	}
}

void gnomesortint(int n, int ar[], int key[]) {
	int i;
	int tmp;
	int kmp;
	i = 1;
	while (i < n) {
		if( i == 0){
		       i++;
		} else if( ar[i-1] <= ar[i]) {
			i++;
		}
		else {
			tmp = ar[i]; ar[i] = ar[i-1]; ar[i-1] = tmp;
			kmp = key[i]; key[i] = key[i-1]; key[i-1] = kmp;
			i--;
		}
	}
}

