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
extern void gnomesortint(int n, int ar[], int key[] ) ;


#define	FILEID_DEFAULT	-1
#define	FILEID_ON	-2
#define	FILEID_OFF	-3
#define	FILEID_LOCATION	-4
#define	FILEID_FORMAT	-5
#define	FILEID_TYPE	-6
#define	FILEID_NAME	-7
#define	FILEID_LIST	-8
#define FILEID_DUMMY	-9
#define FILEID_UR	-10
#define FILEID_UL	-11
#define FILEID_LR	-12
#define FILEID_LL	-13
#define FILEID_UC	-14
#define FILEID_LC	-15
#define FILEID_EQUALS	-16
#define FILEID_COLONS	-17
#define FILEID_NONAMES	-18
#define KZDATE  100
#define KZTIME  101
#define KEVNM   102


struct arghdr fileidarg[] = {
	{FILEID_DEFAULT,"DEFAULT", IHDR, NO, 0, NO, "", 2},
	{FILEID_ON,	"ON"	, IHDR, 0, 0, NO, "",2},
	{FILEID_OFF,	"OFF"	, IHDR, 0, 0, NO, "",2},
	{FILEID_LIST,	"LIST"	, IHDR, 0, 0, NO, "",2},
	{FILEID_FORMAT, "FORMAT", CHDR, 0, 1, NO, "Format EQuals|Colons|NOnames",1},
	{FILEID_TYPE,	"TYPE"	, IHDR, 0, 0, NO, "",2},
	{FILEID_NAME,	"NAME"	, IHDR, 0, 0, NO, "", 1},
	{FILEID_LOCATION,"LOCATION", CHDR, 0, 1, NO, "LOcation UL|UC|UR|LL|LC|LR", 2},
	{ KZDATE, "KZDATE"  , IHDR, YES, 0, YES, "" ,-1},
	{ KZTIME, "KZTIME"  , IHDR, YES, 0, YES, "" ,-1},
	{ KEVNM,  "KEVNM"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_STLA, "STLA"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_STLO, "STLO"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_STEL, "STEL"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_EVLA, "EVLA"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_EVLO, "EVLO"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_EVDP, "EVDP"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_DIST, "DIST"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_AZ  , "AZ"      , IHDR, YES, 0, YES, "" ,-1},
	{ H_BAZ , "BAZ"     , IHDR, YES, 0, YES, "" ,-1},
	{ H_GCARC, "GCARC"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_KSTNM, "KSTNM"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_KCMPNM, "KCMPNM"  , IHDR, YES, 0, YES, "" ,-1},

	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float fileid_real[10];
int   fileid_int [10];
int   fileid_yn;
int   fileid_num;

/* these are prototypes for global variables to be used by the routine */
#define FILEID_MAXLIST 10
int fileid_use = YES;
int fileid_location = FILEID_UR ;
int fileid_format = FILEID_NONAMES;
int fileid_type = FILEID_LIST;
int fileid_n = 2;
int fileid_list[FILEID_MAXLIST] =
	{ 0, 20, -1, -1, -1, -1, -1, -1, -1, -1
	};	/* this is the actual list to be displayed */
char fileid_clist[FILEID_MAXLIST][8] ={
	"KSTNM", "KCMPNM", "", "", "", "", "", "", "", ""
	};	/* this is the list of title strings to be displayed */
int fileid_ar[FILEID_MAXLIST] = 
	{ 0, 1, -1, -1, -1, -1, -1, -1, -1, -1
	};	/* this is used for the sort to ensure the list is in order given */
int fileid_key[FILEID_MAXLIST] =
	{ 0, 1, -1, -1, -1, -1, -1, -1, -1, -1
	}; /* this created by gnomesortint */
char fileid_ostr[FILEID_MAXLIST][30];
	

void gsac_set_param_fileid(int ncmd, char **cmdstr)
{
	int i;
	char instr[100];
	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, fileidarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; fileidarg[i].key[0] != '\0' ; i++){
		if(fileidarg[i].used > 0){
			if(fileidarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, fileidarg[i].key, 
					fileidarg[i].mfit,fileidarg[i].narg, fileid_real);
			} else if(fileidarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, fileidarg[i].key, 
					fileidarg[i].mfit,fileidarg[i].narg, fileid_int );
			} else if(fileidarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, fileidarg[i].key, 
					fileidarg[i].mfit,fileidarg[i].narg, &fileid_yn );
			} else if(fileidarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, fileidarg[i].key, 
					fileidarg[i].mfit,fileidarg[i].narg, &fileid_num );
			} else if(fileidarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, fileidarg[i].key, 
					fileidarg[i].mfit, fileidarg[i].narg, instr );
			}
			switch(fileidarg[i].id){
				case FILEID_LIST:
					fileid_n = 0;
					fileid_type = FILEID_LIST;
				case FILEID_ON:
					fileid_use = YES;
					break;
				case FILEID_OFF:
					fileid_use = NO;
					break;
				case FILEID_DEFAULT:
					fileid_use = YES;
					fileid_location = FILEID_UR ;
					fileid_format = FILEID_NONAMES;
					fileid_type = FILEID_DEFAULT;
					fileid_n = 4;
					fileid_list[0]=H_KSTNM;
					strcpy(fileid_clist[0],"KSTNM");
					fileid_list[1]=H_KCMPNM;
					strcpy(fileid_clist[1],"KCMPNM");
					fileid_list[2]=KZDATE;
					strcpy(fileid_clist[2],"KZDATE");
					fileid_list[3]=KZTIME;
					strcpy(fileid_clist[3],"KZTIME");
					break;
				case FILEID_FORMAT:
					gsac_strupr(instr);
					if(strncmp(instr,"CO",2)==0)
						fileid_format = FILEID_COLONS;
					else if(strncmp(instr,"EQ",2)==0)
						fileid_format = FILEID_EQUALS;
					else if(strncmp(instr,"NO",2)==0)
						fileid_format = FILEID_NONAMES;
					break;
				case FILEID_LOCATION:
					gsac_strupr(instr);
printf("fileid_location %s\n",instr);
					if(strncmp(instr,"UL",2)==0)
						fileid_location = FILEID_UL;
					else if(strncmp(instr,"UR",2)==0)
						fileid_location = FILEID_UR;
					else if(strncmp(instr,"UC",2)==0)
						fileid_location = FILEID_UC;
					else if(strncmp(instr,"LL",2)==0)
						fileid_location = FILEID_LL;
					else if(strncmp(instr,"LR",2)==0)
						fileid_location = FILEID_LR;
					else if(strncmp(instr,"LC",2)==0)
						fileid_location = FILEID_LC;
					break;
				case FILEID_NAME:
					fileid_type = FILEID_NAME;
					break;
				case KZDATE:
				case KZTIME:
				case H_KCMPNM:
				case H_KSTNM:
				case KEVNM:
				case H_GCARC:
				case H_DIST:
				case H_AZ:
				case H_BAZ:
				case H_STLA:
				case H_STLO:
				case H_STEL:
				case H_EVLA:
				case H_EVLO:
				case H_EVDP:
					/* safety so never exceed array dimension  */
					if(fileid_n < FILEID_MAXLIST  ){
					fileid_list[fileid_n]=fileidarg[i].id;;
					/* now annotate the labeling string */
					switch(fileidarg[i].id){
						case KZDATE:
						strcpy(fileid_clist[fileid_n],"KZDATE");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"KZDATE");
						break;
						case KZTIME:
						strcpy(fileid_clist[fileid_n],"KZTIME");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"KZTIME");
						break;
						case H_KCMPNM:
						strcpy(fileid_clist[fileid_n],"KCMPNM");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"KCMPNM");
						break;
						case H_KSTNM:
						strcpy(fileid_clist[fileid_n],"KSTNM");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"KSTNM");
						break;
						case KEVNM:
						strcpy(fileid_clist[fileid_n],"KEVNM");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"KEVNM");
						break;
						case H_GCARC:
						strcpy(fileid_clist[fileid_n],"GCARC");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"GCARC");
						break;
						case H_DIST:
						strcpy(fileid_clist[fileid_n],"DIST");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"DIST");
						break;
						case H_AZ:
						strcpy(fileid_clist[fileid_n],"AZ");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"AZ");
						break;
						case H_BAZ:
						strcpy(fileid_clist[fileid_n],"BAZ");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"BAZ");
						break;
						case H_STLA:
						strcpy(fileid_clist[fileid_n],"STLA");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"STLA");
						break;
						case H_STLO:
						strcpy(fileid_clist[fileid_n],"STLO");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"STLO");
						break;
						case H_STEL:
						strcpy(fileid_clist[fileid_n],"STEL");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"STEL");
						break;
						case H_EVLA:
						strcpy(fileid_clist[fileid_n],"EVLA");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"EVLA");
						break;
						case H_EVLO:
						strcpy(fileid_clist[fileid_n],"EVLO");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"EVLO");
						break;
						case H_EVDP:
						strcpy(fileid_clist[fileid_n],"EVDP");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"EVDP");
						break;
					}
					fileid_key[fileid_n] = fileid_n;
					fileid_n++;
					}
					break;
			}
		}
	}
}

void gsac_exec_fileid(void)
{
}

void gsac_plot_fileid(float x0,float y0,float xlen,float ylen, int k)
{
	/* implement the fileid labeling 
	the bounding box is given by (x0,y0) for the lower left 
		and (x0+xlen,y0+ylen) for the upper right corner
	k is the index of the trace file - note we need a pointer to the
		sacdata structure at the top of this routine
	*/
	int i,j;
	float ht;
	float xp, yp;
	float dy;
	int LCR;	/* Left < 0 , Center = 0 , Right > 0 */
	int clength;	/* maximum length of character string for placement */
	float cw;	/* width of character string*/
	char outstr[30];
	if(fileid_use == YES){
		/* change order of output of used LIST */
		if(fileid_type == FILEID_LIST){
			/* sort the list */
			gnomesortint(fileid_n, fileid_ar, fileid_key);
			ht = MIN(0.4*ylen/fileid_n,0.1);
		} else {
			ht = MIN(0.10*ylen,0.1);
		}
		switch(fileid_location){
			case FILEID_UR:
				xp = x0 + xlen ;
				dy = -1.25*ht;
				yp = y0 + ylen + 2.5*dy;
				LCR = 1;
				break;
			case FILEID_UC:
				xp = x0 + 0.5*xlen;
				dy = -1.25*ht;
				yp = y0 + ylen + 2.5*dy;
				LCR = 0;
				break;
			case FILEID_UL:
				xp = x0 + 2.*ht;
				dy = -1.25*ht;
				yp = y0 + ylen + 2.5*dy;
				LCR = -1;
				break;
			case FILEID_LR:
				xp = x0 + xlen ;
				dy = 1.25*ht;
				yp =  y0 + dy;
				LCR = 1;
				break;
			case FILEID_LC:
				xp = x0 + 0.5*xlen;
				dy = 1.25*ht;
				yp =  y0 + dy;
				LCR = 0;
				break;
			case FILEID_LL:
				xp = x0 + 2.*ht;
				dy = 1.25*ht;
				yp =  y0 + dy;
				LCR = -1;
				break;
		}
		if(fileid_type == FILEID_NAME){
			clength = strlen(sacdata[k].sac_ofile_name);
			cw = clength*ht;
			if(LCR < 0)
				gleft(xp, yp,ht,sacdata[k].sac_ofile_name,0.0);
			else if(LCR == 0)
				gleft(xp-0.5*cw,yp,ht,sacdata[k].sac_ofile_name,0.0);
			else if(LCR > 0)
				gleft(xp -cw-ht-ht,yp,ht,sacdata[k].sac_ofile_name,0.0);
		} else if(fileid_type == FILEID_LIST){
			/* this is complicated by the size of the fields, the
				fact that the fields are numeric or string
				and the formatting - so be patient 
				for SAC compatability the field width depends on the
				numbers - yuck */
			clength = 0;
			for(j = 0 ; j < fileid_n; j++){
				i = fileid_key[j];
				switch(fileid_list[i]){
				case KZDATE:
					printkdatestr(sacdata[k].sachdr.ihdr[0],
                                                        sacdata[k].sachdr.ihdr[1],outstr);
					switch (fileid_format){
						case FILEID_NONAMES:
							sprintf(fileid_ostr[i],"%-18s",outstr);
							break;
						case FILEID_COLONS:
							sprintf(fileid_ostr[i],"%-8s:%-18s",fileid_clist[i],outstr);
							break;
						case FILEID_EQUALS:
							sprintf(fileid_ostr[i],"%8s=%-18s",fileid_clist[i],outstr);
							break;
					}
					break;
				case KZTIME:
					printktimestr(sacdata[k].sachdr.ihdr[2],
                                                        sacdata[k].sachdr.ihdr[3],
                                                        sacdata[k].sachdr.ihdr[4],
                                                        sacdata[k].sachdr.ihdr[5],outstr);
					switch (fileid_format){
						case FILEID_NONAMES:
							sprintf(fileid_ostr[i],"%-12s",outstr);
							break;
						case FILEID_COLONS:
							sprintf(fileid_ostr[i],"%-8s:%-12s",fileid_clist[i],outstr);
							break;
						case FILEID_EQUALS:
							sprintf(fileid_ostr[i],"%8s=%-12s",fileid_clist[i],outstr);
							break;
					}
					break;
				case H_KCMPNM:
				case H_KSTNM:
					switch (fileid_format){
						case FILEID_NONAMES:
							sprintf(fileid_ostr[i],"%-8s",sacdata[k].schdr[fileid_list[i]]);
							break;
						case FILEID_COLONS:
							sprintf(fileid_ostr[i],"%-8s:%-8s",fileid_clist[i],sacdata[k].schdr[fileid_list[i]]);
							break;
						case FILEID_EQUALS:
							sprintf(fileid_ostr[i],"%8s=%-8s",fileid_clist[i],sacdata[k].schdr[fileid_list[i]]);
							break;
					}
					break;
				case KEVNM:
					sprintf(fileid_ostr[i],"%-8s%-8s",
						sacdata[k].schdr[1],
						sacdata[k].schdr[2]);
					break;
				case H_GCARC:
				case H_DIST:
				case H_AZ:
				case H_BAZ:
				case H_STLA:
				case H_STLO:
				case H_STEL:
				case H_EVLA:
				case H_EVLO:
				case H_EVDP:
					switch (fileid_format){
						case FILEID_NONAMES:
							sprintf(fileid_ostr[i],"%-9.3g",sacdata[k].sachdr.rhdr[fileid_list[i]]);
							break;
						case FILEID_COLONS:
							sprintf(fileid_ostr[i],"%-8s:%-9.3g",fileid_clist[i],sacdata[k].sachdr.rhdr[fileid_list[i]]);
							break;
						case FILEID_EQUALS:
							sprintf(fileid_ostr[i],"%8s=%-9.3g",fileid_clist[i],sacdata[k].sachdr.rhdr[fileid_list[i]]);
							break;
					}
					break;
				}
					/* safety so never exceed array dimension  */
				clength=MAX(clength,strlen(fileid_ostr[i]));
			}
			/* now that we have the output strings and the length we
				can properly place them on the plot */
			cw = clength*ht;
			for(j = 0 ; j < fileid_n; j++){
				i = fileid_key[j];
				if(LCR < 0)
					gleft(xp, yp,ht,fileid_ostr[i],0.0);
				else if(LCR == 0)
					gleft(xp-0.5*cw,yp,ht,fileid_ostr[i],0.0);
				else if(LCR > 0)
					gleft(xp -cw-ht-ht,yp,ht,fileid_ostr[i],0.0);
				yp += dy;
			}
		}

	}

}
