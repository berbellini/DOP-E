#include	<stdio.h>
#include	<sys/types.h>
#include	<sys/stat.h>
#include	<unistd.h>
#include	<string.h>

#ifdef WIN32
#include	"glob.h"
#else
#include	<glob.h>
#endif
#include	"gsac.h"

/* safety for old cc compiler suite */
#ifndef GLOB_TILDE
#define  GLOB_TILDE      0x0000  /* Expand tilde names from the passwd file. */
#endif
#include "gsac_docommand.h"
#include "gsac_sac.h"
#include "csstim.h"

struct sacfile_ *sacdata ;
int *sortptr;
float *sortfloat;


glob_t globbuf;
struct stat statbuf;

void gsac_set_param_rh(int ncmd, char **cmdstr)
{
	/* the only parameter to be set is MORE 
	 *
	 */
	int i,j,k;
	int oldmax;
	int DOMORE ;
	char ctmp[1000];
	/* is command is by itself there is nothing to parse
	 * This in fact means that the execution of the command
	 * uses the previous parameters
	 */
	if(ncmd == 1) return;
	/* parse for the phrase MORE */
	DOMORE = 0 ;
	for (i=1 ; i < ncmd ; i++){
		strcpy(ctmp,cmdstr[i]);
		if(strcmp(gsac_strupr(ctmp),"MORE") == 0){
			DOMORE = 1;
			/* we need more arguments to read things in */
			if(ncmd == 2) return;
		}
	}
	oldmax = gsac_control.max_number_traces ;
	if(!DOMORE)
		gsac_control.max_number_traces = 0;
	for(i=1; i < ncmd; i++){
		glob(cmdstr[i],  GLOB_TILDE, NULL, &globbuf);
		for (j = 0; j < globbuf.gl_pathc; j++){
			/* test to see if it is a file - if so proceed */
			stat(globbuf.gl_pathv[j], &statbuf);
			if( S_ISREG(statbuf.st_mode)){
				/* TEST TO SEE IF SAC FILE */
				if (gsac_valid_sacfile(globbuf.gl_pathv[j]) > 0){

				/* allocate a data structure */
				gsac_alloc_trace(oldmax);

				/* now get the data */
				k = gsac_control.max_number_traces -1;
				/* beware of size limitations */
				strcpy(sacdata[k].sac_ifile_name , globbuf.gl_pathv[j]);
				}
			}
		}
	}

}


void gsac_exec_rh(void)
{
	int i,j,iret,k;
	int npts;
	float delta;
	float evla, evlo, stla, stlo;
	float gcarc, az, baz, dist;
	float fmax;

	gsac_control.begmin= 1.0e+30;
	gsac_control.endmax=-1.0e+30;
	gsac_control.fft = NO;
	/* also the trace file is known to be a valid SAC file */
	/* force the state of no traces in memory */
	gsac_control.number_itraces = 0;
	for(k=0, gsac_control.number_iheaders=0, gsac_control.number_otraces=0 ; k < gsac_control.max_number_traces ; k++){
		printf("%s ",sacdata[k].sac_ifile_name);
		strcpy(sacdata[k].sac_ofile_name , sacdata[k].sac_ifile_name);
		/* redefine the output file name */
		sacdata[k].display = NO;
		iret = brsach(sacdata[k].sac_ifile_name,&sacdata[k].sachdr);
		if(iret >=0){
			gsac_control.number_iheaders++;
			/* set epoch time */
			if(sacdata[k].sachdr.ihdr[H_NZYEAR] == -12345 ||
				sacdata[k].sachdr.ihdr[H_NZJDAY] == -12345 ||
				sacdata[k].sachdr.ihdr[H_NZHOUR] == -12345 ||
				sacdata[k].sachdr.ihdr[H_NZMIN] == -12345 ||
				sacdata[k].sachdr.ihdr[H_NZSEC] == -12345 ||
				sacdata[k].sachdr.ihdr[H_NZMSEC] == -12345){
				printf("Bad Reference Time Field: %d %d %d %d %d %d\n"
					,sacdata[k].sachdr.ihdr[H_NZYEAR],
                                sacdata[k].sachdr.ihdr[H_NZJDAY],
                                sacdata[k].sachdr.ihdr[H_NZHOUR],
                                sacdata[k].sachdr.ihdr[H_NZMIN],
                                sacdata[k].sachdr.ihdr[H_NZSEC],
                                sacdata[k].sachdr.ihdr[H_NZMSEC]);
				printf("Using 1970 000 00 00 00 000\n");
				sacdata[k].tzref = 0.0;
			} else {
			htoe1(sacdata[k].sachdr.ihdr[H_NZYEAR], 
				sacdata[k].sachdr.ihdr[H_NZJDAY], 
				sacdata[k].sachdr.ihdr[H_NZHOUR],
				sacdata[k].sachdr.ihdr[H_NZMIN],
				sacdata[k].sachdr.ihdr[H_NZSEC],
				sacdata[k].sachdr.ihdr[H_NZMSEC],
				&sacdata[k].tzref);
			}
			/* make a clean sac header string! */
			for(i=0;i<24;i++){
			for(j=0;j<8;j++){
			strncpy(sacdata[k].schdr[i],sacdata[k].sachdr.chdr[i],8);
			}
			sacdata[k].schdr[i][8] = '\0';
			}
			npts = sacdata[k].sachdr.ihdr[H_NPTS];
			delta = sacdata[k].sachdr.rhdr[H_DELTA];
			sacdata[k].display = YES;
			/* preserve the permint permax */
			sacdata[k].permin = 
				sacdata[k].sachdr.rhdr[H_USER1];
			sacdata[k].permax = 
				sacdata[k].sachdr.rhdr[H_USER2];
			/* fix distance, azimuth, etc if not set 
			 * if LCALDA == TRUE */
			if(sacdata[k].sachdr.ihdr[H_LCALDA] !=0){
				stla = sacdata[k].sachdr.rhdr[H_STLA];
				stlo = sacdata[k].sachdr.rhdr[H_STLO];
				evla = sacdata[k].sachdr.rhdr[H_EVLA];
				evlo = sacdata[k].sachdr.rhdr[H_EVLO];
				az   = sacdata[k].sachdr.rhdr[H_AZ];
				baz  = sacdata[k].sachdr.rhdr[H_BAZ];
				dist  = sacdata[k].sachdr.rhdr[H_DIST];
				gcarc = sacdata[k].sachdr.rhdr[H_GCARC];
				if(stla != -12345. && stlo != -12345.
				&& evla != -12345. && evlo != -12345.
				&& (dist == -12345. || az   == -12345.
				|| baz == -12345. || gcarc == -12345.) ){
				delaz( evla,  evlo,  stla,  stlo,  &gcarc,  &az,  &baz,  &dist);
				sacdata[k].sachdr.rhdr[H_DIST] = dist;
				sacdata[k].sachdr.rhdr[H_AZ] = az;
				sacdata[k].sachdr.rhdr[H_BAZ] = baz;
				sacdata[k].sachdr.rhdr[H_GCARC] = gcarc;
				}
			}
			/* check to see if E is defined */
			if(sacdata[k].sachdr.rhdr[H_E] == -12345.){
				sacdata[k].sachdr.rhdr[H_E] = 
					sacdata[k].sachdr.rhdr[H_B] 
					 + (npts -1)*delta;
			}
			/* ADDED 03 AUG 2005  fixed 03 JUL 2006 */
			/* check to see if E is reasonable */
			if (sacdata[k].sachdr.rhdr[H_E] != 
				sacdata[k].sachdr.rhdr[H_B] +  (npts -1)*delta){
				 sacdata[k].sachdr.rhdr[H_E] =
                                        sacdata[k].sachdr.rhdr[H_B]
                                         + (npts -1)*delta;
			}

			/* set epoch times of first and last sample */
/*
			added 18 APR 2006 to permit picking of frequency
			for IXY plot of spectrum as a time series
*/
			if(sacdata[k].sachdr.ihdr[H_IFTYPE] == ENUM_IXY){
				sacdata[k].tzref = 0.0 ;
				sacdata[k].tzbeg = sacdata[k].sachdr.rhdr[H_B];
				sacdata[k].tzend = sacdata[k].sachdr.rhdr[H_E];
			} else {

				sacdata[k].tzbeg = sacdata[k].tzref + sacdata[k].sachdr.rhdr[H_B];
				sacdata[k].tzend = sacdata[k].tzref + sacdata[k].sachdr.rhdr[H_E];
			}
			sacdata[k].tzbegx = sacdata[k].tzbeg;
			sacdata[k].tzendx = sacdata[k].tzend;
			/* get bounds for absolute plotting */
			if(sacdata[k].tzbeg < gsac_control.begmin)
				gsac_control.begmin = sacdata[k].tzbeg;
			if(sacdata[k].tzend > gsac_control.endmax)
				gsac_control.endmax = sacdata[k].tzend;
			fmax = 0.5/sacdata[k].sachdr.rhdr[H_DELTA];
			if(fmax > gsac_control.fmax)
				gsac_control.fmax = fmax;
			/* same the component name */
			strcpy(sacdata[k].ocmpnm,sacdata[k].schdr[H_KCMPNM]);
		}
	}
	printf("\n");
}

