#ifndef _GSAC_SAC
#define _GSAC_SAC
#include        <stdio.h>
#include        <stdlib.h>
#include        <string.h>
#include        <errno.h>
#include <ctype.h>
struct sachdr_ sachdr;
#ifdef MSDOS
#include <fcntl.h>
#endif


/* this header file contains everything required about the 
 * SAC files in memory. the struct sacfile_ knows about the 
 * external file and also the internal parameters
 */

#define INT int

/* headers */

struct sachdr_  {
	float rhdr[70];
	INT ihdr[40];
	char chdr[24][8];
	}  ;


struct sacfile_ {
	char sac_ifile_name[1000];
	char sac_ofile_name[1000];
	int  sac_file_swab      ;
	struct sachdr_ sachdr   ;
	float *sac_data         ;
	float *sac_spectra      ;
	float df		;	/* df for spectra */
	int npow2		;	/* for spectra and FFT */
	char schdr[24][9]	;
	char ocmpnm[9];
	double tzref		;	/* absolute reference time of trace */
	double tzbeg		;	/* absolute begin time of trace */
	double tzend		;	/* absolute end   time of trace */
	int display		;	/* implementation of a delete */
	double tzbegx		;	/* absolute plot begin time of trace */
	double tzendx		;	/* absolute plot end   time of trace */
					/* these two are required for xlim */
	float permin		;	/* used to set header on a write */
	float permax		;       /* w command so that these are not set
					with a wh command */
	float winmax		;	/* maximum amplitude of a trace segment
						plotted - used in gsac_prs */
	float winmin		;
};



/* function prototypes */

int brsac(char *fname,struct  sachdr_ *sachdr, float **data);
int brsach(char *fname,struct  sachdr_ *sachdrret);
int bwsac(char *fname,struct  sachdr_ sachdr, float *data);
int bwsach(char *fname,struct  sachdr_ sachdr);
int gsac_valid_sacfile(char *name);
void gsac_alloc_trace(int oldmax );


#endif
