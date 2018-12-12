/* Changes:
	20 FEB 2007 - variable zero padding for greater frequency
		resolution
*/
#include	<stdio.h>
#include	<stdlib.h>
#include	"gsac.h"
#include	"gsac_docommand.h"
#include	"gsac_sac.h"
#include        "gsac_arg.h"

extern struct sacfile_ *sacdata;

#define FFT_DFLT 0
#define FFT_LENG 1

static int fft_length = 1;

struct arghdr fftarg[] = {
        {FFT_DFLT, "DEFAULT", IHDR, 0, 0, NO, "",  1},
        {FFT_LENG, "LENGTH" , IHDR, 0, 1, NO, "LENGTH 1|2|4 ", 1},
        {0,     ""              , IHDR, 0, 0, NO, "", -1}
};


/* these are temporary variables only used here */
float fft_real[10];
int   fft_int [10];
int   fft_yn;
int   fft_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_fft(int ncmd, char **cmdstr)
{
	int i;
	/* parsing code here */


	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, fftarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; fftarg[i].key[0] != '\0' ; i++){
		if(fftarg[i].used > 0){
			if(fftarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, fftarg[i].key, 
					fftarg[i].mfit,fftarg[i].narg, fft_int );
			}
			switch(fftarg[i].id){
				case FFT_DFLT:
					fft_length = 1;
					break;
				case FFT_LENG:
					if(fft_int[0] == 1) {
						fft_length = 1;
					} else if(fft_int[0] == 2) {
						fft_length = 2;
					} else if(fft_int[0] == 4) {
						fft_length = 4;
					} else {
						fft_length = 1;
						printf("FFT Usage: Length 1|2|4 - using 1\n");
					}
					break;

			}
		}
	}
	return;
}

void gsac_exec_fft(void)
{
	/* for each trace, determine power of two, 
	 * reallocoate memory for the spectra
	 * compute the fft in place
	 * */
	int i, j, k, ntrc, npts, n2;
	float dt, df;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
	gsac_control.fft = YES;
		
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		dt = sacdata[k].sachdr.rhdr[H_DELTA];
		if(npts > 0){
			n2 = fft_length*npow2(npts);
			/* allocate space for the spectra array */
			sacdata[k].sac_spectra = 
				(float *)realloc(sacdata[k].sac_spectra,
						 (n2+n2)*sizeof(float));
			/* fill the array with the time series */
			/* The even elements are the time series, the odd are 
			 * the imaginary part of the time series which is
			 * zero. The trace is effectively zero filled at 
			 * the end */
			for(i=0, j=0; i < n2 ; i++){
				/* fill in real part of time series */
				if(i<npts)
					sacdata[k].sac_spectra[j++] = 
						sacdata[k].sac_data[i];
				else
					sacdata[k].sac_spectra[j++] = 0.0;
				/* fill in imaginary part of time series */
				sacdata[k].sac_spectra[j++] = 0.0;
			}
			/* get the DFT - note we must define df as something
			 * I define df = 0 to force computation of df from
			 * dt and n2 */
			df = 0.0;
			four(sacdata[k].sac_spectra, n2, -1, &dt, &df);
			sacdata[k].df = df;
			sacdata[k].npow2 = n2;
		}
	}
}


