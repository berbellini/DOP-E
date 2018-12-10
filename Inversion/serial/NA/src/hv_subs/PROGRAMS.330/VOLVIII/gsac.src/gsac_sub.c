#include	<stdio.h>
#include	"gsac.h"
#include	"gsac_sac.h"
#include "gsac_docommand.h"
#include "gsac_arg.h"

#define	SUB_SUB 1

struct arghdr subarg[] = {
	{SUB_SUB, "SUB", RHDR, 0, 1, NO, "SUB value",-1},
	{0,     ""     , IHDR, 0, 0, NO, "",-1}
};

/* these are temporary variables only used here */
float sub_real[1];
float sub_value = 0.0;

extern struct sacfile_ *sacdata;

void gsac_set_param_sub(int ncmd, char **cmdstr)
{
	float tmp;
	if (ncmd == 1)
		return;
	if(isargr(cmdstr[1], &tmp) == 1)
		sub_value = tmp;
}

void gsac_exec_sub(void)
{

	int i, k, ntrc, npts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
		
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		if(npts > 0){
			for(i=0; i < npts ; i++)
				sacdata[k].sac_data[i] -= sub_value;
			getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
			sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
			sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
			sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		}
	}
}


void fillstr1(char *instr, char *str1){
/* fill up str1 with instr but no termination '\0' */
	int ls, i;
	ls=strlen(instr);
	for (i=0 ; i < 8 ; i++)
		str1[i] = ' ';
	for(i=0 ; i < 8 && i < ls ; i ++)
		str1[i] = instr[i];
	str1[8] = '\0';
}
void fillstr12(char *instr, char *str1, char *str2){
/* fill up str1 and str2 with instr but no termination '\0' */
	int ls, i;
	ls=strlen(instr);
	for (i=0 ; i < 8 ; i++){
		str1[i] = ' ';
		str2[i] = ' ';
	}
	for(i=0 ; i < 8 && i < ls ; i ++)
		str1[i] = instr[i];
	for(i=8 ; i < 16 && i < ls ; i ++)
		str2[i-8] = instr[i];
	str1[8] = '\0';
	str2[8] = '\0';
}

int npow2(int n)
{
	int k, nn;
	k = 1; nn = 2;
	while(nn < n){
		nn *=2;
	}
	return(nn);

}
void four(float data[], int nn, int isign, float *dt, float *df){
/*
	subroutine four(data,nn,isign,dt,df)
c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(data(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  data is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     data are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array data, replacing the input data.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c-----
      real*4 data(*)
*/
	int n;
	int i, j, m, mmax, iiii, istep;
	float tempr, tempi;
	float wr, wi;
	float wstpr, wstpi;
	float sinth, theta;
	float dtt, dff;

	dtt = (*dt);
	dff = (*df);

	n = 2 * nn;
	if((dtt) == 0.0) dtt = 1./(nn*(dff)) ;
	if((dff) == 0.0) dff = 1./(nn*(dtt)) ;
	if((dtt) != (nn*(dff))) dff = 1./(nn*(dtt)) ;
	*dt = dtt;
	*df = dff;
	j = 1;
	for (i=1;i<= n ; i+=2){
		if(i < j){
			tempr = data[j-1];
			tempi = data[j  ];
			data[j-1] = data[i-1];
			data[j  ]=data[i  ];
			data[i-1] = tempr;
			data[i  ] = tempi;
		}
		m = n/2;
statement3:		if(j <= m) goto statement4 ;
			j = j-m;
			m = m/2;
			if(m >= 2)goto statement3 ;
statement4:		
		j=j+m;
    	}
	mmax = 2 ;
	while(mmax < n ){
		istep= 2 *mmax;
		theta = 6.283185307/(float)(isign*mmax);
		sinth=sin(theta/2.);
		wstpr=-2.*sinth*sinth;
		wstpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1; m <= mmax ; m +=2){
				for(i = m ; i <= n ; i+=istep){
					j=i+mmax;
					tempr=wr*data[j-1]-wi*data[j  ];
					tempi=wr*data[j  ]+wi*data[j-1];
					data[j-1]=data[i-1]-tempr;
					data[j  ]=data[i  ]-tempi;
					data[i-1]=data[i-1]+tempr;
					data[i  ]=data[i  ]+tempi;
				}
				tempr = wr;
				wr = wr*wstpr-wi*wstpi + wr;
				wi = wi*wstpr+tempr*wstpi + wi;
		}
		mmax = istep;
	}
 /*
  * 	get the correct dimensions for a Fourier Transform 
  * 	from the Discrete Fourier Transform
*/ 
	if(isign > 0){
		/*
		frequency to time domain
		*/
		for (iiii= 0 ; iiii < n ; iiii++){
			data[iiii] = data[iiii] * (*df);
		}
	} else {
		/*
		time to frequency domain
		*/
		for (iiii= 0 ; iiii < n ; iiii++){
			data[iiii] = data[iiii] * (*dt);
		}
 }
}
