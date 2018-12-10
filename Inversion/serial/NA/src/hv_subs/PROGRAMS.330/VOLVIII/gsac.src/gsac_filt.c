#include <stdio.h>
#include <stdlib.h>
#include "gsac.h"
#include "gsac_sac.h"

/* BUG FIXES 
 * 19 JUL 2004 - error in setting USER1 USER2 for hp lp filters
 * 		I had used the frequency not the period
 * 		I also check to ensure that previoous values are not
 * 		unset
 *
 * 11 JAN 2010 - added Chebyshev Type I filter - not that I define a dummy iirap structure
            for the chebyshev type I filter that is filled on setup and never reused
            Digital Filter Designers Handbook with C++ algorithms
            C.  Britton Rorabaugh 2nd Edition
            McGraw Hill, New York, 479 pp 1997  Chapter 5
 
               - also corrected an error in the implementation of the bandreject
                 for a single-pole stage   at line 308 of current code from
                                  } else if (filt_lhp == BANDREJECT){
                                        wc = fwarp(Wc,dt);
                                        wl = fwarp(Wl/ffac,dt);
                                        wu = fwarp(Wu/ffac,dt);
                                        num[0] = wl*wu*dt*dt;
                                        num[1] = 0.0;
                  to
                                  } else if (filt_lhp == BANDREJECT){
                                        wc = fwarp(Wc,dt);
                                        wl = fwarp(Wl/ffac,dt);
                                        wu = fwarp(Wu/ffac,dt);
                                        num[0] = wc*wc*dt*dt;
                                        num[1] = 0.0;

 */


extern struct sacfile_ *sacdata;
void stozc(int n,double fs,double *a);
void getroot(float wl,float wh,float wc,float zeta,float *Wc, float *Zeta, int ndo);
void getchebpole(int n, float eps);


/* for a two pole system, a low pass filter is of the form
 *                      2
 *                     w
 * H(s) =  -------------------------
 *          2                     2
 *         s  + 2 w cos phi  s + w  
 *
 *         where the cos phi = zeta, the damping factor.
 *
 * This form works for Butterworth and Bessel filters. Whereas the
 * Butterworth corner frequency w is the same for all filter stages, the
 * w for each stage of the Bessel filter is different.  The 'f' corner 
 * frequency factor is a number by which the desired 'w' is multiplied in
 * one particular stage. The Butterworth coefficients are computed from
 * first principles; the Bessel coefficients are taken from
 *
 * */
struct iirap {		/* angle and corner frequency of 2 pole */
	int n;
        float gain;     /* filter constant required for Chebyshev-I filter, 1.0 for others 
                           added 11 JAN 2010 */
	float fs  ;	/* corner frequency factor odd order filters
			   note repeated for even order but not used */
	float f[5];	/* corner frequency factor for 2-pole stages */
	float p[5];	/* angle     */
};

/* note that when using the Bessel, the combined filter is adjusted to have
 * a 3 db cutoff at the desired sedign corner of omega = 1. This means that
 * the filter stages are either
 *         1                                        1
 *   ______________             or     __________________________
 *    (s/wb) + 1                       (s/wb)^2 + 2 zeta (s/wb) + 1
 *
 *    The lowpass to lowpass transformation replaces the S == (s/wb)
 *    by (s/wc) whhere wc is the desired corner. or
 *    LP -> LP         s -> (s/wb wc)
 *    LP -> HP         s ->  (wc/s wb)
 *    LP -> BP
 *    */

/* set up a pointer to the coefficient arrays */

/* Bessel coefficients */
struct iirap besfp[] = {
  {0, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{00.000, 00.000, 00.000, 00.000, 00.000 }},
  {1, 1.0, 1.0000, {1.2723, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{30.000, 00.000, 00.000, 00.000, 00.000 }},
  {1, 1.0, 1.3243, {1.4494, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{43.6558, 00.000, 00.000, 00.000, 00.000 }},
  {2, 1.0, 1.3243, {1.6043, 1.4310, 1.0000, 1.0000, 1.0000}, 
	  	{51.6327, 16.6646, 00.000, 00.000, 00.000 }},
  {2, 1.0, 1.5040, {1.7573, 1.5581, 1.0000, 1.0000, 1.0000}, 
	  	{56.9343, 27.4641, 00.000, 00.000, 00.000 }},
  {3, 1.0, 1.5040, {1.9070, 1.6911, 1.6058, 1.0000, 1.0000}, 
	  	{60.7514, 35.1050, 11.5358, 00.000, 00.000 }},
  {3, 1.0, 1.6871, {2.0529, 1.8254, 1.7192, 1.0000, 1.0000}, 
	  	{63.6470, 40.8346, 20.0824, 00.000, 00.000 }},
  {4, 1.0, 1.6871, {2.1914, 1.9556, 1.8343, 1.7806, 1.0000}, 
	  	{65.9270, 45.2996, 26.6836, 8.81062, 00.000 }},
  {4, 1.0, 1.8570, {2.3228, 2.0808, 1.9483, 1.8788, 1.0000}, 
	  	{67.7778, 48.8981, 31.9728, 15.8248, 00.000 }},
  {5, 1.0, 1.8570, {2.4539, 2.2067, 2.0650, 1.9832, 1.9453}, 
	  	{69.3107, 51.8735, 36.3124, 21.5496,  7.1609 }}
};

/* Butterworth coefficients */
struct iirap butfp[] = {
  {0, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{00.000, 00.000, 00.000, 00.000, 00.000 }},
  {1, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{45.000, 00.000, 00.000, 00.000, 00.000 }},
  {1, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{60.000, 00.000, 00.000, 00.000, 00.000 }},
  {2, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{22.500, 67.500, 00.000, 00.000, 00.000 }},
  {2, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{36.000, 72.000, 00.000, 00.000, 00.000 }},
  {3, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{15.000, 45.000, 75.000, 00.000, 00.000 }},
  {3, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{25.714, 51.428, 77.142, 00.000, 00.000 }},
  {4, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{11.250, 33.750, 56.250, 78.750, 00.000 }},
  {4, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{20.000, 40.000, 60.000, 80.000, 00.000 }},
  {5, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{ 9.000, 27.000, 45.000, 63.000, 81.000000 }}
};

/* Chebyshev coefficient place holder */
struct iirap chebIfp[] = {
  {0, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{00.000, 00.000, 00.000, 00.000, 00.000 }},
  {1, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{45.000, 00.000, 00.000, 00.000, 00.000 }},
  {1, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{60.000, 00.000, 00.000, 00.000, 00.000 }},
  {2, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{22.500, 67.500, 00.000, 00.000, 00.000 }},
  {2, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{36.000, 72.000, 00.000, 00.000, 00.000 }},
  {3, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{15.000, 45.000, 75.000, 00.000, 00.000 }},
  {3, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{25.714, 51.428, 77.142, 00.000, 00.000 }},
  {4, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{11.250, 33.750, 56.250, 78.750, 00.000 }},
  {4, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{20.000, 40.000, 60.000, 80.000, 00.000 }},
  {5, 1.0, 1.0000, {1.0000, 1.0000, 1.0000, 1.0000, 1.0000}, 
	  	{ 9.000, 27.000, 45.000, 63.000, 81.000000 }}
};

/* rethink this since I have structure here for only the complex conjugate pairs 
 * For the Bessel I will need to modify the center frequency of the real pole
 * */


#define LOWPASS 3
#define HIGHPASS 2
#define BANDPASS 1
#define BANDREJECT 4

#define BUTTERWORTH 0
#define BESSEL	1
#define CHEBYSHEVI 2

float *x = (float *)NULL;

static void reverse(float *x, int n);
float fwarp(float omega, float dt);

void gsac_filt(float filt_fl, float filt_fh, int filt_np, int filt_p, int filt_lhp, int filt_type, float cheb_eps)
{
	int i, k, m, ntrc, npts;
	int ns;
	float dt;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float tmp;
	int iseven, nstage, nst;
	float wc,wu,wl,wumwl;
	float Wc, Wu, Wl, Zeta;
	float WC;
	float ffac;
	double  xm[3],ym[3];
	float zeta;
	float permin, permax; /*  for USER1 and USER2 header values for
				other programs */
	int nadd; /* additional points to add to temporary based
			 on period of low corner */
	int npass; /* 0 for causal, 1 for zero phase double pass */
	double num[3], den[3]; /* array that defines Laplace transform
		and z-transform recursive filter. The 5 is fixed since we
		only have to be worried about orders 1 and 2 (lowpass, highpass)
		and 3 and 4 (bandpass and bandreject ) stages. */
	int npole;   /* order of polynomial n=2 means s +1, n=3 = s^2 + as +b */
	int ndo; /* for second order stages is one for LP or HP and 2 for BP BR */



	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
	/* special case for Chebyshev Type I - define the coefficients */
	if(filt_type == CHEBYSHEVI){
		getchebpole(filt_np, cheb_eps);
	}


	if( (filt_np%2) == 0)
		iseven = YES;
	else
		iseven = NO ;
	nstage = filt_np/2;
	if(filt_lhp == BANDPASS || filt_lhp == BANDREJECT)
		ndo = 2;
	else
		ndo = 1;
	if(filt_lhp == LOWPASS)
		Wc = 6.2831853*filt_fh ;
	else if(filt_lhp == HIGHPASS)
		Wc = 6.2831853*filt_fl ;
	else if(filt_lhp == BANDPASS){
		Wc = 6.2831853*sqrt(filt_fl * filt_fh) ;
		Wl = 6.2831853*filt_fl;
		Wu = 6.2831853*filt_fh;
		/* safety check */
		if(Wl > Wu){
			tmp = Wl; Wl = Wu ; Wu = tmp;
		}
	} else if(filt_lhp == BANDREJECT){
		Wc = 6.2831853*sqrt(filt_fl * filt_fh) ;
		Wl = 6.2831853*filt_fl;
		Wu = 6.2831853*filt_fh;
		/* safety check */
		if(Wl > Wu){
			tmp = Wl; Wl = Wu ; Wu = tmp;
		}
	}

	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		dt = sacdata[k].sachdr.rhdr[H_DELTA];
		/* define the filter here */
		/* now process */
		if(npts > 0){
			if(filt_lhp == LOWPASS){
				nadd = 2.0/(filt_fh * dt);
			} else if (filt_lhp == HIGHPASS){
				nadd = 2.0/(filt_fl * dt);
			} else if (filt_lhp == BANDPASS){
				nadd = 2.0/(filt_fl * dt);
			} else if (filt_lhp == BANDREJECT){
				nadd = 2.0/(filt_fl * dt);
			}
			x = (float *)realloc(x,(npts+nadd)*sizeof(float));	
			/* fill in the temporary array */
			for(i=0 ; i < npts ; i ++)
				x[i] = sacdata[k].sac_data[i];
			for(i=npts ; i < npts +nadd ; i++)
				x[i] = x[npts-1];
			/* begin two pass for zero phase */
			for(npass=0 ; npass < filt_p ; npass++){

			/* single pole stage  is easy */
			if(iseven == NO){
				/* eventually use a switchable pointer here */
				switch(filt_type){
					case BUTTERWORTH:
						ffac = butfp[filt_np -1].fs;
						break;
					case BESSEL:
						ffac = besfp[filt_np -1].fs;
						break;
					case CHEBYSHEVI:
						ffac = chebIfp[filt_np -1].fs;
						break;
				}
				if(filt_lhp == LOWPASS){
					wc =fwarp(Wc*ffac,dt);
					num[0] = wc*dt;
					num[1] = 0.0;
					den[0] = wc*dt ;
					den[1] = 1.0;
					npole = 1;
					stozc(npole,1.0,num);
					stozc(npole,1.0,den);
					num[2] = 0.0;
					den[2] = 0.0;
				} else if (filt_lhp == HIGHPASS){
					wc =fwarp(Wc/ffac,dt);
					num[0] = 0.0;
					num[1] = 1.0;
					den[0] = wc*dt ;
					den[1] = 1.0;
					npole = 1;
					stozc(npole,1.0,num);
					stozc(npole,1.0,den);
					num[2] = 0.0;
					den[2] = 0.0;
				} else if (filt_lhp == BANDPASS){
					wc =fwarp(Wc,dt);
					wl =fwarp(Wl*ffac,dt);
					wu =fwarp(Wu*ffac,dt);
					num[0] = 0.0;
					num[1] = (wu-wl)*dt;
					num[2] = 0.0;
					den[0] = wc*wc*dt*dt ;
					den[1] = (wu-wl)*dt;
					den[2] = 1.0;
					npole = 2;
					stozc(npole,1.0,num);
					stozc(npole,1.0,den);
				} else if (filt_lhp == BANDREJECT){
					wc = fwarp(Wc,dt);
					wl = fwarp(Wl/ffac,dt);
					wu = fwarp(Wu/ffac,dt);
					num[0] = wc*wc*dt*dt;
					num[1] = 0.0;
					num[2] = 1.0;
					den[0] = wc*wc*dt*dt ;
					den[1] = (wu-wl)*dt;
					den[2] = 1.0;
					npole = 2;
					stozc(npole,1.0,num);
					stozc(npole,1.0,den);
				}
				/* normalize */
				for(i=0 ; i <= npole ; i++){
					num[i] /= den[npole];
					if(i != npole)
						den[i] /= den[npole];
				}
				for(i=0 ; i < 3 ; i++){
					xm[i] = 0;
					ym[i] = 0.0;
				}
					xm[0] = x[0];
				/* recursion */
				/* here xm[j] = xmj, ym[j= - ymj */
				for(i=0; i < npts+nadd ; i++){
					xm[0] = x[i];
					if(npole == 1)
					ym[0] = num[1]*xm[0] + num[0]*xm[1] - den[0]*ym[1];
					else if(npole == 2)
					ym[0] = num[2]*xm[0]+num[1]*xm[1]+num[0]*xm[2] 
						- den[1]*ym[1] - den[0]*ym[2];
					x[i] = ym[0];
					ym[2] = ym[1];
					ym[1] = ym[0];
					xm[2] = xm[1];
					xm[1] = xm[0];
				}
			}
			/*  2-poles stages are mmore work since
			 *  the lp -> BP or LP -> BR transformation of a
			 *  2 poles leads to a 4th order whose bi-linear Z 
			 *  representation is 4 order and noisy. In this case 
			 *  the 4th order is replaced by 2 second order stages
			 *  */

			if(filt_np > 1){
				for(m = 0 ; m < ndo ; m++){
				switch(filt_type){
					case BUTTERWORTH:
					nst = butfp[filt_np -1].n;
						break;
					case BESSEL:
					nst = besfp[filt_np -1].n;
						break;
					case CHEBYSHEVI:
					nst = chebIfp[filt_np -1].n;
				}
				for(ns = 0 ; ns < nst ; ns++){
				switch(filt_type){
					case BUTTERWORTH:
					zeta = cos(butfp[filt_np -1].p[ns]*3.1415927/180.0);
					ffac = butfp[filt_np -1].f[ns];
						break;
					case BESSEL:
					zeta = cos(besfp[filt_np -1].p[ns]*3.1415927/180.0);
					ffac = besfp[filt_np -1].f[ns];
						break;
					case CHEBYSHEVI:
					zeta = cos(chebIfp[filt_np -1].p[ns]*3.1415927/180.0);
					ffac = chebIfp[filt_np -1].f[ns];
						break;
				}
					if(filt_lhp == LOWPASS){
						wc =fwarp(Wc*ffac,dt);
						num[0] = wc*wc*dt*dt;
						num[1] = 0.0;
						num[2] = 0.0;
						den[0] = wc*wc*dt*dt ;
						den[1] = 2.0*zeta*wc*dt;
						den[2] = 1.0;
						npole = 2;
						stozc(npole,1.0,num);
						stozc(npole,1.0,den);
					} else if (filt_lhp == HIGHPASS){
						wc =fwarp(Wc/ffac,dt);
						num[0] = 0.0;
						num[1] = 0.0;
						num[2] = 1.0;
						den[0] = wc*wc*dt*dt ;
						den[1] = 2.0*zeta*wc*dt;
						den[2] = 1.0;
						npole = 2;
						stozc(npole,1.0,num);
						stozc(npole,1.0,den);
					} else if (filt_lhp == BANDPASS){
						wc =fwarp(Wc,dt);
						wl =fwarp(Wl*ffac,dt);
						wu =fwarp(Wu*ffac,dt);
						getroot(wl,wu,wc,zeta,&WC, &Zeta, m);
						wumwl= wu - wl ;
						num[0] = 0.0;
						num[1] = wumwl*dt;
						num[2] = 0.0;
						den[0] = WC*WC*dt*dt ;
						den[1] = 2.0*Zeta*WC*dt;
						den[2] = 1.0;
						npole = 2;
						stozc(npole,1.0,num);
						stozc(npole,1.0,den);
					} else if (filt_lhp == BANDREJECT){
						wc =fwarp(Wc,dt);
						wl =fwarp(Wl/ffac,dt);
						wu =fwarp(Wu/ffac,dt);
						getroot(wl,wu,wc,zeta,&WC, &Zeta, m);
						wumwl= wu - wl ;
						num[0] = wc*wc*dt*dt;
						num[1] = 0.0;
						num[2] = 1.0;
						den[0] = WC*WC*dt*dt ;
						den[1] = 2.0*Zeta*WC*dt;
						den[2] = 1.0;
						npole = 2;
						stozc(npole,1.0,num);
						stozc(npole,1.0,den);
					}
					/* normalize */
					for(i=0 ; i <= npole ; i++){
						num[i] /= den[npole];
						if(i != npole)
							den[i] /= den[npole];
					}
					den[npole] = 1.0;
					for(i=0 ; i < 3 ; i++){
						xm[i] = 0;
						ym[i] = 0.0;
					}
						xm[0] = x[0];
				for(i=0; i < npts +nadd ; i++){
					xm[0] = x[i];
					ym[0] =   num[2]*xm[0]
						+ num[1]*xm[1]
						+ num[0]*xm[2] 
						- den[1]*ym[1] 
						- den[0]*ym[2];
					x[i]  = ym[0];
					ym[2] = ym[1];
					ym[1] = ym[0];
					xm[2] = xm[1];
					xm[1] = xm[0];
				}
					
				}
				}
			}
			if(filt_p == 2)
				reverse(x,npts+nadd);
			}
			/* special case for Chebyshev Type I - apply the gain) */
			if(filt_type == CHEBYSHEVI){
				for(i=0 ; i < npts ; i++)
					x[i] *= chebIfp[filt_np -1].gain;
			}
			/* restore array */
			for(i=0 ; i < npts ; i++)
				sacdata[k].sac_data[i] = x[i];
			getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
			sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;

			sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
			sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
			sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
			/* define permin and permax of filter */
			if(filt_lhp == LOWPASS){
				permax = 100.*(npts * dt);
				permin = 1.0/filt_fh;
			} else if (filt_lhp == HIGHPASS){
				permax = 1.0/filt_fl;
				permin = 2.0*dt;
			} else if (filt_lhp == BANDPASS){
				permax = 1.0/filt_fl;
				permin = 1.0/filt_fh;
			} else if (filt_lhp == BANDREJECT){
				permax = 100.*(npts * dt);
				permin = 2.0*dt;
			}
			/* if filter bound is set at default, finally use */
			if(sacdata[k].permin == -12345.)
				sacdata[k].permin = permin;
			if(sacdata[k].permax == -12345.)
				sacdata[k].permax = permax;
			/* never expand range of periods */
			if(sacdata[k].permin < permin)
				sacdata[k].permin = permin;
			if(sacdata[k].permax > permax)
				sacdata[k].permax = permax;
		}
	}

}


static void reverse(float *x, int n)
{
	/* reverse a time series */
	int i;
	float tmp;
	if(n <= 1)
		return;
	for(i=0; i<= (n-1)/2 ; i++){
		tmp = x[i];
		x[i] = x[n-1-i];
		x[n-1-i]  = tmp;
	}
}


/*
c---------------------------------------------------------------------------c
c                                                                           c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                      c
c      VOLUME I                                                             c
c                                                                           c
c      PROCEDURE: STOZC                                                     c
c                                                                           c
c      COPYRIGHT 1986                                                       c
c      R. B. Herrmann                                                       c
c      Department of Earth and Atmospheric Sciences                         c
c      Saint Louis University                                               c
c      221 North Grand Boulevard                                            c  
c      St. Louis, Missouri 63103                                            c
c      U. S. A.                                                             c
c                                                                           c
c---------------------------------------------------------------------------c
*/
/* This converts an s polynomial to z coefficients */
/* H(s) = sum a(i)s**(i) <==> H(z) = sum a(i)z**(i) */
/* enter first the order of the polynomial, numerator = denominator */
/* enter in sampling rate in samples per second */
/* enter in the n+1 coefficients, starting with a(0) */
/* the output are the z-transform coefficients */

void stozc(int n,double fs,double *a)
	/* int n;	order of polynomial */
	/* double  fs;	sampling rate in Hz */
	/* double  a[];	array of coefficients - note original
            array is destroyed */
{
    double fac,temp;
    int i, k;

    /* introduce sampling frequency */
    fac =1.0;
    for(i=0;i<=n;i++){
        a[i] *= fac;
        fac *= 2.*fs;
    }

    /* synthetic division */
    for(i=0;i<n;i++){
        for(k=n-1;k>=i;k--){
            a[k]=a[k]+a[k+1];
        }
    }

    /* replace by reciprocals */
    k = (n-1)/2 ;
    for(i=0;i<=k;i++){
        temp=a[i];
        a[i]=a[n-i];
        a[n-i]=temp;
    }

    /* scale by -1/2 */
    fac =1.0;
    for(i=n;i>=0;i--){
        a[i] *= fac ;
        fac *= -2. ;
    }

    /* zero shift again */
    for(i=0;i<n;i++){
        for(k=n-1;k>=i;k--){
            a[k]=a[k]+a[k+1];
        }
    }
}

void getroot(float wl,float wh, float wc,float zeta,float *Wc, float *Zeta, int ndo)
{
/*
	we are solving a polynomial equation. The root should
	be between 0 and sqrt(wl*wh). So we start at wl.
	In addition we do not need that many iterations
*/
	float A2, B;
	float C,D,A4;
	float x, dx, fx, dfx, x2;
	float c8, c6, c4, c2, c0;
	float wn1, wn2, zeta1, zeta2;
	int iter;
	A2 = wc*wc;
	B = wh - wl;

	A4 = A2*A2;
	C = 2.*A2 + B*B ;
	D = 4.*A2*zeta*zeta*B*B;

	c8 = 1.0;
	c6 = 2.*A2 -C;
	c4 = D +2.*A2*(A2-C);
	c2 = A4*(2.*A2 - C);
	c0 = A4*A4;
	iter = 0;
	x = wl;
	dx = 1000;
	for( iter = 0 ; iter < 10 && ABS(dx) > 1.0e-5*wl ; iter++){
	dx = 1000;
	x2 = x*x;
	fx = c0 + x2*(c2 + x2*(c4 + x2*(c6 + x2*c8)));
	dfx = x*(2.*c2 + x2*(4.*c4 + x2*(6.*c6 + x2*(8.*c8))));
	dx = - fx/dfx;
	x = x + dx;
	}
 	wn1 = x;
	wn2 = A2 / wn1;
	zeta1 = ( zeta * B) /(wn1 + wn2);
	zeta2 = zeta1;
	if(ndo == 0 ) {
		*Wc = wn1;
		*Zeta = zeta1;
	} else {
		*Wc = wn2;
		*Zeta = zeta2;
	}

}

/* bi-linear frequency warping 
 * wd dt/2 = tan ( wa dt / 2)
 * or wd = 2/dt tan (w dt / 2 )  We have a problem when we get too
 * close to the Nyquist frequency because of the tan() function. So
 * we */
float fwarp(float wc, float dt)
{

	float fac;
	fac = 0.5*wc*dt;
	if(fac > 1.47113)
		return(2.0*10.0/dt);
	else
		return(2.0*sin(fac)/(cos(fac)*dt));
}

void getchebpole(int n, float eps)
{
	float p_real[10], p_imag[10], wn[10];
	float gain;
	float epsi,thetam,ca,sa,sinv,sh,ch;
	int m,i,j,ifortran;
	epsi = 1.0/eps;
 
	/* we only need half of the poles of n is even and half +1 if n is odd */
        m = n/2;
        if(n%2 == 1){
		m++;
	}
        for(i=0 ; i < m ; i++){
		ifortran = i + 1;
		thetam = (2*ifortran-1)*3.1415927/(2.*n);
                ca = cos(thetam);
                sa = sin(thetam);
                /* we could replace the next two lines */
                sinv = (1./n)*log(epsi+sqrt(epsi*epsi+1.));
                sinv = exp(sinv);
                /* by  sinv = (epsi+sqrt(epsi*epsi+1.))**1/n */
                sh = 0.5*(sinv - 1./sinv);
                ch = 0.5*(sinv + 1./sinv);
                p_real[i] = -sh*sa;
		p_imag[i] =  ch*ca;
		wn[i] = sqrt(p_real[i]*p_real[i] + p_imag[i]*p_imag[i]);
	}
	/* now set the chebIfp structure */
        if(n%2 == 0)
        	chebIfp[n-1].gain = 1./sqrt(1.0+eps*eps);
	else
        	chebIfp[n-1].gain = 1.0;
printf("gain: %f ripple: %f\n",chebIfp[n-1].gain,1./sqrt(1.0+eps*eps));
	j = 0;
	for(i=0;i< m;i++){
		ifortran = i + 1;
		if( (2*ifortran-1) == n){
			/* single real pole */
			chebIfp[n-1].fs = wn[i];
		} else {
			chebIfp[n-1].p[j] = acos(- p_real[i]/wn[i])*180.0/3.1415927;
			chebIfp[n-1].f[j] = wn[i];
			j++;
		}
	}
}
