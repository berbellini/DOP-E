#include        <stdio.h>
#include        "gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

/* implement INTERPOLATE 
	CHANGES: 
        26 MAY 2006 routine inter re-written 26 May 2006 INGV Roma
	because logic was bad. Also cubic interpolation works at
	the ends instead of linear.

        08 AUG 2006 - internal x and sx arrays in inter pccp are forced to be double
        there was a problem with delta = 0.005 and one days record of 86400 samples. 
        note that this requires 86400*200=17280000 unique values which gets too close
        to the granularity of single precision (32 bit) floating point numbers

	25 OCT 2010 - special case in gsac_inter now considered. Basically given
	the pair (sx, sy) and a desired x, we estimate y. However if x < sx[1] or 
	x > sx[N-2) then we use linear interpolation. These extremes were not considered
	before.

*/

extern struct sacfile_ *sacdata;

#define	IN_DELTA 0

float in_delta = 0;
static float *sy = (float *)NULL;
static double *sx = (double *)NULL;
static double *x = (double *)NULL;
void gsac_inter(double *sx,float *sy,double *x,float *y,int npts, int *mpts);
void  gsac_pccp(double *x, float *y, double xd, float *yd);



struct arghdr inarg[] = {
	{IN_DELTA, "DELTA", RHDR, 0, 1, YES, "", 1},
	{-10, ""        , CHDR, 0, 0, YES, "" ,-1}
};

/* these are temporary variables only used here */
float in_real[10];
int do_in;


void gsac_set_param_in(int ncmd, char **cmdstr)
{
	int i;
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* note when the testrg routine is used, if the argument is
		NO then you must use internal variables to define the 
		state of the operation - if you use YES, then things are
		not changed until the input is proven correct. An exmple of
		this concept with YES is the following:
		Assume we wish aa LP filter with fc 1 np 2 p 1 
		If we enter  fc 2 np2   there is a syntax error and we
		should not chnge the fc since the np2 is wrong. One way to
		do this in the code would be to do two calls

			if(testarg(ncmd, cmdstr, cmdargs, YES) is OK
			then
				testarc,ncmd, cmdstr, cmdargs, NO)
		*/
	if(testarg(ncmd, cmdstr, inarg, NO, YES))
		return;
	do_in = NO;
	for(i=0 ; inarg[i].key[0] != '\0' ; i++){
		if(inarg[i].used > 0){
			if(inarg[i].ricell == RHDR){
			getargr(ncmd, cmdstr, inarg[i].key, 
				inarg[i].mfit, inarg[i].narg, in_real );
				if(in_real[0] > 0){
					do_in = YES;
					in_delta = in_real[0];
				}
			}

		}
	}
}

void gsac_exec_in(void)
{
	int i, k, ntrc, npts, mpts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float t0,dt;

	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;

	gsac_control.begmin= 1.0e+30;
	gsac_control.endmax=-1.0e+30;

	/* process the traces */
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		/* copy the data into a temporary array */
		if(sx == (double *)NULL)
			sx = (double *)calloc(npts,sizeof(double));
		else
			sx = (double *)realloc(sx,npts*sizeof(double));
		if(sy == (float *)NULL)
			sy = (float *)calloc(npts,sizeof(float));
		else
			sy = (float *)realloc(sy,npts*sizeof(float));
		t0 = sacdata[k].sachdr.rhdr[H_B];
		dt = sacdata[k].sachdr.rhdr[H_DELTA];
		for(i=0;i< npts;i++){
			sx[i] = (double)i*(double)dt;
			sy[i] = sacdata[k].sac_data[i];
		}
		/* define the new number of points */
		mpts = (int)((sacdata[k].sachdr.rhdr[H_E]
			-sacdata[k].sachdr.rhdr[H_B])/in_delta +0.49 ) + 1;
		/* reallocate space for the interpolated array */
		if(x == (double *)NULL)
			x = (double *)calloc(mpts,sizeof(double));
		else
			x = (double *)realloc(x,mpts*sizeof(double));
		sacdata[k].sac_data = (float *)
			realloc(sacdata[k].sac_data,mpts*sizeof(float));
		/* now interpolate */
		for(i=0;i< mpts;i++){
			x[i] =  (double)i*(double)in_delta;
		}
		gsac_inter(sx,sy,x,sacdata[k].sac_data,npts, &mpts);

		/* update the header values */
		getmxmn(sacdata[k].sac_data, mpts,&depmax, &depmin, &depmen,&indmax,&indmin);
		sacdata[k].sachdr.ihdr[H_NPTS] = mpts;
		sacdata[k].sachdr.rhdr[H_DELTA] = in_delta;
		sacdata[k].sachdr.rhdr[H_E] = sacdata[k].sachdr.rhdr[H_B] +
				(mpts -1 )*in_delta;
		sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
		sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
		sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;

		sacdata[k].tzbeg=sacdata[k].tzref+sacdata[k].sachdr.rhdr[H_B];
		sacdata[k].tzend=sacdata[k].tzref+sacdata[k].sachdr.rhdr[H_E];
		/* get bounds for absolute plotting */
		sacdata[k].tzbegx = sacdata[k].tzbeg;
		sacdata[k].tzendx = sacdata[k].tzend;
		if(sacdata[k].tzbeg < gsac_control.begmin)
			gsac_control.begmin = sacdata[k].tzbeg;
		if(sacdata[k].tzend > gsac_control.endmax)
				gsac_control.endmax = sacdata[k].tzend;
	}

}

void gsac_inter(double *sx,float *sy,double *x,float *y,int npts, int *mpts)
{
/* Wiggins, R. A. (1976). Interpolation of digitized curves,
 * Bull. Seism. Soc. Am. 66, 2077-2081.
 *
 * From David Russell's match filter routine. This routine
 * preserves the ability to use unevenly spaced data and then
 * to interpolate that data */
/*
        sx      input array of abscissas        [0,npts-1]
        sy      input array of ordinates        [0,npts-1]

        x       output array of abscissas equally spaced        [0,mpts-1]
        y       output array of ordinates                       [0,mpts-1]
*/
/* 20 AUG 2004 modified FORTRAN code so that never exceed allocation space, and permit possibility the the looping actually produced fewer points
 * */
/* design for the order of interpolation - initially 2 pt for linear 
then 4 point for Wiggins spline */
#define NORD 4
	double xx[NORD];
	float yy[NORD];
	int i, j;
	int k;
	int Mpts;
	double p;
/* initialize - we know the first output point so set up the interpolation
*/
	y[0] = sy[0];
	Mpts = *mpts;
	if(npts < NORD){
		*mpts = 1;
		return;
	} else {
		for(k=0;k < NORD ; k++){
			xx[k] = sx[k];
			yy[k] = sy[k];
		}
	}

	/* special indices - for a cubic interpolater we use a 4
		point interpolation grid. This grid marches along but must never
		sample beyond the ends of the sx array - also when we get near 
		the end we do not use center points but actually implement
		a one sided cubic - this avoids linear interpolation at the
		ends. If the interpolator grid is 0 1 2 3, whenever the
		grid wants to use the 2 3 interval we shift the data to avoid
		this except at the end of the array

		i = index of output array 
		j = index of upper limit of input array that maps into the '2'
			position. hoever to avoid going off the end we
			test so that the j+1 point < npts 

	*/
	/* linear interpolation for sx[0] <=x < sx[1]
           cubic  interpolation for sx[1] <=x <= sx[N-2]
           linear interpolation for sx[N-2] < x <= sx[N-1]
	*/
        i=0;
        while(x[i] < sx[1] && i < Mpts){
		p = (x[i]-sx[0])/(sx[1]-sx[0]);
		y[i]=p*sy[1] + (1.0-p)*sy[0];
		i++;
	}
	for(j=2; i < Mpts && x[i] < sx[npts-1] ; i++){
		/* shift the interpolator if necessary */
		while(x[i] > sx[j] && j < npts -2){
			for(k=0;k < NORD ; k++){
				xx[k] = sx[k+j-1];
				yy[k] = sy[k+j-1];
			}
			j++;
		}
			
		/* do the interpolation */
		gsac_pccp(xx, yy,  x[i], &y[i]);
	}
	while(x[i] <= sx[npts -1] && i < Mpts){
		p = (x[i]-sx[npts-2])/(sx[npts-1]-sx[npts-2]);
		y[i]=p*sy[npts-1] + (1.0-p)*sy[npts-2];
		*mpts = i +1;
		i++;
	}
}
		
void  gsac_pccp(double *x, float *y, double xd, float *yd)
{
/*
  subroutine pccp performs piecewise continuous cubic
  polynomial interpolation following wiggins(1976)
  and akima(1970). the method requires two x and y
  points on each side of position xd where the value of yd is
  determined. weighted averages of the slopes are determined at the
  knots(positions x(2) and x(3)), and the slopes and the values x(2),
  y(2),x(3),y(3) are used to determine yd which corresponds to
  x position xd which falls between x(2) and x(3).
  references
  wiggins,r.a.,bull. seism. soc. am.,v.66,p2077-2081,1976.
  akima,h.,j.assoc.comp.mach.,v.17,p 589-602,1970.

  note the x-values must be distinct!!
*/
double eps, sx2, sx3, sx4, w1, w2,w3, s1, s2, p0, p1, p2, p3;
float xmxd;
	eps=0.001;
	/*
	determine slopes at x(2) and x(3)
	*/
	sx2=(y[1]-y[0])/(x[1]-x[0]);
	sx3=(y[2]-y[1])/(x[2]-x[1]);
	sx4=(y[3]-y[2])/(x[3]-x[2]);
	/*
	weight slopes.
	*/
	w1=1./MAX(ABS(sx2),eps);
	w2=1./MAX(ABS(sx3),eps);
	w3=1./MAX(ABS(sx4),eps);
	s2=(w2*sx3+w3*sx4)/(w2+w3);
	s1=(w1*sx2+w2*sx3)/(w1+w2);
	/*
	evaluate polynomial at x=xd
	with slopes given at x(1) and x(2)
	*/
	p0=y[1];
	p1=s1;
	p2=(3.*(y[2]-y[1])/(x[2]-x[1])-2.*s1-s2)/(x[2]-x[1]);
	p3=(s1+s2-2.*(y[2]-y[1])/(x[2]-x[1]))/((x[2]-x[1])*(x[2]-x[1]));
	xmxd=xd-x[1];
	*yd = p0 + xmxd*(p1 + xmxd*(p2 + xmxd*p3));
}
