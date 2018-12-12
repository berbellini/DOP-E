#ifndef GSAC_PLOT_
#define GSAC_PLOT_

#include "calplot.h"
#include "grphsubc.h"

/* data structure to define contents of a plot window 
 * This will completely describe every aspect of the trace display */
struct plt_ctl {
	int	abstime;	/* YES or NO */
	float	xl; 	/* physical coordinates of the plot frame */
	float	yl;
	float	xh;
	float	yh;
	int	k;		/* index of the trace being plotted */
	int	n;		/* particular plot window top frame = 0 
			 	bottom = numperframe -1 */
	int	npts;
	float	xlen;	/* this is actually xh - xl */
	double	tmin;	/* absolute minimum time */
	double	tmax;	/* absolute maximum time */
	double	tfirst; /* time of first sample plotted */
	double	tlast;  /* time of last  sample plotted */
	int	ifirst; /* array index of first sample plotted */
	int	ilast;  /* array index of last  sample plotted */
	float	delta;	/* sample interval */
	float	ymult;	/* trace amplitude multiplier */
	float	ymin;	/* minimum value of a(t) in window */
	float	ymax;	/* maximum value of a(t) in window */
	float	dyda;	/* mapping of a(t) to y(t) for frame */
	float	ylen;	/* this is actually yh - yl */
	float	uymul;	/* trace multiplier for ppk */
	float depmax;	/* trace maximum amplitude */
	float depmin;	/* trace maximum amplitude */
};

/* define command codes for ylim */
#define YLIM_OFF 0	/* ylim is not controlled */
#define YLIM_ALL 1	/* force all traces to have same scale */
#define YLIM_USR 2	/* for user specified scale */
#endif
