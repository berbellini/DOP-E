/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GINFO                                                 c
c                                                                     c
c      COPYRIGHT (C)  1997 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      3507 Laclede Avenue                                            c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
*/

#include "callib.h"

#include "dosubs.h"
#include "calplot.h"
 
#include <stdio.h>

extern struct Gcplot Gcplt;

void  ginfo(int *HasMouse, float *XminDev, float *YminDev, 
	float *XmaxDev, float *YmaxDev, float *XminClip, 
	float *YminClip, float *XmaxClip, float *YmaxClip, int *Color)
{
	INT gHasMouse, gXminDev, gYminDev, gXmaxDev, gYmaxDev;
	INT gXminClip, gYminClip, gXmaxClip, gYmaxClip, gColor;
	(*do_info)(&gHasMouse, &gXminDev, &gYminDev, 
		&gXmaxDev, &gYmaxDev, &gXminClip, 
		&gYminClip, &gXmaxClip, &gYmaxClip,&gColor);
	*HasMouse = gHasMouse ;
	/* For Device boundaries, use absolute coordinates */
	*XminDev = gXminDev ;
	*YminDev = gYminDev ;
	*XmaxDev = gXmaxDev ;
	*YmaxDev = gYmaxDev ;
	/* return clip coordinates with respect to current origin */
	*XminClip = gXminClip /1000.0 ;
	*YminClip = gYminClip /1000.0 ;
	*XmaxClip = gXmaxClip /1000.0 ;
	*YmaxClip = gYmaxClip /1000.0 ;
	*XminClip = (*XminClip-Gcplt.xold)/Gcplt.xstp;
	*YminClip = (*YminClip-Gcplt.yold)/Gcplt.ystp;
	*XmaxClip = (*XmaxClip-Gcplt.xold)/Gcplt.xstp;
	*YmaxClip = (*YmaxClip-Gcplt.yold)/Gcplt.ystp;
	*Color = gColor;
}
