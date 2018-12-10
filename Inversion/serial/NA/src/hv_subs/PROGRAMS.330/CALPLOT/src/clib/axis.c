/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: AXIS                                                  c
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

#include <math.h>
#define LN10 2.302585093
#include "calplot.h"

void axis(float xpage,float ypage,char *ititl,int nchar,float axlen,
	float angle,float firstv,float deltav)
{
/*
c-----
c       xpage,ypage  coordinates of starting point of axis in
c                    inches relative to current origin
c       ititl   axis title
c       nchar   number of characters in title
c               >0 for tic marks, numbers and title on
c                  counterclockwise side
c               <0 for tic marks, numbers and title on
c                  clockwise side
c       axlen   floating point axis length in inches
c       angle   angle (degrees) that x axis makes with
c               horizontal x-direction
c       firstv  scale value at (xpage,ypage)
c       deltav  change in slace between tic marks per unit
c               length
c               <0 value at xpage,ypage is a max, and values
c               decrease along the axis
c-----
c-----
c       we do not worry about the vagaries of storage of characters
c       since the address is passed on down to subroutine symbol
c-----
*/
	float	HTNUM=0.10,
		SYMHT=0.14,
		TICHT=0.07;
	float a, ex, xval, xdel, ct, st, dx, dy, xn, yn, z, yex, xt, yt;
	int kn, ntic, mtic, i;

        a=1.0;
        kn=nchar;
        if(kn < 0){
                a= -a;
                kn= -kn;
        }
/*
c-----
c       if deltav is too large invoke scientific notation
c-----
*/
        ex = 0.0;
        if(deltav != 0.0){
                yex= log(fabs(deltav))/LN10;
                if(yex < -2.0)
			ex = (int)(yex+0.01) - 1.0;
                else if(yex < -1.0)
			ex = (int)(yex+0.01)      ;
                if(yex >= 2.0) ex = (int)(yex + 0.01);
        }
        xval = firstv*pow(10.0,-ex);
        xdel = deltav*pow(10.0,-ex);
        ct = cos(angle*0.01745329);
        st = sin(angle*0.01745329);
        ntic = axlen + 1.0;
/*
c-----
c       first put in numbers and title
c       adjust offset for numbers
c-----
*/
        dx = - HTNUM;
        dy = 1.5*a*HTNUM - 0.5*HTNUM;
/*
c-----
c       find initial position given rotation
c-----
*/
        xn = xpage + dx*ct - dy*st;
        yn = ypage + dy*ct + dx*st;
        mtic = ntic/2;
	for(i=1;i<=ntic;i++){
                number(xn,yn,HTNUM,xval,angle,2);
                xval=xval+xdel;
                xn=xn+ct;
                yn=yn+st;
/*
c-----
c       halfway down axis put in title
c-----
*/
                if(i == mtic){
                        z=kn;
                        if(ex != 0.0)z=z+7.0;
                        dx = -0.5*SYMHT*z + axlen*0.5;
                        dy = (2.5*a-0.5)*SYMHT;
                        xt=xpage +dx*ct-dy*st;
                        yt=ypage +dy*ct+dx*st;
                        symbol(xt,yt,SYMHT,ititl,angle,kn);
                        if(ex != 0.0){
                                z=kn+2;
                                xt=xt+z*ct*SYMHT;
                                yt=yt+z*st*SYMHT;
                                symbol(xt,yt,SYMHT,"*10",angle,3);
                                xt=xt+(3.0*ct-0.8*st)*SYMHT;
                                yt=yt+(3.0*st+0.8*ct)*SYMHT;
                                number(xt,yt,0.7*SYMHT,ex,angle,-1);
                        }
                }
	}
/*
c-----
c       now put in tic marks
c-----
*/
        plot(xpage+axlen*ct,ypage+axlen*st,3);
        dx = - TICHT*st*a;
        dy =   TICHT*ct*a;
        a = ntic -1;
        xn = xpage + a*ct;
        yn = ypage + a*st;
	for(i=1;i<=ntic;i++){
                plot(xn,yn,2);
                plot(xn+dx,yn+dy,2);
                plot(xn,yn,2);
                xn=xn-ct;
                yn=yn-st;
	}
}
