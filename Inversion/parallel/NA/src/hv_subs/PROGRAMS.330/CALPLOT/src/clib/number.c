/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: NUMBER                                                c
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
#include <math.h>
#include "calplot.h"

extern struct Gcplot Gcplt;
static void fpack(char *num,int *n,float fpv,int mdec,int mwid,int nzflag);

void number (float xpage,float ypage,float height,
		float fpn,float angle,int ndec)
{
/*
c.....     xpage,ypage coordinates of lower left corner of number.
c.....     height   height of plotted number.
c.....     fpn      floating point number to be plotted.
c.....     angle    angle at which number is plotted, in degrees.
c.....     ndec     number of decimal places to be drawn.
c.....     this version of number requires the symbol version with
c.....     999. x, y feature, and  nc = 0 feature.
c.....
c.....     ndec .ge. 1000 invokes an exponential format
c.....     of the Fortran Format Ew.d, where d = ndec - 1000
c.....     and w = 7 + d. The output will look like
c.....     sd.dddddEsdd, where s is the sign and d is a digit
c.....     In all cases, no more than 9 significant figures are
c.....     permitted to the right of the decimal place.
c.....
c.....     ndec .ge. 2000 .le. 2009 invokes exponential but
c.....     plotted in scientific notation, e.g., 1.0 10 ^ exp
*/
	char num[20];
	char cexp[4];
/*
c-----
c    to avoid dregs in memory null first character
c-----
*/
	int i, n, ii, nzflag, mdec, iex, jex;
	int nman, nexp;
	float fp, aex;
	float ang, ca, sa, xl, yl, ht;
	float xx0, yy0, xxx0, yyy0;

	if(xpage < 999.0)Gcplt.x0 = xpage;
	if(ypage < 999.0)Gcplt.y0 = ypage;

	for(i=0;i<20;i++)
		num[i] = '\0';
	ii=0;
	n = 0;
        fp = fpn;
        if(fpn < 0.0){
                num[n]='-';
                n=n+1;
                fp= -fp;
        }
	if(ndec < 1000){
                nzflag = 0;
                mdec = -ndec;
                if(ndec == -1){
                        nzflag = 1;
                        mdec = 1;
                }
                fpack(num,&n,fp,mdec,19,nzflag);
        } else if(ndec >= 1000){
		mdec = 1;
		if( ndec >= 1000 && ndec <= 1009){
                	mdec = ndec -1000;
		} else if ( ndec >= 2000 && ndec <= 2009){
                	mdec = ndec -2000;
		} else {
			mdec = 9;
		}
                mdec = - mdec;
                if(fp == 0.0){
                        iex=0;
                } else {
                        aex = log10(fp);
                        iex = aex;
/*
c----------careful check of integer exponent
*/
                        if(iex < 0 && (float)iex > aex)iex=iex-1;
                }
                fp = fp * pow(10.0, (double)(-iex));
                nzflag = 0;
                fpack(num,&n,fp,mdec,14,nzflag);
/*
c---put in exponent assuming iex < 100
*/
		nman = n;
                num[n++]='E';
                if(iex < 0){
                        num[n++]='-';
                        iex = -iex;
                } else {
                        num[n++]='+';
                }
		nexp = n;
                jex = iex;
                jex = (jex/10) % 10;
                num[n++] = (char)('0' + jex);
                jex = iex % 10;
                num[n++] = (char)('0' + jex);
        }
	if(ndec <= 1009){
        	symbol(Gcplt.x0,Gcplt.y0,height,num,angle,n);
	} else if(ndec >= 2000 && ndec <=2009){
/*
c-----
c       save the exponent, replace the E-EX with 'x10'
c-----
*/
		cexp[0] = num[nexp-1];
		cexp[1] = num[nexp];
		cexp[2] = num[nexp+1];
		cexp[3] = '\0';
		num[nman] = 'x';
		num[nman+1] = '1';
		num[nman+2] = '0';
		num[nman+3] = '\0';
		xx0 = Gcplt.x0;
		yy0 = Gcplt.y0;
        	symbol(Gcplt.x0,Gcplt.y0,height,num,angle,nman+3);
/*
c----
c       get the proper position because of rotation
c-----
*/
		ang = angle*3.1415927/180.0;
		ca = cos(ang);
		sa = sin(ang);
		xl = height*(nman+3);
		yl = 0.6*height;
		xx0 = xx0 + xl*ca ;
		yy0 = yy0 + xl*sa ;
		xxx0 = xx0  - yl*sa;
		yyy0 = yy0  + yl*ca;
		ht = 0.7*height;
		symbol(xxx0,yyy0,ht,cexp,angle,3);
/*
c-----
c		now position at end of proper string
c-----
*/
		xl = 3.0* ht;
		Gcplt.xcur = xx0 + xl*ca ;
		Gcplt.ycur = yy0 + xl*sa ;
		Gcplt.x0 = Gcplt.xcur;
		Gcplt.y0 = Gcplt.ycur;
	}
}

static void fpack(char *num,int *n,float fpv,int mdec,int mwid,int nzflag)
{
	float fp;
	int nzflg, m, maxm, mm, ilw, iup, ipos, k, iex;
        fp = fpv;
        nzflg = nzflag;
/*
c-----since we have a maximum field width of mwid and since
c-----the . takes a position, we can only have mwid -1 figures
*/
	m = mdec ;
        maxm = 9;
        if(m > maxm)m=maxm;
        if(m <  -maxm)m= -maxm;
        mm = m;
        if(m > 0) mm--; 
/*
c----do a simple rounding up
*/
        fp = fp + 0.5 * pow(10.0, (double)mm);
/*
c---- get position of first significant figure
c
c     5 4 3 2 1 0 -1 -2 -3
c
c     d d d d d .  d  d  d
c
c----
*/

	iex = log10(fp);
        if(fp >= 1.0){
                iex++;
        } else {
                iex--;
        }
/*
c----
c----   procede from most significant down
c----   but never less than 10**-9
c----
*/
        
        if(fp <= 1.0e-9){
                fp = 0.0;
                iex = -9;
        }
        ilw = mdec;
        if(fp >= 1.0){
                iup=iex;
        } else {
                iup = 1;
        }
        if(m <= 0){
                ilw= m+1;
        } else {
                ilw = m ;
        }
        if(iex > 9){
                ilw = 1;
                nzflg = 1;
        }
	for(ipos=iup; ipos >= ilw;ipos--){
                k = (int)(fp * pow(10.0,(double)( -ipos +1 )) );
                if(*n < mwid){
                        num[*n] = (char)('0'+k);
			*n = *n + 1;
                }
                fp = fp - (float)(k) * pow(10.0,(double)(ipos-1));
        	if(ipos==1 && nzflg==0 && *n < mwid){
                        num[*n] = '.';
			*n = *n + 1;
                }

	}
}
