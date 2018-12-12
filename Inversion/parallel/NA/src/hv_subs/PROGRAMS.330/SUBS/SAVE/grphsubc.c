#include	"calplot.h"
#include	"grphsubc.h"
#include	<math.h>
#include	<string.h>
#include	<stdio.h>
#include	<stdlib.h>
/*
 * c-----
 * c       changes
 * c
 * c       01 JAN 2004 adjusted position of Y-axis title in doliny dolny etc
 * c               for case when numbers are negative
 * c       08 JAN 2005 "while(c0 != 0 && c != 0) to "while(c0 != 0 && c1!= 0)
 #			in rclip baker@usgs.gov
 * c----
 *
*/


	/* DEFINES */
#define MIN(a,b) ( (b) > (a) ? (a):(b) )
#define MAX(a,b) ( (a) > (b) ? (a):(b) )
#define ABS(a)   ( (a) > (0) ? (a): -(a) )
#define ON	1
#define OFF	0
	/* PROTOTYPES not in grphsubc.h */
static void ffpack(char *num,int *n,float fpv,int mdec,int mwid,int nzflag);
static void linecode(float x,float y,float ymax,float xmin,float ymin,float xmax,int *c) ;
static void mynum (float xpage,float ypage,float height,float fpn,
        float angle,int ndec,char *cmd);


	/* GLOBAL VARIABLES */
static float num10[] =
	{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };


void gbox(float xl,float yl,float xh,float yh)
{
	plot(xl,yl,3);
	plot(xh,yl,2);
	plot(xh,yh,2);
	plot(xl,yh,2);
	plot(xl,yl,2);
}

void gcent(float xx,float yy,float ht,char *string,float angle)
{
	int il;
	float t, ct, st, rl, xxx, yyy; 
	il = strlen(string);
	t = angle*3.1415927/180.0;
	ct = cos(t);
	st = sin(t);
	rl = 0.5*il*ht;
		
	xxx = xx - rl*ct;
	yyy = yy - rl*st;
	symbol(xxx,yyy,ht,string,angle,il);
}

void gright(float xx,float yy,float ht,char *string,float angle)
{
	int il;
	float t, ct, st, rl, xxx, yyy; 
	il = strlen(string);
	t = angle*3.1415927/180.0;
	ct = cos(t);
	st = sin(t);
	rl = 1.0*il*ht;
		
	xxx = xx - rl*ct;
	yyy = yy - rl*st;
	symbol(xxx,yyy,ht,string,angle,il);
}

void gleft(float xx,float yy,float ht,char *string,float angle)
{
	int il;
	float t, ct, st, rl, xxx, yyy; 
	il = strlen(string);
	t = angle*3.1415927/180.0;
	ct = cos(t);
	st = sin(t);
	rl = 0.0*il*ht;
		
	xxx = xx - rl*ct;
	yyy = yy - rl*st;
	symbol(xxx,yyy,ht,string,angle,il);
}

void dology(float x0,float y0,float yleng,float ymax,float ymin,
        float sizey,int ticlft,int lablft,int dopow,
	int ly, char *titley)
{
/*
c	put in tic marks along the y-axis 
c
c	x0	R*4	- position of bottom side of axis
c	y0	R*4	- position of bottom side of axis
c	yleng	R*4	- length of Y axis
c	ymax	R*4	- maximum value of number corresponding to
c				top
c	ymin	R*4	- minimum value of number corresponding to
c				bottom
c	nocy	I*4	- number of cycles along the axis
c	sizey	R*4	- size of numbers, e.g., 10
c	ticlft	L	- .true. tics go left
c			  .false. tics go right
c	lablft	L	- .true. 10**N goes to left of axis
c			  .false. goes to right
c	dopow	L	- .true. put in 10**N labels
c	ly	I*4	length of y-axis title string
c	titley	Ch	y-axis title string
*/
	int nocy; 
	int iy,iiy;
	float ymxlog, ymmin,ymlog,xshft;
	int ja, jpow,ii,jj;
	float xpos, xpos1, xpos2, ypos, ypos1;
	int ifirst;
	float ylow, yscal;
	float xx, yy;
	float tenpow, yval, yv, tlen;
	char *title;
	float ycen;


	/* 	put in the axis	*/
	ycen = y0 + 0.5*yleng ;
	plot(x0,y0,3);
	plot(x0,y0+yleng,2);
	/* set up positioning for 10**N labels */
	/* compute limits for scale */
	ymxlog = log10(ymax);
	ymmin = log10(ymin);
	nocy = ymxlog - ymmin + 1;
	iy = ymmin;
	if(ymmin < 0)iy = iy - 1;

	ymlog = log10(ymax);
	iiy = ymmin;
	if(ymlog < 0)iiy = iiy - 1;

	ja = MAX(ABS(iy),ABS(iiy));
	if(ja<10)
		jpow = 1;
	else if(ja>=10 && ja<100)
		jpow = 2;
	else if(ja>=100 && ja<1000)
		jpow = 3;
	ii = MIN(iy,iiy);
	if(ii<0)jpow = jpow + 1;
	jpow = jpow + 2;
	xshft = (2 + jpow*0.7)*sizey;

	if(lablft){
		if(ticlft){
			xpos = x0 - 0.50*sizey - xshft;
			xpos1 = xpos + 0.7*sizey;
			xpos2 = x0 - 0.50*sizey - xshft -1.2*sizey ;
		} else {
			xpos = x0 + 0.50*sizey - xshft;
			xpos1 = xpos + 0.7*sizey;
			xpos2 = x0 + 0.50*sizey - xshft -1.2*sizey ;
		}
	} else {
		if(ticlft){
			xpos = x0 + 1.00*sizey;
			xpos1 = xpos + 0.7*sizey;
			xpos2 = x0 + 0.00*sizey + 0.2*sizey + xshft;
		} else {
			xpos = x0 + 2.0*sizey;
			xpos1 = xpos + 0.7*sizey;
			xpos2 = x0 + 1.0*sizey + 0.2*sizey + xshft;
		}
	}

	ifirst = 1;
	ylow = y0;
	yscal = yleng/log10(ymax/ymin);
	for(ii=iy ; ii <= iy+nocy+2 ; ii++){
		tenpow = pow(10.0,ii);
		for(jj=1 ; jj <= 9 ; jj++){	
			yval = num10[jj-1]*tenpow;
			/*put in tics where appropriate */
			if(yval >= ymin && yval<=ymax){
				yv = log10(yval/ymin);
				yy = y0 + yv*yscal;
				if(jj==1){
					tlen = 1.5*sizey;
			/*		put in 10**N */
					if(ifirst==1) {
						ypos = yy;
					} else {
						ypos = yy - 0.5*sizey;
						if(yval==ymax)ypos=yy-1.4*sizey;
						if(ypos<y0)ypos = y0;
					}
					if(dopow){
					ypos1 = ypos + 0.7*sizey;
					symbol(xpos,ypos,sizey,"10",0.0,2);
					number(999.,ypos1,0.7*sizey,(float)(ii),0.0,-1);
					}
				} else {
					tlen = sizey;
				}
				ifirst = ifirst + 1;
				plot(x0,yy,3);
				if(ticlft)	
					plot(x0-tlen,yy,2);
				else
					plot(x0+tlen,yy,2);
			}
		}
		plot(x0,y0,3);
	}
/*
c-----
c       put in the title if we put in the numbers
c-----
*/
        if(dopow &&  ly > 0){
		title = (char *)calloc(ly+1,sizeof(char));
		strncpy(title,titley,ly);
		title[ly] = '\0';
                if(lablft){
                        gcent(xpos2,ycen,1.2*sizey,title, 90.0);
                } else {
                        gcent(xpos2,ycen,1.2*sizey,title,-90.0);
                }
		free(title);
        }                 
}

void dologx(float x0,float y0,float xleng,float sxmax,float sxmin,
        float sizex,int ticup,int labtop,int dopow,
	int lx, char *titlex)
{
			/*
c	put in tic marks along the x-axis 
c
c	x0	R*4	- position of left side of axis
c	y0	R*4	- position of left side of axis
c	xleng	R*4	- length of X axis
c	xmax	R*4	- maximum value of number corresponding to
c				far right
c	xmin	R*4	- maximum value of number corresponding to
c				far left
c	sizex	R*4	- size of numbers, e.g., 10
c	ticup	L	- .true. tics go up
c			  .false. tics go down
c	labtop	L	- .true. 10**N goes above axis
c			  .false. goes below
c	dopow	L	- .true. put in 10**N labels
			*/
	int nocx; 
	int ix,iix;
	float xmxlog, xmmin,xmlog,yshft;
	int ja, jpow,ii,jj;
	float xpos, xpos1, ypos, ypos1, ypos2;
	int ifirst;
	float xlow, xscal;
	float xx, yy;
	float tenpow, xval, xv, tlen;
	float xmin, xmax, xshft;
	float xcen;
	char *title;




	xmin = sxmin;
	xmax = sxmax;
			/* 	put in the axis */
	xcen = x0 + 0.5*xleng ;
	plot(x0,y0,3);
	plot(x0+xleng,y0,2);
			/* 	set up positioning for 10**N labels */
	if(labtop){
		if(ticup){
			ypos  = y0   + 2.0 * sizex	;
			ypos1 = ypos + 0.7*sizex	;
			ypos2 = y0   + 3.7 * sizex 	;
		} else {
			ypos  = y0   + sizex		;
			ypos1 = ypos + 0.7*sizex	;
			ypos2 = y0   + 2.7 * sizex 	;
		}
	} else {
		if(ticup){
			ypos  = y0   - 2.0*sizex	;
			ypos1 = ypos + 0.7*sizex	;
			ypos2 = y0   - 3.7*sizex 	;
		} else {
			ypos  = y0   - 3.0*sizex	;
			ypos1 = ypos + 0.7*sizex	;
			ypos2 = y0   - 4.7*sizex 	;
		}
	}
			/* 	compute limits for scale */
	xmxlog = log10(xmax)	;
	xmmin = log10(xmin)	;
	nocx = xmxlog - xmmin + 1;
	ix = xmmin;
	if(xmmin < 0)ix = ix - 1;
	ifirst = 1;
	for(ii=ix ; ii <= ix+nocx+2; ii++){
		tenpow = pow(10.0,ii);
		ja = ABS(ii);
		if(ja<10)
			jpow = 1;
		else if(ja>=10 && ja<100)
			jpow = 2;
		else if(ja>=100 && ja<1000)
			jpow = 3;
		if(ii<0)jpow = jpow + 1;
		xscal = xleng/log10(xmax/xmin);
		for(jj=1 ; jj <= 9 ; jj++){
			xval = num10[jj-1]*tenpow;
			/* 	put in tics where appropriate */
			if(xval >= xmin && xval<=xmax){
				xv = log10(xval/xmin);
				xx = x0 + xv*xscal;
				if(jj==1){
					tlen = 1.5*sizex;
					/* put in 10**N */
					if(ifirst==1){
						xpos = xx;
					} else {
						xshft = (2 + jpow*0.7)*sizex;
						if(xx+xshft > x0+xleng){
							xpos = x0+xleng - xshft;
						} else {
							xpos = xx - xshft/2.0;
						}
						if(xpos<x0)xpos = x0;
					}
					if(dopow){
					symbol(xpos,ypos,sizex,"10",0.0,2);
					number(999.,ypos1,0.7*sizex,(float)(ii),0.0,-1);
					}
				} else {
					tlen = sizex;
					ifirst = ifirst + 1;
				}
				plot(xx,y0,3);
				if(ticup)	
					plot(xx,y0+tlen,2);
				else
					plot(xx,y0-tlen,2);
	
			}
		}
		plot(xx,y0,3);
	}
        if(dopow && lx > 0){
		title = (char *)calloc(lx+1,sizeof(char));
		strncpy(title,titlex,lx);
		title[lx] = '\0';
                gcent(xcen,ypos2,1.2*sizex,title,0.0);
		free(title);
        }  
}

void dolinx(float x0,float y0,float xleng,float xmax,float xmin,
        float sizex,int ticup,int labtop,int dopow,
	int llx, char *titlex)
{
	dolnx(x0,y0,xleng,xmax,xmin,sizex,ticup,labtop,dopow,llx,titlex, ON);
}

void dnlinx(float x0,float y0,float xleng,float xmax,float xmin,
        float sizex,int ticup,int labtop,int dopow,
	int llx, char *titlex)
{
	dolnx(x0,y0,xleng,xmax,xmin,sizex,ticup,labtop,dopow,llx,titlex, OFF);
}

void dolnx(float x0,float y0,float xleng,float xmax,float xmin,
        float sizex,int ticup,int labtop,int dopow,
	int llx, char *titlex, int doasis)
{
			/*
c	put in tic marks along the x-axis 
c
c	x0	R*4	- position of left side of axis
c	y0	R*4	- position of left side of axis
c	xleng	R*4	- length of X axis
c	xmax	R*4	- maximum value of number corresponding to
c				far right
c	xmin	R*4	- maximum value of number corresponding to
c				far left
c	sizex	R*4	- size of numbers, e.g., 10
c	ticup	L	- .true. tics go up
c			  .false. tics go down
c	labtop	L	- .true. number goes above axis
c			  .false. goes below
c	dopow	L	- .true. put in number labels
c	doasis	L	- .true. x-axis numbers increase from left to right
			*/

	int dosci;
	int ndec,nshft;
	float xl, tlen1, tlen2, ypos, ypos1, ypos2, xv, xval, xnorm;
	int nxnorm;
	float xmx, xmn, dxx, dxxx, xc, xe;
	float xb, xx, ticlen, xxx, xxv;
	int n, ifac, ilow, iup, i, j, ndxxx;
	char *title;
	float sizexx;
	float xcen;
	int nscalx;
	int lx;

	char cmd[11];

	lx = llx;
	if(lx < 0)lx = 0;
			/* put in the axis */
	xl = x0 + xleng;
	xcen = x0 + 0.5*xleng ;
	plot(x0,y0,3);
	plot(xl,y0,2);
			/* set up positioning for tics */
	tlen1 = 0.6*sizex;
	tlen2 = 1.2*sizex;
	if(labtop){
		if(ticup){
			ypos  = y0   + 2.0 * sizex;
			ypos1 = ypos + 0.7*sizex;
			ypos2 = y0   + 4.0 * sizex ;
		} else {
			ypos  = y0   + sizex;
			ypos1 = ypos + 0.7*sizex;
			ypos2 = y0   + 3.0 * sizex ;
		}
	} else {
		if(ticup){
			ypos  = y0   - 2.0*sizex;
			ypos1 = ypos + 0.7*sizex;
			ypos2 = y0   - 4.5*sizex ;
		} else {
			ypos  = y0   - 3.0*sizex;
			ypos1 = ypos + 0.7*sizex;
			ypos2 = y0   - 5.5*sizex ;
		}
	}
	/* compute limits for scale */
	xmn = MIN(xmin,xmax);
	xmx = MAX(xmin,xmax);
	xxv = xmn;
	/* safety */
	if(xmn==xmx)xmx = xmn + 1.0;
	/* get some idea of the size of 
		the number to be plotted */
	xval = MAX(ABS(xmn),ABS(xmx));
	/*
c	we set the maximum number of decimal digits as ndig, 
c	if the number is greater than 100 or less than 0.01, use
c	scientific notation
c	get the maximum number of digits for the number
	*/
	if(xval<=0.01 &&xval>0.0 ){
		dosci = ON;
		xnorm = log10(1.0/xval);
		nxnorm = xnorm+1;
		xnorm = pow(10.0,  nxnorm);
		nscalx = -nxnorm;
	}else if( xval>=10000.0){
		dosci = ON;
		xnorm = log10(xval);
		nxnorm = xnorm;
		nscalx =  nxnorm;
		xnorm = pow(10.0, -nxnorm);
	} else {
		dosci = OFF;
		xnorm = 1.0;
		nxnorm = 0.0;
		nscalx = 0;
	}
			/*
c	choose the maximum number of tics somehow on
c	xleng, and size of numbers and necessary decimals
			*/
	dxx = xnorm*(xmx - xmn)/5.0;
	if(dxx == 0.0)dxx = 1.0;
			/*
c	get start limits for search - we will place at most 10 tics
c	here
			*/
	n = log10(dxx);
	if(dxx<1.0)n = n -1;
	ifac = pow(10.0, -n)*dxx;

			/*
c	now get a new minimum value that is some multiple 
c	of ifac 10**n 
			*/
	dxx = ifac * pow(10.0, n);
	ilow = xnorm*xmn/dxx -1;
	iup  = xnorm*xmx/dxx +1;
/* define a resonable range for the numbers */
	if(ifac == 1)
		ndxxx = 5;
	else if(ifac == 2)
		ndxxx = 5;
	else
		ndxxx = ifac;
	dxxx = dxx/ndxxx;

			/*
c	loop for labels
c	xe is the end of the previous label
			*/
	xc = xmn;
	plot(x0,y0,3);
	/* this is the heart of the axis labeling. Everything is done in
	 * terms of the actual physical units of the axis and not in terms 
	 * of the scaled units for the plot
	 *
	 * A particular tic mark corresponds to the value
	 *
	 * xxv = i*dxx + j*dxxx , where
	 * i = [ilow,iup]
	 * j = [0, ndxxx ]
	 * The major tic mark and number label corresponds in the i index,
	 * The Semi tic mark corresponds to the j index
	 * The minor tic mark corresponds to the k index
	 * For example, it we have something that looks like
	 *
	 * |                   |
	 * |   |   |   |   |   |      j=[0,5]
	 *
	 *
	 * x0 is physical position of left axis end in plot space
	 * xl is physical position of right axis termination in plot space
	 * xx is the physical position in plot space
	 * xb is position of number to be plotted
	 * xe is position of end last number to be plotted
	 */
	xe = x0;
	for(i=ilow ; i<=iup ; i++){
		for(j=0 ; j <= ndxxx-1 ; j++){
		
			xxv = i*dxx + j*dxxx ;
			xxv = xxv/xnorm;
			if(xxv>=xmn && xxv <=xmx){
				xx = x0 + (xxv - xmn)*xleng/(xmx - xmn);
				plot(xx,y0,3);
				if(j==0 ){
			/*
c	estimate the number of digits for proper labeling
c	only do this once for the plot
c		dosci 	= .false.    +0.000
c			= .true.     +0.000 10^+00
			*/
					/* get width of number string */
					if(doasis==ON)
						xv = xxv*xnorm;
					else
						xv = - xxv*xnorm;
					nshft = 1;
					if(xv < 0.0)nshft++;
					if(dosci){
						nshft += 1;
						ndec = 2;
						nshft += ndec;
					} else {
						if(ABS(dxx)> 1){
							ndec = -1;
						} else {
							if(ABS(dxx)< 0.1){
								ndec = 3;
							} else {
								ndec = 2;
							}
							nshft += ndec;
						}
						if(ABS(xv) >= 10.0)
							nshft++;
						if(ABS(xv) >= 100.0)
							nshft++;
						if(ABS(xv) >= 1000.0)
							nshft++;
					}

					xxx = 0.5*(nshft -0.5)*sizex;
					xb = xx - xxx;
					/* number must fit */
					if((xb - xe) > 2.*sizex){
					if(xb>=x0 && (xx+xxx)<=xl ){
					if(dopow){
					number(xb,ypos,sizex,xv,0.0,ndec);
					}
					xe = xb + xxx + sizex;
					}
						ticlen = tlen2;
					} else {
						/* put in intermediate tic */
						ticlen = 0.8*tlen2;
					}
					
				} else {
					ticlen = tlen1;
				}
				plot(xx,y0,3);
				if(ticup)	
					plot(xx,y0+ticlen,2);
				else
					plot(xx,y0-ticlen,2);
			}
		}
	}
/*
c-----
c       put in the title if we put in the numbers
c-----
*/
        if(dopow ){
        	sizexx = 1.2*sizex;
                if(dosci){
			title = (char *)calloc(lx+6,sizeof(char));
			strncpy(title,titlex,lx);
			title[lx]='\0';
			strcat(title," *10 ");
                        gcent(xcen,ypos2,sizexx,title,0.0);
                        xe = xcen + 0.5*sizexx*(lx+4);
                        number(xe,ypos2+0.7*sizexx,0.7*sizexx,
                               (float)(nscalx),0.0,-1);
                } else {
			title = (char *)calloc(lx+1,sizeof(char));
			title[lx]='\0';
			strncpy(title,titlex,lx);
                        gcent(xcen,ypos2,sizexx,title,0.0);
                }
		free(title);
        }                  
}


void doliny(float x0,float y0,float yleng,float ymax,float ymin,
        float sizey,int ticlft,int lablft,int dopow,
	int lly, char *titley)
{
	dolny(x0,y0,yleng,ymax,ymin,sizey,ticlft,lablft,dopow,lly, titley,ON);
}

void dnliny(float x0,float y0,float yleng,float ymax,float ymin,
        float sizey,int ticlft,int lablft,int dopow,
	int lly, char *titley)
{
	dolny(x0,y0,yleng,ymax,ymin,sizey,ticlft,lablft,dopow,lly, titley,OFF);
}

void dolny(float x0,float y0,float yleng,float ymax,float ymin,
        float sizey,int ticlft,int lablft,int dopow,
	int lly, char *titley, int doasis)
{
/*
 * c-----
 * c       put in tic marks along the y-axis
 * c
 * c       x0      R*4     - position of bottom side of axis
 * c       y0      R*4     - position of bottom side of axis
 * c       yleng   R*4     - length of Y axis
 * c       ymax    R*4     - maximum value of number corresponding to
 * c                               top
 * c       ymax    R*4     - maximum value of number corresponding to
 * c                               bottom
 * c       sizey   R*4     - size of numbers, e.g., 10
 * c       ticlft  L       - .true. tics go left
 * c                         .false. tics go right
 * c       lablft  L       - .true. numbers goes to left of axis
 * c                         .false. goes to right
 * c       dopow   L       - .true. put in numbers labels
 * c       ly      I*4     length of y-axis title string
 * c       titley  Ch      y-axis title string
 * c       doasis  L       .true. plot values
 * c                       .false. annotate with negative (useful for depth plot)
 * c-----
*/

float tlen1, tlen2, xpos, xpos2, ymn, ymx, yv, yyv, yval, ynorm;
float dyy;
int ndec, nshft, nshftmx, dosci,n,j ,i,nynorm;
float yc, yl, yb, ye, yyy, ticlen, yy, yp, dyyy;
int ilow, iup, ndyyy, ifac;
char *title;
float ycen;
float sizeyy;
int naxdig;
int nscaly;
int ly;
int maxdig;

	char cmd[11];


	ndec = 2;
	ly = lly;
	if(ly < 0)ly = 0;
			/* 	put in the axis */
	yl = y0 + yleng;
	ycen = y0 + 0.5*yleng  ;
	plot(x0,y0,3); 
	plot(x0,yl,2) ;
			/* set up positioning for tics */
	tlen1 = 0.6*sizey;
	tlen2 = 1.2*sizey;
	if(lablft){
		if(ticlft){
			xpos = x0 - 2*sizey ;
		} else {
			xpos = x0 - sizey ;
		}
		strcpy(cmd,"RIGHT");
	} else {
		if(ticlft){
			xpos = x0 + 1.00*sizey ;
		} else {
			xpos = x0 + 2.0*sizey ;
		}
		strcpy(cmd,"RIGHT");
	}
			/* compute limits for scale */
	ymn = MIN(ymin,ymax);
	ymx = MAX(ymin,ymax);
	yyv = ymn;
			/* safety */
	if(ymn==ymx)ymx = ymn + 1.0;
			/*get some idea of the size of 
				the number to be plotted */
	yval = MAX(ABS(ymn),ABS(ymx));
			/*
c	we set the maximum number of decimal digits as 3, 
c	if the number is greater than 100 or less than 0.01, use
c	scientific notation
c	get the maximum number of digits for the number
			*/
	if(yval<=0.01 && yval>0.0 ){
		dosci = ON;
		ynorm = log10(1.0/yval);
		nynorm = ynorm+1;
		ynorm = pow(10.0, nynorm);
		nscaly = -nynorm;
		maxdig = 4;
	}else if( yval>=10000.0){
		dosci = ON;
		ynorm = log10(yval);
		nynorm = ynorm;
		ynorm = pow(10.0,-nynorm);
		nscaly =  nynorm;
		maxdig = 4;
	} else {
		dosci = OFF;
		ynorm = 1.0;
		nynorm = 0.0;
		nscaly = 0;
		if(yval >= 1. && yval < 10.)
			maxdig = 4;
		else if(yval >= 10.0 && yval < 100.)
			maxdig = 2;
		else if(yval >= 100.0 && yval < 1000.)
			maxdig = 2;
		else
			maxdig = 4;
	}
	if(ymn < 0.0 || ymx < 0.0)
		maxdig = maxdig + 1 ;
	if(!lablft)
		xpos += maxdig*sizey;

			/*
c	choose the maximum number of tics somehow on
c	yleng, and size of numbers and necessary decimals
			*/
	dyy = ynorm*(ymx - ymn)/5.0;
	if(dyy == 0.0)dyy = 1.0;
			/*
c	get start limits for search - we will place at most 10 tics
c	here
			*/
	n = log10(dyy);
	if(dyy<1.0)n = n -1;
	ifac = (int)(pow(10.0, -n)*dyy);
	
/*
c	now get a new minimum value that is some multiple 
c	of ifac 10**n 
*/
	dyy = ifac * pow(10.0, n);
	ilow = ynorm*ymn/dyy -1;
	iup  = ynorm*ymx/dyy +1;
	if(ifac == 1)
		ndyyy = 5;
	else if(ifac == 2)
		ndyyy = 5;
	else
		ndyyy = ifac;

	dyyy = dyy/ndyyy;

			/*
c	loop for labels
c	ye is the end of the previous label
			*/
	yc = ymn;
	plot(x0,y0,3);
			/*
c	here
c	yy is the position of the plot corresponding to a value of yyv
			*/
	ye = y0 ;
	nshftmx = 0;
	for(i=ilow ; i <= iup ; i++){
		for(j=0 ; j <= ndyyy -1 ; j++){
		
		yyv = i*dyy + j*dyyy;

		yyv = yyv/ynorm;
		if(yyv>=ymn && yyv <=ymx){
			yy = y0 + (yyv - ymn)*yleng/(ymx - ymn);
			plot(x0,yy,3);
			if(j==0 ){
				if(doasis == ON)
					yv = yyv * ynorm;
				else
					yv = - yyv * ynorm;
				nshft = 1;
				if(yv< 0.0)nshft++;
				if(dosci){
					nshft +=1;
					ndec = 2;
					nshft += 2;
				} else {
					if(ABS(dyy) > 1){
						ndec = -1;
					} else {
						if(ABS(dyy)< 0.1){
							ndec = 3;
						} else {
							ndec = 2;
						}
						nshft += ndec;
					}
					if(ABS(yv) >= 10.0)
						nshft++;
					if(ABS(yv) >= 100.0)
						nshft++;
					if(ABS(yv) >= 1000.0)
						nshft++;

				}
				if(nshft > nshftmx)nshftmx = nshft;
				yp = yy - 0.5*sizey;
				yb = yy;
				if((yb - ye ) > 2.*sizey){
				if(yy>=MIN(ye,yl) && yy<=MAX(ye,yl)){
				if(dopow){
				mynum(xpos,yp,sizey,yv,0.0,ndec,cmd);
				}
					ye = yb ;
				}
					ticlen = tlen2;
				} else {
					/* put in intermediate tic */
					ticlen = 0.8*tlen2;
				}
			} else {
				ticlen = tlen1;
			}
			plot(x0,yy,3);
			if(ticlft)	
				plot(x0-ticlen,yy,2);
			else
				plot(x0+ticlen,yy,2);
		}
		}
	}
/*
c-----
c       put in the title if we put in the numbers
c-----
*/
	if(lablft){
		if(ticlft){
			xpos = x0 - 2*sizey ;
			naxdig = -(nshftmx+1+2) ;
			xpos2  = x0 + (naxdig)* sizey -1.2*sizey ;
		} else {
			xpos = x0 - sizey ;
			naxdig = -(nshftmx+0+2) ;
			xpos2  = x0 + (naxdig)* sizey -1.2*sizey ;
		}
	} else {
		if(ticlft){
			xpos = x0 + 1.00*sizey ;
			naxdig =  (nshftmx+0+2) ;
			xpos2  = x0 + (naxdig)* sizey +1.0*sizey ;
		} else {
			xpos = x0 + 2.0*sizey ;
			naxdig =  (nshftmx+1+2) ;
			xpos2  = x0 + (naxdig)* sizey +1.0*sizey ;
		}
	}
        if(dopow && ly > 0){
        	sizeyy = 1.2*sizey;
                if(dosci){
			title = (char *)calloc(ly+6,sizeof(char));
			strncpy(title,titley,ly);
			title[ly]='\0';
			strcat(title," *10 ");
                        if(lablft){
                                gcent(xpos2,ycen,sizeyy,title, 90.0);
                                ye = ycen + 0.5*sizeyy*(ly+4);
                                number(xpos2-0.7*sizey,ye,0.7*sizeyy,
   					(float)(nscaly), 90.0,-1);
			} else {
                                gcent(xpos2,ycen,sizeyy,title,-90.0);
                                ye = ycen - 0.5*sizeyy*(ly+4);
                                number(xpos2+0.7*sizeyy,ye,0.7*sizeyy,
                                       (float)(nscaly),-90.0,-1);
                        }
			free(title);
                } else {
			title = (char *)calloc(ly+1,sizeof(char));
			strncpy(title,titley,ly);
			title[ly]='\0';
                        if(lablft){
                                gcent(xpos2,ycen,sizeyy,title, 90.0);
                        } else {
                                gcent(xpos2,ycen,sizeyy,title,-90.0);
                        }
			free(title);
                }
        }                 
}


void  fillit(char *cmd,float rad,float x0,float y0)
{
	int ipatx, ipaty;
	float xlen, ylen, r2, dx, dy, x1, x2, x3, y1, y2 ,y3;
	float ang, sa, ca, xnew, ynew;
	float xval[370], yval[370];
	int i, jj;
/*
c-----
c	fill in a solid symbol
c-----
*/
	ipatx = 0;
	ipaty = 0;
	xlen = 0.01;
	ylen = 0.01;
	r2 = rad ;
	if(strncmp(cmd,"SQ",2)==0){
		dx = r2;
		jj = 0;
		for(i=45;i<=405;i+=90){
			ang = 3.1415927*i/180.0;
			sa = sin(ang);
			ca = cos(ang);
			xnew = x0 + r2*ca;
			ynew = y0 + r2*sa;
			xval[jj] = xnew;
			yval[jj] = ynew;
			jj++;
		}
		shadep(jj,xval,yval);
	} else if(strncmp(cmd,"TR",2)==0){ 
		dx = r2;
		jj = 0;
		for(i=90;i<=450;i+=120){
			ang = 3.1415927*i/180.0;
			sa = sin(ang);
			ca = cos(ang);
			xnew = x0 + r2*ca;
			ynew = y0 + r2*sa;
			xval[jj] = xnew;
			yval[jj] = ynew;
			jj++;
		}
		shadep(jj,xval,yval);
	} else if(strncmp(cmd,"HX",2)==0){ 
		x1 = x0 - 0.866*r2;
		x2 = x0 + 0.866*r2;
		x3 = x0;
		y1 = y0 - 0.500*r2;
		y2 = y0 - 0.500*r2;
		y3 = y0 + r2;
		shadet(x1,y1,x2,y2,x3,y3,ipatx,ipaty,xlen,ylen);
		y1 = y0 + 0.500*r2;
		y2 = y0 + 0.500*r2;
		y3 = y0 - r2;
		shadet(x1,y1,x2,y2,x3,y3,ipatx,ipaty,xlen,ylen);
	} else if(strncmp(cmd,"DI",2)==0){ 
		dx = r2;
		jj = 0;
		for(i=0;i<=360;i+=90){
			ang = 3.1415927*i/180.0;
			sa = sin(ang);
			ca = cos(ang);
			xnew = x0 + r2*ca;
			ynew = y0 + r2*sa;
			xval[jj] = xnew;
			yval[jj] = ynew;
			jj++;
		}
		shadep(jj,xval,yval);
	} else if(strncmp(cmd,"CI",2)==0){ 
		dx = r2;
		jj = 0;
		for(i=0;i<=360;i+=10){
			ang = 3.1415927*i/180.0;
			sa = sin(ang);
			ca = cos(ang);
			xnew = x0 + r2*ca;
			ynew = y0 + r2*sa;
			xval[jj] = xnew;
			yval[jj] = ynew;
			jj++;
		}
		shadep(jj,xval,yval);
	}
}

void  curvit(char *cmd,float rad,float x0,float y0)
{
	int ipatx, ipaty;
	float xlen, ylen, r2, dx, dy, x1, x2, x3, y1, y2 ,y3;
	float ang, sa, ca, xnew, ynew;
	float xval[370], yval[370];
	int i, jj;
/*
c-----
c	fill in a solid symbol
c-----
*/
	ipatx = 0;
	ipaty = 0;
	xlen = 0.01;
	ylen = 0.01;
	r2 = rad ;
	if(strncmp(cmd,"SQ",2)==0){
		dx = r2;
		jj = 0;
		for(i=45;i<=405;i+=90){
			ang = 3.1415927*i/180.0;
			sa = sin(ang);
			ca = cos(ang);
			xnew = x0 + r2*ca;
			ynew = y0 + r2*sa;
			xval[jj] = xnew;
			yval[jj] = ynew;
			jj++;
		}
		drawcv(jj,xval,yval);
	} else if(strncmp(cmd,"TR",2)==0){ 
		dx = r2;
		jj = 0;
		for(i=90;i<=450;i+=120){
			ang = 3.1415927*i/180.0;
			sa = sin(ang);
			ca = cos(ang);
			xnew = x0 + r2*ca;
			ynew = y0 + r2*sa;
			xval[jj] = xnew;
			yval[jj] = ynew;
			jj++;
		}
		drawcv(jj,xval,yval);
/*
	} else if(strncmp(cmd,"HX",2)==0){ 
		x1 = x0 - 0.866*r2;
		x2 = x0 + 0.866*r2;
		x3 = x0;
		y1 = y0 - 0.500*r2;
		y2 = y0 - 0.500*r2;
		y3 = y0 + r2;
		shadet(x1,y1,x2,y2,x3,y3,ipatx,ipaty,xlen,ylen);
		y1 = y0 + 0.500*r2;
		y2 = y0 + 0.500*r2;
		y3 = y0 - r2;
		shadet(x1,y1,x2,y2,x3,y3,ipatx,ipaty,xlen,ylen);
*/
	} else if(strncmp(cmd,"DI",2)==0){ 
		dx = r2;
		jj = 0;
		for(i=0;i<=360;i+=90){
			ang = 3.1415927*i/180.0;
			sa = sin(ang);
			ca = cos(ang);
			xnew = x0 + r2*ca;
			ynew = y0 + r2*sa;
			xval[jj] = xnew;
			yval[jj] = ynew;
			jj++;
		}
		drawcv(jj,xval,yval);
	} else if(strncmp(cmd,"CI",2)==0){ 
		dx = r2;
		jj = 0;
		for(i=0;i<=360;i+=10){
			ang = 3.1415927*i/180.0;
			sa = sin(ang);
			ca = cos(ang);
			xnew = x0 + r2*ca;
			ynew = y0 + r2*sa;
			xval[jj] = xnew;
			yval[jj] = ynew;
			jj++;
		}
		drawcv(jj,xval,yval);
	}
}

void drawcv(int jj, float xval[], float yval[])
{
	int i;
	plot(xval[0], yval[0], 3);
	for(i=0 ; i < jj ; i++){
		plot(xval[i], yval[i], 2);
	}
	plot(xval[jj], yval[jj], 3);
}
		


#define      Xmin    01
#define      Xmax    02
#define      Ymax   010
#define      Ymin    04
#define      None    00      

void rclip(float xmin,float xmax,float ymin,float ymax,
		float *xc1,float *yc1,float *xc2,float *yc2,
		float x0,float yy0,float x1,float yy1,int *iplt)
{
/*
c-----
c	perform a line clipping for the plot region
c-----
c	xmin	R*4	- minimum value of X
c	xmax	R*4	- maximum value of X
c	ymin	R*4	- minimum value of Y
c	ymax	R*4	- maximum value of Y
c	xc1	R*4	- coordinate #1 of plotted line segment
c	yc1	R*4
c	xc2	R*4	- coordinate #2 of plotted line segment
c	yc2	R*4
c	x0	R*4	- first coordinate
c	yy0	R*4
c	x1	R*4	- second coordinate
c	yy1	R*4
c	iplt	I*4	- > 0 plot the segment, otherwise do not
c-----
c-----
c	C implementation of Cohen and Sutherland Line Clipping Routine with floats
c-----
*/
	int c0, c1, c;
	float xx0, yz0, xx1, yz1;
	float slope, slopeinv;
	float x, y;

	*iplt = 0;
	xx0 = x0;
	yz0 = yy0;
	xx1 = x1;
	yz1 = yy1;
	linecode(xx0,yz0,ymax,xmin,ymin,xmax,&c0);
	linecode(xx1,yz1,ymax,xmin,ymin,xmax,&c1);
	if( xx0 != xx1)
		slope = (yz1-yz0)/(xx1-xx0);
	if( yz0 != yz1) 
		slopeinv = (xx1-xx0)/(yz1-yz0);
	while(c0 != 0 && c1 != 0){
		if((c0&c1) !=0)
			return;
		if(c0 == c1)
			return;
                if(c0 == 0 ){
                        c = c1;
                } else { 
                        c = c0;
		}
		if(c & Xmin){ 
                        y = yz0 + slope*(xmin-xx0);
                        x = xmin;
                }
                if( c & Xmax){ 
			y = yz0 + slope*(xmax-xx0);
                        x = xmax;
                }
                if( c & Ymax){ 
                        x = slopeinv*(ymax-yz0) + xx0;
                        y = ymax;
                }
                if( c & Ymin){  
                        x = slopeinv*(ymin-yz0) + xx0;
                        y = ymin;
                }
                if(c == c0){ 
                        xx0 = x;
                        yz0 = y;
			linecode(xx0,yz0,ymax,xmin,ymin,xmax,&c0);
                } else { 
                        xx1 = x;
                        yz1 = y;
			linecode(xx1,yz1,ymax,xmin,ymin,xmax,&c1);
                }
	}
        *xc1 = xx0;
        *xc2 = xx1;
        *yc1 = yz0;
        *yc2 = yz1;
	*iplt = 1;
}

static void linecode(float x,float y,float ymax,float xmin,float ymin,float xmax,int *c)
{
	
	*c = None;
	if(x < xmin) {
		*c = Xmin;
	} else if(x  >  xmax){
		*c = Xmax;
	}
	if(y < ymin) {
		*c |= Ymin;
	} else if(y > ymax){
		*c |= Ymax;
	}
}

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


/* reproduce here the callib.h from CALPLOT beware of incompatibilities */
/* include file to get control parameters from CALPLOT C procedures */
#ifndef _CALLIB_
#define _CALLIB_
struct Gcplot {
	float owidth;
	float xstp;
	float ystp;
	float xold;
	float yold;
	float xcur;
	float ycur;
	int iunit;
	int ifont;
	float xs;
	float ys;
	float x0;
	float y0;
	int is;		/* index to bit is pattern ipat - plotd */
	int ipatld;	/*  previous pattern used to init - plotd */
	float xsv;	/* offset within current bin in xlen unit - plotd */
	} ;
#endif
extern struct Gcplot Gcplt;

static void mynum (float xpage,float ypage,float height,
		float fpn,float angle,int ndec,char *cmd)
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
	char cexp[4];
	char num[21];
/*
c-----
c    to avoid dregs in memory null first character
c-----
*/
	int i, n, ii, nzflag, mdec, iex, jex;
	int nman, nexp;
	float fp, aex;
	float ang, ca, sa, xl, yl, x0, y0,ht;
	float xx0, yy0, xxx0, yyy0;


	if(xpage < 999.0)x0 = xpage;
	if(ypage < 999.0)y0 = ypage;

	/* to avoid dregs in memory null all characters */
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
		if( ndec < 0){
			mdec = mdec -1;
		}
                if(ndec == -1){
                        nzflag = 1;
                        mdec = 1;
                }
                ffpack(num,&n,fp,mdec,19,nzflag);
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
                ffpack(num,&n,fp,mdec,14,nzflag);
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
		if(strcmp(cmd,"LEFT") == 0){
                	symbol(x0,y0,height,num,angle,n);
		} else if(strcmp(cmd,"CENTER")==0){
                	symbol(x0-0.5*n*height,y0,height,num, angle,n);
		} else if(strcmp(cmd,"RIGHT")==0){
                	symbol(x0-n*height,y0,height,num,angle,n);
		}
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
		xx0 = x0;
		yy0 = y0;
        	symbol(x0,y0,height,num,angle,nman+3);
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

static void ffpack(char *num,int *n,float fpv,int mdec,int mwid,int nzflag)
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

void msupsc (float xpage,float ypage,float height,char *str,int npow,float angle)
{
/*
c.....     xpage,ypage coordinates of lower left corner of number.
c.....     height   height of plotted number.
c.....     str	Character string
c.....     npow power
c.....     fpn      floating point number to be plotted.
c.....     angle    angle at which number is plotted, in degrees.
*/

	float x0, y0, xx0, yy0, ca, sa, xl, yl, ang;
	float xcur, ycur, xxx0, yyy0, ht;
	int il;
	if(xpage < 999.0) {
		x0 = xpage;
	} else {
		x0 = Gcplt.x0;
	}
	if(ypage < 999.0) {
		y0 = ypage;
	} else {
		y0 = Gcplt.y0;
	}
	xx0 = x0;
	yy0 = y0;
	il = strlen(str);
        symbol(x0,y0,height,&str[0], angle,il);
	/* get the proper position because of rotation */
                ang = angle*3.1415927/180.0;
                ca = cos(ang);
                sa = sin(ang);
                xl = height*(il);
                yl = 0.6*height;
		xx0 = xx0 + xl*ca;
		yy0 = yy0 + xl*sa;
                xxx0 = xx0  - yl*sa;
                yyy0 = yy0  + yl*ca;
                ht = 0.7*height;
                number(xxx0,yyy0,ht,(float)npow,angle,-1);
		/* now position at end of proper string */
		xl =  3.0*ht;
		xcur = xx0 + xl*ca;
		ycur = yy0 + xl*sa;
		x0 = xcur;
		y0 = ycur;
}

void gsubsc(float x,float y,float ht,char *s1,int n1,char *s2,int n2)
{
	int infont;
	int nn, n;
	infont = Gcplt.ifont;
	if(n1 < 0) {
		gfont(4);
		nn = -n1;
	} else {
		nn = n1;
		gfont(infont);
	}
	if(nn > 0)symbol(x,y,ht,s1,0.0,nn);
	if(n2 < 0){
		gfont(4) ;
		n = -n2 ;
	} else {
		n = n2 ;
		gfont(infont) ;
	}
	if(n >0)symbol(x+nn*ht+0.4*ht,y-0.4*ht,0.7*ht,s2,0.0,n);
	if(n1 < 0 || n2 < 0)gfont(infont);
}

void gsupsc(float x, float y,float ht,char *s1,int n1,char *s2,int n2)
{
	int infont;
	int nn, n;
	infont = Gcplt.ifont;
	if(n1 < 0){
		gfont(4);
		nn = -n1;
	} else {
		nn = n1;
		gfont(Gcplt.ifont);
	}
	if(nn > 0)
		symbol(x,y,ht,s1,0.0,nn);
	if(n2 < 0){
		gfont(4);
		n = -n2;
	} else {
		n = n2;
		gfont(infont);
	}
	if(n > 0)symbol(x+nn*ht+0.4*ht,y+0.4*ht,0.7*ht,s2,0.0,n);
	if(n1 < 0  ||  n2 < 0)gfont(infont);
}

void gsubsup(float x,float y,float ht,char *s1,int n1,char *s2,int n2,char *s3,int n3)
{
	int infont;
	int nn, n;
	infont = Gcplt.ifont;
	if(n1 < 0){
		gfont(4);
		nn = -n1;
	} else {
		nn = n1;
		gfont(infont);
	}
	if(nn > 0)
		symbol(x,y,ht,s1,0.0,nn);
	if(n2 < 0){
		gfont(4);
		n = -n2;
	} else {
		n = n2;
		gfont(infont);
	}
	if(n > 0)
		symbol(x+nn*ht+0.4*ht,y-0.4*ht,0.7*ht,s2,0.0,n);
	if(n3 < 0){
		gfont(4);
		n = -n3;
	} else {
		n = n3;
		gfont(infont);
	}
	if(n > 0)symbol(x+nn*ht+0.4*ht,y+0.4*ht,0.7*ht,s3,0.0,n);
	if(n1 < 0  ||  n2 < 0 || n3 < 0)gfont(infont);
}
