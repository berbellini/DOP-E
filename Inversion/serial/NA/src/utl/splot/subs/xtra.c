/* xtra.c - typset simulation using OpenWindows fonts */

/* X include files */

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include <X11/X.h>
#include <X11/keysym.h>
#include <stdio.h>
#include <math.h> 
#include <sys/types.h>
#include <sys/uio.h>
#include <strings.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <utmp.h>

/* Program macros and parameter constants */

#define ABS(x) (((x) < 0) ? -x : (x)) 

/* Fortran Common Block equivalents (defined in xpakw.c) */
extern struct {int lplot,irot,il34; float a,b,c,d,asp,theta;} p00000_;
extern struct {int xorig,yorig;} p00001_;
extern struct {float a1,a2,b1,b2,c1,c2,d1,d2;} p00002_;
extern struct {float psca; int ixo, iyo, iox, ioy;} a00000_;

/* External variables (defined in xpakw.c) */
 
extern Display        *mydisplay;
extern Window         mywindow;
extern GC             mygc;       
extern Pixmap         mypixmap;
extern XSizeHints     myhint;
extern Colormap       cmap; 
extern XColor         color[];
extern unsigned long  pixels[];
extern unsigned long  plane_masks;
extern int            ncol;
extern unsigned int   ncolors;
extern XEvent         myevent;
extern KeySym         mykey;
extern Cursor         mycursor;
extern int            myscreen;
extern unsigned int   depth, pixheight, pixwidth;
extern int            mylinestyle;
extern int            mylinewidth;
         
/* 'LOCAL' GLOBAL VARIABLES */
Font           tfont1,tfont2,tfont3,tfont4,tfont5;
char           *fntnm11,*fntnm12,*fntnm13,*fntnm14,*fntnm15;
char           *fntnm21,*fntnm22,*fntnm23,*fntnm24,*fntnm25;
char           *fntnm31,*fntnm32,*fntnm33,*fntnm34,*fntnm35;
char           *fntnm41,*fntnm42,*fntnm43,*fntnm44,*fntnm45;
char           *fntnm51,*fntnm52,*fntnm53,*fntnm54,*fntnm55;
char           *fntnm61,*fntnm62,*fntnm63,*fntnm64,*fntnm65;
char           *fntnm71,*fntnm72,*fntnm73,*fntnm74,*fntnm75;
char           *fntnm81,*fntnm82,*fntnm83,*fntnm84,*fntnm85;
char           *fntnm91,*fntnm92,*fntnm93,*fntnm94,*fntnm95;
char           *fntnmA1,*fntnmA2,*fntnmA3,*fntnmA4,*fntnmA5; 
  
/* Subroutine library */
                                                 
typset_(xa,ya) 

float *xa, *ya;

/* set up font names from Openwin library */
                              
{
 
 fntnm11 = "lucidasans-10";
 fntnm12 = "lucidasans-12";
 fntnm13 = "lucidasans-14";
 fntnm14 = "lucidasans-18";
 fntnm15 = "lucidasans-24";

 fntnm21 = "lucidasans-bold-10";
 fntnm22 = "lucidasans-bold-12";
 fntnm23 = "lucidasans-bold-14";
 fntnm24 = "lucidasans-bold-18";
 fntnm25 = "lucidasans-bold-24";
  
 fntnm31 = "palatino-roman-10";
 fntnm32 = "palatino-roman-12";
 fntnm33 = "palatino-roman-14";
 fntnm34 = "palatino-roman-18";
 fntnm35 = "palatino-roman-24";

 fntnm41 = "times-roman-10";
 fntnm42 = "times-roman-12";
 fntnm43 = "times-roman-14";
 fntnm44 = "times-roman-18";
 fntnm45 = "times-roman-24";

 fntnm51 = "lucidasans-italic-10";
 fntnm52 = "lucidasans-italic-12";
 fntnm53 = "lucidasans-italic-14";
 fntnm54 = "lucidasans-italic-18";
 fntnm55 = "lucidasans-italic-24";

 fntnm61 = "palatino-italic-10";
 fntnm62 = "palatino-italic-12";
 fntnm63 = "palatino-italic-14";
 fntnm64 = "palatino-italic-18";
 fntnm65 = "palatino-italic-24";

 fntnm71 = "times-italic-10";
 fntnm72 = "times-italic-12";
 fntnm73 = "times-italic-14";
 fntnm74 = "times-italic-18";
 fntnm75 = "times-italic-24";

 fntnm81 = "palatino-bold-10";
 fntnm82 = "palatino-bold-12";
 fntnm83 = "palatino-bold-14";
 fntnm84 = "palatino-bold-18";
 fntnm85 = "palatino-bold-24";

 fntnm91 = "palatino-bolditalic-10";
 fntnm92 = "palatino-bolditalic-12";
 fntnm93 = "palatino-bolditalic-14";
 fntnm94 = "palatino-bolditalic-18";
 fntnm95 = "palatino-bolditalic-24";

 fntnmA1 = "symbol-10";
 fntnmA2 = "symbol-12";
 fntnmA3 = "symbol-14";
 fntnmA4 = "symbol-18";
 fntnmA5 = "symbol-24";

/* set up default font no 1 */
 
 tfont1 = XLoadFont(mydisplay,fntnm11);
 tfont2 = XLoadFont(mydisplay,fntnm12);
 tfont3 = XLoadFont(mydisplay,fntnm13);
 tfont4 = XLoadFont(mydisplay,fntnm14);
 tfont5 = XLoadFont(mydisplay,fntnm15);

 XSetFont(mydisplay,mygc,tfont1);  
 XSetFont(mydisplay,mygc,tfont2);  
 XSetFont(mydisplay,mygc,tfont3);  
 XSetFont(mydisplay,mygc,tfont4);  
 XSetFont(mydisplay,mygc,tfont5);  

}
/*------------------------------*/
zpick_(ifnm,iall,istati)

/* assign requisite font name and sizes into Xfonts */

int *ifnm, *iall, *istati;

{
int ifn2;


XUnloadFont(mydisplay,tfont1);
XUnloadFont(mydisplay,tfont2);
XUnloadFont(mydisplay,tfont3);
XUnloadFont(mydisplay,tfont4);
XUnloadFont(mydisplay,tfont5);

ifn2 = *ifnm;

switch (ifn2){

 case (0) :
     tfont1 = XLoadFont(mydisplay,fntnm11);
     tfont2 = XLoadFont(mydisplay,fntnm12);
     tfont3 = XLoadFont(mydisplay,fntnm13);
     tfont4 = XLoadFont(mydisplay,fntnm14);
     tfont5 = XLoadFont(mydisplay,fntnm15);
 break;
 case (1) :
     tfont1 = XLoadFont(mydisplay,fntnm11);
     tfont2 = XLoadFont(mydisplay,fntnm12);
     tfont3 = XLoadFont(mydisplay,fntnm13);
     tfont4 = XLoadFont(mydisplay,fntnm14);
     tfont5 = XLoadFont(mydisplay,fntnm15);
 break;
 case (2) :
     tfont1 = XLoadFont(mydisplay,fntnm21);
     tfont2 = XLoadFont(mydisplay,fntnm22);
     tfont3 = XLoadFont(mydisplay,fntnm23);
     tfont4 = XLoadFont(mydisplay,fntnm24);
     tfont5 = XLoadFont(mydisplay,fntnm25);
 break;
 case (3) :
     tfont1 = XLoadFont(mydisplay,fntnm31);
     tfont2 = XLoadFont(mydisplay,fntnm32);
     tfont3 = XLoadFont(mydisplay,fntnm33);
     tfont4 = XLoadFont(mydisplay,fntnm34);
     tfont5 = XLoadFont(mydisplay,fntnm35);
 break;
 case (4) :
     tfont1 = XLoadFont(mydisplay,fntnm41);
     tfont2 = XLoadFont(mydisplay,fntnm42);
     tfont3 = XLoadFont(mydisplay,fntnm43);
     tfont4 = XLoadFont(mydisplay,fntnm44);
     tfont5 = XLoadFont(mydisplay,fntnm45);
 break;
 case (5) :
     tfont1 = XLoadFont(mydisplay,fntnm51);
     tfont2 = XLoadFont(mydisplay,fntnm52);
     tfont3 = XLoadFont(mydisplay,fntnm53);
     tfont4 = XLoadFont(mydisplay,fntnm54);
     tfont5 = XLoadFont(mydisplay,fntnm55);
 break;
 case (6) :
     tfont1 = XLoadFont(mydisplay,fntnm61);
     tfont2 = XLoadFont(mydisplay,fntnm62);
     tfont3 = XLoadFont(mydisplay,fntnm63);
     tfont4 = XLoadFont(mydisplay,fntnm64);
     tfont5 = XLoadFont(mydisplay,fntnm65);
 break;
 case (7) :
     tfont1 = XLoadFont(mydisplay,fntnm71);
     tfont2 = XLoadFont(mydisplay,fntnm72);
     tfont3 = XLoadFont(mydisplay,fntnm73);
     tfont4 = XLoadFont(mydisplay,fntnm74);
     tfont5 = XLoadFont(mydisplay,fntnm75);
 break;
 case (8) :
     tfont1 = XLoadFont(mydisplay,fntnm81);
     tfont2 = XLoadFont(mydisplay,fntnm82);
     tfont3 = XLoadFont(mydisplay,fntnm83);
     tfont4 = XLoadFont(mydisplay,fntnm84);
     tfont5 = XLoadFont(mydisplay,fntnm85);
 break;
 case (9) :
     tfont1 = XLoadFont(mydisplay,fntnm91);
     tfont2 = XLoadFont(mydisplay,fntnm92);
     tfont3 = XLoadFont(mydisplay,fntnm93);
     tfont4 = XLoadFont(mydisplay,fntnm94);
     tfont5 = XLoadFont(mydisplay,fntnm95);
 break;
 case (10) :
     tfont1 = XLoadFont(mydisplay,fntnmA1);
     tfont2 = XLoadFont(mydisplay,fntnmA2);
     tfont3 = XLoadFont(mydisplay,fntnmA3);
     tfont4 = XLoadFont(mydisplay,fntnmA4);
     tfont5 = XLoadFont(mydisplay,fntnmA5);
 break;
} /*switch*/
}
/*---------------------------*/
typstr_(x,y,size,iword,angl,nchar,dum)

/* write a string on the plot (plotter units) in font
   selected by zpick*/

float         *x,*y,*size,*angl;
char          *iword;
int           *nchar;
int           dum;
                     
{
static        float xsiz;
float         xp, yp;
short         i, ixv, iyv;
XGCValues     values;
unsigned long valuemask;

xp = *x;
yp = *y;
if(p00000_.irot != 0) {
  yp = *x;
  xp = -*y;
  if(p00000_.il34 == 0) {
    if(*size > 0) xp = 27.2 - *y ;}
  else if(p00000_.il34 == 1) {
    if(*size > 0) xp = 40.1 - *y ;} }

ixv = (int) (xp * a00000_.psca);
iyv = (int) (yp * a00000_.psca);

if (*size < 0.0) {
  ixv = a00000_.iox + ixv;
  iyv = a00000_.ioy - iyv;} 
else {
  ixv = p00001_.xorig + ixv;
  iyv = a00000_.iyo - p00001_.yorig - iyv;}

if(ABS(*size) != xsiz) xsiz=ABS(*size); 
  if(xsiz <= 0.25) XSetFont(mydisplay,mygc,tfont1);
  else if(xsiz <= 0.35) XSetFont(mydisplay,mygc,tfont2);  
  else if(xsiz <= 0.45) XSetFont(mydisplay,mygc,tfont3);
  else if(xsiz <= 0.55) XSetFont(mydisplay,mygc,tfont4);
  else XSetFont(mydisplay,mygc,tfont5); 

XDrawString(mydisplay,mywindow,mygc,ixv,iyv,iword,*nchar);

}

/*-------------------------------------*/
typnum_(x,y,size1,fnum1,angle,ndec1,dum)
      
/* write a number on the plot (plotter units) in font
   selected by zpick
   format if nsf > 0 , nsf is number of dec places
          if nsf < 0 , rn is fixed to an integer*/

float *x, *y, *size1, *fnum1, *angle;
int *ndec1;
int dum;

{                   
  float size, rn, pow;
  int ndec, iadd, np, nsf;
  int itot, idpl, ipt, ir;
  char iform[7];
  char iword[30];
  char *iptr;
                    
  ndec = *ndec1;
  size = ABS(*size1);
  rn = *fnum1;
  iadd = 1;
  if (ndec >= 0) iadd = 2;
  if (rn < 0.0) iadd = iadd +1;
  if (ABS(rn) < 1.0) iadd = iadd + 1;
  if (rn == 0.0) {
    if (ndec < 0) np = 1;
    if (ndec >= 0) np = 2;} 
   else {
    pow = log10(ABS(rn));
    np = (int) (pow + 0.01) + iadd;}
  if (ndec < 0) nsf = -np;
  if (ndec >= 0) nsf = 10 * (np + ndec) +ndec;

  if (nsf >= 0) {
    itot = nsf / 10;
    idpl = nsf - itot * 10;
    sprintf(iform,"%%%d.%df",itot,idpl);
    sprintf(iword,iform,rn);}
  else {              
    itot = -nsf;
    sprintf(iform,"%%%dd",itot);
    ir = (int) rn;
    sprintf(iword,iform,ir);}
       
  iptr = &iword[0];

  typstr_(x,y,&size,iptr,angle,&itot,dum);
  
}               
/*-------------------------------------*/
symnum_(x,y,size,rn,angl,nsf,dum)

/* write a number on the plot (plotter units) - standard font
   format if nsf > 0 , nsf is number of dec places
          if nsf < 0 , rn is fixed to an integer*/

float *x,*y,*size,*rn,*angl;
int *nsf;
int dum;

{
int itot, idpl, ipt, ir;
int *it;
char iform[7];
char iword[30];
char *iptr;
 
if (*nsf > 0 ) { 
  idpl = *nsf;
  ipt  = (int) log10(*rn) + 1;
  if(ipt > 0) {
    itot = ipt + idpl + 1 ;}
  else {
    itot = idpl + 2 ;}
  sprintf(iform,"%%%d.%df",itot,idpl);
  sprintf(iword,iform,*rn);}
else {
  ipt = (int) log10(*rn) + 1;
  ir = (int) *rn;           
  sprintf(iword,"%d",ir);}
       
iptr = &iword[0];

symbol_(x,y,size,iptr,angl,it,dum);
}
