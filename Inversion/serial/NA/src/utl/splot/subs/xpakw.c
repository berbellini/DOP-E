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

#define ABS(x) (((x) < 0) ? -(x) :(x))
#define PIXELS    230
#define MAXDOT    500

/* Fortran Common Block equivalents */

struct {int lplot,irot,il34; float a,b,c,d,asp,thet;} p00000_;
struct {int xorig,yorig;} p00001_;
struct {float a1,a2,b1,b2,c1,c2,d1,d2;} p00002_;
struct {float psca; int ixo, iyo, iox, ioy;} a00000_;

/* Global variables */

Display          *mydisplay;
Window	 	 mywindow;
GC               mygc;                     
Pixmap           mypixmap;
XSizeHints       myhint;
Font             font1,font2,font3,font4,font5;
Colormap         cmap;
XColor           color[PIXELS];
XColor           colorgrey[PIXELS];
unsigned long    pixels[PIXELS];
unsigned long    plane_masks;
int              ncol;
short            xdot[MAXDOT][PIXELS],ydot[MAXDOT][PIXELS];
unsigned int     ncolors;
XEvent           myevent;
KeySym           mykey;
Cursor           mycursor;
int              myscreen;
unsigned int     depth, pixheight, pixwidth;
int              mylinestyle;
int              mylinewidth;
unsigned int     i2max();
unsigned int     i2min();
char    	 mywindowname[256];

static char      dash1[]={3,3};
static char      dash2[]={6,3};
static char      dash3[]={9,3};
static char      dash4[]={9,3,3,3};
static char      dash5[]={9,3,6,3};
static char      dash6[]={9,3,6,3,6,3};

static int       red[] =   {255,  0,  0,255, 63, 31,255,255,
                            191,255,159,255,255,159,159,159,
                             15, 31, 47, 63, 79, 95,111,127,
                            143,159,175,191,207,223,239,255};
static int       green[] = {255,  0,  0, 31, 63,225,255,159,
                            127,191, 63,159,159,159,255,159,
                             15, 31, 47, 63, 79, 95,111,127,
                            143,159,175,191,207,223,239,255};                          static int       blue[] =  {255,  0, 63, 31,255, 31, 31, 31,
                             63, 97,159,159,255,255,159,159,
                             15, 31, 47, 63, 79, 95,111,127,
                            143,159,175,191,207,223,239,255};

/* Subroutine Library */
italic_(theta)
float *theta;
{return;}

/*-----------------------------------------------------------------*/
hplots_(ion,iro,lpl,ils)

/* Hplots initializes the graphics process:
        if ils = 0 - a4 window
        if ils = 1 - a3 window
        if ils = 2 - a3 window
   and initializes handshaking characteristics (ion=1,2) or terminates
   plotting (ion=0). */ 

int *ion,*iro,*lpl,*ils;

{ 

/* declarations */
unsigned long    myforeground, mybackground;
int              i;
char             text[10],dspnm[40];
int              done;
int              sync;

int count;  
int maxnames;
char **list;
char pattern[80];

if (*ion == 0) goto finish;

/* Initialize common block variables */
            
p00000_.lplot = *lpl;
p00000_.irot = *iro;
p00000_.il34 = *ils;
p00000_.a = 1.0;
p00000_.b = 0.0;
p00000_.c = 1.0;
p00000_.d = 0.0;
p00000_.asp = 0.66666;

p00001_.xorig = 0;
p00001_.yorig = 0;
ncol = 0;

/* Establish the current environment */

i =40; 
dispnm_(dspnm,&i,i);

/* Fix bug with opening a window */ 
/* on an x-terminal */ 
/* mydisplay=XOpenDisplay(dspnm); */

mydisplay=XOpenDisplay(NULL);

myscreen=DefaultScreen(mydisplay);
mybackground=XWhitePixel(mydisplay,myscreen);
myforeground=XBlackPixel(mydisplay,myscreen);  
mylinewidth = 1;                
depth=XDefaultDepth(mydisplay,myscreen);
/*mylinestyle=LineSolid;
XSetLineAttributes (mydisplay,mygc,mylinewidth,mylinestyle,
                      CapButt,JoinBevel);*/                            
sync = 1;
/*XSynchronize(mydisplay,sync);*/

/* Load fonts */
/*maxnames = 1000;
strcpy(pattern,"*");
list = XListFonts(mydisplay,pattern,maxnames,&count);
for (i = 0; i < count; i++) printf("%s, %d\n", list[i],i);*/
/*font1=XLoadFont(mydisplay,"6x10");
font2=XLoadFont(mydisplay,"6x10");
font3=XLoadFont(mydisplay,"6x13");
font4=XLoadFont(mydisplay,"6x13");
font5=XLoadFont(mydisplay,"9x15");
font1=XLoadFont(mydisplay,"lucidasans-10");
font2=XLoadFont(mydisplay,"lucidasans-12");
font3=XLoadFont(mydisplay,"lucidasans-14");
font4=XLoadFont(mydisplay,"lucidasans-18");
font5=XLoadFont(mydisplay,"lucidasans-24");*/
 
/* Perform appropriate scaling */

                                 
if (*ils==0) {
  a00000_.psca = 32.0;
  myhint.width = (int) (a00000_.psca * 27.2);                                 
  myhint.height = (int) (a00000_.psca * 19.1);}
else {
  a00000_.psca = 25.0;
  myhint.width = (int) (a00000_.psca * 40.1);                                 
  myhint.height = (int) (a00000_.psca * 27.4);}
                                     
myhint.x = (1024 - myhint.width) / 2;                                 
myhint.y = (800 - myhint.height) / 2;       
                          
a00000_.ixo = 0;                                 
a00000_.iyo = myhint.height;
a00000_.iox = 0;                                 
a00000_.ioy = 0;                                 

myhint.flags=PPosition|PSize;

/* Create colour map */

CreateColourMap();

/* Create window and set default properties */
mywindow=XCreateSimpleWindow (mydisplay, DefaultRootWindow (mydisplay), 
                              myhint.x,myhint.y,myhint.width,myhint.height, 
                              2,color[3].pixel,mybackground);
   
XSetStandardProperties(mydisplay,mywindow,mywindowname,mywindowname,
  None,0,0,&myhint);

mycursor=XCreateFontCursor(mydisplay,XC_center_ptr);
XDefineCursor(mydisplay,mywindow,mycursor);
mygc=XCreateGC(mydisplay,mywindow,0,0);
XSetBackground(mydisplay,mygc,mybackground);
XSetForeground(mydisplay,mygc,myforeground);

/* input event selection */
XSelectInput(mydisplay,mywindow,
  ButtonPressMask|KeyPressMask|ExposureMask);

/* Map the window to the screen */
XMapRaised(mydisplay,mywindow);
 
myevent.type = NoExpose;
while (myevent.type != Expose)   /* Wait for initial expose event to occur */
                             XNextEvent (mydisplay,&myevent);                  
return;

finish:

/* Send message to screen */

printf ("Click on window to exit\n");

/* printf(" got to finish ion = %d \n",*ion); */

/* Copy window to pixmap after plotting */

depth = XDefaultDepth(mydisplay,myscreen);
pixheight = (unsigned) myhint.height;
pixwidth = (unsigned) myhint.width;
mypixmap=XCreatePixmap(mydisplay,mywindow,
                       pixwidth,pixheight,depth);
XCopyArea(mydisplay,mywindow,mypixmap,mygc,0,0,
          myhint.width,myhint.height,0,0);

/* Redraw window if obscured */

done=0;
while(done==0)
{
/* read the next event */
  XNextEvent(mydisplay,&myevent);
  switch(myevent.type)
  {
    /* repaint window on expose events */
    case Expose:
      if (myevent.xexpose.count==0)
        XCopyArea(mydisplay,mypixmap,mywindow,mygc,0,0,
                  myhint.width,myhint.height,0,0);
/*      if (myevent.xexpose.count==0)
      XDrawImageString(
        myevent.xexpose.display,myevent.xexpose.window,mygc,
        50,50,
        hello,strlen(hello));     */
      break;

    /* process keyboard mapping changes */
    case MappingNotify:
      XRefreshKeyboardMapping(&myevent);
      break;

    case ButtonPress:
      done=1; break;

    /* process keyboard input 
    case KeyPress:
      i=XLookupString(&myevent,text,10,&mykey,0);
      if (i==1 && text[0]=='q') done=1;
      break;*/

    } /* switch(myevent.type) */
} /* while (done==0) */
/* termination */             
ncolors=PIXELS;
XFreePixmap(mydisplay,mypixmap);
if (ncol != 0 ) {
if (depth <9) XFreeColors (mydisplay,cmap,pixels,ncol,plane_masks);
}
XFreeGC(mydisplay,mygc);
XFreeColormap(mydisplay,cmap);
XDestroyWindow(mydisplay,mywindow);
XCloseDisplay(mydisplay);
exit(0);
}   /* hplots */ 


/* ------------------------------------------------- */
scale_(xmin,xmax,px1,px2,ymin,ymax,py1,py2)

/* Set up scale factors used to go between user and plotter coordinates */ 

float     *xmin,*xmax,*ymin,*ymax;
float     *px1,*px2,*py1,*py2;

{
  p00000_.a = (*px2 - *px1) / (*xmax - *xmin);
  p00000_.b = *px1 - p00000_.a * (*xmin);
  p00000_.c = (*py2 - *py1) / (*ymax - *ymin);
  p00000_.d = *py1 - p00000_.c * (*ymin);

  p00002_.a1 = *xmin;
  p00002_.a2 = *xmax;
  p00002_.b1 = *px1;
  p00002_.b2 = *px2;
  p00002_.c1 = *ymin;
  p00002_.c2 = *ymax;
  p00002_.d1 = *py1;
  p00002_.d2 = *py2;

} 
/*------------------------------------------------------*/
plot_(x,y,i)

/* Raises (i=3) or lowers (i=2) pen and moves to coordinates (x,y)
   if i > 0 or to current position plus (x,y) if i < 0 */

float     *x,*y;
int       *i;

{
float            xp, yp;
short            ixv, iyv;

xp = *x;
yp = *y;
  
if(p00000_.irot != 0){
  yp = *x;
  xp = -*y;
  if(p00000_.il34 == 0){
    if(*i > 0) xp = 27.2 - *y;}
  else if(p00000_.il34 == 1){
    if(*i > 0) xp = 40.1 - *y;}
}
ixv = (int) (xp *  a00000_.psca);
iyv = (int) (yp *  a00000_.psca);

if (*i < 0) {
  ixv = a00000_.iox + ixv;
  iyv = a00000_.ioy - iyv;}
else {
  ixv = ixv + p00001_.xorig;
  iyv = a00000_.iyo - p00001_.yorig - iyv;}
 
if(ABS(*i) == 2) {
  XDrawLine (mydisplay,mywindow,mygc,
             a00000_.iox,a00000_.ioy,ixv,iyv);}

/* Update the current position */

a00000_.iox = ixv;
a00000_.ioy = iyv;
}
/*------------------------------------------------------*/
plotu_(x,y,i)

/* Scales user coordinates to plotter coordinates and plots as in plot */

float     *x,*y;
int       *i;

{
float            xp, yp;

xp = *x * p00000_.a;
yp = *y * p00000_.c;

if(*i > 0) {
  xp = xp + p00000_.b;
  yp = yp + p00000_.d;}

plot_(&xp,&yp,i);
}
/*------------------------------------------------------*/
dashln_(ldash,lpat)

int              *ldash, *lpat;

{        
if(*ldash < 0 || *ldash >= 12){
  mylinestyle=LineSolid;
  XSetLineAttributes (mydisplay,mygc,mylinewidth,mylinestyle,
                      CapButt,JoinBevel);}
  
switch (*ldash){

  case (0):
    XSetDashes (mydisplay,mygc,3,dash1,2);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,mygc,mylinewidth,mylinestyle,
                        CapButt,JoinBevel);
  break;
  case (1):

    XSetDashes (mydisplay,mygc,0,dash1,2);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,mygc,mylinewidth,mylinestyle,
                        CapButt,JoinBevel);
  break;
  case (2):
    XSetDashes (mydisplay,mygc,0,dash2,2);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,mygc,mylinewidth,mylinestyle,
                        CapNotLast,JoinBevel);
  break;
  case (3):
    XSetDashes (mydisplay,mygc,0,dash3,2);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,mygc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
  break;
  case (4):
    XSetDashes (mydisplay,mygc,0,dash4,4);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,mygc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
  break;
  case (5):
    XSetDashes (mydisplay,mygc,0,dash5,4);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,mygc,mylinewidth,mylinestyle,
                        CapButt,JoinBevel);
  break;
  case (6):
    XSetDashes (mydisplay,mygc,0,dash6,6);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,mygc,mylinewidth,mylinestyle,
                        CapButt,JoinBevel);
  break;
}/* switch */
}
/*------------------------------------------------------*/
symbol_(x,y,size,iword,angl,nchar,dum)

/* Writes a Hollerwith string on the plot (plotter units) */ 

float            *x,*y,*size,*angl;
char             *iword;    
int              *nchar;
int              dum;

{
static           float xsiz;
float            xp, yp;
short            i, ixv, iyv;

xp = *x;
yp = *y;
if(p00000_.irot != 0){
  yp = *x;
  xp = -*y;
  if(p00000_.il34 == 0){
    if(*size > 0) xp = 27.2 - *y;}
  else if(p00000_.il34 == 1){
    if(*size > 0) xp = 40.1 - *y;}}
    
ixv = (int) (xp *  a00000_.psca);
iyv = (int) (yp *  a00000_.psca);

if (*size < 0.0) {
  ixv = a00000_.iox + ixv;
  iyv = a00000_.ioy - iyv;}
else {
  ixv = p00001_.xorig + ixv;
  iyv = a00000_.iyo - p00001_.yorig - iyv;}

if(ABS(*size) != xsiz){
  xsiz=ABS(*size);
  if(xsiz <=  0.25) XSetFont(mydisplay,mygc,font1);
  else if(xsiz <= 0.35) XSetFont(mydisplay,mygc,font2);
  else if(xsiz <= 0.45) XSetFont(mydisplay,mygc,font3);
  else if(xsiz <= 0.55) XSetFont(mydisplay,mygc,font4);
  else {
    XSetFont(mydisplay,mygc,font5);}}


XDrawString(mydisplay,mywindow,mygc,
            ixv,iyv,iword,*nchar);

}
/*------------------------------------------------------*/
symbu_(x,y,size,iword,angl,nchar,dum)

/* Writes a Hollerwith string on the plot (user units) */ 

float            *x,*y,*size,*angl;
char             *iword;      
int              *nchar;
int              dum;

{
/* move pen to symbol location */
static          float xsiz;
float           xp, yp;
short           i, ixv, iyv;

xp = *x;
yp = *y;
if(p00000_.irot != 0){
  yp = *x;
  xp = -*y;
  if(p00000_.il34 == 0){
    if(*size > 0) xp = 27.2 - *y;}
  else if(p00000_.il34 == 1){
    if(*size > 0) xp = 40.1 - *y;}}
    
xp = xp * p00000_.a;
yp = yp * p00000_.c;

if(*size > 0) {
  xp = xp + p00000_.b;
  yp = yp + p00000_.d;}

ixv = (int) (xp *  a00000_.psca);
iyv = (int) (yp *  a00000_.psca);

if (*size < 0.0) {
  ixv = a00000_.iox + ixv;
  iyv = a00000_.ioy - iyv;}
else {
  ixv = p00001_.xorig + ixv;
  iyv = a00000_.iyo - p00001_.yorig - iyv;}

if(ABS(*size) != xsiz){
  xsiz=ABS(*size);
  if(xsiz <=  0.25) XSetFont(mydisplay,mygc,font1);
  else if(xsiz <= 0.35) XSetFont(mydisplay,mygc,font2);
  else if(xsiz <= 0.45) XSetFont(mydisplay,mygc,font3);
  else if(xsiz <= 0.55) XSetFont(mydisplay,mygc,font4);
  else XSetFont(mydisplay,mygc,font5);}

XDrawString(mydisplay,mywindow,mygc,
            ixv,iyv,iword,*nchar);

}
/*-------------------------------------------------------------*/
number_(x,y,size,rn,angl,nsf,dum)

/* Writes a number on the plot: if nsf=klm, format is fkl.m (FORTRAN),
                                if nsf = -lm, rn is fixed to an integer 
                                   and format is ilm (FORTRAN) */
float *x,*y,*size,*rn,*angl;
int *nsf;                       
int dum;
            
{
int itot, idpl, ir;
int *it;                         
char iform[7];
char iword[30];
char *iptr;                             

if(*nsf > 0) {
  itot = *nsf / 10;
  it = &itot;
  idpl = *nsf - (*nsf /10) * 10;
  sprintf(iform,"%%%d.%df",itot,idpl);
  sprintf(iword,iform,*rn);}
else {
  ir = (int) *rn;
  it = nsf;
  sprintf(iword,"%d",ir);}

iptr = &iword[0];

symbol_(x,y,size,iptr,angl,it,dum);
}
 
/*---------------------------------------------------*/
csymbl_(x,y,ip,size,it)
 
float *x, *y, *size;
int *ip, *it;

{              
float radius, xp, yp;                               
float angl;
int il;
int dum;

static char isym[] = "ox*+#$@8hz";
static int icir[] = { 0,3,4,5,6,7,8,9,10,11,
                     100,103,104,105,106,107,108,109,110,111};

plotu_(x,y,ip);
angl = 0.0;
if(*it <= 10) {
  xp = *x - *size / 2.0;
  yp = *y - *size / 2.0;               
  il = 1;
  symbu_(&xp,&yp,size,&isym[*it-1],&angl,&il,dum);}
 else {
   radius = *size * 0.75;
   circle_(&radius, &icir[*it-11]);}
} 
/*------------------------------------------------------*/
pcirclef_(x,y,size)

float *x, *y, *size;

{
int *ip, *it;
*ip = 3;
*it = 21;
csymbl_(x,y,ip,size,it);
}
/*------------------------------------------------------*/
circle_(radius,nsides)

float *radius;
int *nsides;

{
XPoint *points;
float sta, ang, rpa;
float xnar[72], ynar[72];
int i, n, icx, icy, ixv, iyv, icol;                  
int shape, mode;
Bool fflag;

fflag = False;
n = *nsides;
if (*nsides >= 100) {
  fflag = True;
  n = *nsides-100;}
if(n == 0) n = 72;

if (fflag == True) {points = (XPoint *) malloc(n * sizeof(XPoint));}
icx = a00000_.iox;
icy = a00000_.ioy;
sta = 0.0;
ang = 6.283185308 / ((float) n);
rpa = *radius * a00000_.psca;
ixv = (int) rpa + icx;
iyv = icy;
a00000_.iox = ixv;
a00000_.ioy = iyv;
for (i = 0; i < n; i = i+1){
  sta = sta + ang;
  if (fflag == False){
    ixv = (int)(rpa * cos(sta)) + icx;
    iyv = (int)(rpa * sin(sta)) + icy;
    XDrawLine (mydisplay,mywindow,mygc,
             a00000_.iox,a00000_.ioy,ixv,iyv);
    a00000_.iox = ixv;
    a00000_.ioy = iyv;} 
  else if (fflag == True){
    ixv = (int)(rpa * cos(sta)) + icx;
    iyv = (int)(rpa * sin(sta)) + icy;
    points[i].x = ixv;
    points[i].y = iyv;}
}
 
if (fflag == True) {
  shape = Complex;
  mode = CoordModeOrigin;
  XFillPolygon(mydisplay,mywindow,mygc,points,n,shape,mode);
  free( (void *) points);}

a00000_.iox = icx;
a00000_.ioy = icy;
}
/*------------------------------------------------------------*/
CreateColourMap()

{
Bool contig;
unsigned int nplanes;
Status result;
int i;

contig=False;
nplanes=0;
                        
cmap=XDefaultColormap(mydisplay,myscreen);                           

for (i = 0;i < 32; i++)
  {
  color[i].red  =  red[i] * 256;
  color[i].green=  green[i] * 256;
  color[i].blue =  blue[i] * 256;
  color[i].pixel=  pixels[i];
  result=XAllocColor(mydisplay,cmap,&color[i]);
  }
} 
/*--------------------------------------------------*/
pen_(ipen,ithk)

int *ipen, *ithk;

{
unsigned long valuemask;
XGCValues values;

valuemask = GCForeground;
values.foreground = color[*ipen].pixel;
if(*ithk != 0) mylinewidth = *ithk;
if(*ithk == -1) mylinewidth = 0.5;
if(*ithk == -2) mylinewidth = 0.35;
if(*ithk == -3) mylinewidth = 0.25;
XSetLineAttributes (mydisplay,mygc,mylinewidth,mylinestyle,
                     CapButt,JoinBevel);

XChangeGC(mydisplay,mygc,valuemask,&values);
}
/*---------------------------------------------------*/
fillpoly_(xf,yf,nf)
 
float *xf, *yf;
int *nf;

{                  
XPoint *points;
int ixv, iyv;
int i, shape, mode;

points = (XPoint *) malloc(*nf * sizeof(XPoint));

for(i = 0; i < *nf; i++){
  ixv = (int)(a00000_.psca * *(xf+i)) + p00001_.xorig;
  iyv = a00000_.iyo - (int)(a00000_.psca * *(yf+i)) - p00001_.yorig;
  points[i].x = ixv;
  points[i].y = iyv;}

shape = Complex;
mode = CoordModeOrigin;

XFillPolygon(mydisplay,mywindow,mygc,points,*nf,shape,mode);

free( (void *) points);
}
/*--------------------------------------------------*/
origin_(x,y,iorig)
      
float *x, *y;
int *iorig;

{
if(*iorig == 0){
  p00001_.xorig = a00000_.psca * (int) *x;
  p00001_.yorig = a00000_.psca * (int) *y;}
else if(*iorig > 0){
  p00001_.xorig = p00001_.xorig + a00000_.psca * (int) *x;
  p00001_.yorig = p00001_.yorig + a00000_.psca * (int) *y;}
else if(*iorig < 0){  
  if(a00000_.psca == 0){ 
    printf("ORIGIN error: zero scale\n");}
  else{
    *x = p00001_.xorig / a00000_.psca;
    *y = p00001_.yorig / a00000_.psca;}}

}
/*--------------------------------------------------*/
where_(x,y,rfact)
      
float *x, *y, *rfact;

/* Check this */

{
float r;
     
r = 1.0;
rfact = &r;
*x = (float) a00000_.iox * a00000_.psca;
*y = (float) a00000_.ioy * a00000_.psca;
factor_(rfact);
}
/*--------------------------------------------------*/
factor_(fact)

float *fact;

{    
static int first = 0;
static float factl = 1.0;
float osca;

if(first == 1){
  osca = a00000_.psca;
  first = 1;}
if(*fact > 0){
  a00000_.psca = osca * *fact;
  factl = *fact;}
else{
  *fact = factl;}

} 
/*---------------------------------------------------*/
shadrt_(xi,yi)
 
float *xi, *yi;

{             
unsigned int dx,dy;

dx =  (int) a00000_.psca * *xi;   
dy =  (int) a00000_.psca * *yi;   

XFillRectangle(mydisplay,mywindow,mygc,a00000_.iox,
               a00000_.ioy-dy,dx,dy);

}
/*---------------------------------------------------*/
edgert_(xi,yi)
 
float *xi, *yi;

{             
unsigned int dx,dy;

dx =  (int) a00000_.psca * *xi;   
dy =  (int) a00000_.psca * *yi;   

XDrawRectangle(mydisplay,mywindow,mygc,a00000_.iox,
               a00000_.ioy-dy,dx,dy);
}

  
/*---------------------------------------------------*/
dispnm_(dspnm,nchr,dum)
char *dspnm;
int *nchr,dum;

/* Routine to find the display name on which client is working*/

{
int m,i,j,k;
struct utmp ut;

m = sizeof ut;
k = ttyslot();				/* Get offset of my user's entry in utmp */
j = open("/etc/utmp","rb",0);		/* Open utmp and */
lseek(j,k*m,L_SET); 			/*    seek to the appropriate entry. */
read(j,&ut,sizeof ut);			/* Read the entry. */  
close(j);

if (strlen(ut.ut_line) == 0)		/* If no remote host, use the local one */
      gethostname(dspnm,*nchr);
else
      strcpy(dspnm,ut.ut_line);

strcat(dspnm,":0");			/* Append the :0 - probably not needed */

return;
}

/*---------------------------------------------------*/
/*
      	This version of ldcolr replaces old version which produced an 
        error with X11 server but not the xnews server. 
							MS 9/4/96  
*/
ldcolrnew_(lunit)

int *lunit;
       
{
Bool contig;
unsigned int nplanes;
/*extern cfopen(int &, int &, int *, int *, int *, int *);*/
extern cfopen();
Status result;
int lu;
int id[PIXELS],ir[PIXELS],ig[PIXELS],ib[PIXELS];
int i;
unsigned int j;

contig=False;
nplanes=0;
            
lu = *lunit;       

/* printf ("Into ldcolr \n");  */

cmap=XDefaultColormap(mydisplay,myscreen);                           

      
cfopen_(&lu, &ncol, id, ir, ig, ib);
/* printf ("%4d %4d %4d %4d %4d\n",ncol,*(id),*(ir),*(ig),*(ib)); */

result=XAllocColorCells (mydisplay,cmap,contig,
                         plane_masks,nplanes,pixels,ncol+32);
if (result==0) printf ("1 Error in allocating the color cells...\n");
/* printf("ncol is, %3d \n",ncol); */

for (i = 0;i < 32; i++)
  {
  color[i].red  =  red[i] * 256;
  color[i].green=  green[i] * 256;
  color[i].blue =  blue[i] * 256;
  color[i].flags=  DoRed | DoGreen | DoBlue;
  color[i].pixel=  pixels[i];

  XStoreColor(mydisplay,cmap,&color[i]);
  }
for (i = 0;i < ncol; i++)
  {     
/*  printf("%3d %3d %3d %3d\n",id[i],ir[i],ig[i],ib[i]);*/
  color[id[i]].red  =  ir[i] * 256;
  color[id[i]].green=  ig[i] * 256;
  color[id[i]].blue =  ib[i] * 256;
  color[id[i]].flags=  DoRed | DoGreen | DoBlue;
  color[id[i]].pixel=  pixels[i+32];
  XStoreColor (mydisplay,cmap,&color[id[i]]);
  }




/* XStoreColors (mydisplay,cmap,color,ncolors); */

/*printf (" Xstorecolours \n");*/                        
}
/*---------------------------------------------------*/
/*
      	This version of ldcolr produced an error with X11 server but not 
	the xnews server. A new version has been written for the X11 server
							MS 9/4/96  
*/
ldcolr_(lunit)

int *lunit;
       
{
Bool contig;
unsigned int nplanes;
/*extern cfopen(int &, int &, int *, int *, int *, int *);*/
extern cfopen();
Status result;
int lu;
int id[PIXELS],ir[PIXELS],ig[PIXELS],ib[PIXELS];
int i;
unsigned int j;
unsigned int displaydepth;

displaydepth=XDefaultDepth(mydisplay,myscreen);

contig=False;
nplanes=0;
            
lu = *lunit;       

/* printf ("Into ldcolr \n");  */

cmap=XDefaultColormap(mydisplay,myscreen);                           

/*printf ("Out of defcm \n");*/                       
      
cfopen_(&lu, &ncol, id, ir, ig, ib);
/* printf ("%4d %4d %4d %4d %4d\n",ncol,*(id),*(ir),*(ig),*(ib)); */

/*   printf(" start XAllocColorCells\n"); */
if (displaydepth<9) 
{result=XAllocColorCells (mydisplay,cmap,contig,
                         plane_masks,nplanes,pixels,ncol+32);
if (result==0) printf ("1 Error in allocating the color cells...\n");}

for (i = 0;i < 32; i++)
  {
  color[i].red  =  red[i] * 256;
  color[i].green=  green[i] * 256;
  color[i].blue =  blue[i] * 256;
  color[i].flags=  DoRed | DoGreen | DoBlue;
  color[i].pixel=  pixels[i];

  if (displaydepth<9) {XStoreColor(mydisplay,cmap,&color[i]);}
  else {color[i].pixel=0; XAllocColor (mydisplay,cmap, &color[i]);}
  }
for (i = 0;i < ncol; i++)
  {     
/*  printf("%3d %3d %3d %3d\n",id[i],ir[i],ig[i],ib[i]);*/
  color[id[i]].red  =  ir[i] * 256;
  color[id[i]].green=  ig[i] * 256;
  color[id[i]].blue =  ib[i] * 256;
  color[id[i]].flags=  DoRed | DoGreen | DoBlue;
  color[id[i]].pixel=  pixels[i+32];
  if (displaydepth<9) {XStoreColor(mydisplay,cmap,&color[id[i]]);}
  else {color[id[i]].pixel=0; XAllocColor (mydisplay,cmap, &color[id[i]]);}
  }

/* XStoreColors (mydisplay,cmap,color,ncolors); */

/*printf (" Xstorecolours \n");*/                        
}
/*---------------------------------------------------*/
pimag4_(xori,yori,xxl,yyl,nsxx,nsyy,arr,nth11,nth22,sth11,sth22)

/* POSTPAK pimag4 variables nth11,nth22 (ranging between 1 and 16) are
mapped into a 1,256 range and incorporated into a call to pimag8 (see
below)
*/

float *xori,*yori,*xxl,*yyl,*sth11,*sth22;
float *arr;
int *nsxx,*nsyy,*nth11,*nth22;

{
int n1,n2;
 
n1 = (*nth11 - 1) * 8 + 4;
n2 = (*nth22 - 1) * 8 + 4;
pimag8_(xori,yori,xxl,yyl,nsxx,nsyy,arr,&n1,&n2,sth11,sth22);
}
/*---------------------------------------------------*/
picol_(xori,yori,xxl,yyl,nsxx,nsyy,arr,nth11,nth22,sth11,sth22,imap)

/* Routine to simulate POSTPAK picol. Variables are defined
as follows:
      Passed variables:
        xori,yori - coordinates of bottom left corner of image.
        xxl,yyl - length of x and y sides (plotter units as above). 
        nsxx,nsyy - no of individual colour cells along x and y borders.
        arr - the array containing the values to be contoured. 
        nth11,nth22 - range of indices in POSTPAK grey shades (1:256).
        sth11,sth22 - corresponding array maximum and minimum values.
        imap - switch to colour outside bounds.

      Work variables:
        npx,npy - no of screen pixels along x and y borders.
        ixo,iyo - identification of plot origin pixels.
        ix,iy - identification of current pixel being considered. (0:npx-1)
        nx,ny - identification of cell in which current pixel belongs. (1:nsxx) 
     
   N.B. The maximum number of permitted color indices in a colour map 
        is PIXELS (i.e. colour indices in the range 0 --> PIXELS-1). If 
	either nth11 or nth22 is greater than PIXELS-1 it is set equal
	to PIXELS-1.

	If either the pen range or the parameter range are swapped around
	then the colour mapping is reversed. All array values which fall
	beyond the parameter range (s1,s2) are mapped onto the end colours
	if imap = 0 and are not plotted if imap = 1.
*/
float *xori,*yori,*xxl,*yyl,*sth11,*sth22;
float arr[];
int *nsxx,*nsyy,*nth11,*nth22,*imap;

{ 
int           npx,npy;
int           ix,iy,nx,ny;
unsigned long valuemask;
XGCValues     values;
XPoint        points[MAXDOT];
int           count[PIXELS];
int           pix,ixo,iyo,i,k,ndiff,id;          
float         a1,sdiff;
unsigned int iwd,nth1,nth2,n1,n2;
unsigned int iht; 
 
sdiff = *sth22 - *sth11;
nth1 = *nth11;
nth2 = *nth22;
n1 = i2min(nth1,nth2);
n2 = i2max(nth1,nth2);
ndiff = nth2 - nth1;
npx = (int) (*xxl *  a00000_.psca);
npy = (int) (*yyl *  a00000_.psca);
ixo = (int) (*xori * a00000_.psca) + p00001_.xorig;
iyo = a00000_.iyo - p00001_.yorig - (int) (*yori * a00000_.psca);


for (i = 0 ; i < PIXELS; i++) count[i]=0;

for (iy = 0; iy < npy; iy++){
  ny = 1 + (int) ((float) (iy) * (float) (*nsyy) / (float) (npy)); 
  for (ix = 0; ix < npx; ix++ ){
    nx = 1 + (int) ((float) (ix) * (float) (*nsxx) / (float) (npx)); 
    id = nx + (ny - 1) * (*nsxx) - 1;
    a1 = *(arr+id);
    pix = (int) ((a1 - *sth11) * (float) (ndiff) / sdiff) + nth1;
    if (*imap > 0){
      if (pix <= n2){
      if (pix >= n1){
      xdot[count[pix]][pix] = ix + ixo;
      ydot[count[pix]][pix] = iyo - iy;
      count[pix]++;
        
/* Colour in groups of individual pixels when the number of pixels with
a given pixel value pix exceeds MAXDOT. Note that this is efficient when
dealing with large arrays (arr).  */

      if (count[pix] == MAXDOT){
        for(k = 0; k <= MAXDOT; k++){
          points[k].x=xdot[k][pix];
          points[k].y=ydot[k][pix];
        }
        valuemask = GCForeground;
        values.foreground = color[pix].pixel;
        XChangeGC(mydisplay,mygc,valuemask,&values);
        XDrawPoints (mydisplay,mywindow,mygc,points,MAXDOT,CoordModeOrigin);
        count[pix]=0;
      }
    }}}
    else {
      if (pix > n2) pix = n2;
      if (pix < n1) pix = n1;
      xdot[count[pix]][pix] = ix + ixo;
      ydot[count[pix]][pix] = iyo - iy;
      count[pix]++;

/* Colour in groups of individual pixels when the number of pixels with
a given pixel value pix exceeds MAXDOT. Note that this is efficient when
dealing with large arrays (arr).  */

      if (count[pix] == MAXDOT){
        for(k = 0; k <= MAXDOT; k++){
          points[k].x=xdot[k][pix];
          points[k].y=ydot[k][pix];
        }
        valuemask = GCForeground;
        values.foreground = color[pix].pixel;
        XChangeGC(mydisplay,mygc,valuemask,&values);
        XDrawPoints (mydisplay,mywindow,mygc,points,MAXDOT,CoordModeOrigin);
        count[pix]=0;
      }
    }
  }
}

/* Colour in remaining values of pix where count[pix] does not exceed MAXDOT.*/
 
for(pix = n1; pix <= n2; pix++){
  if(count[pix] != 0){
    for(k=0; k <= count[pix]; k++){
      points[k].x=xdot[k][pix];
      points[k].y=ydot[k][pix];
    }
     valuemask = GCForeground;
     values.foreground = color[pix].pixel;
     XChangeGC(mydisplay,mygc,valuemask,&values);
     XDrawPoints(mydisplay,mywindow,mygc,points,count[pix],CoordModeOrigin);
  }
}
}
/*---------------------------------------------------*/
pimag8_(xori,yori,xxl,yyl,nsxx,nsyy,arr,nth11,nth22,sth11,sth22)

/* 
Routine to simulate Postscript library routine pimag8. 
This is an obsolete routine and has been superseeded by picol which 
does the same thing in colour. 
Picol uses the current colour table for the pixels (a greyscale can
be achieved by reading in a greyscale using the routine ldcolr,
e.g. the one in /usr/server/plot/pak/Examples.

*/
float *xori,*yori,*xxl,*yyl,*sth11,*sth22;
float arr[];
int *nsxx,*nsyy,*nth11,*nth22;

{ 
/*    write out error message */ 

printf("Error routine pimag8 is obsolete in the xpak library \n");
printf("use subroutine picol with a greyscale colour table \n");
printf("read in with routine ldcolr\n");
}

/*---------------------------------------------------*/
unsigned int i2max(n1,n2)
/* program to calculate the greater of two integers*/

int n1,n2;
{
int	n;
	n = n1;
	if(n2 > n1) n = n2;
	return (n);
}
/*---------------------------------------------------*/
unsigned int i2min(n1,n2)
/* program to calculate the lesser of two integers*/

int n1,n2;
{
int	n;
	n = n1;
	if(n2 < n1) n = n2;
	return (n);
}
/* -------------------------------------------------- */
xname_(fnme,nlen,nfnme)

/* to give a name to the window and icon.
   note: to have any effect it must be called before hplots*/

char     *fnme;
int      *nfnme;
int      *nlen;

{

int     i;

for ( i=0 ; i< *nlen ; i++) mywindowname[i]=fnme[i];
mywindowname[ *nlen]='\0';

}

/* -------------------------------------------------- */
plottype_(n)

/* Tells user if the X-pak or postscript library has been compiled */

int      *n;

{
*n=1;
}

