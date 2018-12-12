/* debug subroutines */
#include	<stdio.h>

#include <stdio.h>
#define INT int
#include	<stdlib.h>
int dc_hasmouse = 0;
#ifndef INT
#ifdef MSDOS
#define INT long
#else
#define INT int
#endif
#endif
/* Global Variables		*/
INT     NumX    =0 ;
INT     NumY    =0 ;

void onintr()
{
	exit( 1 );
}

void di_gcmdln(argc,argv)
int argc;
char **argv;
{
}

#define W(x)	fprintf( stdout, "%s", x )
#define WN(x)	fprintf( stdout, "%s\n", x )

void di_openpl(void)
{
	WN( "Open plot." );
}

void di_closepl(int mode)
{
	fprintf(stdout, "Close plot Mode %d\n",mode);
}


void di_move( INT x, INT y)
{
	fprintf( stdout, "Move to (%d,%d).\n", x, y );
}

void di_linec( INT X0, INT Y0, INT X1, INT Y1)
{
	fprintf( stdout, "Line from (%d,%d) to (%d,%d).\n",
		X0, Y0, X1, Y1 );
}

void di_erase( INT mode )
{
	fprintf(stdout,"Erase(%d).\n",mode );
}


void di_point( INT x, INT y)
{
	fprintf( stdout, "Point at (%d,%d).\n", x, y );
}

void di_cont( INT x, INT y)
{
	fprintf( stdout, "Continue line to (%d,%d).\n", x, y );
}

void di_control( int type, int i1, int i2, int i3, int i4)
{
	fprintf( stdout, "Control type %d args %d,%d,%d,%d.\n", type, i1,i2,i3,i4 );
}


void di_space( INT X0, INT Y0, INT X1, INT Y1)
{
	fprintf( stdout, "Space  (%d,%d) to (%d,%d)\n", X0, Y0, X1, Y1 );
}

void di_arc( INT Xi, INT Yi, INT X0, INT Y0, INT X1, INT Y1 )
{
	fprintf( stdout, "Arc  (%d,%d)  (%d,%d)  (%d,%d)\n",
		Xi, Yi, X0, Y0, X1, Y1 );
}

void di_circle( INT x, INT y, INT r)
{
	fprintf( stdout, "Circle at (%d,%d) with radius %d.\n", x, y, r );
}

void di_linemod( char *s )
{
	fprintf( stdout, "Line type: \"%s\"\n", s );
}

void di_pen( INT n )
{
	fprintf( stdout, "Pen %d.\n", n );
}

void di_fillt( INT X0,INT Y0,INT X1,INT Y1,INT X2,INT y2,INT patx,INT paty,INT lenx,INT leny)
{
	fprintf(stdout,"Fill Triangle(%d,%d),(%d,%d),(%d,%d) patx=%d paty=%d lenx=%d leny=%d\n",
		X0,Y0,X1,Y1,X2,y2,patx,paty,lenx,leny);
}

void di_fillr( INT X0,INT Y0,INT X1,INT Y1,INT patx,INT paty,INT lenx,INT leny)
{
	fprintf(stdout, "Fill Rectangle (%d,%d) to (%d,%d) patx=%d paty=%d lenx=%d leny=%d\n", 
		X0,Y0,X1,Y1,patx,paty,lenx,leny);
}

void di_fills( INT X0,INT Y0,INT ixy,INT istnd,INT iplmn)
{
	fprintf(stdout, "Variable Area Trace (%d,%d) ixy=%d istnd=%d iplmn=%d\n"		,X0,Y0,ixy,istnd,iplmn);
}

void di_fillp( INT n , INT *x, INT *y)
{
	int i;
	fprintf(stdout, "Polygon Fill n=%d\n",n);
	for(i=0; i<n; i++)
		fprintf(stdout, "(%d, %d)\n",x[i],y[i]);

	
}



void di_label( char *s )
{
	fprintf( stdout, "Label: \"%s\".\n", s );
}

void di_gottxt( char *s )
{
	fprintf( stdout, "Gottxt: \"%s\".\n", s );
}

void di_font(INT fontvalue)
{
	fprintf( stdout, "Font: %d\n",fontvalue);
}

void di_gwid(INT width)
{
	fprintf( stdout, "Width: %d\n",width);
}

void di_cursor(INT curcmd)
{
	if(curcmd == 1)
		fprintf(stdout,"Gcursor: Arrow\n");
	else if(curcmd == 2)
		fprintf(stdout,"Gcursor: Crosshairs\n");
	else if(curcmd == 3)
		fprintf(stdout,"Gcursor: Plus\n");
}

void di_gsymb(INT x,INT y,INT ht,INT ang,INT nchar,char *str)
{
	int sval;
	fprintf( stdout, "Gsymb: at (%d,%d) ht=%d, angle=%d nochar=%d\n"
			,x,y,ht,ang,nchar);
	if(nchar>0)
		fprintf(stdout, "          string=%s\n",str);
	else {
		sval = (int)str[0];
		fprintf(stdout, "          symbol=%d\n",sval);
	}
}

void di_clip(INT reset, INT X0, INT Y0, INT X1, INT Y1 )
{
	if(reset == 1)
	fprintf( stdout, "Set Clip  (%d,%d) to (%d,%d)\n", X0, Y0, X1, Y1 );
	else
		fprintf( stdout, "Reset Clip  \n" );
}


/* get text string from terminal */
void di_gintxt(int cnt,char *s)
{
	s[0]='\0';
}
void di_cross(INT *ix, INT *iy, char *c)
{
	*ix = 0;
	*iy = 0;
	c[0] = '\0';
}

void di_info(INT *HasMouse, INT *XminDev, INT *YminDev, 
	INT *XmaxDev, INT *YmaxDev, INT *XminClip, 
	INT *YminClip, INT *XmaxClip, INT *YmaxClip, INT *Color)
{
	*HasMouse = 0;
	*XminDev  = 1;
	*YminDev  = 2;
	*XmaxDev  = 3;
	*YmaxDev  = 4;
	*XminClip = 5;
	*YminClip = 6;
	*XmaxClip = 7;
	*YmaxClip = 8;
	*Color    = 2;
}
void dv_mesg(char *mesg)
{
	fprintf(stdout,"gmesg: %s\n",mesg);
}
