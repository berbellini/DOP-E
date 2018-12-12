#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "dosubs.h"

#define LSTR 100
static char str1[LSTR];
static char str2[LSTR];

void clospl(int mode );
void openpl(int ls,char *s, int lc, char *c);



#ifdef MSDOS
void dfclip(INT *cmd, INT *X0,INT *Y0,INT *X1,INT *Y1)
#else
void dfclip_(INT *cmd, INT *X0,INT *Y0,INT *X1,INT *Y1)
#endif
{
	(*do_clip)( *cmd,  *X0, *Y0, *X1, *Y1);
}

#ifdef MSDOS
void dfcont(INT *X0,INT *Y0)
#else
void dfcont_(INT *X0,INT *Y0)
#endif
{
	(*do_cont)( *X0, *Y0);
}

#ifdef MSDOS
void dfcontrol(int *type, int *i1, int *i2, int *i3, int *i4)
#else
void dfcontrol_(int *type, int *i1, int *i2, int *i3, int *i4)
#endif
{
	(*do_control)( *type, *i1, *i2, *i3, *i4);
}

#ifdef MSDOS
void dfcros(int *X0,int *Y0, char *c)
#else
void dfcros_(int *X0,int *Y0, char *c, long lc)
#endif
{
	(*do_cross)( X0, Y0, c);
}

#ifdef MSDOS
void dfcurs(INT *curstyp)
#else
void dfcurs_(INT *curstyp)
#endif
{
	(*do_cursor)( *curstyp);
}


#ifdef MSDOS
void dferas(INT *mode)
#else
void dferas_(INT *mode)
#endif
{
	INT pode;
	pode = *mode;
	(*do_erase)(pode);
}

#ifdef MSDOS
void dffilr(INT *X0,INT *Y0,INT *X1,INT *Y1,
	INT *patx,INT *paty,INT *lenx,INT *leny)
#else
void dffilr_(INT *X0,INT *Y0,INT *X1,INT *Y1,
	INT *patx,INT *paty,INT *lenx,INT *leny)
#endif
{
	(*do_fillr)( *X0, *Y0, *X1, *Y1,
		*patx, *paty, *lenx, *leny);
}

#ifdef MSDOS
void dffils(INT *X0,INT *Y0,INT *ixy,INT *istnd,INT *iplmn)
#else
void dffils_(INT *X0,INT *Y0,INT *ixy,INT *istnd,INT *iplmn)
#endif
{
	(*do_fills)( *X0, *Y0, *ixy, *istnd, *iplmn);
}

#ifdef MSDOS
void dffilt(INT *X0,INT *Y0,INT *X1,INT *Y1,INT *X2,INT *Y2,
	INT *patx,INT *paty,INT *lenx,INT *leny)
#else
void dffilt_(INT *X0,INT *Y0,INT *X1,INT *Y1,INT *X2,INT *Y2,
	INT *patx,INT *paty,INT *lenx,INT *leny)
#endif
{
	(*do_fillt)( *X0, *Y0, *X1, *Y1, *X2, *Y2,
		*patx, *paty, *lenx, *leny);
}

#ifdef MSDOS
void dffont(INT *Xi)
#else
void dffont_(INT *Xi)
#endif
{
	(*do_font)( *Xi);
}


#ifdef MSDOS
void dfgsym(INT *X0,INT *Y0,INT *X1,INT *Y1,INT *n,char *s)
#else
void dfgsym_(INT *X0,INT *Y0,INT *X1,INT *Y1,INT *n,char *s, long lstr)
#endif
{
	(*do_gsymb)( *X0, *Y0, *X1, *Y1, *n, s);
}

#ifdef MSDOS
void dfgwid(INT *wid)
#else
void dfgwid_(INT *wid)
#endif
{
	(*do_gwid)( *wid);
}

#ifdef MSDOS
void dfinfo(INT *HasMouse,INT  *XminDev,INT  *YminDev, 
	 INT *XmaxDev,INT  *YmaxDev,INT  *XminClip, 
	 INT *YminClip,INT  *XmaxClip,INT  *YmaxClip, INT *Color)
#else
void dfinfo_(INT *HasMouse,INT  *XminDev,INT  *YminDev, 
	 INT *XmaxDev,INT  *YmaxDev,INT  *XminClip, 
	 INT *YminClip,INT  *XmaxClip,INT  *YmaxClip, INT *Color)
#endif
{
	(*do_info)(HasMouse,  XminDev,  YminDev, 
	 	XmaxDev,  YmaxDev,  XminClip, 
	 	YminClip,  XmaxClip,  YmaxClip, Color);
}


#ifdef MSDOS
void dfmove(INT *X0,INT *Y0)
#else
void dfmove_(INT *X0,INT *Y0)
#endif
{
	(*do_move)( *X0, *Y0);
}

#ifdef MSDOS
void dfpenn(INT *Xi)
#else
void dfpenn_(INT *Xi)
#endif
{
	(*do_pen)( *Xi);
}

#ifdef MSDOS
void dfpont(INT *X0,INT *Y0)
#else
void dfpont_(INT *X0,INT *Y0)
#endif
{
	(*do_point)( *X0, *Y0);
}

#ifdef MSDOS
void dfspce(INT *X0,INT *Y0,INT *X1,INT *Y1)
#else
void dfspce_(INT *X0,INT *Y0,INT *X1,INT *Y1)
#endif
{
	(*do_space)( *X0, *Y0, *X1, *Y1);
}

void gomesg(int cnt, char *mesg);
#ifdef MSDOS
void dfmesg(INT *lmesg,char *mesg)
#else
void dfmesg_(INT *lmesg,char *mesg, long lstr)
#endif
{
	int lst;
	lst = *lmesg;
	gomesg(*lmesg, mesg);
}


#ifdef MSDOS
void dffilp(INT *narr,float *xarr,float *yarr, float * xold, float *yold,
	float *xcur, float* ycur, float *xstp, float *ystp)
#else
void dffilp_(INT *narr,float *xarr,float *yarr, float * xold, float *yold,
	float *xcur, float* ycur, float *xstp, float *ystp)
#endif
{
	float xx, yy;
	INT i;
	INT ilw, jlw;
	INT *x, *y;
	INT n;
	n = *narr;
	if((x = (INT *)calloc(n, sizeof(INT))) == NULL)
		return;
	if((y = (INT *)calloc(n, sizeof(INT))) == NULL)
		return;
	n = *narr;
	for(i = 0 ; i < n ; i++){
       		xx = 1000. * (xarr[i]* (*xstp) + *xold);
       		yy = 1000. * (yarr[i]* (*ystp) + *yold);
		if(xx > 1000000000.0)
			ilw = 1000000000;
		else if(xx < -1000000000.0)
			ilw = -1000000000;
		else
			ilw = xx ;
		if(yy > 1000000000.0)
			jlw = 1000000000;
		else if(yy < -1000000000.0)
			jlw = -1000000000;
		else
			jlw = yy ;
		x[i] = ilw;
		y[i] = jlw;
	}
	(*do_fillp)(n,x,y);
	free(x);
	free(y);
}

#ifdef MSDOS
void dfopen(char *lfname, char *lcname, INT *lstr, INT *lcon)
{
	int lf, lc;
	lf = *lstr;
	lc = *lcon;
	openpl(lf,lfname,lc,lcname);
}
#else
void dfopen_(char *lfname, char *lcname, INT *lstr, INT *lcon,
		long lstr1, long lstr2)
{
	int i;
	int lf, lc;
	lf = *lstr;
	lc = *lcon;
	for(i=0;i<lf;i++)
		str1[i] = lfname[i];
	str1[lf]='\0';
	for(i=0;i<lc;i++)
		str2[i] = lcname[i];
	str2[lc]='\0';
	
	openpl(lf,lfname,lc,lcname);
}
#endif

#ifdef MSDOS
void dfclos(int mode)
#else
void dfclos_(int mode)
#endif
{
	clospl(mode);
}

#ifdef MSDOS
void dfgott(int *cnt, char *text)
#else
void dfgott_(int *cnt, char *text, long lstr)
#endif
{
	int ls;
	ls = *cnt;
	gottxt(ls,text);
}


#ifdef MSDOS
void dfgint(long *cnt, char *s)
#else
void dfgint_(long *cnt, char *s, long lstr)
#endif
{
	int ls;
	ls = *cnt;
	gintxt(ls,s);
	ls = strlen(s);
#ifdef MSDOS
	if(ls < *cnt)
		s[ls++] = ' ';
#else
	s[ls] = ' ';
#endif
}


void dv_gread(int lf,char *fname, INT NumX, INT NumY, INT LowX, INT LowY,
        INT HighX, INT HighY, INT Num, INT Sclx, INT Scly);
#ifdef MSDOS
void dfread(INT *lfname,char *fname, INT *NumX, INT *NumY, INT *LowX,
	INT* LowY, INT * HighX, INT *HighY, INT *Num, INT *Sclx, INT *Scly)
#else
void dfread_(INT *lfname,char *fname, INT *NumX, INT *NumY, INT *LowX,
	INT* LowY, INT * HighX, INT *HighY, INT *Num, INT *Sclx, INT *Scly, long lstr)
#endif
{
	int lf;
	lf = *lfname;
	dv_gread(lf, fname, *NumX, *NumY, *LowX, *LowY,
        	*HighX, *HighY, *Num, *Sclx, *Scly);
}
