#ifndef INT
#define INT int
#endif

extern void dfarcc_(INT *Xi,INT *Yi,INT *X0,INT *Y0,INT *X1,INT *Y1);
extern void dfcirc_(INT *Xi,INT *Yi,INT *r);
extern void dfclip_(INT *cmd, INT *X0,INT *Y0,INT *X1,INT *Y1);
extern void dfcont_(INT *X0,INT *Y0);
extern void dfcontrol_(int *type, int *i1, int *i2, int *i3, int *i4);
extern void dfcros_(int *X0,int *Y0, char *c, long lc);
extern void dfcurs_(INT *curstyp);
extern void dferas_(void);
extern void dffilp_(INT *n,float *x,float *y, float * xold, float *yold,
	float *xcur, float* ycur, float *xstp, float *ystp);
extern void dffilr_(INT *X0,INT *Y0,INT *X1,INT *Y1,
	INT *patx,INT *paty,INT *lenx,INT *leny);
extern void dffils_(INT *X0,INT *Y0,INT *ixy,INT *istnd,INT *iplmn);
extern void dffilt_(INT *X0,INT *Y0,INT *X1,INT *Y1,INT *X2,INT *Y2,
	INT *patx,INT *paty,INT *lenx,INT *leny);
extern void dffont_(INT *Xi);
extern void dfgint_(int *cnt, char *s, long lstr);
extern void dfgott_(int *cnt, char *s, long lstr);
extern void dfgsym_(INT *X0,INT *Y0,INT *X1,INT *Y1,INT n*,char *s,long lstr);
extern void dfgwid_(INT *wid);
extern void dfinfo_(INT *HasMouse, INT *XminDev, INT *YminDev, 
	INT *XmaxDev, INT *YmaxDev, INT *XminClip, 
	INT *YminClip, INT *XmaxClip, INT *YmaxClip);
extern void dfmove_(INT *X0,INT *Y0);
extern void dfpenn_(INT *Xi);
extern void dfpont_(INT *X0,INT *Y0);
extern void dfspce_(INT *X0,INT *Y0,INT *X1,INT *Y1);
extern void dfopen_(int ls, char *fname, long lstr)
extern void dfclos_(void);
