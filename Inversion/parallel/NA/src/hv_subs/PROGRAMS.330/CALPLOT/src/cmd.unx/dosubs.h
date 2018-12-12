#ifndef INT
#define INT int
#endif
extern void (*do_arc)(INT Xi,INT Yi,INT X0,INT Y0,INT X1,INT Y1);
extern void (*do_circle)(INT Xi,INT Yi,INT r);
extern void (*do_clip)(INT cmd, INT X0,INT Y0,INT X1,INT Y1);
extern void (*do_cont)(INT X0,INT Y0);
extern void (*do_control)(int type, int i1, int i2, int i3, int i4);
extern void (*do_cross)(int *X0,int *Y0, char *c);
extern void (*do_cursor)(INT curstyp);
extern void (*do_erase)(INT mode);
extern void (*do_fillp)(INT n,INT *x,INT *y);
extern void (*do_fillr)(INT X0,INT Y0,INT X1,INT Y1,
	INT patx,INT paty,INT lenx,INT leny);
extern void (*do_fills)(INT X0,INT Y0,INT ixy,INT istnd,INT iplmn);
extern void (*do_fillt)(INT X0,INT Y0,INT X1,INT Y1,INT X2,INT Y2,
	INT patx,INT paty,INT lenx,INT leny);
extern void (*do_font)(INT Xi);
extern void (*do_gintxt)(int cnt, char *s);
extern void (*do_gottxt)(char *s);
extern void (*do_gsymb)(INT X0,INT Y0,INT X1,INT Y1,INT n,char *s);
extern void (*do_gwid)(INT wid);
extern void (*do_info)(INT *HasMouse, INT *XminDev, INT *YminDev, 
	INT *XmaxDev, INT *YmaxDev, INT *XminClip, 
	INT *YminClip, INT *XmaxClip, INT *YmaxClip, INT *Color);
extern void (*do_label)(char *s);
extern void (*do_linec)(INT X0,INT Y0,INT X1,INT Y1);
extern void (*do_linemod)(char *s);
extern void (*do_move)(INT X0,INT Y0);
extern void (*do_pen)(INT Xi);
extern void (*do_point)(INT X0,INT Y0);
extern void (*do_space)(INT X0,INT Y0,INT X1,INT Y1);
