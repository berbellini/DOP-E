/* Prototype for Device Dependent Functions */
#ifndef INT
#define INT int
#endif
extern void di_arc(INT Xi,INT Yi,INT X0,INT Y0,INT X1,INT Y1);
extern void di_circle(INT Xi,INT Yi,INT r);
extern void di_clip(INT cmd, INT X0,INT Y0,INT X1,INT Y1);
extern void di_closepl(int mode);
extern void di_cont(INT X0,INT Y0);
extern void di_control(int type, int i1, int i2, int i3, int i4);
extern void di_cross(int *X0,int *Y0, char *c);
extern void di_cursor(INT curstyp);
extern void di_erase(INT mode);
extern void di_fillp(INT n,INT *x,INT *y);
extern void di_fillr(INT X0,INT Y0,INT X1,INT Y1,
	INT patx,INT paty,INT lenx,INT leny);
extern void di_fills(INT X0,INT Y0,INT ixy,INT istnd,INT iplmn);
extern void di_fillt(INT X0,INT Y0,INT X1,INT Y1,INT X2,INT Y2,
	INT patx,INT paty,INT lenx,INT leny);
extern void di_font(INT Xi);
extern void di_gintxt(int cnt, char *s);
extern void di_gottxt(char *s);
extern void di_gsymb(INT X0,INT Y0,INT X1,INT Y1,INT n,char *s);
extern void di_gwid(INT wid);
extern void di_info(INT *HasMouse, INT *XminDev, INT *YminDev, 
	INT *XmaxDev, INT *YmaxDev, INT *XminClip, 
	INT *YminClip, INT *XmaxClip, INT *YmaxClip, INT *Color);
extern void di_label(char *s);
extern void di_linec(INT X0,INT Y0,INT X1,INT Y1);
extern void di_linemod(char *s);
extern void di_move(INT X0,INT Y0);
extern void di_openpl(int showmenu);
extern void di_pen(INT Xi);
extern void di_point(INT X0,INT Y0);
extern void di_space(INT X0,INT Y0,INT X1,INT Y1);
