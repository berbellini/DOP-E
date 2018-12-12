#ifndef INT
#define INT int
#endif
void putpi(INT a);
void dp_arc(INT Xi,INT Yi,INT X0,INT Y0,INT X1,INT Y1);
void dp_circle(INT Xi,INT Yi,INT r);
void dp_clip(INT cmd, INT X0,INT Y0,INT X1,INT Y1);
void dp_cont(INT X0,INT Y0);
void dp_control(int type, int i1, int i2, int i3, int i4);
void dp_cross(int *X0,int *Y0, char *c);
void dp_cursor(INT curstyp);
void dp_erase(INT mode);
void dp_fillp(INT n,INT *x,INT *y);
void dp_fillr(INT X0,INT Y0,INT X1,INT Y1,
	INT patx,INT paty,INT lenx,INT leny);
void dp_fillt(INT X0,INT Y0,INT X1,INT Y1,INT X2,INT Y2,
	INT patx,INT paty,INT lenx,INT leny);
void dp_fills(INT X0,INT Y0,INT ixy,INT istnd,INT iplmn);
void dp_font(INT Xi);
void dp_gintxt(int cnt, char *s);
void dp_gottxt(char *s);
void dp_gsymb(INT X0,INT Y0,INT X1,INT Y1,INT n,char *s);
void dp_gwid(INT wid);
void dp_info(int *HasMouse, INT *XminDev, INT *YminDev, 
	INT *XmaxDev, INT *YmaxDev, INT *XminClip, 
	INT *YminClip, INT *XmaxClip, INT *YmaxClip, INT *Color);
void dp_label(char *s);
void dp_linec(INT X0,INT Y0,INT X1,INT Y1);
void dp_linemod(char *s);
void dp_move(INT X0,INT Y0);
void dp_pen(INT Xi);
void dp_point(INT X0,INT Y0);
void dp_space(INT X0,INT Y0,INT X1,INT Y1);
