void algaxe(float xaxlen,float yaxlen,int nocx,int nocy,char *ttlx,char *ttly,int mtx,int mty,float x1,float y1,float deltax,float deltay);
void axis(float xpage,float ypage,char *ititl,int nchar,float axlen, float angle,float firstv,float deltav);
void curuxy(float *xx, float *yy, float x1, float y1, float deltax, float deltay,int nocx,int nocy,char *ic);
void factor (float fact) ;
void frame(void );
void gend(int mode );
void gframe(int mode );
void gclip(char *cmd, float xlow,float ylow,float xhigh,float yhigh );
void gcontrol(int type, float p1, float p2, float p3, float p4);
void gcursor(char *ctyp);
void gfont(int ifont);
void ginfo(int *HasMouse, float *XminDev, float *YminDev, 
	float *XmaxDev, float *YmaxDev, float *XminClip, 
	float *YminClip, float *XmaxClip, float *YmaxClip,int *Color);
void ginitf(char *string, char *iconstr);
void gmesg(char *mesg);
void gread(char *fname, float X0, float Y0, float Xlow, float Ylow, 
	float Xhigh, float Yhigh, int PageNum, float scalx, float scaly);
void gwrtxt(float xx,float yy,char *text,int flgbln);
void grdtxt(char *text,int lstr);
void gunit(char *str);
void gwidth(float width);
void line(float x[],float y[],int n,int inc,int lintyp,int inteq);
void lined(float x[],float y[],int n,int inc,int lintyp, int inteq,int ipat,float xlen);
void newpen(int j);
void number(float xpage,float ypage,float height, float fpn,float angle,int ndec);
void pend(void );
void pinit(void);
void pinitf(char *string);
void plot(float x,float y,int ipen);
void plotd(float xx,float yy,int ipat,float xlen);
void plots(int ibuf,int nbuf,int ldev);
void pltlgd(float *x,float *y,int n,float x1,float y1, float deltax,float deltay,int lintyp,int inteq, float ht,int nocx,int nocy,int ipat,float xlen);
void pltlog(float x[],float y[],int n,float x1,float y1, float deltax,float deltay,int lintyp,int inteq, float ht,int nocx,int nocy);
void pltscl(float x[],float axlen,int n,float *x1,float *deltax,int nocx);
void gscale(float x[],float xaxlen,int n,int inc);
void shadep(int narr, float *xarr, float *yarr);
void shader(float x1,float y1,float x2,float y2,int ipatx, int ipaty,float xlen,float ylen);
void shadet(float x1,float y1,float x2,float y2, float x3,float y3,int ipatx,int ipaty,float xlen,float ylen);
void shdsei(float x1,float y1,int ixy,int istnd,int iplmn);
void symbol (float xloc,float yloc,float height, char *inbuf,float angle,int nocar);
void where(float *xpag, float *ypag, float *fct);

void curixy(int *ix,int *iy,char *ic);
void curaxy(float *xx,float *yy,char *ic);
void currxy(float *xx,float *yy,char *ic);
void cross(int *ix,int *iy, char *s);
void gcont(float x, float y);
void gmove(float x, float y);
