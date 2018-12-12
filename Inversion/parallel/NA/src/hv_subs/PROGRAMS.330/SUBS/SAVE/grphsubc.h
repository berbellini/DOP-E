#ifndef _GRAPHSUB_C
void gbox(float xl,float yl,float xh,float yh);
void gcent(float xx,float yy,float ht,char *string,float angle);
void gright(float xx,float yy,float ht,char *string,float angle);
void gleft(float xx,float yy,float ht,char *string,float angle);
void dology(float x0,float y0,float yleng,float ymax,float ymin,
	float sizey,int ticlft,int lablft,int dopow, int ly, char *str);
void dologx(float x0,float y0,float xleng,float sxmax,float sxmin,
	float sizex,int ticup,int labtop,int dopow, int lx, char *str);
void dnlinx(float x0,float y0,float xleng,float xmax,float xmin,
	float sizex,int ticup,int labtop,int dopow, int lx, char *str);
void dolinx(float x0,float y0,float xleng,float xmax,float xmin,
	float sizex,int ticup,int labtop,int dopow, int lx, char *str);
void dolnx(float x0,float y0,float xleng,float xmax,float xmin,
        float sizex,int ticup,int labtop,int dopow,
	int llx, char *titlex, int doasis);
void dnliny(float x0,float y0,float yleng,float ymax,float ymin,
	float sizey,int ticlft,int lablft,int dopow, int ly, char *str);
void doliny(float x0,float y0,float yleng,float ymax,float ymin,
	float sizey,int ticlft,int lablft,int dopow, int ly, char *str);
void dolny(float x0,float y0,float yleng,float ymax,float ymin,
        float sizey,int ticlft,int lablft,int dopow,
	int lly, char *titley, int doasis);
void fillit(char *cmd, float rad, float x0, float y0);
void curvit(char *cmd, float rad, float x0, float y0);
void drawcv(int jj, float xval[], float yval[]);
void rclip(float xmin,float xmax,float ymin,float ymax,float *xc1,
	float *yc1,float *xc2,float *yc2, float x0,float yy0,float x1,float yy1,int *iplt)  ;
void gsubsc(float x,float y,float ht,char *s1,int n1,char *s2,int n2);
void gsupsc(float x, float y,float ht,char *s1,int n1,char *s2,int n2);
void gsubsup(float x,float y,float ht,char *s1,int n1,char *s2,int n2,char *s3,int n3);

#define _GRAPHSUB_C
#endif /* grphsubc.h */
