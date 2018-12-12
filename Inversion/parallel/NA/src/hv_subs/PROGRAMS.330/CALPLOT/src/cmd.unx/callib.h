#ifndef _CALLIB_
#define _CALLIB_
struct Gcplot {
	float owidth;
	float xstp;
	float ystp;
	float xold;
	float yold;
	float xcur;
	float ycur;
	int iunit;
	int ifont;
	float xs;
	float ys;
	float x0;
	float y0;
	int is;		/* index to bit is pattern ipat - plotd */
	int ipatld;	/*  previous pattern used to init - plotd */
	float xsv;	/* offset within current bin in xlen unit - plotd */
	} ;
#endif
