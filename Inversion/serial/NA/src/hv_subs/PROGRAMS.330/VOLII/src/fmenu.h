#ifndef _FMENU_H
#define _FMENU_H
struct fmenu { 
	float xl;
	float yl;
	float xh ;
	float yh ;
	char *str ;
	int action;
	int lstrmx;
	int type;
	int line;
	int fsize;
	int nsamp;
	char kstnm[9];
	char kcmpnm[9];
	char datetime[24];
	int page;
	float dist;
	float az;
	float baz;
	int used;
};
#endif
