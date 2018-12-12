#ifndef _NFMENU_H
#define _NFMENU_H
typedef struct Fmenu { 
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
	int used;
	float dist;
	float az;
	float baz;
	struct Fmenu *next;
} fmenu;

void *getmemory(int h);
fmenu *getnode(void);
void deletenode(fmenu *p);
void insertnode(fmenu *p, float xl, float yl, float xh, float yh,
	char *str, int action, int lstrmx, int type, int line, int fsize,
	int nsamp, char *kstnm, char *kcmpnm, char *datetime,
	int page, int used, float dist, float az, float baz);
void appendnode(float xl, float yl, float xh, float yh,
	char *str, int action, int lstrmx, int type, int line, int fsize,
	int nsamp, char *kstnm, char *kcmpnm, char *datetime,
	int page, int used, float dist, float az, float baz);
#endif
