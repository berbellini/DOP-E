#ifndef _NMENU_H
#define _NMENU_H
struct menu { 
	float xl;
	float yl;
	float xh ;
	float yh ;
	char *str ;
	int action;
	int lstrmx;
	int visible;	/* 0 not visible */
	int fileptr;
};
#endif
