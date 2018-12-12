
#include	"nfmenu.h"
#include	"nmenu.h"
#include	"calplot.h"
#include	<stdio.h>
#include	<string.h>

void show_menu (float x0, float y0, struct menu *men, int size, int *nm);
int inside(float xv, float yv, 
	float xlb, float ylb, float xhb, float yhb);
void clearregion(float xl, float yl, float xh, float yh);
int do_page2(char *fname);

#define		MENU_FILE_QUIT	-1
#define		MENU_FILE_NEXT	-2
#define		MENU_FILE_PREV	-3
#define		MENU_FILE_NUMB	10

static char strout[80];
extern fmenu **file_menu;
fmenu *q;
extern int ndfiles;
char ostrs[80]; 
char ostrr[80]; 

extern int HasMouse; 
extern float XminDev, YminDev, 
	XmaxDev, YmaxDev, XminClip, 
	YminClip, XmaxClip, YmaxClip;
extern int Color;

struct menu menu_p1[] = {
	{  -1.0, -1.0, -1.0, -1.0, "Quit\0" , MENU_FILE_QUIT, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, " <  \0" , MENU_FILE_PREV, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, " >  \0" , MENU_FILE_NEXT, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1}
};

int do_page1(int npage, int *curpage)
{
	int i, cmd, nmd ;
	char c[2];
	float xv, yv;
	float xl, yl, xh, yh;
	int cur_page;
	int ret;
	/* output a mode menu for mode selection */
	cur_page = *curpage;  
	show_menu(1.0, 0.5, menu_p1,sizeof(menu_p1),&nmd);
	/* place current mode at top of page */
	cmd = -1;
	sprintf(strout,"Page %d of %d",cur_page+1, npage+1);
	gmesg(strout);
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nmd ; i++){
			xl = menu_p1[i].xl;
			yl = menu_p1[i].yl;
			xh = menu_p1[i].xh;
			yh = menu_p1[i].yh;
			if(inside(xv,yv,xl,yl,xh,yh)) {
				cmd = menu_p1[i].action;
				if(cmd == MENU_FILE_QUIT){
					gframe(1);
					return(-1);
				} else if(cmd == MENU_FILE_PREV){
					cur_page--;
					if(cur_page < 0)cur_page=npage ;
					*curpage = cur_page;
					gframe(2);
					ginfo(&HasMouse, &XminDev, &YminDev, 
						&XmaxDev, &YmaxDev, &XminClip, 
						&YminClip, &XmaxClip,&YmaxClip,&Color);
					return (1);
				} else if(cmd == MENU_FILE_NEXT){
					cur_page++;
					if(cur_page >npage )cur_page= 0;
					*curpage = cur_page;
					gframe(2);
					ginfo(&HasMouse, &XminDev, &YminDev, 
						&XmaxDev, &YmaxDev, &XminClip, 
						&YminClip, &XmaxClip,&YmaxClip,&Color);
					return (1);
				} else if(cmd >= 1 && cmd <= 10){
					*curpage = cur_page;
					gframe(2);
					ginfo(&HasMouse, &XminDev, &YminDev, 
						&XmaxDev, &YmaxDev, &XminClip, 
						&YminClip, &XmaxClip,&YmaxClip,&Color);
					q = file_menu[i];
					ret=do_page2(q->str); 
					q->used = ret;
	if(ret > 1){
	/* insert nodes for the matched and residual files */
	strcpy(ostrs,file_menu[i]->str);
	strcpy(ostrr,file_menu[i]->str);
	strcat(ostrs,"s");   
	strcat(ostrr,"r");   
	insertnode(q, xl, yl, xh, yh,
        	ostrs,  q->action,  q->lstrmx,
		q->type,  q->line,  q->fsize,
        	q->nsamp, q->kstnm, q->kcmpnm, 
		q->datetime, q->page, q->used, q->dist, q->az, q-> baz) ;
	insertnode(q, xl, yl, xh, yh,
        	ostrr,  q->action,  q->lstrmx,
		q->type,  q->line,  q->fsize,
        	q->nsamp, q->kstnm, q->kcmpnm, 
		q->datetime, q->page, q->used, q->dist, q->az, q-> baz) ;
	}
					gframe(1);
					return(1);
				}
				break;
			}
		}
	}
	/* never get here but we must have a return */
	return (-1);
}
