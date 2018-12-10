
#include	"nfmenu.h"
#include	"nmenu.h"
#include	"calplot.h"
#include	<stdio.h>

void show_menu (float x0, float y0, struct menu *men, int size, int *nm);
int inside(float xv, float yv, 
	float xlb, float ylb, float xhb, float yhb);
void clearregion(float xl, float yl, float xh, float yh);
int do_page3(void);
int create_cmdfil(void );
void reset_used(void);
void setall_used(void);
static int select_flag = 0;
void mgwrtxt(float x0, float y0, char *str, int cmd, int color); 

#define		MENU_FILE_QUIT	-1
#define		MENU_FILE_NEXT	-2
#define		MENU_FILE_PREV	-3
#define		MENU_FILE_DMIN	-4
#define		MENU_FILE_DMAX	-5
#define		MENU_FILE_SLCT	-6
#define		MENU_FILE_DELT	-7
#define		MENU_FILE_REST	-8
#define		MENU_FILE_DPOM	-9
#define		MENU_FILE_SLCA	-10
#define		MENU_FILE_NUMB	10

static char strout[80];
extern fmenu **file_menu;
fmenu *q;
extern int ndfiles;
char ostrs[80]; 
char ostrr[80]; 

int page_entry_max = 9;       /* first non-zero in menu_p1 */
int menu_p1_entry = 20;
struct menu menu_p1[] = {
	{  -1.0, -1.0, -1.0, -1.0, "   Quit  \0"   , MENU_FILE_QUIT, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    <    \0"   , MENU_FILE_PREV, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    >    \0"   , MENU_FILE_NEXT, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    Dmin \0"   , MENU_FILE_DMIN, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    Dmax \0"   , MENU_FILE_DMAX, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "SelectALL\0"   , MENU_FILE_SLCA, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "  Select \0"   , MENU_FILE_SLCT, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "  Reject \0"   , MENU_FILE_DELT, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "  Reset  \0"   , MENU_FILE_REST, -1, 1, -1},
	{   1.0,  7.5,  2.0,  8.0, "  Do POM \0"   , MENU_FILE_DPOM, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "         \0"   , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "         \0"   , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "         \0"   , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "         \0"   , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "         \0"   , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "         \0"   , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "         \0"   , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "         \0"   , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "         \0"   , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "         \0"   , 10, -1, 0, -1}
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
	show_menu(0.5, 1.0, menu_p1,sizeof(menu_p1),&nmd);
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
					gframe(1);
					return (1);
				} else if(cmd == MENU_FILE_NEXT){
					cur_page++;
					if(cur_page >npage )cur_page= 0;
					*curpage = cur_page;
					gframe(1);
					return (1);
				} else if(cmd == MENU_FILE_DPOM){
					if(create_cmdfil() > 0){
						gframe(1);
						ret = do_page3();
						gframe(1);
						return (1);
					}
				} else if(cmd == MENU_FILE_REST){
					reset_used();
					gframe(1);
					return (1);
				} else if(cmd == MENU_FILE_DELT){
					if(select_flag == 0){
						select_flag = -1;
					} else if(select_flag == 1){
						select_flag = -1;
						gframe(1);
						return (1);
					} else {
						select_flag =  0;
					}
				} else if(cmd == MENU_FILE_SLCT){
					if(select_flag == 0){
						select_flag =  1;
					} else if(select_flag == -1){
						select_flag = 1;
						gframe(1);
						return (1);
					} else {
						select_flag =  0;
					}
				} else if(cmd == MENU_FILE_SLCA){
					setall_used();
					gframe(1);
					return (1);
				} else if(cmd >= 1 && cmd <= page_entry_max+8){
					*curpage = cur_page;
					q = file_menu[i];
					if(select_flag !=  0){
						q = file_menu[i];
						q->used = select_flag;
					}
				}
				break;
			}
		}
	}
}


extern fmenu *Start, *End;  
fmenu *qp;
FILE *cmdfil;
int create_cmdfil(void )
{
	/* create the file cmdfil, placing one SAC file name per line
		if at least one return > 0 else < 0
	*/
	int ret = -1;
	qp = Start;
	while ( qp != End){
		if(qp->used > 0){
			ret++;
			if(ret == 0){
				cmdfil = fopen("cmdfil","w+");
			}
			fprintf(cmdfil,"%s\n",qp->str);
		}
		qp = qp->next;
	}
	if(ret >= 0)
		fclose(cmdfil);
	if(ret < 0 )
		gmesg("No SAC files selected");
	return (ret);
}
void reset_used(void)
{
	/* reset the used parameter field to zero */
	qp = Start;
	while ( qp != End){
		qp->used =  0;
	 	qp = qp->next;
	}
}
void setall_used(void)
{
	/* reset the used parameter field to zero */
	qp = Start;
	while ( qp != End){
		qp->used =  1;
	 	qp = qp->next;
	}
}
