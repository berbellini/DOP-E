/*
 * 	CHANGES
 * 	
  02 NOV 2003 - cleaned up WIN32/MSDOS interface
 */
#include	"nmenu.h"
#include	"calplot.h"
#include	<stdio.h>
#include	<stdlib.h>
#include	<unistd.h>
#include	<string.h>



/* defaults activated when program runs 
	if reset, then the reset values will be used
*/
float permin = 4.0;
float permax = 100.0;
int   nray   = 20 ;
int   nlength = 1;
float vmin   = 2.0;	 
float vmax   = 5.0;	 
int	doshade=0;
int	dotype=0;
int	XaxisPeriod = 1;
int	XaxisLog = 1;
int nper;


void show_menu (float x0, float y0, struct menu *men, int size, int *nm);
int inside(float xv, float yv, 
	float xlb, float ylb, float xhb, float yhb);
void clearregion(float xl, float yl, float xh, float yh);
int do_page4(char *type, int dotype);
float do_perpow(void);
int do_perfac(void);
int do_nryfac(void);
int do_nlnfac(void);
void mgwrtxt(float x0, float y0, char *str, int cmd, int color);
char *my_pathfind( char *path, const char *name, const char *mode);
void gwrstr(float x0, float y0, int color, char *fmt, char *str);
void gwrflt(float x0, float y0, int color, char *fmt, float fval);
void gwrint(float x0, float y0, int color, char *fmt, int   ival);

#define	MENU_PER_p0	 0
#define	MENU_PER_p1	 1
#define	MENU_PER_p2	 2
#define	MENU_PER_p3	 3
#define	MENU_PER_p4	 4
#define	MENU_PER_p5	 5
#define	MENU_PER_p6	 6

float pfac[] = { 0.00001, 0.0001, 0.001, 0.010, 0.100, 1.00, 10.0};

#define	MENU_POMCTL_PMN	0
#define	MENU_POMCTL_PMX	1
#define	MENU_POMCTL_NRY	2
#define	MENU_POMCTL_SHD	3
#define	MENU_POMCTL_POM	4
#define	MENU_POMCTL_RTN	5
#define	MENU_POMCTL_VMN	6
#define	MENU_POMCTL_VMX	7
#define	MENU_POMCTL_XAX	8
#define	MENU_POMCTL_XLN	9
#define	MENU_POMCTL_TYP	10
#define	MENU_POMCTL_EXT	11

static struct menu menu_p3[] = {
	{  -1.0, -1.00, -1.0, -1.0, "Return\0" , MENU_POMCTL_RTN, -1, 1, -1},
	{  -1.0,  -1.0, -1.0, -1.0, "Do POM\0" , MENU_POMCTL_POM, -1, 1, -1},
	{   0.1,  5.50,  4.0, 6.00, "MinPer\0" , MENU_POMCTL_PMN, -1, 1, -1},
	{   0.1,  5.00,  4.0, 5.50, "MaxPer\0" , MENU_POMCTL_PMX, -1, 1, -1},
	{   0.1,  4.50,  4.0, 5.00, "Nray  \0" , MENU_POMCTL_NRY, -1, 1, -1},
	{   0.1,  4.00,  4.0, 4.50, "Shade \0" , MENU_POMCTL_SHD, -1, 1, -1},
	{   0.1,  3.50,  4.0, 4.00, "Type  \0" , MENU_POMCTL_TYP, -1, 1, -1},
	{   5.1,  5.50,  9.0, 6.00, "Vmin  \0" , MENU_POMCTL_VMN, -1, 1, -1},
	{   5.1,  5.00,  9.0, 5.50, "Vmax  \0" , MENU_POMCTL_VMX, -1, 1, -1},
	{   5.1,  4.50,  9.0, 5.00, "X-Axis\0" , MENU_POMCTL_XAX, -1, 1, -1},
	{   5.1,  4.00,  9.0, 4.50, "X-Axis\0" , MENU_POMCTL_XLN, -1, 1, -1},
	{   5.1,  3.50,  9.0, 4.50, "Length\0" , MENU_POMCTL_EXT, -1, 1, -1}
};

static struct menu menu_pow[] = {
	{  -1.0, -1.00, -1.0, -1.0, ".0001\0" , MENU_PER_p0, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, ".001\0"  , MENU_PER_p1, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "0.01\0"  , MENU_PER_p2, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "0.10\0"  , MENU_PER_p3, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "1.00\0"  , MENU_PER_p4, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "10.0\0"  , MENU_PER_p5, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "100.\0"  , MENU_PER_p6, -1, 1, -1}
};

float do_perpow(void)
{
	int i, cmd, nmd;
	char c[2];
	float xv, yv;
	float xl, yl, xh, yh;
	show_menu(0.5, 2.0, menu_pow,sizeof(menu_pow),&nmd);
	cmd = -1;
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nmd ; i++){
			xl = menu_pow[i].xl;
			yl = menu_pow[i].yl;
			xh = menu_pow[i].xh;
			yh = menu_pow[i].yh;
			if(inside(xv,yv,xl,yl,xh,yh)) {
				cmd = menu_pow[i].action;
				clearregion(0.0,0.0,10.0,2.5);
				return (pfac[cmd]);
			}
		}
	}
}


#define	MENU_PER_10	10
#define	MENU_PER_11	11
#define	MENU_PER_12	12
#define	MENU_PER_13	13
#define	MENU_PER_14	14
#define	MENU_PER_15	15
#define	MENU_PER_16	16
#define	MENU_PER_17	17
#define	MENU_PER_18	18
#define	MENU_PER_19	19

#define MENU_PER_20	20
#define MENU_PER_21	21
#define MENU_PER_22	22
#define MENU_PER_23	23
#define MENU_PER_24	24
#define MENU_PER_25	25
#define MENU_PER_27	27
#define MENU_PER_26	26
#define MENU_PER_28	28
#define MENU_PER_29	29

#define MENU_PER_30	30
#define MENU_PER_32	32
#define MENU_PER_34	34
#define MENU_PER_36	36
#define MENU_PER_38	38
#define MENU_PER_40	40
#define MENU_PER_42	42
#define MENU_PER_44	44
#define MENU_PER_46	46
#define MENU_PER_48	48

#define MENU_PER_50	50
#define MENU_PER_55	55
#define MENU_PER_60	60
#define MENU_PER_65	65
#define MENU_PER_70	70
#define MENU_PER_75	75
#define MENU_PER_80	80
#define MENU_PER_85	85
#define MENU_PER_90	90
#define MENU_PER_95	95

static struct menu menu_fac[] = {
	{  -1.0, -1.00, -1.0, -1.0, "1.0\0" , MENU_PER_10, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "1.1\0" , MENU_PER_11, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "1.2\0" , MENU_PER_12, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "1.3\0" , MENU_PER_13, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "1.4\0" , MENU_PER_14, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "1.5\0" , MENU_PER_15, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "1.6\0" , MENU_PER_16, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "1.7\0" , MENU_PER_17, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "1.8\0" , MENU_PER_18, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "1.9\0" , MENU_PER_19, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "2.0\0" , MENU_PER_20, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "2.1\0" , MENU_PER_21, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "2.2\0" , MENU_PER_22, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "2.3\0" , MENU_PER_23, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "2.4\0" , MENU_PER_24, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "2.5\0" , MENU_PER_25, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "2.6\0" , MENU_PER_26, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "2.7\0" , MENU_PER_27, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "2.8\0" , MENU_PER_28, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "2.9\0" , MENU_PER_29, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "3.0\0" , MENU_PER_30, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "3.2\0" , MENU_PER_32, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "3.4\0" , MENU_PER_34, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "3.6\0" , MENU_PER_36, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "3.8\0" , MENU_PER_38, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "4.0\0" , MENU_PER_40, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "4.2\0" , MENU_PER_42, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "4.4\0" , MENU_PER_44, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "4.6\0" , MENU_PER_46, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "4.8\0" , MENU_PER_48, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "5.0\0" , MENU_PER_50, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "5.5\0" , MENU_PER_55, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "6.0\0" , MENU_PER_60, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "6.5\0" , MENU_PER_65, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "7.0\0" , MENU_PER_70, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "7.5\0" , MENU_PER_75, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "8.0\0" , MENU_PER_80, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "8.5\0" , MENU_PER_85, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "9.0\0" , MENU_PER_90, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "9.5\0" , MENU_PER_95, -1, 1, -1}
};

int do_perfac(void)
{
	int i, cmd, nmd;
	char c[2];
	float xv, yv;
	float xl, yl, xh, yh;
	show_menu(0.1, 2.0, menu_fac,sizeof(menu_fac),&nmd);
	cmd = -1;
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nmd ; i++){
			xl = menu_fac[i].xl;
			xh = menu_fac[i].xh;
			yl = menu_fac[i].yl;
			yh = menu_fac[i].yh;
			if(inside(xv,yv,xl,yl,xh,yh)) {
				cmd = menu_fac[i].action;
				clearregion(0.0,0.0,10.0,2.5);
				return (cmd);
			}
		}
	}
}

char *shdstr[] ={"FALSE", " TRUE"};
char *typstr[] = {" UNKNOWN", "    LOVE", "RAYLEIGH"};

/* RBH   not used 
	cmdfil is the list of sac files 
*/

char *pathname;
/* display PER file header and permit changes, rewrite of SAC file */
int do_page3(void)
{
	int i, cmd, nmd ;
	int ret;
	int fac_pmin,fac_pmax;
	float pow_pmin, pow_pmax;
	char c[2];
	char outstr[512];
	char outstr1[512];
	float xv, yv;
	float xl, yl, xh, yh;
	show_menu(1.0, 7.5, menu_p3,sizeof(menu_p3),&nmd);
	for(i=0 ; i < nmd ; i++){
		cmd = menu_p3[i].action;
		xl = menu_p3[i].xl;
		yl = menu_p3[i].yl;
		switch (cmd){
			case MENU_POMCTL_PMN:
				gwrflt(xl+1.0,yl, 1, "%10.4f", permin);
				break;
			case MENU_POMCTL_PMX:
				gwrflt(xl+1.0,yl, 1, "%10.4f", permax);
				break;
			case MENU_POMCTL_NRY:
				gwrint(xl+1.0,yl, 1, "%10d  ", nray );
				break;
			case MENU_POMCTL_SHD:
				gwrstr(xl+1.0,yl, 1, "%10s" , shdstr[doshade]);
				break;
			case MENU_POMCTL_TYP:
				gwrstr(xl+1.0,yl, 1, "%10s" , typstr[dotype]);
				break;
			case MENU_POMCTL_VMN:
				gwrflt(xl+1.0,yl, 1, "%10.3f", vmin );
				break;
			case MENU_POMCTL_VMX:
				gwrflt(xl+1.0,yl, 1, "%10.3f", vmax );
				break;
			case MENU_POMCTL_XAX:
				if(XaxisPeriod)
					gwrstr(xl+1.0,yl,1," %-9s","   Period");
				else
					gwrstr(xl+1.0,yl,1," %-9s","Frequency");
				break;
			case MENU_POMCTL_XLN:
				if(XaxisLog)
					gwrstr(xl+1.0,yl,1," %-9s","      Log");
				else
					gwrstr(xl+1.0,yl,1," %-9s","      Lin");
				break;
			case MENU_POMCTL_EXT:
				gwrint(xl+1.0,yl, 1, "%10d  ", nlength );
				break;
			default:
				break;
		}
	}
	/* response to mouse selections */
	cmd = -1;
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		for(i=0 ; i < nmd ; i++){
			xl = menu_p3[i].xl;
			yl = menu_p3[i].yl;
			xh = menu_p3[i].xh;
			yh = menu_p3[i].yh;
			if(inside(xv,yv,xl,yl,xh,yh))
			{
				cmd = menu_p3[i].action;
				switch(cmd ){
					case MENU_POMCTL_PMN:
						pow_pmin = do_perpow();
						fac_pmin = do_perfac();
						permin = pow_pmin*(float)fac_pmin;
					gwrflt(xl+1.0,yl, 1, "%10.3f", permin);
						cmd = -1;
						break;
					case MENU_POMCTL_PMX:
						pow_pmax = do_perpow();
						fac_pmax = do_perfac();
						permax = pow_pmax*(float)fac_pmax;
					gwrflt(xl+1.0,yl, 1, "%10.3f", permax);
						cmd = -1;
						break;
					case MENU_POMCTL_NRY:
						nray = do_nryfac();
					gwrint(xl+1.0,yl, 1, "%10d", nray);
						cmd = -1;
						break;
					case MENU_POMCTL_EXT:
						nlength = do_nlnfac();
					gwrint(xl+1.0,yl, 1, "%10d", nlength);
						cmd = -1;
						break;
					case MENU_POMCTL_VMN:
						pow_pmax = do_perpow();
						fac_pmax = do_perfac();
						vmin = pow_pmax*(float)fac_pmax;
					gwrflt(xl+1.0,yl, 1, "%10.3f", vmin);
						cmd = -1;
						break;
					case MENU_POMCTL_VMX:
						pow_pmax = do_perpow();
						fac_pmax = do_perfac();
						vmax = pow_pmax*(float)fac_pmax;
					gwrflt(xl+1.0,yl, 1, "%10.3f", vmax);
						cmd = -1;
						break;
					case MENU_POMCTL_XAX:
						XaxisPeriod = 1 - XaxisPeriod;
						if(XaxisPeriod)
					gwrstr(xl+1.0,yl,1," %-9s","   Period");
						else
					gwrstr(xl+1.0,yl,1," %-9s","Frequency");
						cmd = -1;
						break;
					case MENU_POMCTL_XLN:
						XaxisLog = 1 - XaxisLog;
				if(XaxisLog)
					gwrstr(xl+1.0,yl,1," %-9s","      Log");
				else
					gwrstr(xl+1.0,yl,1," %-9s","      Lin");
						cmd = -1;
						break;
					case MENU_POMCTL_SHD:
						doshade = 1 - doshade;
				gwrstr(xl+1.0,yl, 1, "%10s" , shdstr[doshade]);
						cmd = -1;
						break;
					case MENU_POMCTL_TYP:
						dotype++;
						dotype %= 3;
				gwrstr(xl+1.0,yl, 1, "%10s" , typstr[dotype]);
						cmd = -1;
						break;
					case MENU_POMCTL_RTN:
						return(0);
						break;
					case MENU_POMCTL_POM:
fprintf(stderr,"pathname\n");
#ifdef MSDOS
if(nlength == 1)
sprintf(outstr1,"sacpom96 -C cmdfil -PMIN %f -PMAX %f -nray %d -A -VMIN %f -VMAX %f ", permin,permax,nray,vmin,vmax);
else if(nlength == 2)
sprintf(outstr1,"sacpom96 -C cmdfil -PMIN %f -PMAX %f -nray %d -A -VMIN %f -VMAX %f -2", permin,permax,nray,vmin,vmax);
else if(nlength == 4)
sprintf(outstr1,"sacpom96 -C cmdfil -PMIN %f -PMAX %f -nray %d -A -VMIN %f -VMAX %f -4", permin,permax,nray,vmin,vmax);
else if(nlength == 8)
sprintf(outstr1,"sacpom96 -C cmdfil -PMIN %f -PMAX %f -nray %d -A -VMIN %f -VMAX %f -8", permin,permax,nray,vmin,vmax);
#else

        pathname = (char *)my_pathfind(getenv("PATH"), "sacpom96", "rx");
					/* make up the command line */
fprintf(stderr,"%s\n",pathname);
               sprintf(outstr1," -C cmdfil -PMIN %f -PMAX %f -nray %d -A -VMIN %f -VMAX %f ",
                        permin,permax,nray,vmin,vmax);
fprintf(stderr,"%s\n",outstr1);
#endif
					if(dotype == 0)
						strcat(outstr1,"   ");
					else if(dotype == 1)
						strcat(outstr1," -L");
					else if(dotype == 2)
						strcat(outstr1," -R");
					if(doshade)strcat(outstr1," -S");
					if(!XaxisPeriod)strcat(outstr1," -FREQ");
					if(!XaxisLog)strcat(outstr1," -XLIN");
					if(nlength==2)strcat(outstr1," -2");
					if(nlength==4)strcat(outstr1," -4");
					if(nlength==8)strcat(outstr1," -8");
#ifdef MSDOS
					strcpy(outstr,outstr1);
#else
					strcpy(outstr,pathname);
					strcat(outstr,outstr1);
#endif
					printf("%s\n",outstr);
					/* put up message */
					strcpy(outstr1,"Executing: ");
					strncat(outstr1,outstr,60);
					outstr1[70]='\0';
					gmesg(outstr1);
						ret = system(outstr);
						if(ret >= 0){
							gframe(1);
							ret = do_page4(typstr[dotype],dotype);
							return(ret);
						} else {
							return(-1);
						}
					default:
						cmd = -1;
				}
			}
		}
	}
}

#define	MENU_NRY_020	0
#define	MENU_NRY_040	1
#define	MENU_NRY_080	2
#define	MENU_NRY_100	3
#define	MENU_NRY_250	4
#define	MENU_NRY_500	5

int nryval[] = {  20, 40, 80, 100, 250, 500 };


static struct menu menu_nry[] = {
	{  -1.0, -1.0, -1.0, -1.0, "  20 \0" , MENU_NRY_020   , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "  40 \0" , MENU_NRY_040   , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "  80 \0" , MENU_NRY_080   , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, " 100 \0" , MENU_NRY_100   , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, " 250 \0" , MENU_NRY_250   , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, " 500 \0" , MENU_NRY_500   , -1, 1, -1}
};


int do_nryfac(void)
{
	int i, cmd, nmd;
	char c[2];
	float xv, yv;
	float xl, yl, xh, yh;
	show_menu(0.1, 2.0, menu_nry,sizeof(menu_nry),&nmd);
	cmd = -1;
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nmd ; i++){
			xl = menu_nry[i].xl;
			yl = menu_nry[i].yl;
			xh = menu_nry[i].xh;
			yh = menu_nry[i].yh;
			if(inside(xv,yv,xl,yl,xh,yh)) {
				cmd = menu_nry[i].action;
				clearregion(0.0,0.0,10.0,2.5);
				return (nryval[cmd]);
			}
		}
	}
}

#define	MENU_NLN_001	0
#define	MENU_NLN_002	1
#define	MENU_NLN_004	2
#define	MENU_NLN_008	3

int nlnval[] = {  1, 2, 4, 8 };


static struct menu menu_nln[] = {
	{  -1.0, -1.0, -1.0, -1.0, "   1 \0" , MENU_NLN_001   , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "   2 \0" , MENU_NLN_002   , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "   4 \0" , MENU_NLN_004   , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "   8 \0" , MENU_NLN_008   , -1, 1, -1}
};


int do_nlnfac(void)
{
	int i, cmd, nmd;
	char c[2];
	float xv, yv;
	float xl, yl, xh, yh;
	show_menu(0.1, 2.0, menu_nln,sizeof(menu_nln),&nmd);
	cmd = -1;
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nmd ; i++){
			xl = menu_nln[i].xl;
			yl = menu_nln[i].yl;
			xh = menu_nln[i].xh;
			yh = menu_nln[i].yh;
			if(inside(xv,yv,xl,yl,xh,yh)) {
				cmd = menu_nln[i].action;
				clearregion(0.0,0.0,10.0,2.5);
				return (nlnval[cmd]);
			}
		}
	}
}

char *tpath = (char *)NULL;
char *my_pathfind( char *path, const char *name, const char *mode)
{

	char *fpath;
	char *token;
	tpath = (char *)realloc(tpath,(strlen(path)+2)*sizeof(char));
	strcpy(tpath,path);
/*
fprintf(stderr,"strlen(tpath) %d\n",strlen(tpath));
fprintf(stderr,"path  :%s\n",path);
fprintf(stderr,"strlen(path) %d\n",strlen(path));
fprintf(stderr,"name  :%s\n",name);
fprintf(stderr,"mode  :%s\n",mode);
fprintf(stderr,"tpath :%s\n",tpath);
*/
	token = (char *)strtok(tpath,":");
	while(token != NULL){
		fpath = calloc(strlen(token)+strlen(name)+2,sizeof(char));
		strcpy(fpath,token);
		strcat(fpath,"/");
		strcat(fpath,name);
		if(access(fpath, R_OK)==0){
printf("fpath  :%s\n",fpath);
/*
			free(tpath);
	*/
			return(fpath);
		}
		free(fpath);
		token = strtok(NULL,":");
	}
	/*
	free(tpath);
	*/
	fprintf(stderr,"Cannot find the execuable %s \n",name);
	fprintf(stderr,"TERMINATING!!!!\n");
	exit(1);

}
