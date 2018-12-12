
#include	"nmenu.h"
#include	"calplot.h"
#include	<stdio.h>
#include "csstime.h"

extern int Units;

extern struct date_time dt_begin, dt_origin, dt_ptime, dt_stime;
extern struct date_time dt_refer;
extern float evla, evlo, evdp, stla, stlo, delta, dist, baz, az, gcarc, b, e, o, a, t0;
extern int npts, nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec;
extern char kstnm[9], kcmpnm[9], kevnm[9], kevnmc[9];

void show_menu (float x0, float y0, struct menu *men, int size, int *nm);
int inside(float xv, float yv, 
	float xlb, float ylb, float xhb, float yhb);
void clearregion(float xl, float yl, float xh, float yh);
int do_page2(char *fname);
int do_page3(char *fname, char *unitstr);
int type_file(char *cstr,int *nsamp, int *fsize, char *datetime, char *kstnm, char *kcmpnm);
void gwrstr(float x0, float y0, int color, char *fmt, char *str);
void gwrflt(float x0, float y0, int color, char *fmt, float fval);
void gwrint(float x0, float y0, int color, char *fmt, int   ival);
int do_units(void);

#define		FILE_SAC_BINARY	1
#define		FILE_SAC_ASCII 	2
#define		FILE_UNKNOWN	0

#define		MENU_SAC_REJECT  -1
#define		MENU_SAC_RETURN   0
#define		MENU_SAC_DOMFT	  1
#define		MENU_SAC_UNITS	  2
#define		MENU_SAC_STNM	  3
#define		MENU_SAC_CMPNM	  4
#define		MENU_SAC_EVLA	  5
#define		MENU_SAC_EVLO	  6
#define		MENU_SAC_STLA	  7
#define		MENU_SAC_STLO	  8
#define		MENU_SAC_BAZ	  9
#define		MENU_SAC_AZ	  10
#define		MENU_SAC_DELTA	  11
#define		MENU_SAC_DIST	  12
#define		MENU_SAC_GCARC	  13
#define		MENU_SAC_NPTS	  14
#define		MENU_SAC_OT	  15
#define		MENU_SAC_T0	  16
#define		MENU_SAC_TP	  17
#define		MENU_SAC_TS	  18

static struct menu menu_p2[] = {
	{  -1.0, -1.00, -1.0, -1.0, "Return\0" , MENU_SAC_RETURN, -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "Do MFT\0" , MENU_SAC_DOMFT , -1, 1, -1},
	{  -1.0, -1.00, -1.0, -1.0, "Reject\0" , MENU_SAC_REJECT, -1, 1, -1},
	{   0.1,  6.50,  4.0, 7.00, "Units \0" , MENU_SAC_UNITS , -1, 1, -1},
	{   0.1,  6.00,  4.0, 6.50, "Stanam\0" , MENU_SAC_STNM  , -1, 1, -1},
	{   0.1,  5.50,  4.0, 6.00, "CmpNam\0" , MENU_SAC_CMPNM , -1, 1, -1},
	{   0.1,  5.00,  4.0, 5.50, "EvtLat\0" , MENU_SAC_EVLA  , -1, 1, -1},
	{   0.1,  4.50,  4.0, 5.00, "Evtlon\0" , MENU_SAC_EVLO  , -1, 1, -1},
	{   0.1,  4.00,  4.0, 4.50, "StaLat\0" , MENU_SAC_STLA  , -1, 1, -1},
	{   0.1,  3.50,  4.0, 4.00, "StaLon\0" , MENU_SAC_STLO  , -1, 1, -1},
	{   0.1,  3.00,  4.0, 3.50, "  Az  \0" , MENU_SAC_AZ    , -1, 1, -1},
	{   0.1,  2.50,  4.0, 3.00, " Baz  \0" , MENU_SAC_BAZ   , -1, 1, -1},
	{   0.1,  2.00,  4.0, 2.50, "  DT  \0" , MENU_SAC_DELTA , -1, 1, -1},
	{   0.1,  1.50,  4.0, 2.00, " Dist \0" , MENU_SAC_DIST  , -1, 1, -1},
	{   0.1,  1.00,  4.0, 1.50, " Gcarc\0" , MENU_SAC_GCARC , -1, 1, -1},
	{   4.1,  6.00,  9.0, 6.50, " NPTS \0" , MENU_SAC_NPTS  , -1, 1, -1},
	{   4.1,  5.50,  9.0, 6.00, "  OT  \0" , MENU_SAC_OT    , -1, 1, -1},
	{   4.1,  5.00,  9.0, 5.50, "  T0  \0" , MENU_SAC_T0    , -1, 1, -1},
	{   4.1,  4.50,  9.0, 5.00, "  TP  \0" , MENU_SAC_TP    , -1, 1, -1},
	{   4.1,  4.00,  9.0, 4.50, "  TS  \0" , MENU_SAC_TS    , -1, 1, -1},
};

#define	MENU_UNITS_UNK		0
#define	MENU_UNITS_COUNTS	1
#define	MENU_UNITS_M		2
#define	MENU_UNITS_MPS		3
#define	MENU_UNITS_MPSPS	4
#define	MENU_UNITS_CM		5
#define	MENU_UNITS_CMPS		6
#define	MENU_UNITS_CMPSPS	7
#define	MENU_UNITS_NM		8
#define	MENU_UNITS_NMPS		9
#define	MENU_UNITS_NMPSPS	10


struct menu menu_units[] = {
	{  -1.0, -1.0, -1.0, -1.0, "Unknown\0" , MENU_UNITS_UNK    , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "Counts \0" , MENU_UNITS_COUNTS , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "m      \0" , MENU_UNITS_M      , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "m/s    \0" , MENU_UNITS_MPS    , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "m/s/s  \0" , MENU_UNITS_MPSPS  , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "cm     \0" , MENU_UNITS_CM     , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "cm/s   \0" , MENU_UNITS_CMPS   , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "cm/s/s \0" , MENU_UNITS_CMPSPS , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "nm     \0" , MENU_UNITS_NM     , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "nm/s   \0" , MENU_UNITS_NMPS   , -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "nm/s/s \0" , MENU_UNITS_NMPSPS , -1, 1, -1},
};


/* display SAC file header and permit changes, rewrite of SAC file */
int do_page2(char *fname)
{
	int i, cmd, nmd ;
	char c[2];
	float xv, yv;
	float xl, yl, xh, yh;
	int retval;
	char kstnm[9], kcmpnm[9], datetime[24]; 
	char ostr[80];
	int nsamp, fsize; 
	/* obtain information about the file 	*/
	retval = type_file(fname,&nsamp,&fsize,datetime,kstnm,kcmpnm);
	/* */
	if(retval < 1 || retval > 2)return(MENU_SAC_REJECT);
	/* put up the menus			*/
	gmesg(fname);
/*
	mgwrtxt(1.0,7.25,fname,0,1000);
*/
	show_menu(1.0, 7.5, menu_p2,sizeof(menu_p2),&nmd);
	
	/* display the header values */
	for(i=0 ; i < nmd ; i++){
		yl = menu_p2[i].yl;
		xl = menu_p2[i].xl;
		yh = menu_p2[i].yh;
		xh = menu_p2[i].xh;
		cmd = menu_p2[i].action;
		switch(cmd ){
			case MENU_SAC_UNITS:
				gwrstr(xl+1.0,yl,1,"    %8s",menu_units[Units].str);
				break;
			case MENU_SAC_STNM:
				gwrstr(xl+1.0,yl,1,"    %8s",kstnm);
				break;
			case MENU_SAC_CMPNM:
				gwrstr(xl+1.0,yl,1,"    %8s",kcmpnm);
				break;
			case MENU_SAC_EVLA:
				gwrflt(xl+1.0,yl,1,"%12.5f",evla);
				break;
			case MENU_SAC_EVLO:
				gwrflt(xl+1.0,yl,1,"%12.5f",evlo);
				break;
			case MENU_SAC_STLA:
				gwrflt(xl+1.0,yl,1,"%12.5f",stla);
				break;
			case MENU_SAC_STLO:
				gwrflt(xl+1.0,yl,1,"%12.5f",stlo);
				break;
			case MENU_SAC_BAZ:
				gwrflt(xl+1.0,yl,1,"%12.5f",baz);
				break;
			case MENU_SAC_AZ:
				gwrflt(xl+1.0,yl,1,"%12.5f",az);
				break;
			case MENU_SAC_DELTA:
				gwrflt(xl+1.0,yl,1,"%12.5f",delta);
				break;
			case MENU_SAC_DIST:
				gwrflt(xl+1.0,yl,1,"%12.5f",dist);
				break;
			case MENU_SAC_GCARC:
				gwrflt(xl+1.0,yl,1,"%12.5f",gcarc);
				break;
			case MENU_SAC_NPTS:
				gwrint(xl+1.0,yl,1,"%12d",npts);
				break;
			case MENU_SAC_OT:
				timeprintstr(&dt_origin,ostr); 
				gwrstr(xl+1.0,yl,1,"%32s",ostr);
				break;
			case MENU_SAC_T0:
				timeprintstr(&dt_begin,ostr); 
				gwrstr(xl+1.0,yl,1,"%32s",ostr);
				break;
			case MENU_SAC_TP:
				if(a != -12345.){ 
					timeprintstr(&dt_ptime,ostr); 
					gwrstr(xl+1.0,yl,1,"%32s",ostr);
				}
				break;
			case MENU_SAC_TS:
				if(t0 != -12345.){ 
					timeprintstr(&dt_stime,ostr); 
					gwrstr(xl+1.0,yl,1,"%32s",ostr);
				}
				break;
		}
	}
	/* response to mouse selections */
	cmd = -1;
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		for(i=0 ; i < nmd ; i++){
			yl = menu_p2[i].yl;
			xl = menu_p2[i].xl;
			yh = menu_p2[i].yh;
			xh = menu_p2[i].xh;
			if(inside(xv,yv,xl,yl,xh,yh)) {
				cmd = menu_p2[i].action;
				switch(cmd ){
					case MENU_SAC_REJECT:
						gframe(1);
						return(MENU_SAC_REJECT);
					case MENU_SAC_RETURN:
						return(MENU_SAC_RETURN);
					case MENU_SAC_DOMFT:
						gframe(1);
						return (do_page3(fname,menu_units[Units].str));
					case MENU_SAC_UNITS:
						Units = do_units();
						clearregion(0.0,0.0,10.0,0.9);
				gwrstr(xl+1.0,yl,1,"    %8s",menu_units[Units].str);
						cmd = -1;
						break;
					case MENU_SAC_CMPNM:
						printf("MENU_SAC_CMPNM   \n");
						cmd = -1;
						break;
					case MENU_SAC_EVLA:
					 	printf("MENU_SAC_EVLA    \n");
						cmd = -1;
						break;
					case MENU_SAC_EVLO:
						printf("MENU_SAC_EVLO    \n");
						cmd = -1;
						break;
					case MENU_SAC_STLA:
						printf("MENU_SAC_STLA    \n");
						cmd = -1;
						break;
					case MENU_SAC_STLO:
					 	printf("MENU_SAC_STLO    \n");
					 	cmd = -1;
						break;
					case MENU_SAC_BAZ:
					 	printf("MENU_SAC_BAZ     \n");
					 	cmd = -1;
						break;
					case MENU_SAC_AZ:
					 	printf("MENU_SAC_AZ      \n");
					 	cmd = -1;
						break;
					case MENU_SAC_DELTA:
					 	printf("MENU_SAC_DELTA   \n");
					 	cmd = -1;
						break;
					case MENU_SAC_DIST:
					 	printf("MENU_SAC_DIST    \n");
					 	cmd = -1;
					case MENU_SAC_GCARC:
					 	printf("MENU_SAC_GCARC\n");
					 	cmd = -1;
						break;
					case MENU_SAC_NPTS:
					 	printf("MENU_SAC_NPTS    \n");
					 	cmd = -1;
						break;
					default:
						cmd = -1;
				}
			}
		}
	}
	/* control should never get this far but we need a dummy return */
	return(MENU_SAC_RETURN);
}

int do_units(void)
{
	int i, cmd, nmd;
	char c[2];
	float xv, yv;
	float xl, yl, xh, yh;
	show_menu(0.5, 0.5, menu_units,sizeof(menu_units),&nmd);
	cmd = -1;
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nmd ; i++){
			xl = menu_units[i].xl;
			yl = menu_units[i].yl;
			xh = menu_units[i].xh;
			yh = menu_units[i].yh;
			if(inside(xv,yv,xl,yl,xh,yh)) {
			cmd = menu_units[i].action;
			switch(cmd){
				case MENU_UNITS_UNK:
					return(MENU_UNITS_UNK);
				case MENU_UNITS_COUNTS:;
					return(MENU_UNITS_COUNTS);
				case MENU_UNITS_M:;
					return(MENU_UNITS_M);
				case MENU_UNITS_MPS:;
					return(MENU_UNITS_MPS);
				case MENU_UNITS_MPSPS:;
					return(MENU_UNITS_MPSPS);
				case MENU_UNITS_CM:;
					return(MENU_UNITS_CM);
				case MENU_UNITS_CMPS:;
					return(MENU_UNITS_CMPS);
				case MENU_UNITS_CMPSPS:;
					return(MENU_UNITS_CMPSPS);
				case MENU_UNITS_NM:;
					return(MENU_UNITS_NM);
				case MENU_UNITS_NMPS:;
					return(MENU_UNITS_NMPS);
				case MENU_UNITS_NMPSPS:;
					return(MENU_UNITS_NMPSPS);
				default:
					cmd = -1;
			}
			}
		}
	}
	/* control should never get this far - this is for safety */
	return(MENU_UNITS_UNK);
}
