/*	Changes
 *	30 OCT 2002 extended to very short periods for use in exploration
 *	23 JUL 2003 Following Email from Meijian An, 
 *	Department of Geophysics, Institute of Astronomics, 
 *	Geophysical and Atmospheric Sciences, *	University of Sao Paulo
 *	Sao Paulo, Brazil, carefully free memory
 *
 *	27 JAN 2004 - add option to increase frequency resolution by adding 
 *		zeros to time series
 *	12 FEB 2004 - put in error return and termination if the program
 *		sacpom96 cannot be found
 *		23 OCT 2004 - minor clean up in clearscreen in do_pom4.c
 */
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include <sys/stat.h>
#include <unistd.h>

#define 	PNAME 		"do_pom"
#include	"nfmenu.h"
#include	"nmenu.h"

/* DEFINES */

#define		FILE_SAC_BINARY	1
#define		FILE_SAC_ASCII 	2
#define		FILE_UNKNOWN	0

#define		MENU_FILE_QUIT	-1
#define		MENU_FILE_NEXT	-2
#define		MENU_FILE_PREV	-3
#define		MENU_FILE_NUMB	10
 
extern struct menu menu_p1[] ;

#define ON      1
#define OFF     0
/* plot window for POM96 graphs 
	   The plot will be placed between (XL,YL) and (XH,YH) 
	   The POM96.PLT will be shifted by adding (XOFF,YOFF) */
#define XL      0.0
#define XH      10.00
#define YL      1.5
#define YH      8.0
#define XOFF	0.0
#define YOFF	-0.5;
#define MIN(a,b) ( (b) > (a) ? (a):(b) )
#define MAX(a,b) ( (a) > (b) ? (a):(b) )
#define ABS(a)   ( (a) > (0) ? (a): -(a) )
#define LIN 0
#define LOG 1
#include "sacsubc.h"
#include "csstime.h"
#define MAXSACARR 100000

char *ftype[] = { "UNK", "BIN", "ASC"  }; 
char Strout[100], str1out[80];
#include "calplot.h"
#include "grphsubc.h"

/* GLOBAL VARIABLE DECLARATIONS */
	/* display information */
int HasMouse; 
float XminDev, YminDev, 
	XmaxDev, YmaxDev, XminClip, 
	YminClip, XmaxClip, YmaxClip;
int Color;

int black = 1;	/* 1 white background, 0 black background */
int kolor = 1;	/* 1 gray scale, 2 color */

fmenu **file_menu;
int ndfiles;

struct date_time dt_begin, dt_origin, dt_ptime, dt_stime;
struct date_time dt_refer;
float evla, evlo, evdp, stla, stlo, delta, 
	dist, baz, az, gcarc, b, e, o, a, t0;
int npts, nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec;
char kstnm[9], kcmpnm[9], kevnm[9], kevnmc[9];

/* prototypes */
int main(int argc, char **argv , char **arge);
void Inquire_file(int argc, char **cp, int *npage);
int type_file(char *cstr,int *nsamp, int *fsize, char *datetime, char *kstnm, char *kcmpnm);
void mgwrtxt(float x0, float y0, char *str, int cmd, int color);
int Pick_file(int npage, int argc, int *curpage);
	/* ret =  */
int do_page1(int npage, int *curpage);
	/*  ret = */
int do_page2(char *fname);
	/*  ret = */
int do_units(void);
int do_page3(void );
	/*  ret = */
int do_page4(char *type, int dotype);
	/*  ret = */
void do_check(float xl, float yl, float xh, float yh);
void do_reject(float xl, float yl, float xh, float yh);
void show_menu (float x0, float y0, struct menu *men, int size, int *nm);
int inside(float xv, float yv, 
	float xlb, float ylb, float xhb, float yhb);
void clearregion(float xl, float yl, float xh, float yh);

int Units = 0;

char **cp = NULL;

fmenu *p;
extern fmenu *Start, *End;
extern int menu_p1_entry ; 

int main(int argc, char **argv, char **arge)
{
	int i;
	int npage;
	int aargc;
	int curpage = 0;
	int ret, line, page;
	float xl, yl, xh, yh;
	aargc = argc;
	/* make private copy of file names to be examined */
	cp = (char **) calloc(aargc, sizeof (char *));
	/* argc[0] is the name of the executable */

#ifdef DEBUG
		printf("argc: %d aargc: %d\n",argc,aargc);
#endif

	for (i=1 ; i < aargc; i++){
		cp[i-1] = (char *)calloc(strlen(argv[i])+1,
			sizeof (char));
		strcpy(cp[i-1], argv[i]);
		strcat(cp[i-1],"\0");
	}
	/* open graphics - note the calxvig window size can be
		fixed by adding simething like this to
		.Xdefaults:
		plotxvig.calxvig.nm*geometry: 800x600+00+00 
		where nm is temporarily the PNAME
	*/
	ginitf("INTER", PNAME);
	/* get information about the display */
	ginfo(&HasMouse, &XminDev, &YminDev, 
		&XmaxDev, &YmaxDev, &XminClip, 
		&YminClip, &XmaxClip, &YmaxClip,&Color);
		if(Color >= 4)
			black = 0;
		else
			black = 1;
		kolor = Color%4;

	gmesg("Examining Files ");
	/* initialize the nodes for the file menu */
		Start = End = getnode();
	/* file_menu relates to those displayed by page */
		file_menu = (fmenu **)calloc(menu_p1_entry,sizeof(fmenu *));
	/* grab information about all files on the command line */
		Inquire_file(aargc-1,cp,&npage);
	/* OK choose the file for use */
		ret = 0;
		while(ret >= 0){
			/* recompute the page and line entries */
			p = Start;
			ndfiles = 0;
			npage = 0;
			while(p != End){
				ndfiles++;
				line = (ndfiles-1)%10; 
				page = (ndfiles-1)/10; 
				if(page > npage)
					npage = page;
				xl = 0.1;
				yl = 6.0 -0.5*line;
				xh = 9.9;
				yh = yl + 0.5;
				p->xl = xl;
				p->yl = yl;
				p->xh = xh;
				p->yh = yh;
				p->page = page;
#ifdef DEBUG
printf("%d %d %d %d %s\n",p->lstrmx, p->page, p->line, p->type,p->str);
#endif
				p = p->next;
			}
			/* select an ooption */
			ret = Pick_file(npage, aargc, &curpage );
			/* if ret > 0 append */
			if(ret > 0 ){
			}
		}
	/* terminate the program */
	/* eventually invoke the signal mechanism to cause a
		clean up of all memory allocated

		also this may be a place to save the processing
		summary
		so that the command will be nnn Dir
		or nnn will look for the .xxxx file int he current directory
	*/
	p = Start;
	while(p != End)
		deletenode(p);
	gend(1);
	/*
	 * Remove free the cp memory
	 */
	for (i=1 ; i < aargc; i++){
		free(cp[i-1]);
	}
	free (cp);
	return (0);
}

extern int page_entry_max;
int Pick_file(int npage, int aargc, int *curpage)
{
	int i, ls;
	int page_entry;
	i=0;
	gframe(1);
	newpen(1);
	gwidth(0.02);
	gcent(5.0,7.2,0.2,"POM96 File Selection",0.0);
	gwidth(0.00);
	mgwrtxt(0.2,6.8,
		"Type File",
		1,   1);
	mgwrtxt(0.2,6.6,
	"  Stnm    Cmpnm          Npts    Bytes    First Sample Time          Dist Proc",
		0,   4);
	page_entry = page_entry_max;
	p = Start;
	while(p != End){
#ifdef DEBUG
		printf("%d %s %s %d %d %d %s %s %s %d %f\n",
			p->type, ftype[p->type], 
			p->str, p->line, 
			p->fsize, p->nsamp, 
			p->kstnm, p->kcmpnm, 
			p->datetime, p->page,p->dist);
#endif
		if(p->xl > 0.0 && *curpage == p->page){
			sprintf(Strout ,"%3s %-85s", ftype[p->type], p->str);
			ls = strlen(Strout);
			if(ls > 79)ls = 79;
			Strout[ls] = '\0';
			/* special fix for exploration processing */
			if(p->dist < 1.0)
 				sprintf(str1out,
				"   %-8s %-8s %8d %8d %23s %10.5f",
					p->kstnm, p->kcmpnm, p->nsamp,
					p->fsize, p->datetime, p->dist);
			else
 				sprintf(str1out,
				"   %-8s %-8s %8d %8d %23s %10.3f",
					p->kstnm, p->kcmpnm, p->nsamp,
					p->fsize, p->datetime, p->dist);
			ls = strlen(str1out);
			if(ls > 79)ls = 79;
			str1out[ls] = '\0';
			mgwrtxt(p->xl+0.1,p->yl+0.25,Strout, 0,  1);
			mgwrtxt(p->xl+0.1,p->yl+0.05,str1out,0,  4);
			gbox(0.1,p->yl,9.9,p->yh);
			if(p->used > 0){
				do_check(p->xl,p->yl, p->xh,p->yh);
			} else if(p->used < 0){
				do_reject(p->xl,p->yl, p->xh,p->yh);
			}
			page_entry++;
			menu_p1[page_entry].fileptr = i;
			menu_p1[page_entry].action = page_entry -2;
			menu_p1[page_entry].xl = p->xl;
			menu_p1[page_entry].yl = p->yl;
			menu_p1[page_entry].xh = p->xh;
			menu_p1[page_entry].yh = p->yh;
			file_menu[page_entry] = p;
		
		}
		p = p->next;
	}
	/* set menu entries with no data */
	for( i= page_entry +1; i < MENU_FILE_NUMB+3; i++){
		
		menu_p1[i].xl = -1;
		menu_p1[i].yl = -1;
		menu_p1[i].xh = -1;
		menu_p1[i].yh = -1;
		menu_p1[i].action = 0;
	}

	/* now display the menu */
	
	return(do_page1(npage, curpage));
}



/* put in a check mark */
void do_check(float xl, float yl, float xh, float yh)
{
	float x1,y1,x2,y2,x3,y3,yc,xc,ht,hht;
	newpen(2);
	hht = yh - yl;
	ht = 0.6*hht;
	yc = yl + 0.5*hht - 0.5*ht;
	xc = xh - hht;
	x1 = xc - 0.2*ht; y1 = yc + 0.2*ht;
	x2 = xc; y2 = yc;
	x3 = xc + 0.2*ht; y3 = yc + 0.8*ht;
	plot(x1,y1,3);plot(x2,y2,2);plot(x3,y3,2);plot(x3,y3,3);
	newpen(1);
}

/* put in a rejection mark */
void do_reject(float xl, float yl, float xh, float yh)
{
	float x1,y1,x2,y2,x3,y3,x4,y4,yc,xc,ht,hht;
	newpen(1);
	hht = yh - yl;
	ht = 0.4*hht;
	yc = yl + 0.5*hht ;
	xc = xh - hht;
	x1 = xc - 0.5*ht; y1 = yc - 0.5*ht;
	x2 = xc + 0.5*ht; y2 = yc + 0.5*ht;
	x3 = xc - 0.5*ht; y3 = yc + 0.5*ht;
	x4 = xc + 0.5*ht; y4 = yc - 0.5*ht;
	plot(x1,y1,3);plot(x2,y2,2);plot(x3,y3,3);plot(x4,y4,2);plot(x4,y4,3);
	newpen(1);
}


/* select file to be processed from list valid files 
	this takes time since we read the SAC files */
void Inquire_file(int aargc, char **cp, int *npage)
{
	int i,retval;
	int l, lstr;
	struct stat statbuf;
	float xl, yl, xh, yh;
	int action, lstrmx, type, line, fsize, nsamp, page, used;
	char kstnm[9], kcmpnm[9], datetime[24];
	/* get the maximum string length in this menu category */
	action = -1;
	lstr = 0;
	for(i=0 ; i < aargc ; i++){
		l = strlen(cp[i]);
		if(l > lstr)lstr = l;
	}
	ndfiles = 0;
	*npage = 0;
	for(i=0;i< aargc; i++){
		retval = type_file(cp[i],&nsamp,&fsize,datetime,kstnm,kcmpnm);
		if(retval == 1 || retval == 2){
			if( stat(cp[i],&statbuf)==0)
				fsize = statbuf.st_size;
			else
				fsize = -1 ;
			/* get rid of dregs */
			used = 0;
			lstrmx = lstr;
			type = retval;
			ndfiles++;
			line = (ndfiles-1)%10; 
			page = (ndfiles-1)/10; 
			if(page > *npage)
				*npage = page;
			xl = 0.1;
			yl = 6.0 -0.5*line;
			xh = 9.9;
			yh = yl + 0.5;
		
			appendnode(xl, yl, xh, yh,
				cp[i], action, lstrmx, type, line, fsize,
				nsamp, kstnm, kcmpnm, datetime,
				page, used, dist, az, baz);
		}
	}
}


int type_file(char *cp,int *nsamp, int *fsize, char *datetime, char *ksnm, char *kcmnm)
{
	int nberr, naerr, nerr;
	float *data;
	char ostr[80];
	
/*
	struct date_time dt_begin, dt_origin, dt_ptime, dt_stime;
	struct date_time dt_refer;
	float evla, evlo, evdp, stla, stlo, delta, 
		dist, baz, az, gcarc, b, e, o, a, t0;
	int npts, nzyear, nzjday, nzhour, nzmin, 
		nzsec, nzmsec;
	char kstnm[9], kcmpnm[9], kevnm[9], kevnmc[9];
*/


	/* attempt to open as a SAC file */
	nberr = -100;
	naerr = -100;
	brsac(MAXSACARR,cp,  &data, &nberr);
	if(nberr >= 0){
		free(data);
		
	}
	if(nberr < 0){
		arsac(MAXSACARR,cp,  &data, &naerr);
		if(naerr >= 0){
		free(data);
		}
	}
	if(naerr >= 0 || nberr >=0){
		getfhv("EVLA",&evla,&nerr);
		getfhv("EVLO",&evlo,&nerr);
		getfhv("EVDP",&evdp,&nerr);
		getfhv("STLA",&stla,&nerr);
		getfhv("STLO",&stlo,&nerr);
		getfhv("BAZ",&baz,&nerr);
		getfhv("AZ",&az,&nerr);
		getfhv("DELTA",&delta,&nerr);
		getfhv("DIST",&dist,&nerr);
		getfhv("GCARC",&gcarc,&nerr);
		getfhv("B",&b,&nerr);
		getfhv("O",&o,&nerr);
		getfhv("A",&a,&nerr);
		getfhv("T0",&t0,&nerr);
		getnhv("NPTS",&npts,&nerr);
		getnhv("NZYEAR",&nzyear,&nerr);
		getnhv("NZJDAY",&nzjday,&nerr);
		getnhv("NZHOUR",&nzhour,&nerr);
		getnhv("NZMIN",&nzmin,&nerr);
		getnhv("NZSEC",&nzsec,&nerr);
		getnhv("NZMSEC",&nzmsec,&nerr);
		getkhv("KSTNM",kstnm,&nerr);

		getkhv("KCMPNM",kcmpnm,&nerr);
		getkhv("KEVNM",kevnm,&nerr);
		getkhv("KEVNMC",kevnmc,&nerr);


		/* convert the reference to epoch time */
		dt_refer.date = 1000L*nzyear + nzjday;
		dt_refer.hour = nzhour;
		dt_refer.minute = nzmin;
		dt_refer.second = (float)nzsec + (float)nzmsec/1000.0;
		/* convert to epoch */
		htoe(&dt_refer);
		etoh(&dt_refer);
		timeprintstr(&dt_refer,ostr);  
		/* create an entry for origin time */
		dt_origin.epoch = dt_refer.epoch + o;
		etoh(&dt_origin);
		/* now create one for first sample time */
		dt_begin.epoch = dt_refer.epoch + b;
		etoh(&dt_begin);
		/* now create one for first P time */
		if(a != -12345.){
			dt_ptime.epoch = dt_refer.epoch + a;
			etoh(&dt_ptime);
		}
		/* now create one for first S time */
		if(t0 != -12345.){
			dt_stime.epoch = dt_refer.epoch + t0;
			etoh(&dt_stime);
		}
		timestr(&dt_begin,ostr);  
		*nsamp = npts;
		strcpy(datetime, ostr);
		strcpy(ksnm,kstnm);
		strcpy(kcmnm,kcmpnm);
	}
	/* safety check for bad file */
	if(npts <= 0)
		return(FILE_UNKNOWN);
	if(nberr >= 0){
		return(FILE_SAC_BINARY);
	} else if(naerr >= 0 ){
		return(FILE_SAC_ASCII);
	} else {
		return(FILE_UNKNOWN);
	}
}

