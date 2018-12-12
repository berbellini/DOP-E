/* NOTES
 * before ploting traces in a group determine the extremes of the absolute
 * B and E so that the time scale can be defined
 */

/* Changes:
	30 MAR 2011 - added a q to the y/n/b query 
*/
#include	<stdio.h>
#include	<string.h>
#include	"gsac.h"

#include	"gsac_plot.h"
#include	"gsac_sac.h"
#include	"gsac_arg.h"
#include	"gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;
extern int  sfgetline(FILE *fp, char s[], int lim);
extern void gsac_show_plotpk(float x0, float y0, float xlen, float ylen, int ns, int ne, float ts, float te, float dy, int ntrc, double twin,int numperframe, int pabsolute, float* uy, float* ly, int setup, int ylimctrl, float ylimlow , float ylimhigh, int p1overlay);
extern void gsac_exec_doxlim(void);
extern void gsac_exec_hold(void);


#define	P1_PERPLOT	0
#define	P1_OVERLAY	1
#define	P1_ABSOLUTE	2
#define	P1_RELATIVE	3


struct arghdr p1arg[] = {
	{P1_PERPLOT,  "PERPLOT"	, NHDR, 0, 1, NO, "PERLOT [n|OFF]", 1},
	{P1_OVERLAY,  "OVERLAY"	, YHDR, 0, 1, NO, "OVERLAY [ON|OFF] ", 1},
	{P1_ABSOLUTE, "ABSOLUTE", RHDR, 0, 0, NO, "ABSOLUTE", 1},
	{P1_RELATIVE, "RELATIVE", RHDR, 0, 0, NO, "RELATIVE", 1},
	{0,	""		, IHDR, 0, 0, NO, "",-1}
};


/* these are temporary variables only used here */
float p1_real[10];
int   p1_int [10];
int   p1_yn;
int   p1_num;

static int p1perplot = -1;
static int p1overlay = NO;
static int p1absolute = YES;

static int p1inp2mode = NO;

extern struct plt_ctl *pmap ;


void gsac_set_param_plot(int ncmd, char **cmdstr)
{
	int i;
	int HasMouse;
	float XminDev, YminDev,
	        XmaxDev, YmaxDev, XminClip,
	        YminClip, XmaxClip, YmaxClip;
	int Color;

	if(p1inp2mode == YES){
		/* force a reset */
		p1overlay = NO ;
		p1inp2mode = NO;
	}
	/* initialize graphics */
	if(gsac_control.plotinit == NO){
		if(gsac_control.plotdevice==WIN){
			ginitf("INTEM","GSAC");
			printf("Initializing Interactive Graphics\n");
			gmesg("Initializing Interactive Graphics");
			gsac_control.everinteractive = YES;
			gsac_control.plotinit = YES;
			gsac_control.plotchange = NO;
			gsac_control.XmaxDev = XmaxDev;

			ginfo(&HasMouse, &XminDev, &YminDev, 
				&XmaxDev, &YmaxDev, &XminClip, 
				&YminClip, &XmaxClip, &YmaxClip,&Color);
			gsac_control.XmaxDev = XmaxDev;
			gsac_control.YmaxDev = YmaxDev;
			if(Color >= 4)
				gsac_control.black = 0;
			else
				gsac_control.black = 1;
			gsac_control.kolor = Color%4;

		}
	}
	/* AUGUST 16, 2007 
		if the command is PLOT1  or P1 force multitrace
		if the command is PLOT2  or P2 force overlay
		Requested by Ghassan Al-Eqabi */
	if(strncmp(gsac_strupr(cmdstr[0]),"PLOT1",5) == 0 || 
			strncmp(gsac_strupr(cmdstr[0]),"P1",2) == 0 ){
		p1perplot = -1;	
	printf("recognized PLOT1\n");
		p1overlay = NO;
	}
	if(strncmp(gsac_strupr(cmdstr[0]),"PLOT2",5) == 0 || 
			strncmp(gsac_strupr(cmdstr[0]),"P2",2) == 0 ){
		p1overlay = YES;
	printf("recognized PLOT2\n");
	p1inp2mode = YES;
	}
	/* if no option on the command line, continue */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, p1arg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; p1arg[i].key[0] != '\0' ; i++){
		if(p1arg[i].used > 0){
			if(p1arg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, p1arg[i].key, 
					p1arg[i].mfit,p1arg[i].narg, p1_real);
			} else if(p1arg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, p1arg[i].key, 
					p1arg[i].mfit,p1arg[i].narg, p1_int );
			} else if(p1arg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, p1arg[i].key, 
					p1arg[i].mfit,p1arg[i].narg, &p1_yn );
			} else if(p1arg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, p1arg[i].key, 
					p1arg[i].mfit,p1arg[i].narg, &p1_num );
			}
			switch(p1arg[i].id){
				case P1_PERPLOT:
					p1perplot = p1_num;
					break;
				case P1_ABSOLUTE:
					p1absolute = YES;
					break;
				case P1_RELATIVE:
					p1absolute = NO;
					break;
				case P1_OVERLAY:
					if(p1_yn == NO)
						p1overlay = NO;
					else if(p1_yn == YES)
						p1overlay = YES;
					break;

			}
		}
	}

}


void gsac_exec_plot(void)
{
	int k, kkk;
	int ntrc;
	float x0, y0, xlen, ylen, dy;
	char pltname[10];
	int numperframe;
	char instr[100];

	double ts, te;
	double twin;

	float uy, ly;	/* not used here indicates clip region */


	/* if there are no traces return */
	ntrc = gsac_control.number_otraces;
	gsac_exec_hold();
	if(ntrc < 1)
		return;
	gsac_exec_doxlim();
	/* initialize */
	if(gsac_control.plotdevice==WIN){
		if(gsac_control.hold == NO){
			gframe(2);
		} else if(gsac_control.hold == 1){
			gframe(2);
			gsac_control.hold++ ;
		}
	}
	if(pmap == (struct plt_ctl *)NULL)
		pmap = (struct plt_ctl *)calloc(ntrc, sizeof(struct plt_ctl));
	else
		pmap = (struct plt_ctl *) realloc(pmap,ntrc*sizeof(struct plt_ctl));
	gsac_control.uxmul = 1.0;
	gsac_control.uxcen = 0.5;
	/* beginning of trace plot 
	 * For absolute time plot this is just 
	 * 	gsac_control.tzend - gsac_control.tzbeg
	 * For relative, we need the maximum time window of
	 * 	a trace
	 * Note then when we build in plot limits this will change
	 * But here we introduce the concept time for the X-axis */
	if(p1absolute == YES) {
		twin = (gsac_control.endmaxx - gsac_control.begminx);
		ts = 0.0;
		te = twin;
	} else {
		twin = 0.0;
		for ( k=0 ; k < ntrc ; k++){
			ts = sacdata[k].tzbegx - sacdata[k].tzref ;
			te = sacdata[k].tzendx - sacdata[k].tzref ;
			if((te-ts) > twin )
				twin = te - ts ;
		}

	}

	if(p1perplot > 0)
		if(p1perplot > ntrc)
			numperframe = ntrc;
		else
			numperframe = p1perplot;
	else
		numperframe = ntrc;

	/* we must define the time window for the plot */

	xlen = gsac_control.xlen ;
	ylen = gsac_control.ylen ;
	x0   = gsac_control.x0 ;
	y0   = gsac_control.y0 ;

	if(p1overlay == YES){
		dy = ylen ;
	} else {
		dy = ylen / numperframe;
	}
	for ( kkk=0 ; kkk < ntrc  ; kkk+=numperframe){
		if(gsac_control.plotdevice==PLT){
			if(gsac_control.hold == NO || 
				(gsac_control.hold != NO && 
				 gsac_control.inpltmode == NO)){
			sprintf(pltname,"P%3.3d.PLT",gsac_control.plotcount++);
			printf("Initializing %s\n",pltname);
			ginitf(pltname,"GSAC");
			gsac_control.inpltmode = YES ;
			}
		} else {
			if(gsac_control.hold == NO)gframe(2);
		}
		gclip("off", x0, y0, x0+xlen, y0+ylen);
		/* put up the traces */
		if(p1overlay == YES)
			if(gsac_control.ylim_ctrl != YLIM_USR)
		gsac_show_plotpk( x0,  y0,  xlen,  ylen,  kkk,  
			MIN(kkk+numperframe,ntrc),  ts,  te,  dy, ntrc, 
			twin, numperframe, p1absolute, &uy, &ly, YES, 
			YLIM_ALL, gsac_control.ylim_low, gsac_control.ylim_high, YES);
			else
		gsac_show_plotpk( x0,  y0,  xlen,  ylen,  kkk,  
			MIN(kkk+numperframe,ntrc),  ts,  te,  dy, ntrc, 
			twin, numperframe, p1absolute, &uy, &ly, YES, 
			YLIM_USR, gsac_control.ylim_low, gsac_control.ylim_high, YES);
		else
		gsac_show_plotpk( x0,  y0,  xlen,  ylen,  kkk,  
			MIN(kkk+numperframe,ntrc),  ts,  te,  dy, ntrc, 
			twin, numperframe, p1absolute, &uy, &ly, YES, 
			gsac_control.ylim_ctrl, gsac_control.ylim_low, 
				gsac_control.ylim_high, NO);
		if(gsac_control.plotdevice!=PLT){
			if( p1perplot > 0){
				printf("More? y/n/b/q\n");
				sfgetline(stdin, instr, 100);
				switch (instr[0]){
				case 'n' :
				case 'N' :
				case 'q' :
				case 'Q' :
					goto jump;
				case 'b' :
				case 'B' :
					kkk = MAX(-numperframe, 
						kkk-numperframe-numperframe);
					break;
				default:
					break;
				}
			}
		}
		if(gsac_control.plotdevice==PLT && gsac_control.hold == NO){
			/* force new Pnnnn.PLT on next call */
			gsac_control.plotchange = NO; 
			gend(0);
			gsac_control.inpltmode = NO;
		}
	}

	gclip("off", x0, y0, x0+xlen, y0+dy);
	/* clean up */
jump:
	if(gsac_control.hold == NO){
		if(gsac_control.everinteractive == YES){
			ginitf("INTEM","GSAC");
		}
	}



}
