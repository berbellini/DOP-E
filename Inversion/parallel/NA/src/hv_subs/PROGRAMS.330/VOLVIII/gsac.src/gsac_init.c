/* CHANGES
	11 JUL 2010 - add support for the ~ screendump
*/
/* initialize control parameters */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "gsac.h" 
#include "gsac_version.h"
extern int markt_on;

void gsac_init(void)
{
	/* set up buffering for line input output 
	 * This is required so that terminal input is 
	 * responsive under WINDOWS XP
	 * */
	setvbuf(stdin , (char *)NULL, _IONBF, 0);
	setvbuf(stdout, (char *)NULL, _IONBF, 0);
	setvbuf(stderr, (char *)NULL, _IONBF, 0);
	/* set the initial control values */
	gsac_control.number_itraces=0;
	gsac_control.max_number_traces=0;
	gsac_control.more=0;
	gsac_control.plotdevice=WIN;
	gsac_control.plotinit=NO;
	gsac_control.plotchange=NO;
	gsac_control.plotcount=1;
	gsac_control.plotcount_prs=1;
	gsac_control.plotcount_refr=1;
	gsac_control.plotcount_dump=1;
	gsac_control.begmin= 1.0e+30;
	gsac_control.endmax=-1.0e+30;
	gsac_control.fmax = 0.0;
	gsac_control.docut = NO;
	gsac_control.doxlim = NO;
	gsac_control.cutint[0] = H_B;	/* initially refer to start of trace */
	gsac_control.cutoff[0] = 0.0;	
	gsac_control.cutint[1] = H_E;	/* initially refer to  end  of trace */
	gsac_control.cutoff[1] = 0.0;	
	gsac_control.cutepoch[0] = 0.0;
	gsac_control.cutepoch[1] = 0.0;
	gsac_control.xlimint[0] = 5;	/* initially refer to start of trace */
	gsac_control.xlimoff[0] = 0.0;	
	gsac_control.xlimint[1] = 6;	/* initially refer to  end  of trace */
	gsac_control.xlimoff[1] = 0.0;	
	strcpy(gsac_control.xlimkey[0],"B");
	strcpy(gsac_control.xlimkey[1],"E");
	gsac_control.qdp = -1;
	gsac_control.x0 = 1.25;
	gsac_control.y0 = 1.0;
	gsac_control.xlen = 8.0;
	gsac_control.ylen = 6.0;
	gsac_control.grid = NO;		/* underlying grid */
	gsac_control.local = YES;	/* force conversion to local order */
	gsac_control.fft = NO;
	gsac_control.hold = NO;
	gsac_control.inpltmode = NO;
	gsac_control.xvigenv = NO;
	gsac_control.uxcen = 0.5;
	gsac_control.uxmul = 1.0;
	gsac_control.ylim_low =  -1.0;;
	gsac_control.ylim_high = 1.0;
	gsac_control.ylim_ctrl = 0 ;	/* 0 = YLIMOFF */
	gsac_control.ylim_rsc = NO ;	/* user rescale in PPK */
	gsac_control.plotlinx = YES;
	gsac_control.plotliny = YES;
	printf("%s\n",gsac_version);
	markt_on = NO;
	gsac_control.prs = 0;	/* 0=prs 1=prs-refraction 2=prs-reflection */
	gsac_control.prshist = NULL;
	gsac_control.refrpick = NULL;
	gsac_control.xgrid = NO;
	gsac_control.xgrid_type = 2;
	gsac_control.xgrid_color = 1030;
	gsac_control.xgrid_minor = YES;
	gsac_control.ygrid = NO;
	gsac_control.ygrid_type = 2;
	gsac_control.ygrid_color = 1030;
	gsac_control.ygrid_minor = NO;
	gsac_control.background = NO;
	gsac_control.background_color = 0;


}

