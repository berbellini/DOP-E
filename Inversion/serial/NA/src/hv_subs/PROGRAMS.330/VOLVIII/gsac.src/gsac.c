#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#define MAIN
#include	"gsac.h"
#include	"calplot.h"
#undef MAIN

/* Changes:
	28 MAR 2009 - corrected xlim off command
	30 MAR 2009 - added new command to rotate3, e.g.,
		rotate3 to uvw  to recoverr individual sensor
        21 MAY 2009 - add FAULT to combine Grene's functions for DC 
		ugly since this smacks of accessing a data base but I needit
*/


int main(int argc, char **argv)
{
	/* initialize */
	gsac_init();
	
	/* Do the work  */
	while(gsac_parse_command())
		;
	printf("\n");
	if(gsac_control.plotinit == YES)
		/* if in PLT mode close the device */
		if(gsac_control.plotdevice == PLT)
			gend(1);
		/* if we use the plot mode and close that the
		 * interactive window may not get closed so
		 * ensure that we can access interactive and then close it */ 
		if(gsac_control.everinteractive == YES){
			ginitf("INTEM","GSAC");
			gend(1);
		}
	if(gsac_control.prshist != NULL)
		fclose(gsac_control.prshist);
	if(gsac_control.refrpick != NULL){
		fclose(gsac_control.refrpick);
	}
	return (0);
}
