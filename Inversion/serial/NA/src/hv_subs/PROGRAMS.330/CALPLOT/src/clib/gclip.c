#include "callib.h"
#include "dosubs.h"
#include "calplot.h"
extern struct Gcplot Gcplt;

void gclip(char *cmd, float xlow,float ylow,float xhigh,float yhigh )
{
	float xx, yy;
	INT ilw, jlw, iup, jup;
	INT icmd;
       	xx = 1000. * (xlow*Gcplt.xstp + Gcplt.xold);
       	yy = 1000. * (ylow*Gcplt.ystp + Gcplt.yold);
	if(xx > 1000000000.0)
		ilw = 1000000000;
	else if(xx < -1000000000.0)
		ilw = -1000000000;
	else
		ilw = xx ;
	if(yy > 1000000000.0)
		jlw = 1000000000;
	else if(yy < -1000000000.0)
		jlw = -1000000000;
	else
		jlw = yy ;
       	xx = 1000. * (xhigh*Gcplt.xstp + Gcplt.xold);
       	yy = 1000. * (yhigh*Gcplt.ystp + Gcplt.yold);
	if(xx > 1000000000.0)
		iup = 1000000000;
	else if(xx < -1000000000.0)
		iup = -1000000000;
	else
		iup = xx ;
	if(yy > 1000000000.0)
		jup = 1000000000;
	else if(yy < -1000000000.0)
		jup = -1000000000;
	else
		jup = yy ;
	if(cmd[0] == 'o' || cmd[0] == 'O'){
		if(cmd[1] == 'n' || cmd[1] == 'N')
			icmd = 1;
		else
			icmd = 0;
	} else
		icmd = 0;
	(*do_clip)(icmd, ilw, jlw, iup, jup);
}
