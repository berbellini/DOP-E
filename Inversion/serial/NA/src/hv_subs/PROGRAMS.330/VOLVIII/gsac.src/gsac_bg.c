#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	"gsac.h"
#include	"gsac_docommand.h"
#include	"gsac_arg.h"
#include	"calplot.h"

#define BG_WIN   0
#define BG_PLT   1
#define BG_GRAY  2
#define BG_COLOR 3
#define BG_GEOM  4
#define BG_NORM  5
#define BG_REV   6
#define BG_DEF   7

struct arghdr bgarg[] = {
	{BG_PLT, "PLT"		, IHDR, 0, 0, NO, "", 1},
	{BG_WIN, "X"		, IHDR, 0, 0, NO, "", 1},
	{BG_WIN, "W"		, IHDR, 0, 0, NO, "", 1},
	{BG_GRAY,"GRAY"		, IHDR, 0, 0, NO, "", 2},
	{BG_COLOR,"COLOR"	, IHDR, 0, 0, NO, "", 1},
	{BG_GEOM,"GEOM"		, IHDR, 0, 2, NO, "GEOM width height", 2},
	{BG_REV ,"REVERSE"	, IHDR, 0, 0, NO, "", 1},
	{0,	""		, IHDR, 0, 0, NO, "",-1}
};

static int bg_int[10];
static int bg_wid=800;
static int bg_hgt=640;
static int bg_color = YES;  /* YES = color else gray */
static int bg_revvideo = NO; /* YES reversed video */
static int bg_putenv = NO;
static char bg_envstr[100];



void gsac_set_param_bg(int ncmd, char **cmdstr)
{
	int i;

	if(ncmd == 1)
		return;
	/* reinitialize used since I do not want to preserve the value */
	for(i=0 ; bgarg[i].key[0] != '\0' ; i++){
		bgarg[i].used = 0 ;
	}
	if(testarg(ncmd, cmdstr, bgarg, NO, YES))
	       	return	;
	/* systematically go through the arguments to see if they have
	 * been correctly invoked. If the have, then get the values */
	for(i=0 ; bgarg[i].key[0] != '\0' ; i++){
		if(bgarg[i].used > 0){
			if(bgarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, bgarg[i].key, 
					bgarg[i].mfit,bgarg[i].narg, bg_int );
			}
			switch(bgarg[i].id){
				case BG_PLT:
					if(gsac_control.plotdevice != PLT)
						gsac_control.plotchange = YES;
					gsac_control.plotdevice = PLT;
					break;
				case BG_WIN:
					if(gsac_control.inpltmode ==YES){
						gend(0);
						gsac_control.inpltmode = NO;
					}
					if(gsac_control.plotdevice != WIN)
						gsac_control.plotchange = YES;
					gsac_control.plotdevice = WIN;
					break;
				case BG_GRAY:
					bg_color = NO;
					bg_putenv = YES;
					break;
				case BG_COLOR:
					bg_color = YES;
					bg_putenv = YES;
					break;
				case BG_REV:
					bg_revvideo = YES;
					bg_putenv = YES;
					break;
				case BG_NORM:
					bg_revvideo = NO;
					bg_putenv = YES;
					break;
				case BG_GEOM:
					bg_wid = bg_int[0];
					bg_hgt = bg_int[1];
					bg_putenv = YES;
					break;

			}
		}
	}
	if(bg_putenv == YES && gsac_control.xvigenv == NO &&
			gsac_control.everinteractive == NO){
		sprintf(bg_envstr,"PLOTXVIG=:-g:%dx%d:",bg_wid,bg_hgt);
		if(bg_revvideo == YES)
			strcat(bg_envstr,"-I:");
		if(bg_color == YES)
			strcat(bg_envstr,"-K:");
		else
			strcat(bg_envstr,"-G:");
		gsac_control.xvigenv = YES;
		putenv(bg_envstr);
	}
}

void gsac_exec_bg(void)
{
}


