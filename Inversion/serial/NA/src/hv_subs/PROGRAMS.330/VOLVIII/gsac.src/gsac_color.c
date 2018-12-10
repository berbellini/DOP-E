#include "gsac_docommand.h"
#include	<stdio.h>
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include	<string.h>

#define	COLOR_ON	0
#define	COLOR_OFF	1
#define	COLOR_DEFAULT	2
#define	COLOR_RAINBOW	3
#define	COLOR_LIST	4
#define COLOR_DUMMY	-1

static int do_color = NO;
static int color_model = COLOR_DEFAULT;

#define COLOR_WHITE  0
#define COLOR_BLACK  1
#define COLOR_RED    2
#define COLOR_GREEN  3
#define COLOR_BLUE   4
#define COLOR_ORANGE 5
#define COLOR_CYAN   6
#define COLOR_YELLOW 7

static int color_default[] = {
	COLOR_WHITE, COLOR_BLACK, COLOR_RED, COLOR_GREEN, COLOR_BLUE, COLOR_ORANGE, COLOR_CYAN,  COLOR_YELLOW};
#define COLOR_DEF_NUM 8
#define COLOR_LIST_NUM 80
static int color_list[COLOR_LIST_NUM];
static int list_num ;

struct colortable {
	int kolor;
	char *name;
};

struct colortable colors[] = {
	{ COLOR_WHITE,	"WHITE"},
	{ COLOR_BLACK,	"BLACK"},
	{ COLOR_RED,	"RED"  },
	{ COLOR_GREEN,	"GREEN"},
	{ COLOR_BLUE,	"BLUE" },
	{ COLOR_ORANGE,	"ORANGE"},
	{ COLOR_CYAN,	"CYAN"},
	{ COLOR_YELLOW,	"YELLOW"}
};

/* we make some dummy vlues for the colors */

struct arghdr colorarg[] = {
	{COLOR_ON,  "ON"	 , YHDR, 0, 0, NO, "",2},
	{COLOR_OFF,   "OFF"	 , YHDR, 0, 0, NO, "",2},
	{COLOR_DEFAULT, "DEFAULT", YHDR, 0, 0, NO, "",1},
	{COLOR_RAINBOW, "RAINBOW", YHDR, 0, 0, NO, "",2},
	{COLOR_LIST, "LIST", YHDR, 0, 0, NO, "",1},
	{COLOR_DUMMY, "BLACK", YHDR, 0, 0, NO, "",3},
	{COLOR_DUMMY, "RED", YHDR, 0, 0, NO, "",2},
	{COLOR_DUMMY, "GREEN", YHDR, 0, 0, NO, "", 1},
	{COLOR_DUMMY, "BLUE", YHDR, 0, 0, NO, "",3},
	{COLOR_DUMMY, "ORANGE", YHDR, 0, 0, NO, "", 2},
	{COLOR_DUMMY, "CYAN", YHDR, 0, 0, NO, "", 1},
	{COLOR_DUMMY, "YELLOW", YHDR, 0, 0, NO, "", 1},
	{COLOR_DUMMY, "WHITE", YHDR, 0, 0, NO, "", 1},
	{0,	""		 , IHDR, 0, 0, NO, "",-1}
};

static char tmpstr[1000];
extern struct sacfile_ *sacdata;

void gsac_set_param_color(int ncmd, char **cmdstr)
{
	int i, k, n ;
	int docmpcolor;
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, colorarg, NO, YES))
	       	return	;
	/* parse the commnds */

	docmpcolor = NO;
	for(i=0 ; colorarg[i].key[0] != '\0' ; i++){
		if(colorarg[i].used > 0){
			switch(colorarg[i].id){
				case COLOR_ON:
					do_color = YES;
					break;
				case COLOR_OFF:
					do_color = NO;
					break;
				case COLOR_DEFAULT:
					color_model = COLOR_DEFAULT;
					do_color = YES;
					break;
				case COLOR_RAINBOW:
					color_model = COLOR_RAINBOW;
					do_color = YES;

					break;
				case COLOR_LIST:
					/* cycle through cmdstr until we
					 * see the word LIST then get the 
					 * colors */
					list_num = 0;
					for(k=0 ; k < ncmd; k++){
						strcpy(tmpstr,cmdstr[k]);
						gsac_strupr(tmpstr);
						if(docmpcolor){
							for(n=0;n<COLOR_DEF_NUM;n++){

								if(strcmp(tmpstr,colors[n].name) ==0){
					if(list_num < COLOR_LIST_NUM -1)
									color_list[list_num++] = colors[n].kolor;						
								}
							}
						}
						if(strcmp(tmpstr,"LIST") ==0){
							docmpcolor = YES;
						}
					}
					break;
			}
		}
	}
	/* safety check */
	if(docmpcolor  && list_num > 0 ){
		color_model = COLOR_LIST;
		do_color = YES;
	}
}

void gsac_exec_color(void)
{
}

void gsac_setcolor(int onoff, int k, int ntrc)
{
	int kolor, kn;
	if(onoff){
		if(do_color){
			switch(color_model){
				case COLOR_RAINBOW:
					if(ntrc == 1)
						kolor = 1000 ;
					else
						kolor = 1000 + 100*(k)/(ntrc-1) ;
					/* safety checks */
					if(kolor < 1000)
						kolor = 1000;
					if(kolor > 1100)
						kolor = 1100;
					newpen(kolor);
					break;
				case COLOR_DEFAULT:
					kn = k%COLOR_DEF_NUM;
					kolor = color_default[kn];
					newpen(kolor);
					break;
				case COLOR_LIST:
					kn = k%list_num;
					kolor = color_list[kn];
					newpen(kolor);
					break;
			}
		} else {
			newpen(1);
		}
	} else {
		newpen(1);
	}
}
