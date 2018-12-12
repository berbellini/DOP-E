/* File>>> color.c
--
-- %M% -- version %I% (IMEC)            last updated: %E%
--
-- Copyright (c) 1993
-- IMEC vzw
-- Kapeldreef 75
-- B-3001 LEUVEN
-- BELGIUM
--
-- Author   : A. Demaree
--
-- Date     : February 1, 1993
--
-- Function :
--
-- Comment  :
--
-- Review   :
--
*/

/* Revision history
	09 01 99 - put in the second pass at getting colors for
		depth = 8 displays - required for Laptop running Gnome
        12 08 00 - put in a check for 24 bit color
*/


/*------------------------------------------------------------------------------
-- Global include files
------------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <X11/Xlib.h>

/*------------------------------------------------------------------------------
-- Local include files
------------------------------------------------------------------------------*/

#include "xviglocal.h"

/*------------------------------------------------------------------------------
-- Static variable declarations
------------------------------------------------------------------------------*/

static Colormap colormap;
static XColor *color;
static int nr_of_colors;
static int def_colmap[XviG_NR_OF_COLORS][3] =
                              { {   0,   0,   0 },     /* black        */
                                { 255, 255, 255 },     /* white        */
                                { 255,   0,   0 },     /* red          */
                                {   0, 255,   0 },     /* green        */
                                {   0,   0, 255 },     /* blue         */
                                {   0, 255, 255 },     /* cyan         */
                                { 255,   0, 255 },     /* magenta      */
                                { 255, 255,   0 },     /* yellow       */
                                { 255, 127,   0 },     /* orange       */
                                { 127, 255,   0 },     /* green-yellow */
                                {   0, 255, 127 },     /* green-cyan   */
                                {   0, 127, 255 },     /* blue-cyan    */
                                { 127,   0, 255 },     /* blue-magenta */
                                { 255,   0, 127 },     /* red-magenta  */
                                {  95,  95,  95 },     /* dark-gray    */
                                { 175, 175, 175 } };   /* light-gray   */

/*------------------------------------------------------------------------------
-- Local function declarations
------------------------------------------------------------------------------*/
/* RBH REQUIREMENTS FOR XVIG */
static int dodither = 0;
struct mydith {
	int patt;
	int fore;
	int back;
	int stip;
};
static int actualbg ;
#ifndef RBH_NR_OF_DITHERS
#define RBH_NR_OF_DITHERS  22
#endif
static struct mydith dith[35];
/* order of allocating CALPLOT colors to permit
	successful dithering later */
static int calcol[35] =
	{0, 34, 1, 33, 17, 11, 23, 6, 27,
	12, 3, 8, 20, 30, 14, 25, 32,
	2, 4, 5, 7, 10, 16, 19, 22, 29,
	9, 13, 15, 18, 21, 24, 26, 28, 31};  
static void resetdither();
static void Fill01(int colmap[][3], int nr);
static void rgb_to_hsv(float r, float g, float b , float *h, float *s, float *v);
#define ABS(a)   ( (a) > 0   ? (a):-(a) )
  



/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_SetColor(int nr)
{
  int local_nr;
  int k,l,m,n;
  unsigned long fg_color;
  unsigned long bg_color;
  unsigned long f_color;
  unsigned long b_color;
  int nr_tmp;

  nr_tmp = nr % nr_of_colors;

  /* safety for Xoring and color map stuff */
  if(RBH_xoring_state == ON){
	if(color[Abs(nr_tmp)].pixel == 0L){    /* anything XOR'd with 0 will not change */
		local_nr = 1;			/* so you will not see an XOR line or color */
	} else {				/* we do not want invisibility of cursors */
		local_nr = nr_tmp;
	}
  } else {
	local_nr = nr_tmp;
  }

if(!dodither){

	if (XviG_cursor_mode){
		fg_color = color[Abs(local_nr)].pixel;
		bg_color = color[Abs(0)].pixel;
	} else {
		fg_color = color[Abs(local_nr)].pixel;
		bg_color = color[Abs(0)].pixel;
	}

	XSetForeground(XviG_display, XviG_gc, fg_color);
	XSetForeground(XviG_display, XviG_gc_fill, fg_color);

	if(actualbg == 0){
		/* black background - we need a white Xhair */
	XviG_xhair_color = color[Abs(local_nr)].pixel
		^ BlackPixel(XviG_display, XviG_screen_nr);
	} else {
	XviG_xhair_color = color[Abs(local_nr)].pixel
		^ WhitePixel(XviG_display, XviG_screen_nr);
	}
/*
fprintf(stderr,"fg_color: %ld\n",fg_color);
fprintf(stderr,"bg_color: %ld\n",bg_color);
fprintf(stderr,"   local_nr %d, color[Abs(local_nr)].pixel %ld\n",local_nr,color[Abs(local_nr)].pixel);
fprintf(stderr,"   XviG_xhair_color: %ld\n",XviG_xhair_color);
fprintf(stderr,"   BlackPixel(XviG_display, XviG_screen_nr) %ld\n",BlackPixel(XviG_display, XviG_screen_nr));
*/

} else {  /* do dither */

    XSetLineAttributes(XviG_display, XviG_gc,
                       1, LineSolid, CapButt, JoinMiter);
		k = dith[Abs(local_nr)].fore;
		l = dith[Abs(local_nr)].back;
		f_color = color[k].pixel;
		b_color = color[l].pixel;
		if (XviG_cursor_mode){
			fg_color = f_color
				^ BlackPixel(XviG_display, XviG_screen_nr);
			bg_color = b_color
				^ BlackPixel(XviG_display, XviG_screen_nr);
		} else {
			fg_color = f_color;
			bg_color = b_color;
		}
		if(dith[Abs(local_nr)].stip == 0){
			/* use original color */
			RBH_SetDither(0);
		} else {
			/* use dither */
			RBH_SetDither(dith[Abs(local_nr)].patt);
		}
		XSetForeground(XviG_display, XviG_gc, fg_color);
		XSetForeground(XviG_display, XviG_gc_fill, fg_color);
		XSetBackground(XviG_display, XviG_gc, bg_color);
		XSetBackground(XviG_display, XviG_gc_fill, bg_color);

	if(actualbg == 0){
		/* black background - we need a white Xhair */
	XviG_xhair_color = color[Abs(local_nr)].pixel
		^ BlackPixel(XviG_display, XviG_screen_nr);
	} else {
	XviG_xhair_color = color[Abs(local_nr)].pixel
		^ WhitePixel(XviG_display, XviG_screen_nr);
	}
	}
	/*
	-- Save this color to set it back after an XviG_ClearWindow call
	*/

	XviG_save_color = fg_color;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

int XviG_CreateColors(int colmap[][3],
                       int nr, int depth)
{
  int gotcolors = 0;
	float h, s, v;
	float r,g,b;
	float hh, ss, vv;
	float rr,gg,bb;
	int found = 0;
	int k,nr2;
  XColor kolor;

  if (nr == 0)
  {
    colmap = def_colmap;
    nr_of_colors = XviG_NR_OF_COLORS;
  }
  else
    nr_of_colors = nr;

  color = (XColor *) malloc(nr_of_colors * sizeof(XColor));

  /*
  -- Create a color map and allocate the colors in the color map
  */

  colormap = DefaultColormap(XviG_display, XviG_screen_nr);

  {
    int i,ii;
    int j;


	actualbg = colmap[0][0];
    for(ii = 0; ii < nr_of_colors; ii++)
    {
	i = calcol[ii];
	color[i].red   = 257*colmap[i][0];
	color[i].green = 257*colmap[i][1];
	color[i].blue  = 257*colmap[i][2];
	color[i].flags = DoRed | DoGreen | DoBlue;

	/* initialize the dither matrix stuff */
	dith[ i].patt =  0 ; 
	dith[ i].fore = -1 ; 
	dith[ i].back = -1 ; 
	dith[ i].stip = -1 ; 
      if (XAllocColor(XviG_display, colormap, &color[i])){
        gotcolors++;
	dith[i].fore = i;
	dith[i].back = i;
	dith[i].stip = 0;
/*
printf("XviG_CreateColors: ii %d i %d %ld %d %d %d \n",ii,i,color[i].pixel, color[i].red,
	color[i].green,color[i].blue);
*/
	}
    }
/* second pass -- XallocColor actually tries to make a perfect match. We will permit
   a less than perfect match - note the program 'xv' also goes through these shenanigans.
   The source of the difficulty is that the X11 colormap is a finite quantity on 8-bit
   systems (I do not have much experience with difficulties with 15, 16 or 24 bit depths).
   For a depth = 8, other software can grab the colormap cells and we must interpolate.
   Some window manager systems, e.g., the new 'gnome' define all 256 entries 
*/
  if((depth == 8 || depth == 24) && gotcolors != nr_of_colors){
	/* flag the colors that are not yet defined */
	for(ii = 0; ii < nr_of_colors; ii++){
		i = calcol[ii];
		if(dith[i].stip == -1){
		color[i].red   = 257*colmap[i][0];
		color[i].green = 257*colmap[i][1];
		color[i].blue  = 257*colmap[i][2];
		color[i].flags = DoRed | DoGreen | DoBlue;
		/* XColor specifies that r,b,g are unsigned short (0,65535) */
 		r = (float)color[i].red  /65535.;
 		g = (float)color[i].green/65535.;
 		b = (float)color[i].blue /65535.;
		rgb_to_hsv(r, g, b , &h, &s, &v);
		/* now search through the color base - only do for 256 */
  		found = 0;
  		for(k=0 ; k < 256 && found == 0 ; k++){
  			kolor.pixel = k;
  			XQueryColor(XviG_display, colormap, &kolor);
  			/* XColor specifies that r,b,g are unsigned short (0,65535) */
  			rr = (float)kolor.red  /65535.;
  			gg = (float)kolor.green/65535.;
  			bb = (float)kolor.blue /65535.;
  			rgb_to_hsv(rr, gg, bb , &hh, &ss, &vv);
			if(s < 0.1){
				/* gray */
			if( ABS(s -ss ) < 0.05 && ABS(v-vv) < 0.05 && found == 0){
				found = k;
				color[i].red = kolor.red ;
				color[i].green = kolor.green ;
				color[i].blue = kolor.blue ;
				color[i].flags = DoRed | DoGreen | DoBlue;
      				if (XAllocColor(XviG_display, colormap, &color[i])){
        				gotcolors++;
					dith[i].fore = i;
					dith[i].back = i;
					dith[i].stip = 0;
				}
			}
			} else {
				/* color */
			if(ABS(h-hh) < 5 && ABS(s-ss) < 0.1 && ABS(v-vv) < 0.10 && found == 0){
				found = k;
				color[i].red = kolor.red ;
				color[i].green = kolor.green ;
				color[i].blue = kolor.blue ;
				color[i].flags = DoRed | DoGreen | DoBlue;
      				if (XAllocColor(XviG_display, colormap, &color[i])){
        				gotcolors++;
					dith[i].fore = i;
					dith[i].back = i;
					dith[i].stip = 0;
				}
			}
			}
	
   		}



		}
	}
	/* now that we have filled up some of the entries in the color mapping,
		undo some to permit dithering of neighboring colors 
		Als0 so check that the important entrie are filled
		e.g.,
		33 = BLUE
		01 = RED
        */
        /* get rid of duplicate, but do not change 00 01 33 34 */
	for(i=2; i< nr_of_colors-2 ; i++){
		if(color[i].red == color[i+1].red )
			if(color[i].green == color[i+1].green )
				if(color[i].blue == color[i+1].blue ){
					dith[ i].patt =  0 ; 
					dith[i].fore = -1;
					dith[i].back = -1;
					dith[i].stip = -1;
					gotcolors--;
				}
	}

  }
/*
	for(i=0;i<nr_of_colors;i++)
printf("%2d %2d %2d %2d %2d (%6ld,%6ld,%6ld)\n",
	i,
	dith[i].patt, dith[i].fore, dith[i].back, dith[i].stip,
		color[i].red, color[i].green, color[i].blue);
*/

    if(depth == 1)gotcolors = 2;
    if (gotcolors < nr_of_colors){
  	printf("calxvig: Allocating only %d of  %d CALPLOT colors: Using dither.",
		gotcolors,nr_of_colors);
		/* RBH 01/07/1998 */
		/* hard wire for plotxvig */
		/* manually get the basic colors knowing that
			plotxvig uses 35 unique colors 
			IN THIS TEST ALLOCATE ONLY THE FIRST 5*/
			if(gotcolors < 5){
				Fill01(colmap, 1);
				printf(" Using  monochrome\n");
			} else {
				dodither = 1;
				resetdither();
				printf("\n");
			}

  } else {
  	printf("calxvig: Allocated all %d CALPLOT colors\n",nr_of_colors);

  }
}
  return (gotcolors);
}

void Fill01(int colmap[][3], int nr)
{
	/* fill with a 16 element color map 
		a 64 bit map would be better in terms of
		gradation but would loose sharpness 
	*/
	int i, gotcolors = 0, j;
     		color[0].pixel = BlackPixel(XviG_display, XviG_screen_nr);
     		color[1].pixel = WhitePixel(XviG_display, XviG_screen_nr);
		dodither = 1;
		/* define the dithering according to the Create_Filldither */
	dith[ 0].patt = 5;dith[ 0].fore= 1;dith[ 0].back= 0;dith[ 0].stip=0;
	dith[ 1].patt = 6;dith[ 1].fore= 0;dith[ 1].back= 1;dith[ 1].stip=1;
	dith[ 2].patt = 6;dith[ 2].fore= 0;dith[ 2].back= 1;dith[ 2].stip=1;
	dith[ 3].patt = 7;dith[ 3].fore= 0;dith[ 3].back= 1;dith[ 3].stip=1;
	dith[ 4].patt = 7;dith[ 4].fore= 0;dith[ 4].back= 1;dith[ 4].stip=1;
	dith[ 5].patt = 7;dith[ 5].fore= 0;dith[ 5].back= 1;dith[ 5].stip=1;
	dith[ 6].patt = 8;dith[ 6].fore= 0;dith[ 6].back= 1;dith[ 6].stip=1;
	dith[ 7].patt = 8;dith[ 7].fore= 0;dith[ 7].back= 1;dith[ 7].stip=1;
	dith[ 8].patt = 8;dith[ 8].fore= 0;dith[ 8].back= 1;dith[ 8].stip=1;
	dith[ 9].patt = 9;dith[ 9].fore= 0;dith[ 9].back= 1;dith[ 9].stip=1;
	dith[10].patt = 9;dith[10].fore= 0;dith[10].back= 1;dith[10].stip=1;
	dith[11].patt =10;dith[11].fore= 0;dith[11].back= 1;dith[11].stip=1;
	dith[12].patt =10;dith[12].fore= 0;dith[12].back= 1;dith[12].stip=1;
	dith[13].patt =11;dith[13].fore= 0;dith[13].back= 1;dith[13].stip=1;
	dith[14].patt =11;dith[14].fore= 0;dith[14].back= 1;dith[14].stip=1;
	dith[15].patt =12;dith[15].fore= 0;dith[15].back= 1;dith[15].stip=1;
	dith[16].patt =12;dith[16].fore= 0;dith[16].back= 1;dith[16].stip=1;
	dith[17].patt =13;dith[17].fore= 0;dith[17].back= 1;dith[17].stip=1;
	dith[18].patt =13;dith[18].fore= 0;dith[18].back= 1;dith[18].stip=1;
	dith[19].patt =14;dith[19].fore= 0;dith[19].back= 1;dith[19].stip=1;
	dith[20].patt =14;dith[20].fore= 0;dith[20].back= 1;dith[20].stip=1;
	dith[21].patt =15;dith[21].fore= 0;dith[21].back= 1;dith[21].stip=1;
	dith[22].patt =15;dith[22].fore= 0;dith[22].back= 1;dith[22].stip=1;
	dith[23].patt =16;dith[23].fore= 0;dith[23].back= 1;dith[23].stip=1;
	dith[24].patt =16;dith[24].fore= 0;dith[24].back= 1;dith[24].stip=1;
	dith[25].patt =17;dith[25].fore= 0;dith[25].back= 1;dith[25].stip=1;
	dith[26].patt =17;dith[26].fore= 0;dith[26].back= 1;dith[26].stip=1;
	dith[27].patt =18;dith[27].fore= 0;dith[27].back= 1;dith[27].stip=1;
	dith[28].patt =18;dith[28].fore= 0;dith[28].back= 1;dith[28].stip=1;
	dith[29].patt =19;dith[29].fore= 0;dith[29].back= 1;dith[29].stip=1;
	dith[30].patt =19;dith[30].fore= 0;dith[30].back= 1;dith[30].stip=1;
	dith[31].patt =19;dith[31].fore= 0;dith[31].back= 1;dith[31].stip=1;
	dith[32].patt =20;dith[32].fore= 0;dith[32].back= 1;dith[32].stip=1;
	dith[33].patt =20;dith[33].fore= 0;dith[33].back= 1;dith[33].stip=1;
	dith[34].patt =21;dith[34].fore= 0;dith[34].back= 1;dith[34].stip=0;
}


void resetdither()
{
	/* look at the dither structure  which is fixed at 35 entires
		note this is manually set here and in ../../cmd/plotxvig */
	int i, j, k, jj, kk, jjkk;
	float frac;
	int lastcolor, nextcolor;
	/* do a linear search, which is not efficient but which will
		work. Exclude the first [0] and last [34] entries */
	lastcolor = dith[1].fore ;
	if(lastcolor < 0)lastcolor = dith[34].fore;
	/* test that dith[1].fore is set, else rethink everything */
	for(i=2; i < 33; i++){
		if(dith[i].fore < 0){
			j= i+1;
			/* search for nextcolor */
			jj = j;
			for (k=j, kk=0 ; k < 34 && kk == 0 ; k++){
				if(dith[k].fore >= 0){
					nextcolor = dith[k].fore;
					kk = k;
			}
			}
			/* now interpolate */
			jjkk = nextcolor - lastcolor;
			dith[i].fore = nextcolor;
			dith[i].back = lastcolor;
			dith[i].stip = 1;
			frac = (float)(i-lastcolor)/(float)(nextcolor - lastcolor);
			if(jjkk > 4){
				/* use 16 dither */
				dith[i].patt = 5 + (int)(16.0*frac);
			} else {
				/* use 4 dither */
				dith[i].patt = (int)(4.0*frac);
			}

		} else {
			lastcolor = dith[i].fore  ;	
		}

	}
}

/* RBH for second pass */

#define MAX(a,b) ( (a) > (b) ? (a):(b) )
#define MIN(a,b) ( (b) > (a) ? (a):(b) )

/* from Foley and Van Dam */
static void rgb_to_hsv(float r, float g, float b , float *h, float *s, float *v)
{
	float max, min;
	float rc, gc, bc;
	max = MAX(r,g);
	max = MAX(max,b);
	min = MIN(r,g);
	min = MIN(min,b);
	*v = max;
	if(max != 0.0){
		*s = ( max - min ) / max ;
	} else {
		*s = 0;
	}
	if( *s == 0.0){
		*h = 0.0;
	} else {
		rc = ( max - r)/(max - min);
		gc = ( max - g)/(max - min);
		bc = ( max - b)/(max - min);
		if( r == max){
			*h = bc - gc ;
		} else if(g == max){
			*h = 2 + rc - bc;
		} else if(b == max){
			*h = 4 + gc - rc ;
		}
		*h = *h * 60.0;
		if(*h < 0.0) *h = *h + 360.0;
	}
}
