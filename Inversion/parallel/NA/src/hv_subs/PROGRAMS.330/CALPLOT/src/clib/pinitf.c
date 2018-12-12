/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PINITF                                                c
c                                                                     c
c      COPYRIGHT (C)  1997 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      3507 Laclede Avenue                                            c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
*/

#include "sysinit.h"
#include "callib.h"
#include "dosubs.h"
#include "calplot.h"
void openpl(int ls, char *s, int lis, char *iconstr);


struct Gcplot Gcplt = {
	0.0, 	/* owidth */
	1.0, 	/* xstp */
	1.0, 	/* ystp */
	0.0, 	/* float xold */
	0.0, 	/* float yold */
	0.0, 	/* float xcur */
	0.0, 	/* float ycur */
	0,	/* iunit */
	0,	/* ifont */
	0.0, 	/* float xs */
	0.0, 	/* float ys */
	0.0, 	/* float x0 */
	0.0, 	/* float y0 */
	0,	/* int is */	
	-1,	/* int ipatld */	
	0.0	/* int xsv */	
	} ;

void pinitf(char *string)
{
	char iconstr[1];
	INT x0 = 0;
	INT y0 = 0;
	INT x1 = 7620;
	INT y1 = 7620;
	iconstr[0] = '\0';
	openpl(strlen(string),string,strlen(iconstr),iconstr);
/*
c-----
c map 11" to max tek screen
c this is really done in plot where
c we let 11" = 11000 points
c-----
*/
	(*do_space)(x0,y0,x1,y1);
	(*do_clip)(0,-1,-1,-1,-1);
	gwidth(0.0);
	factor(1.0);
	plot(0.0,0.0,-3);
}

void ginitf(char *string, char *iconstr)
{
	INT x0 = 0;
	INT y0 = 0;
	INT x1 = 7620;
	INT y1 = 7620;
	openpl(strlen(string),string,strlen(iconstr),iconstr);
/*
c-----
c map 11" to max tek screen
c this is really done in plot where
c we let 11" = 11000 points
c-----
*/
	(*do_space)(x0,y0,x1,y1);
	(*do_clip)(0,-1,-1,-1,-1);
	gwidth(0.0);
	factor(1.0);
	plot(0.0,0.0,-3);
}
