/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GFONT                                                 c
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
Revision History:
	27 DEC 2000 -- saved the font in the Gcplt structure
*/
#include "dosubs.h"
#include "calplot.h"
#include "callib.h"

extern struct Gcplot Gcplt;

void gfont(int ifont)
{
	INT i2font;
	i2font = ifont;
	Gcplt.ifont = ifont;
	if(i2font < 0)i2font = 0;
	(*do_font)(i2font);
}
