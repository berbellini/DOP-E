/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLOTS                                                 c
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

#include "callib.h"
#include "calplot.h"

extern struct Gcplot Gcplt;

void plots(int ibuf,int nbuf,int ldev)
{
/*
c----- arguments are not used
c----- kept for compatibility with
c----- older uses of call plots
c----- 
c----- USE call pinit() instead
c-----
*/
      pinitf("plot");
}
