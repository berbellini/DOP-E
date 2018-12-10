/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GUNIT                                                 c
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
#include "sysinit.h"
#include "calplot.h"

extern struct Gcplot Gcplt;

void gunit(char *str)
{
	int junit;
	if(strncmp(str,"cm",2)==0 || strncmp(str,"CM",2)==0 )
			junit = 1;
		else
			junit = 0;
	if(junit != Gcplt.iunit){
			if(junit == 1){
				Gcplt.xstp = Gcplt.xstp * 0.3937;
				Gcplt.ystp = Gcplt.ystp * 0.3937;
			} else {
				Gcplt.xstp = Gcplt.xstp / 0.3937;
				Gcplt.ystp = Gcplt.ystp / 0.3937;
			}
		}
	Gcplt.iunit = junit;
}
