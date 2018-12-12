/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GOTTXT                                                c
c                                                                     c
c      COPYRIGHT (C)  1986, 1989 R. B. Herrmann                       c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
*/
/* implementation of gottxt for devices using internal label routine */

#ifdef MSDOS
#define INT long
#else
#define INT int
#endif

extern void di_gwid(INT x);
extern void di_label(char *s);

void di_gottxt(char *s)
{
	di_gwid(0);
	di_label(s);
}
