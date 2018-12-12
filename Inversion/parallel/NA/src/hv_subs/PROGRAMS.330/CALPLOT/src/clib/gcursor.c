/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GCURSOR                                               c
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

#include "dosubs.h"


void gcursor(char *cmd)
{
	INT icmd;
/*	gcursor
	Arrow
	XORArrow
	Plus
	Box
	RubberBand
	Off

	Crosshairs
*/
	if(cmd[0] == 'A' || cmd[0] == 'a')
		/* Arrow */
		icmd = 0;
	else if(cmd[0] == 'X' || cmd[0] == 'x')
		/* XORArrow */
		icmd = 1;
	else if(cmd[0] == 'C' || cmd[0] == 'c')
		/* Crosshairs */
		icmd = 2;
	else if(cmd[0] == 'P' || cmd[0] == 'p')
		/* Plus */
		icmd = 3;
	else if(cmd[0] == 'B' || cmd[0] == 'b')
		/* Plus */
		icmd = 4;
	else if(cmd[0] == 'R' || cmd[0] == 'r')
		/* Plus */
		icmd = 5;
	else if(cmd[0] == 'O' || cmd[0] == 'o')
		/* Plus */
		icmd = 6;
	else if(cmd[0] == 'H' || cmd[0] == 'h')
		/* Hyperbola */
		icmd = 7;
	else 
		/* none */
		icmd = 0;
	(*do_cursor)(icmd);
}
