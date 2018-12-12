#include <stdio.h>
#include <string.h>
#include "calplot.h"

main()
{
	char txt[20];
	char txt1[10];
	char tetx[80];
	char ic[2];
	float xx, yy;

	pinitf("INTER");
	gunit("in");
	newpen(1);
	gwrtxt(1.0,2.0,"HELLO WORLD   ",0);
	gwrtxt(1.0,1.5,"HELLO WORLD   ",1);
	gwrtxt(1.0,1.0,"HELLO WORLD   ",2);
	strcpy(txt,"TODAY    :");
	newpen(2);
	gwrtxt(3.0,5.0,txt,0);
	gwrtxt(3.0,4.5,txt,1);
	gwrtxt(3.0,4.0,txt,2);
	grdtxt(txt1,10);
	newpen(4);
	gwrtxt(0.0,7.0,txt1,1);
/*
c-----
c	put a cursor on the screen, get coordinates when a key is
c	hit. Then plot this at the position. Note that the (0,0)
c	of the character is its lower left corner
c-----
*/
	gcursor("CROSS");
	currxy(&xx,&yy,ic);
/*
c-----
c	echo character input -- beware when doing this with curuxy
c-----
*/
	if(ic[0] <  32){
		/* special case */
		if(ic[0] == 1)
			gwrtxt(xx,yy,"Left Mouse button",0);
		else if(ic[0] == 2)
			gwrtxt(xx,yy,"Right Mouse button",0);
		else if(ic[0] == 3)
			gwrtxt(xx,yy,"Center Mouse button",0);
		else {
	       		ic[0] = ' ';
			gwrtxt(xx,yy,ic,0);
		}
	} else if(ic[0] > 126){
	       	ic[0] = ' ';
		gwrtxt(xx,yy,ic,0);
	} else {
		gwrtxt(xx,yy,ic,0);
	}
/*
c------
c	safety on input character
c-----
*/
/*
c-----
c	put legend at the bottom
c-----
*/
	sprintf(tetx,"For input xx=%10.2f yy=%10.2f character= %c ",xx,yy,ic[0]);
	gwrtxt(0.0,0.1,tetx,1);
	gwrtxt(8.0,0.1,"CR TO CONTINUE",2);
	grdtxt(txt1,2);
	frame();
	pend();
}
	
