#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	"sacsubc.h"
#include	"csstime.h"

/*
	Program to examine sac header
	Usage:
	
	saclhdr -ISSAC sacfile
		returns 0 for success -1 for no success 
	saclhdr -anyheadervalue sacfile
		returns actual header value

CHANGES:
	07 JUN 2002 - the output formats for NZYEAR NZJDAY NZHOUR
		NZMIN and NZMSEC have fixed width with
		0 fill (%4.4d %3.3d %2.2d %2.2d %2.2d %3.3d)
		to permit returned values to define unique file/
		directory names
	07 MAR 2003 - wrong test in for loop to remove trailing
		blanks from string
	11 MAR 2003 - stupidities repaired
	16 JAN 2005 - use brsach instead of brsac - do not need data
		but want header of potentially large trace file
	15 NOV 2005 - permit more than one entry per line, also implement
			-NL to put a new line at the end
	18 MAR 2006 - added -NZMON -NZDAY to return Month and Day this
		complements the -NZJDAY - the if(donl) to avoid the space
		is bad - could put if(donl)printf(" ") for each loop
	31 MAR 2007 - modified format for print real from %f to %g
		to handle -1e+38 correctly
*/
void gcmdln(int argc, char **argv,  int *ncmd, int *donl, char *sacfile);
void usage(void);
#define CMAX 100
#define MAXSACARR 100
#define NUMSTR 50
char command_string[NUMSTR][CMAX];
char sacfile[CMAX];
struct date_time t;
char str[100];

main(int argc, char **argv)
{
	int nberr, naerr, nerr;
	float *data;	/* array of trace values not needed */
	float fval;
	int doy, sec, msec;
	int ival;
	int donl = 0;	/* output new line at the end */
	int i, j, k, nc, ncmd;
	char cval[9];

	/* parse command line */
	gcmdln(argc, argv,  &ncmd, &donl, sacfile);
	/* now we respond to the command_string 
		to return the desired information */
	/* attempt to open the SAC file */
	
	/* attempt to open as a SAC file */
	nberr = -100;
	naerr = -100;
	brsach(sacfile,  &nberr);
	if(nberr < 0){
		arsach(sacfile,  &naerr);
	}
	if(naerr >= 0 || nberr >=0){
		/* check if valid SAC file by looking at number of points */
		getnhv("NPTS    ",&ival,&nerr);
		if(ival <= 0 ) {
			printf("-1");exit(0);
		}
		/* this is a valid SAC file */
		/* get information for time */
		getnhv("NZYEAR",&t.year, &nerr);
		getnhv("NZJDAY",&doy, &nerr);
		getnhv("NZHOUR",&t.hour, &nerr);
		getnhv("NZMIN",&t.minute, &nerr);
		getnhv("NZSEC",&sec, &nerr);
		getnhv("NZMSEC",&msec, &nerr);
		t.doy = (long)doy;
		t.second = (float)sec + 0.001*(float)msec;
		/* convert from human to epoch */
		month_day(&t);
		mdtodate(&t);
		sec = (int)t.second;
		/* parse the commands */
		for(nc=0; nc < ncmd; nc++){
		if(strncmp(command_string[nc],"ISSAC",5) == 0 ){
			printf("1");exit(0);
		}
		
		/* since I do not know what the command is, check each type */
			/* float header */
			getfhv(command_string[nc],&fval,&nerr);
			if(nerr==0){
				if(donl) {
					printf("%g ",fval);
				} else {
					printf("%g",fval);
				}
			}
			/* special test for Month Day time value */
				/* special for month and day */
			if ( strncmp(command_string[nc],"NZMON",5)==0){
				if(donl)
					printf("%2.2d ",t.month);
				else
					printf("%2.2d",t.month);
			} else if ( strncmp(command_string[nc],"NZDAY",5)==0){
				if(donl)
					printf("%2.2d ",t.day);
				else
					printf("%2.2d",t.day);
			}
			/* integer header */
			getnhv(command_string[nc],&ival,&nerr);
			if(nerr==0){
			/* special output for some time/date fields */
			if ( strncmp(command_string[nc],"NZYEAR",6)==0
				&& ival >= 0 && ival <= 9999) {
				if(donl)
					printf("%4.4d ",ival);
				else
					printf("%4.4d",ival);
			} else if ( strncmp(command_string[nc],"NZJDAY",6)==0
				&& ival >= 0 && ival <= 366) {
				if(donl)
					printf("%3.3d ",ival);
				else
					printf("%3.3d",ival);
			} else if ( strncmp(command_string[nc],"NZHOUR",6)==0
				&& ival >= 0 && ival < 24) {
				if(donl)
					printf("%2.2d ",ival);
				else
					printf("%2.2d",ival);
			} else if ( strncmp(command_string[nc],"NZMIN",5)==0
				&& ival >= 0 && ival < 60) {
				if(donl)
					printf("%2.2d ",ival);
				else
					printf("%2.2d",ival);
			} else if ( strncmp(command_string[nc],"NZSEC",5)==0
				&& ival >= 0 && ival < 60) {
				if(donl)
					printf("%2.2d ",ival);
				else
					printf("%2.2d",ival);
			} else if ( strncmp(command_string[nc],"NZMSEC",6)==0
				&& ival >= 0 && ival < 1000) {
				if(donl)
					printf("%3.3d ",ival);
				else
					printf("%3.3d",ival);
			} else {
				if(donl)
					printf("%d ",ival);
				else
					printf("%d",ival);
			}
			}
			/* string header */
			getkhv(command_string[nc],cval,&nerr);
			if(nerr==0){
			/* for string remove trailing but not leading blanks */
			/* note getkhv uses strncpy but forces cval[8]='\0' */
				j = strlen(cval);
				for(i=j-1,k= -1;i>=0 ; i--){
					if(cval[i] !=' ')
						k = i;
					if(k == -1)
						cval[i] = '\0';
				}
				if(donl)
					printf("%s ",cval);
				else
					printf("%s",cval);
				}
		}
			if(donl)
				printf("\n");
			exit(0);
	} else {
		printf("-1");exit(0);
	}

}

void gcmdln(int argc, char **argv,  int *ncmd, int *donl, char *sacfile)
{
	char *cp;
	int nc;
	nc = 0;
	*donl = 0;
	*ncmd = 0;
	strcpy(command_string[nc]," ");
	strcpy(sacfile       ," ");
	if(argc ==1)usage();
	while(argc-- > 1 ) {
		if(*argv[1] == '-'){
			cp = argv[1];
			if(strncmp(cp,"-h",2) == 0)usage();
			if(strncmp(cp,"-?",2) == 0)usage();
			cp++;
			if(strncmp(cp,"NL",2) == 0){
				*donl = 1;
			} else {
				strcpy(command_string[nc],cp);
				nc++ ;
				strcpy(command_string[nc]," ");
			}
		} else {
			cp = argv[1];
			strcpy(sacfile,cp);
		}
		argv++;
	}
	*ncmd = nc;
}

void usage(void)
{	
	fprintf(stderr,"saclhdr [-?] [-h] -Cmd[s] -NL sacfile]\n");
	fprintf(stderr," Return SAC header value\n");
	fprintf(stderr," -Cmd\tSAC header entry, e.g., -DIST -AZ\n");
	fprintf(stderr," -NL\t (default false) force newline - do not use this for setting shell variables\n");
	fprintf(stderr," sacfile\tSAC binary or alpha file name\n");
	fprintf(stderr," \fIn addition -NZMON -NZDAY return month day\n");
	fprintf(stderr," -?\tThis help screen\n");
	fprintf(stderr," -h\tThis help screen\n");
	exit(0);
}
