#include <stdio.h>
#include <stdlib.h>
#include "csstim.h"

/*
	this is a simple program to change the date time group
	by the addition or subtraction of an offset.
	The purpose of this is to specify the origin time, but to
	acquire the data 60 seconds before the origin time. This
	is required for data at short distance so that the deconvolution and
	band pass filtering is not strongly affected by noise

 Usage: redodate YEAR MONTH DAY HOUR MINUTE SEC MILLISECOND OFFSET
	redodate 2008 12 22 01 02 03 456  -60

	yields
		 2008 12 22 01 01 03 456
*/

/* protytypes */
void gcmdln(int argc, char **argv, int *nzyear, int *nzmonth, int *nzday,
              int *nzhour, int *nzmin, int *nzsec, int *nzmsec, int *offset);
void usage ();


main(int argc, char **argv)
{
	double epoch;
	int nzyear, nzmonth, nzday, nzhour, nzmin, nzsec, nzmsec;
	int offset;
	int nzdoy;

	/* parse the command line */
	gcmdln(argc, argv, &nzyear, &nzmonth, &nzday,
              &nzhour, &nzmin, &nzsec, &nzmsec, &offset);
/*
	printf("%4d %2.2d %2.2d %2.2d %2.2d %2.2d %3.3d\n",
            nzyear,  nzmonth,  nzday,  nzhour,  nzmin, nzsec,  nzmsec);
*/
	/* convert from human to epoch */
	htoe2(nzyear,  nzmonth,  nzday,  nzhour,  nzmin,  
              nzsec,  nzmsec, &epoch);
	/* apply the offset */
	epoch += offset;
	/* convert from epoch to human */
	etoh( epoch, &nzyear, &nzdoy, &nzmonth, &nzday, 
              &nzhour, &nzmin, &nzsec, &nzmsec) ;
        /* output the new time group */
	printf("%4d %2.2d %2.2d %2.2d %2.2d %2.2d %3.3d\n",
            nzyear,  nzmonth,  nzday,  nzhour,  nzmin, nzsec,  nzmsec);

	
}

void gcmdln(int argc, char **argv, int *nzyear, int *nzmonth, int *nzday,
              int *nzhour, int *nzmin, int *nzsec, int *nzmsec, int *offset)
{
	if(argc != 9)usage();
	*nzyear  = atoi( argv[1] );
	*nzmonth = atoi( argv[2] );
	*nzday   = atoi( argv[3] );
	*nzhour  = atoi( argv[4] );
	*nzmin   = atoi( argv[5] );
	*nzsec   = atoi( argv[6] );
	*nzmsec  = atoi( argv[7] );
	*offset  = atoi( argv[8] );
}

void usage()
{
fprintf(stderr," Usage: redodate YEAR MONTH DAY HOUR MINUTE SEC MILLISECOND OFFSET\n");
fprintf(stderr,"        redodate 2008 12 22 01 02 03 456  -60\n");
exit (0);
}
