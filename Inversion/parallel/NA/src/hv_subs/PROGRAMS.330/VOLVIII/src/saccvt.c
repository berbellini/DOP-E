/* convert SAC binary SUN to SAC binary Intel */
/* Robert B. Herrmann
   May 8, 2001
   Seoul National University

	Changes:

	27 APRIL 2003 - build in intelligence to determine if byte order
		should be reversed. the criteria is that 
		o at least one integer header = -12345
		o and NO integer header values are < -99999
		  note we could also try to look at the year doy field
		  for time series (but may not work for spectra other stuff
*/

#include	<stdio.h>
#include	<stdlib.h>
#ifdef MSDOS
#include	<fcntl.h>
#endif
#include	<strings.h>

#define BUFFLOAT	280
#define BUFINTEGER	160
#define BUFCHAR		192
#define BUF		8192

#define	DO_FLOAT	1
#define	DO_INTEGER	2

char	arr[BUF];
int	iarr[BUFINTEGER/4];

void ieee2intel( char *arr, int nbytes, int type);
void gcmdln(int argc, char **argv, int *do_intelligence);
void usage();
int do_sac_swap(void);

int aargc;
char **aargv;



main(int argc, char**argv)
{
int nread;
int do_intelligence;
int do_convert;
	
	/* parse fommand line arguments */
	aargc = argc;
	aargv = argv;
	gcmdln(aargc, aargv,  &do_intelligence);

	/* apply test that tries to guess the order */
	if(do_intelligence != 0){
		rewind(stdin);
		do_convert = do_sac_swap();
		rewind(stdin);
	} else {
		do_convert = 1;
	}
fprintf(stderr,"do_intelligence %d do_convert %d\n",do_intelligence, do_convert);


#ifdef MSDOS
	setmode(fileno(stdin),O_BINARY);
	setmode(fileno(stdout),O_BINARY);
#endif
	if ((nread=fread(arr, sizeof(char), BUFFLOAT, stdin)) == BUFFLOAT){
		if(do_convert==1)ieee2intel( arr, BUFFLOAT, DO_FLOAT );
		fwrite(arr,sizeof(char), BUFFLOAT, stdout);
	} else {
		exit (1);
	}

	if ((nread=fread(arr, sizeof(char), BUFINTEGER, stdin)) == BUFINTEGER){
		if(do_convert==1)ieee2intel( arr, BUFINTEGER, DO_FLOAT );
		fwrite(arr,sizeof(char), BUFINTEGER, stdout);
	} else {
		exit (1);
	}

	if ((nread=fread(arr, sizeof(char), BUFCHAR, stdin)) == BUFCHAR){
		fwrite(arr,sizeof(char), BUFCHAR, stdout);
	} else {
		exit (1);
	}

	while((nread=fread(arr, sizeof(char), BUF, stdin))> 0) {
		if(do_convert==1)ieee2intel( arr, nread/4, DO_FLOAT );
		fwrite(arr, sizeof(char), nread, stdout) ; 
	}
	
}

void ieee2intel( char *arr, int nbytes, int type)
/*	arr	array of bytes to be converted in place 
	nbytes	number of  bytes to be converted
	type	0 float
*/
{
	char a, b, c, d;
	char *cp;
	int i;
	for (i=0; i < nbytes*4 ; i+=4 ){
		a= arr[i]  ;
		b= arr[i+1];
		c= arr[i+2];
		d= arr[i+3];
		arr[i]   = d;
		arr[i+1] = c;
		arr[i+2] = b;
		arr[i+3] = a;
	}
}

void gcmdln(int argc, char **argv, int *do_intelligence){

	/* parse the command line flags */
	*do_intelligence = 0;
	while(argc-- > 1){
		if(*argv[1] == '-'){
			switch(argv[1][1]){ 
			case 'I':
				*do_intelligence = 1;
				break;
			case '?':
                                usage();
                                break;
			case 'h':
                                usage();
                                break;
                        default:
                                break;                                                                  }
                }
                argv++;
        }                                 
}

void usage()
{
fprintf(stderr,	"Usage: saccvt [-I] [-h] [-?]\n");
fprintf(stderr,	"   Convert SAC binary  IEEE to INTEL format \n");
fprintf(stderr,	"   Convert SAC binary  INTEL to IEEE format \n");
fprintf(stderr,	"   All 4 byte integers and floats (a,b,c,d) are\n");
fprintf(stderr,	"   transposed to (d,c,b,a)\n");
fprintf(stderr, " Example: saccvt < SAC_BINARY > tmp ; mv tmp SAC_BINARY\n");
fprintf(stderr,	" -I   (default none)  intelligently guess whether to convert\n");
fprintf(stderr,	" -h   (default none)  this help message\n");
fprintf(stderr,	" -?   (default none)  this help message\n");
exit(0);
}

int do_sac_swap(void)
{
	/* the SAC file consists of a
	 * FLOAT HEADER		280 bytes
	 * INTEGER HEADER	160 bytes
	 * CHARACTER HEADER	192 bytes
	 * FLOAT DATA
	 */
	/* skip to integer header */ 
	int neg12345 = 0;
	int havebigneg = 0;
	int i;
	int nread;

	fseek(stdin, (long) BUFFLOAT, SEEK_SET);
	/* read integer header */
	if ((nread=fread(iarr,sizeof(int),BUFINTEGER/4,stdin))==BUFINTEGER/4){
		/* now apply the rules */
		for(i=0;i < BUFINTEGER/4;i++){
			if(iarr[i] == -12345)neg12345++;
			if(iarr[i] < -99999)havebigneg++;
		}
	} else {
		exit (1);
	}
	if(neg12345 ==0 && havebigneg > 0 )
		return(1);
	else
		return(0);

}
