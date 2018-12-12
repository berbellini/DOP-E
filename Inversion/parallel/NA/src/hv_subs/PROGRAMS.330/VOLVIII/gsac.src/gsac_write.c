/* Change History
18 AUG 2007 
	At the request of Ghassan Al'Eqabi of Washington University
	Implement the options APPEND suffix PREPEND prefix
	If there are N traces in memory, the only permissible
	command sequences are

	WRITE
	WRITE file1 .... fileN
	APPEND suffix
	PREPEND prefix
	APPEND suffix PREPEND prefix

	NOTE char filename[1000] is not safe - later
	use alloc/calloc to avoid overflow

20 JUL 2009 
        Removed the testarg - which means that the syntax of prepend and append is not checked 
        but that the write file1 file2 .. filentrc will work
        I should really do a strcmp her by hand
21 JUL 2009
	put in the proper logic - the real problem is the lack of
	logic in the syntax which forces a lot of special exceptions in the
	parsing
*/

#include        <string.h>
#include "gsac.h"
#include "gsac_docommand.h"
#include "gsac_sac.h"
#include "gsac_sachdr.h"
#include "gsac_arg.h"
#include "csstim.h"
#include	<libgen.h>

extern struct sacfile_ *sacdata;
extern int *sortptr;
static int wrongnumber;
static int overwriteoldtrace;
static char **pchr;
void write_syntax(int ncmd, char **cmdstr);

#define WRITE_APPEND	1
#define WRITE_PREPEND	2

static char  write_append[80];
static char  write_prepend[80];
static int   write_doappend = NO ;
static int   write_doprepend= NO ;
static int   write_ntrace = 0;


struct arghdr writearg[] = {
        {WRITE_APPEND   , "APPEND"  , CHDR, 0, 1, NO, "APPEND suffix ", -1},
        {WRITE_PREPEND  , "PREPEND" , CHDR, 0, 1, NO, "PREPEND prefix ", -1},
        {0,     ""              , IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float write_real[10];
int   write_int [10];
int   write_yn;
int   write_num;


void gsac_set_param_write(int ncmd, char **cmdstr)
{
	int i;
	char instr[80];

	int ntrc;
	/*  systematically go through the options */
	ntrc = gsac_control.number_otraces;

	/* set the default behavior */
	overwriteoldtrace = YES;
	write_doprepend   = NO ;
	write_doappend    = NO ;
	wrongnumber       = NO ;
	write_ntrace = 0;


	
	if(ncmd == 1){
		overwriteoldtrace = YES;
		return;
	}
	/* now test the arguments */
	if(testarg(ncmd, cmdstr, writearg, NO, NO))
		return;
	/* now parse the commands noting if the APPEND or PREPEND are set */
	for(i=0 ; writearg[i].key[0] != '\0' ; i++){
		if(writearg[i].used > 0){
			if(writearg[i].ricell == CHDR){
                                getargs(ncmd, cmdstr, writearg[i].key,
                                        writearg[i].mfit, writearg[i].narg, instr );
			}

			switch(writearg[i].id){
				case WRITE_APPEND:
					if(strlen(instr) < 80){
						strcpy(write_append,instr);
						write_doappend = YES;
						overwriteoldtrace = NO ;
					}
					break;
				case WRITE_PREPEND:
					if(strlen(instr) < 80){
						strcpy(write_prepend,instr);
						write_doprepend = YES;
						overwriteoldtrace = NO ;
					}
					break;
			}
		}
	 }
	/* final safety checks */
	if(write_doappend == YES){
		if(write_doprepend == YES){
			if( (ncmd-1) != 4) {
                        printf("command syntax  error - no write: APPEND suffix PREPEND prefix\n");
				wrongnumber = YES;
			}
		} else {
			if( (ncmd-1) != 2) {
				wrongnumber = YES;
                        	printf("command syntax  error - no write: APPEND suffix \n");
			}
		}
	} else if(write_doappend == NO) {
		if(write_doprepend == YES){
			if( (ncmd-1) != 2) {
                        	printf("command syntax  error - no write: PREPEND prefix\n");
				wrongnumber = YES;
			}
		} else {
			if( (ncmd-1) != ntrc) {
                        	printf("command syntax  error - no write: filelist does not have %d entries  \n",ntrc);
				 wrongnumber = YES;
				overwriteoldtrace = NO;
			} else {
				pchr = cmdstr ;
				overwriteoldtrace = NO;
			}
		}
	}
	if(wrongnumber == YES){
		write_syntax(ncmd, cmdstr);
		return;
	}
}

void gsac_exec_write(void)
{
	int k, ntrc;
	char *tname;
	char *fname ;
	char filename[1000];
	ntrc = gsac_control.number_otraces;

/*
printf("overwriteoldtrace %d\n",overwriteoldtrace);
printf("write_doappend    %d\n",write_doappend   );
if(write_doappend == YES)
	printf("     %s\n",write_append);
printf("write_doprepend   %d\n",write_doprepend   );
if(write_doprepend == YES)
	printf("     %s\n",write_prepend);
printf("wrongnumber       %d\n",wrongnumber   );
printf("write_ntrace      %d\n",write_ntrace   );
*/


	if(ntrc < 1)
		return;
	/* return if we are writing the wrong number of traces */
	if(wrongnumber){
		return;
	}
	/* do we overwrite old trace files */
	if(overwriteoldtrace){
                printf("overwriting traces:\n");
		for ( k=0 ; k < ntrc ; k ++){
			sacdata[k].sachdr.rhdr[H_USER1] = sacdata[k].permin;
			sacdata[k].sachdr.rhdr[H_USER2] = sacdata[k].permax;
			bwsac(sacdata[k].sac_ofile_name,sacdata[k].sachdr,
				sacdata[k].sac_data);
			printf("%s ",sacdata[k].sac_ofile_name);
		}
	} else {
		for ( k=0 ; k < ntrc ; k ++){
			sacdata[k].sachdr.rhdr[H_USER1] = sacdata[k].permin;
			sacdata[k].sachdr.rhdr[H_USER2] = sacdata[k].permax;
			if(write_doprepend == YES || write_doappend == YES){
				fname = strdup(sacdata[k].sac_ofile_name);
				tname = basename(fname);
				if(write_doprepend == YES && write_doappend == NO){
					strcpy(filename,write_prepend);
					strcat(filename,tname);
				} else if(write_doprepend == YES && write_doappend == YES){
					strcpy(filename,write_prepend);
					strcat(filename,tname);
					strcat(filename,write_append);
				} else if(write_doprepend == NO && write_doappend == YES){
					strcpy(filename,tname);
					strcat(filename,write_append);
				}
				/* replace any non-printing characters and space by an underscore */
				cleanstring(filename);
				bwsac(filename,sacdata[k].sachdr,
					sacdata[k].sac_data);
				printf("%s ",filename);
			} else {
				/* use the new names from the command line */
				/* replace any non-printing characters and space by an underscore */
				cleanstring(pchr[k+1]);
				bwsac(pchr[k+1],sacdata[k].sachdr,sacdata[k].sac_data);
				printf("%s ",pchr[k+1]);
			}
		}
	}
	printf("\n");
}

void write_syntax(int ncmd, char **cmdstr)
{
	int i;
	for(i=0;i<ncmd;i++)
		printf("%s ",cmdstr[i]);
	printf("\nThere are %d files in memory\n",ncmd);
	printf("Error: the correct syntax is\n");
	printf("     Write\n");
        printf("          to overwrite traces in memory\n");
	printf("     Write file1 ... fileN\n");
        printf("          to create new files for the N traces in memory\n");
	printf("     Write Append Suffix \n");
	printf("     Write Prepend Prefix \n");
	printf("     Write Append Suffix Prepend Prefix \n");
}
