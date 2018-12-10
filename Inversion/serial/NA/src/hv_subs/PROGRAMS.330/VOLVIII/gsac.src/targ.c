#include <stdio.h>
#include <string.h>
#include "gsac_arg.h"

#define NO  0
#define YES 1

#define WRITE_DFLT      0
#define WRITE_APPEND    1
#define WRITE_PREPEND   2
#define WRITE_NUMBER   3

static char  write_append[80];
static char  write_prepend[80];
static int   write_doappend = NO ;
static int   write_doprepend= NO ;
static int   write_ntrace = 0;

int ttestarg(int ncmd, char **cmdstr, struct arghdr *tstarg, int keepused);
char *gsac_strupr(char *s);
static int mystrcmp(char *argstr, char *tmpstr, int mfit);



struct arghdr writearg[] = {
        {WRITE_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
        {WRITE_APPEND   , "APPEND"  , CHDR, 0, 1, NO, "APPEND string ", -1},
        {WRITE_PREPEND  , "PREPEND" , CHDR, 0, 1, NO, "PREPEND string ", -1},
        {WRITE_NUMBER  , "NUMBER" , IHDR, 0, 1, NO, "APPEND string ", 1},
        {0,     ""              , IHDR, NO, 0, NO, "", -1}
};

static char tmpstr[1000];


char *cmdstr[] = { "write" , "append" , "suf", "file2", "file1", "prepend", "prefix"};
int ncmd = 7;




main()
{
	int i;
	int ret;
	for(i=0;i<ncmd;i++)
		printf("%d: %s\n",i,cmdstr[i]);
	ret = ttestarg(ncmd, cmdstr, writearg, NO );
	printf("ret %d\n",ret);
}
int ttestarg(int ncmd, char **cmdstr, struct arghdr *tstarg, int keepused)
{
	int i, j, foundit, jval;
	int jcnt;
	int ival;
	float rval;
		/* no arguments means a perfect syntax and usually to
		 * repeat the last use of this command */
	if(ncmd == 1)
		return (0);
	if(keepused != YES){
		for(j=0; tstarg[j].key[0] != '\0'; j++){
			tstarg[j].used = NO;
		}
	}
	for(i=1; i < ncmd ; i++) {
		strcpy(tmpstr,cmdstr[i]);
		gsac_strupr(tmpstr);
		/* scroll through the option list */
		/* also reset the .used field */
		foundit = NO;
		for(j=0; foundit == NO && tstarg[j].key[0] != '\0'; j++){
			jval = j;
			/* check for abbreviated definition */
			if(mystrcmp(tmpstr,tstarg[j].key, 
				tstarg[j].mfit) == 0)
					foundit = YES;
		}
/*
		if(foundit == NO){
			argerror(ncmd, cmdstr, "incorrect option",i);
			return (i);
		}
*/
		/* kludge */
if(foundit == YES){
		if(tstarg[jval].narg == 0){
			tstarg[jval].used += YES;
		} else {
			for(jcnt = 0 ; jcnt < tstarg[jval].narg ; jcnt++){
				if(i < ncmd -1 )
					i++;
				else {
					argerror(ncmd,cmdstr," ran out of input",i);
					return(i);
				}
				/* for each header type determine if the syntax is
				 * correct and then note that this value can be set */
				if(tstarg[jval].ricell == IHDR){
					if(isargi(cmdstr[i],&ival)  != 1){
						argerror(ncmd, cmdstr, tstarg[jval].errormessage,i);
/*
						printf("Looking for IHDR jval %d i %d cmdstr %s narg %d\n",jval,i,cmdstr[i],tstarg[jval].narg);
*/
						return (i);
					} else {
						tstarg[jval].used += YES;
					}
				} else if(tstarg[jval].ricell == LHDR){
					if(isargl(cmdstr[i],&ival)  != 1){
						argerror(ncmd, cmdstr, tstarg[jval].errormessage,i);
/*
						printf("Looking for LHDR jval %d i %d cmdstr %s narg %d\n",jval,i,cmdstr[i],tstarg[jval].narg);
*/
						return (i);
					} else {
						tstarg[jval].used += YES;
					}
				} else if(tstarg[jval].ricell == CHDR){
						/* by definition we always 
						 * have a string */
						tstarg[jval].used += YES;
/*
printf("testarg %d CHDR\n",jval);
*/
				} else if(tstarg[jval].ricell == CHDL){
						/* by definition we always 
						 * have a string */
						tstarg[jval].used += YES;
/*
printf("testarg %d CHDL\n",jval);
*/
				} else if(tstarg[jval].ricell == RHDR){
					if(isargr(cmdstr[i],&rval) != 1 ){
						argerror(ncmd, cmdstr, tstarg[jval].errormessage,i);
/*
						printf("Looking for RHDR jval %d i %d cmdstr %s narg %d %d\n",jval,i,cmdstr[i],tstarg[jval].narg,tstarg[jval].used);
*/
						return (i);
					} else {
						tstarg[jval].used += YES;
					}
				} else if(tstarg[jval].ricell == YHDR){
						tstarg[jval].used += YES;
				} else if(tstarg[jval].ricell == NHDR){
						tstarg[jval].used += YES;
				}
			}
		}
	}
	}
	return(0);
}

char *gsac_strupr(char *s)
{
        if (s != NULL )
        {
                char *p;

                for( p=s ; *p; p++)
                        *p = toupper(*p);
        }
        return s;
}

/* version of strcmp and strcmp that uses the minimum fit
 * parameter to define the minumum length fo the string to be fit
 * */
static int mystrcmp(char *argstr, char *tmpstr, int mfit)
{
        if(mfit > 0){
                return(strncmp(argstr,tmpstr, mfit)) ;
        } else {
                /* check for complete fit */
                return(strcmp(argstr,tmpstr));
        }
}

