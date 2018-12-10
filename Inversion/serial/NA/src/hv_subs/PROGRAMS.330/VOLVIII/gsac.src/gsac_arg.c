/* test strings to determine if they are compatible with
 * expected type
 */
/*
	isargr    - is argument real - YES/NO
	isargi    - is argument int  - YES/NO
	isargl    - is argument logical (T/F)  - YES/NO
	argerrr   - print command error with offending syntax
	testarg   - look at all arguments to see if they fit the grammar
	getargr   - get the real argument(s)
	getargi   - get the int  argument(s)
	getargl   - get the logical  argument
	getargs   - get the string   argument(s)
	getargs2  - get the two string   argument(s)
	getargyn  - get the look for ON or OFF or YES or NO 
	getargn   - is the argument OFF or NO
*/

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	"gsac.h"
#include	"gsac_arg.h"

/* clean this up later */
char tmpstr[1000];
static int mystrcmp(char *argstr, char *tmpstr, int mfit);


/* Is the string a FLOAT   - if YES the float is returned as rval */
int isargr(char *arg, float *rval)
{
	char *v;
	*rval = strtod(arg, &v);
	if(*v == '\0')
		return YES;
	else
		return NO;
}

/* Is the string a INT   - if YES the integer is returned as rval */
int isargi(char *arg, int   *ival)
{
	char *v;
	*ival = strtol(arg, &v, 10);
	if(*v == '\0')
		return YES;
	else
		return NO;
}

/* Is the string a YES or NO or TRUE or FALSE  - if YES the YES is returned  */
int isargl(char *arg, int   *lval)
{
	if(arg[0] == 'T' || arg[0] == 't' ){
		*lval = YES;
		return YES;
	} else if(arg[0] == 'F' || arg[0] == 'f'){
		*lval = YES;
		return YES;
	} else
		return NO;
}

void argerror(int ncmd, char **cmdstr, char *mesg, int pos);
void outchar(int n, int c);

/* this is the first important test of the command string
 * before the command string is actually taken apart we determine
 * whether the ENTIRE string is valid 
 *
 * consider the entire argument string, looking for inconsistencies 
 * return (0) is all is OK
 * return (N > 0 ) for offending systax. Note that the N = 0 is just
 * the command which we know exists
 * keepused - for listheader we may only display a subset - if we reset
 * the show that we lose this memory of the subset. In other cases 
 * the subset is preserved by external variables and we can always start with a clean slate
 *
 * ncmd	   - number of entries on the command line , for example
 * 	GSAC> bandpass           has ncmd = 1
 * 	GSAC> bandpass c 0.1 2.0  has ncmd = 4
 * cmdstr  - actual strings on the command line. For the second case these are
 * 		"bandpass", "c", "0.1" and "2.0"
 * tstarg  - this is the list of acceptable command arguments -  this list
 * 		defines the type and number of subarguments. The definition
 * 		for corner frequency states that TWO real numbers are expected
 * keepused - this is part of a command memory mechanism when the subcommand
 * 		does not have any arguments. If YES previous use is preserved
 * 16 AUG 2007 - add testfoundit so that we test syntax and ignore the rest, 
 *     such as write append suffice file1 file2
 */

int testarg(int ncmd, char **cmdstr, struct arghdr *tstarg, int keepused, int testfoundit)
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
		if(foundit == NO && testfoundit){
			argerror(ncmd, cmdstr, "incorrect option",i);
			return (i);
		}
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

/* if there is an error in command syntax, stop processing but also 
 * indicate where the error occurs */
void argerror(int ncmd, char **cmdstr, char *mesg, int pos)
{
	int i;
	printf("Error in %s: %s\n",cmdstr[0],mesg);
	outchar(4,' ');
	for(i=0 ; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	outchar(4,' ');
	for (i=0 ; i < pos  ; i++){
		outchar(strlen(cmdstr[i])+1,' ');
	}
	outchar(strlen(cmdstr[pos]),'^');
	printf("\n");
}

/* simple character output */
void outchar(int n, int c)
{
        int i;
        for(i=0 ; i < n ; i++)
                printf("%c",c);
}

/* get real or float argrments 
 * This is used after the entire command line is shown to have correct syntax
 * Now we searth for a pattern and then */
int getargr(int ncmd, char **cmdstr,  char *argstr, int mfit, int narg, float *rarr)
{
	int i, j;
	float rval;
	for(i=1; i < ncmd ; i++) {
		strcpy(tmpstr,cmdstr[i]);
		gsac_strupr(tmpstr);
		if(mystrcmp(argstr,tmpstr,mfit) ==0){
			for(j=0 ; j < narg ; j++){
				i++;
				isargr(cmdstr[i],&rval);
				rarr[j] = rval;
			}
			return i;
		}
	}
	return -1;
}

/* get an integer */
int getargi(int ncmd, char **cmdstr,  char *argstr, int mfit, int narg, int *iarr)
{
	int i, j;
	int ival;
	for(i=1; i < ncmd ; i++) {
		strcpy(tmpstr,cmdstr[i]);
		gsac_strupr(tmpstr);
		if(mystrcmp(argstr,tmpstr,mfit) ==0){
			for(j=0 ; j < narg ; j++){
				i++;
				isargi(cmdstr[i],&ival);
				iarr[j] = ival;
			}
			return i;
		}
	}
	return -1;
}

/* get logical TRUE or FALSE arguments */
int getargl(int ncmd, char **cmdstr, char *argstr, int mfit, int narg, int   *lval)
{
	/* we look at the first character of the argument which is
	 * either T or F the value is returned as 
	 * 1 for T, 0 for F  */
	int i, j;
	char c;
	for(i=1; i < ncmd ; i++) {
		strcpy(tmpstr,cmdstr[i]);
		gsac_strupr(tmpstr);
		if(mystrcmp(argstr,tmpstr,mfit) ==0){
			for(j=0 ; j < narg ; j++){
				i++;
				c = cmdstr[i][0] ;
				if(c == 'T' || c == 't' )
					*lval = 1;
				else if(c == 'F' || c == 'f')
					*lval = 0;
				else 
					*lval = -1;
			}
			return i;
		}
	}
	return -1;
}

/* get a string */
int getargs(int ncmd, char **cmdstr, char *argstr, int mfit, int narg, char *sval)
{
	int i;
	for(i=1; i < ncmd ; i++) {
		strcpy(tmpstr,cmdstr[i]);
		gsac_strupr(tmpstr);
		if(mystrcmp(argstr,tmpstr,mfit) ==0){
				i++;
				strcpy(sval, cmdstr[i]);
			return i;
		}
	}
	return -1;
}

/* get two strings */
int getargs2(int ncmd, char **cmdstr, char *argstr, int mfit, int narg, char *sval, char *tval)
{
	int i;
	for(i=1; i < ncmd ; i++) {
		strcpy(tmpstr,cmdstr[i]);
		gsac_strupr(tmpstr);
		if(mystrcmp(argstr,tmpstr,mfit) ==0){
				i++;
				strcpy(sval, cmdstr[i]);
				i++;
				strcpy(tval, cmdstr[i]);
			return i;
		}
	}
	return -1;
}

int getargyn(int ncmd, char **cmdstr, char *argstr, int mfit, int narg, int   *lval)
{
	/* look for a ON or OFF */
	int i;
	for(i=1; i < ncmd ; i++) {
		strcpy(tmpstr,cmdstr[i]);
		gsac_strupr(tmpstr);
		if(mystrcmp(argstr,tmpstr,mfit) ==0){
			i++;
			strcpy(tmpstr,cmdstr[i]);
			gsac_strupr(tmpstr);
			if(strcmp(tmpstr,"ON") == 0 || strcmp(tmpstr,"YES") == 0 )
				*lval = YES;
			else if(strcmp(tmpstr, "OFF") == 0 || strcmp(tmpstr, "NO") )
				*lval = NO;
			else
				*lval = -1;
			return i;
		}
	}
	return -1;
}


/* get ON OFF or INTEGER */
int getargn(int ncmd, char **cmdstr, char *argstr, int mfit, int narg, int   *lval)
{
	/* look for a ON or OFF */
	int i;
	for(i=1; i < ncmd ; i++) {
		strcpy(tmpstr,cmdstr[i]);
		gsac_strupr(tmpstr);
		if(mystrcmp(argstr,tmpstr,mfit) ==0){
			i++;
			strcpy(tmpstr,cmdstr[i]);
			gsac_strupr(tmpstr);
			if(strcmp(tmpstr,"OFF") == 0){
				*lval = 0;
			} else {
				if(isargi(tmpstr ,   lval)==NO)
					*lval = 0;
			}
			return i;
		}
	}
	return -1;
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

/* identify the postion of argstr in cmdstr 
	return -1 if not found */
int whicharg(int ncmd, char **cmdstr, char *argstr)
{
	int i;
	for(i = 0 ; i < ncmd; i++){
		strcpy(tmpstr,cmdstr[i]);
		gsac_strupr(tmpstr);
		if(strcmp(tmpstr,argstr)==0)
			return i;
	}
	return -1;
}


