#ifndef _GSACARGHDR
#define _GSACARGHDR


#define XHDR	-1
#define	RHDR	0	/* floats */
#define IHDR	1	/* integers */
#define CHDR	2	/* 8 character string */
#define EHDR	3	/* extended */
#define LHDR	4	/* logical - true or false */
#define YHDR	5	/* ON or OFF */
#define NHDR	6	/* NO or n as in perplot 1 of perplot OFF */
#define CHDL	7	/* 16 character string , only KEVNM */

struct arghdr {
	int id;
	char *key;
	int ricell;
	int used;	/* originally yes or no, now can be > 1, e.g.,
				cut o -10 o 150
			*/
	int narg;
	int show;
	char *errormessage;
	int mfit;	/* mimimum part of key to be fit */
};

int isargi(char *arg, int   *ival);
int isargl(char *arg, int   *lval);
int isargr(char *arg, float *rval);
int isargs(char *arg);

/* we return the index of the last part of the command string that
 * was parsed. this is because we must consider the possibility of
 * cut o -10 o 100 in addition to the simpler cut a -10 t0 +20 
 * */
int testarg(int ncmd, char **cmdstr, struct arghdr *tstarg ,int keepused, int testfoundit);
int getargr(int ncmd, char **cmdstr, char *argstr, int mfit, int narg, float *rarr);
int getargi(int ncmd, char **cmdstr, char *argstr, int mfit, int narg, int   *iarr);
int getargl(int ncmd, char **cmdstr, char *argstr, int mfit, int narg, int   *lval);
int getargs(int ncmd, char **cmdstr, char *argstr, int mfit, int narg, char   *sval);
int getargs2(int ncmd, char **cmdstr, char *argstr, int mfit, int narg, char   *sval, char *tval);
int getargyn(int ncmd, char **cmdstr, char *argstr, int mfit, int narg, int   *lval);
int getargn(int ncmd, char **cmdstr, char *argstr, int mfit, int narg, int   *lval);
int whicharg(int ncmd, char **cmdstr, char *argstr);
/* added 21 JUL 2009 for use in gsac_write - test arguments but do nto complain */
int wtestarg(int ncmd, char **cmdstr, struct arghdr *tstarg ,int keepused, int testfoundit);

#endif
