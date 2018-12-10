/* changes 
 * 01 APR 2004 - forced wsac1 to set depmax depmin depmen
 */
#ifndef _GSACHDR
#define _GSACHDR

#define XHDR	-1
#define	RHDR	0	/* floats */
#define IHDR	1	/* integers */
#define CHDR	2	/* 8 character string */
#define EHDR	3	/* extended */
#define LHDR	4	/* logical - true or false */
#define YHDR	5	/* ON or OFF */
#define NHDR	6	/* NO or n as in perplot 1 of perplot OFF */
#define CHDL	7	/* 16 character string , only KEVNM */

struct _lhdr {
	char key[9];
	int  ricel ;
	int  pos   ;
	int  show  ;
	int  change;
} ;



/* special keywords for gsac_set_param_lh */

#define H_NOP  0
#define H_REAL 1
#define H_INT  2
#define H_STR  3
#define H_LOG  4

struct _lh_param   {
	char *cmd ;
	int  type;
	int  narg;
	};


#endif
