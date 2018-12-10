#include        <stdio.h>
#include        "gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"



void checkalloc(char **cmdstr, char *tcmdstr);
static char tmpstr[100];
int desac(int mcmd, char **tcmdstr, int *ncmd, char *cmdstr[5])
{

/* the objective is to convert the SAC abbreviated syntax for xlim and
 * cut to the complete syntax used by GSAC 
 *
 * To do this, we first define the current syntax
 * Then we construct the correct sequence
 * */

int i,j,m;
int ival;
float v;
int incmd;
int nwind;	/* must output two time windows */

/* FOOL PROOF CHECK 
	We can only have the following syntax for CUT of XLIM

 3:	CUT 0 20
 4:	CUT a 10 t0 
 5:	CUT a 10 t0 20
 8:	CUT a GMT 2004 123 12 22 33 444
10:	CUT a 20 GMT 2004 123 12 22 33 444
10:	CUT a CAL 2004 05 02 12 22 33 444 
11:	CUT a 20 CAL 2004 05 02 12 22 33 444 
15:	CUT GMT 2004 123 12 22 33 444 GMT 2004 123 13 22 33 444
16:	CUT CAL 2004 05 02 12 22 33 444 GMT 2004 123 13 22 33 444
17:	CUT CAL 2004 05 02 12 22 33 444 CAL 2004 05 02 13 22 33 444
*/
	if(mcmd > 17)
		return (NO);

	incmd = NO;
	nwind = 0;
	/* CHECK THIS IT WAS INCORRECTLY
	for (i=0, j=0; i < mcmd ; j < 17  ){
	*/
	for (i=0, j=0; i < mcmd && j < 17 ; ){
		if(i == 0){
			/* get command string */
			checkalloc(&cmdstr[j], tcmdstr[i]);
			strcpy(cmdstr[j++],tcmdstr[i]);
			i++;
		} else {
			/* test to see if we see GMT or CAL */
			strcpy(tmpstr, tcmdstr[i]);
			gsac_strupr(tmpstr);
			if(strcmp(tmpstr,"GMT") == 0){
					if(incmd == YES){
					checkalloc(&cmdstr[j], "0.0");
					strcpy(cmdstr[j++],"0.0");
					nwind++;
					}
					incmd = NO;
				for(m=0;m<7 && i < mcmd ; m++){
					if(m> 0){
						if(isargi(tcmdstr[i],&ival) 
							!= YES){
							return(NO);
						}
					}
					checkalloc(&cmdstr[j], tcmdstr[i]);
					strcpy(cmdstr[j++],tcmdstr[i]);
					i++ ;
					if(i > mcmd)return(NO);
				}
				nwind++;
			} else if(strcmp(tmpstr,"CAL") == 0){
					if(incmd == YES){
					checkalloc(&cmdstr[j], "0.0");
					strcpy(cmdstr[j++],"0.0");
					nwind++;
					}
					incmd = NO;
				for(m=0;m<8 && i < mcmd ; m++){
					if(m> 0){
						if(isargi(tcmdstr[i],&ival) 
							!= YES){
							return(NO);
						}
					}
					checkalloc(&cmdstr[j], tcmdstr[i]);
					strcpy(cmdstr[j++],tcmdstr[i]);
					i++ ;
					if(i > mcmd)return(NO);
				}
				nwind++;
			} else if(strncmp(tmpstr,"ON",2) == 0){
					checkalloc(&cmdstr[j], tcmdstr[i]);
					strcpy(cmdstr[j++],tcmdstr[i]);
					i++ ;
			} else if(strncmp(tmpstr,"OFF",2) == 0){
					checkalloc(&cmdstr[j], tcmdstr[i]);
					strcpy(cmdstr[j++],tcmdstr[i]);
					i++ ;
				/* abbreviated syntax now check for the */
			} else {
				if(isargr(tcmdstr[i],&v) == NO){
					if(incmd == YES){
					checkalloc(&cmdstr[j], "0.0");
					strcpy(cmdstr[j++],"0.0");
					nwind++;
					}
					incmd = YES;
					checkalloc(&cmdstr[j], tcmdstr[i]);
					strcpy(cmdstr[j++],tcmdstr[i]);
					i++;
					if(i > mcmd)return(NO);
				} else {
					if(incmd == NO){
						checkalloc(&cmdstr[j], "B");
						strcpy(cmdstr[j++],"B");
					}
					checkalloc(&cmdstr[j], tcmdstr[i]);
					strcpy(cmdstr[j++],tcmdstr[i]);
					i++;
					if(i > mcmd)return(NO);
					incmd = NO ;
					nwind++;
				}
			}
	
		}
	}
	if(incmd ==YES){
		checkalloc(&cmdstr[j], "0.0");
		strcpy(cmdstr[j++],"0.0");
		nwind++;
	}
	if(nwind < 2 && nwind > 0){
		if(j<17){
			checkalloc(&cmdstr[j], "B");
			strcpy(cmdstr[j++],"B");
		}
		if(j<17){
			checkalloc(&cmdstr[j], "0.0");
			strcpy(cmdstr[j++],"0.0");
		}
		nwind++;
	}
	*ncmd = j;
	return(YES);
}

void checkalloc(char **cmdstr, char *tcmdstr)
{
	if(*cmdstr == (char *)NULL){
		*cmdstr = (char *)calloc( strlen(tcmdstr)+1, sizeof(char));
	} else {
		*cmdstr = (char *)realloc(*cmdstr,
					(strlen(tcmdstr)+1)*sizeof(char));
	}
}




