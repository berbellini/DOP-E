
#include	<stdio.h>
#include	"sisdef.h"

/* UNIX (Sun ) SIS FORMAT on CRAY this is DISK FILE FORMAT */
#define	SIZEOFCHAR	1
#define	SIZEOFREAL	4
#define	SIZEOFINT	4
#define	SZTRHD		256


extern int	xargc;
extern char	**xargv;
extern	int iform;


double atof();
int atoi();
int strlen();
int strncmp();
char *strcpy();

argi4(key,value,defval,noval)
char *key;
int *value;
int defval, noval;
{
	int i,keylen;
	char *cp;
	char **bp;
	keylen = strlen(key);
	bp = xargv;
	for(i=0;i<xargc-1;i++){
		cp = *++bp;
		if(strncmp(key,cp,keylen) == 0){
			while(keylen > 0){
				cp++;
				keylen--;
			}
			if(strlen(cp) > 0){
				*value = atoi(cp);
				return(atoi(cp));
			}
			else {
				*value = defval;
				return(defval);
			}
		}
	}
	*value = noval;
	return(noval);
}

double argr4(key,value,defval,noval)
char *key;
double *value;
double defval, noval;
{
	int i,keylen;
	char *cp;
	char **bp;
	keylen = strlen(key);
	bp = xargv;
	for(i=0;i<xargc-1;i++){
		cp = *++bp;
		if(strncmp(key,cp,keylen) == 0){
			while(keylen > 0){
				cp++;
				keylen--;
			}
			if(strlen(cp) > 0){
				*value = atof(cp);
				return(atof(cp));
			}
			else {
				*value = defval;
				return(defval);
			}
		}
	}
	*value = noval;
	return(noval);
}

int argis(key)
char *key;
{
	int i,keylen;
	char *cp;
	char **bp;
	keylen = strlen(key);
	bp = xargv;
	for(i=0;i<xargc-1;i++){
		cp = *++bp;
		if(strncmp(key,cp,keylen) == 0)
			return(1);
	}
	return(0);
}

int argstr(key,value,defval,noval)
char *key;
char *value;
char *defval, *noval;
{
	int i,keylen;
	char *cp;
	char **bp;
	keylen = strlen(key);
	bp = xargv;
	for(i=0;i<xargc-1;i++){
		cp = *++bp;
		if(strncmp(key,cp,keylen) == 0){
			while(keylen > 0){
				cp++;
				keylen--;
			}
			if(strlen(cp) > 0){
				strcpy(value,cp);
				return(1);
			}
			else {
				strcpy(value,defval);
				return(1);
			}
		}
	}
	strcpy(value,noval);
	return(1);
}

getln(luin,ntap,mode,defval)
int *luin;
char *ntap;
char *mode;
int defval;
{
	if(strcmp(ntap," ") != 0)
		lbopen(luin,ntap,mode);
	else
		*luin = defval;
}

cmdchk(ns,ne,rs,re,ntrc,nrec)
int 	*ns,
	*ne,
	*rs,
	*re,
	ntrc,
	nrec;
{
	if(*ns <=0)*ns = 1;
	if(*ne <=0)*ne = ntrc;
	if(*ne > ntrc)*ne = ntrc;
	if(*ns > *ne)*ns = *ne;
	if(*rs <= 0)*rs = 1;
	if(*re <= 0)*re = nrec;
	if(*re > nrec)*re = nrec;
	if(*rs > *re)*rs = *re;
}

/* for the corrent record, skip traces ns to ne */
skptrc(rec,ns,ne,luin,nsamp,ntrc,itr,lbytes,nbytes,iform)
int 	rec,
	ns,
	ne,
	luin,
	nsamp,
	ntrc,
	itr[],
	lbytes,
	*nbytes,
	iform;
{
	int lbyts4, nbyts4, ntrskp, nbyskp, k;
	int szsamp;
	if(luin == 0){
		if(ns > ne)return;
		for(k=ns; k <= ne ; k++)
			rtape(luin,itr,nbytes);
	}
	else {
		if(iform == 3)
			szsamp = SIZEOFREAL;
		else if(iform == 1000)
			szsamp = SIZEOFCHAR;
		*nbytes = SZTRHD + szsamp  * nsamp;
		if(rec < 0)rec = 1;
		lbyts4 = lbytes + SIZEOFINT;
		nbyts4 = *nbytes + SIZEOFINT;
		if(ne < 0)ne = 0;
		if(ne > ntrc)ne = ntrc;
		ntrskp = ne ;
		nbyskp = lbyts4 + ((rec -1 )*ntrc+ntrskp)*
			nbyts4;
		seekt(luin,nbyskp);
	}
}

/* skip to end of record re */
/* this is only used for sequential reads, skptrc takes care of 
   direct access for file input	*/
skprec(rs,re,luin,ntrc,itr)
int	rs,
	re,
	luin,
	ntrc,
	itr[];
{
	int jj, kk, nbytes;
	if(luin == 0){
		if(rs > re)return(0);
		for(jj=rs;jj <= re;jj++){
			for(kk=0;kk<ntrc;kk++){
				nbytes = 0;
				rtape(luin,itr,&nbytes);
			}
		}
	}
}
			
			

double fabs();

maxsn(lx,x,xmax,jndex)
int lx;
float x[];
float *xmax;
int *jndex;
{
	int i;
	*jndex = 1;
	for(i=0;i < lx; i++)
		if(fabs(x[*jndex]) < fabs(x[i]) ) *jndex = i;
	*xmax = fabs(x[*jndex]);
}



help(hlp)
struct help *hlp;
{
	struct help *hp;
	for(hp = hlp; hp->item != (char *)0;hp++)
		fprintf(stderr,"%s\n",hp->item);
	exit(0);
}
