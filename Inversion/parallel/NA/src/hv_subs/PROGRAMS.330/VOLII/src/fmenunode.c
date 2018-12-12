#include "nfmenu.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
fmenu *start, *end;

/* linked list code
	Ammeraal, L,
	Programs and Data Structures in C, Second Edition,
		Based on ANSI C and C++,
	John Wiley & Sons, Chichester
	1992.
	272 pp

	pages 111-112
*/

void *getmemory(int n)
{
	void *p=malloc(n);
	if(p == NULL){
		fprintf(stderr,"Not enough memory\n");
		exit(1);
	}
	return p;
}

fmenu *getnode(void)
{
	return (fmenu *)getmemory(sizeof (fmenu));
}

void deletenode(fmenu *p)
{
	fmenu *q;
	free(p->str);
	q = p->next;
	if(q == end)
		end = p;
	else
		*p = *q;
	free(q);
}

void insertnode(fmenu *p, float xl, float yl, float xh, float yh,
	char *str, int action, int lstrmx, int type, int line, int fsize,
	int nsamp, char *kstnm, char *kcmpnm, char *datetime,
	int page, int used, float dist, float az, float baz)
{
	fmenu *q;
	q = getnode();
	if(p == end)
		end = q;
	else
		*q = *p;
	p->next = q;
	p->xl = xl;
	p->yl = yl;
	p->xh = xh;
	p->yh = yh;
	p->str = (char *)getmemory(strlen(str)+1);
	strcpy(p->str, str);
	p->action = action;
	p->lstrmx = lstrmx;
	p->type   = type  ;
	p->line   = line  ;
	p->fsize  = fsize ;
	strcpy(p->kstnm, kstnm);
	strcpy(p->kcmpnm, kcmpnm);
	strcpy(p->datetime, datetime);
	p->nsamp  = nsamp ;
	p->page   = page  ;
	p->used   = used  ;
	p->dist = dist;
	p->az = az;
	p->baz = baz;
}


void appendnode(float xl, float yl, float xh, float yh,
	char *str, int action, int lstrmx, int type, int line, int fsize,
	int nsamp, char *kstnm, char *kcmpnm, char *datetime,
	int page, int used, float dist, float az, float baz)
{
	fmenu *p=end;
	end = getnode();
	p->next = end;
	p->xl = xl;
	p->yl = yl;
	p->xh = xh;
	p->yh = yh;
	p->str = (char *)getmemory(strlen(str)+1);
	strcpy(p->str, str);
	p->action = action;
	p->lstrmx = lstrmx;
	p->type   = type  ;
	p->line   = line  ;
	p->fsize  = fsize ;
	p->nsamp  = nsamp ;
	strcpy(p->kstnm, kstnm);
	strcpy(p->kcmpnm, kcmpnm);
	strcpy(p->datetime, datetime);
	p->page  = page ;
	p->used   = used  ;
	p->dist = dist;
	p->az = az;
	p->baz = baz;
}



