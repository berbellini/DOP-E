#include <stdio.h>

#include	<string.h>
#include "csstime.h"

static int days_in_month[] = {31,28,31,30,31,30,31,31,30,31,30,31,31};

void htoe(struct date_time *dt)	/* convert from human to epoch */
{
	double dtoepoch();
	dt->epoch = 
		dtoepoch(dt->date) + 
		dt->hour * 3600. + 
		dt->minute * 60. +
		dt->second;
}

void timeprintstr(struct date_time *dt,char *str)
{
	sprintf(str,"%8ld %s %2d,%4d %2d:%02d:%06.3f",
	dt->date,
	dt->mname,
	dt->day,
	dt->year,
	dt->hour,
	dt->minute,
	dt->second);
}

void timestr(struct date_time *dt, char *str)
{
	sprintf(str,"%04d %02d %02d %02d:%02d:%06.3f",
	dt->year,
	dt->month,
	dt->day,
	dt->hour,
	dt->minute,
	dt->second);
	str[23]='\0';
}

void mdtodate(struct date_time *dt)
{
	int i,dim;
	dt->doy = 0L;
	for( i = 0 ; i < dt->month - 1 ; i++ ){
		dim = days_in_month[i];
		if( i == 1 && ISLEAP(dt->year) ) dim++;
		dt->doy += (long)dim;
	}
	dt->doy += (long)dt->day;
	dt->date = (1000L * dt->year) + dt->doy;
}


static char *month_name[] =
{"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

#define mod(a,b)	(a) - ((int)((a)/(b))) * (b)
void etoh(struct date_time *dt)
{
	int diy;
	double secleft;

	dt->doy = (long)(dt->epoch / 86400.);
	secleft = mod(dt->epoch,86400.0);
	dt->hour = dt->minute = dt->second = 0;

	if(secleft) {			/* compute hours minutes seconds */
		if(secleft < 0) {	/* before 1970 */
			dt->doy = dt->doy - 1L;		/* subtract a day */
			secleft += 86400;	/* add a day */
		}
		dt->hour = secleft/3600;
		secleft = mod(secleft,3600.0);
		dt->minute = secleft/60;
		dt->second = mod(secleft,60.0);
	}

	if(dt->doy >= 0L){
		for( dt->year = 1970 ; ; dt->year++ ){
			diy = ISLEAP(dt->year) ? 366:365;
			if( dt->doy < (long)diy ) break;
			dt->doy -= (long)diy;
		}
	}
	else{
		for( dt->year = 1969 ; ; dt->year-- ){
			diy = ISLEAP(dt->year) ? 366:365;
			dt->doy += (long)diy;
			if( dt->doy >= 0L ) break;
		}
	}
	dt->doy = dt->doy + 1L;
	dt->date = (long)dt->year * 1000L + (long)dt->doy;
	month_day(dt);
}
void month_day(struct date_time *dt)
{
	int i,dim,leap;

	leap = ISLEAP(dt->year);
	dt->day = dt->doy;
	for( i = 0 ; i < 12 ; i ++ ){
		dim = days_in_month[i];
		if( leap && i == 1 ) dim++;
		if( dt->day <= dim ) break;
		dt->day -= dim;
	}
	dt->month = i + 1;
	strcpy(dt->mname,month_name[i]);
}


/*
 * convert julian date to epoch time
 */
double dtoepoch(long date)
{
	register long	cnt;
	long	days;

	cnt = (long)(date / 1000);
	days = 0L;
	if (cnt > 1970L)
		while (--cnt >= 1970L)
			days += ISLEAP(cnt) ? 366 : 365;
	else if (cnt < 1970L)
		while (cnt < 1970L) {
			days -= ISLEAP(cnt) ? 366 : 365;
			cnt++;
		}
	return( (days + ((date - 1) % 1000)) * 86400.);
}
