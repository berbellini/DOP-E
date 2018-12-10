#include	<stdio.h>
#include	<string.h>
#include	"csstim.h"
#include	<math.h>

/* CHANGES
	15 JAN 2005
		in etoh - time not correct for YEAR < 1970
		changed to for( *year = 1969 ; ; (*year)-- ){
		from    *year--
	10 AUG 2005 msec not initialized in etoh 
	02 JAN 2009 - added roundoff to compute msec in etoh()
*/

static int days_in_month[] = {31,28,31,30,31,30,31,31,30,31,30,31,31};

void htoe1(int nzyear, int nzjday, int nzhour, int nzmin, int nzsec, int nzmsec, double *epoch)	/* convert from human to epoch */
{
	*epoch = 
		dtoepoch(nzyear, nzjday) + 
		nzhour * 3600. + 
		nzmin * 60. +
		nzsec + 0.001*nzmsec;
}

void htoe2(int nzyear, int nzmonth, int nzday, int nzhour, int nzmin, int nzsec, int nzmsec, double *epoch)	/* convert from human to epoch */
{
	int nzjday;
	nzjday = mdtodate(nzyear,nzmonth, nzday);
	*epoch = 
		dtoepoch(nzyear, nzjday) + 
		nzhour * 3600. + 
		nzmin * 60. +
		nzsec + 0.001*nzmsec;
}



int mdtodate(int nzyear, int nzmonth, int nzday)
{
	int i,dim;
	int doy;
	doy = 0L;
	for( i = 0 ; i < nzmonth - 1 ; i++ ){
		dim = days_in_month[i];
		if( i == 1 && ISLEAP(nzyear) ) dim++;
		doy += (long)dim;
	}
	doy += nzday;
	return(doy);
}


static char *month_name[] =
{"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

#define mod(a,b)	(a) - ((int)((a)/(b))) * (b)
void etoh(double epoch,int *year,int *doy,int *month,int *day,int *hour,int *minute,int *second,int *msec)
{
	int diy;
	double secleft;

	*doy = (epoch / 86400.);
	secleft = mod(epoch,86400.0);
	*hour = *minute = *second = *msec = 0;

	if(secleft != 0.0) {		/* compute hours minutes seconds */
		if(secleft < 0) {	/* before 1970 */
			*doy = *doy - 1;		/* subtract a day */
			secleft += 86400;	/* add a day */
		}
		*hour = secleft/3600;
		secleft = fmod(secleft,3600.0);
		*minute = secleft/60;
		*second = fmod(secleft,60.0);
		secleft = fmod(secleft,60.0) - (*second);
		*msec = 1000*(secleft+0.00049);
	}


	if(*doy >= 0){
		for( *year = 1970 ; ; (*year)++ ){
			diy = ISLEAP(*year) ? 366:365;
			if( *doy < (long)diy ) break;
			(*doy) -= (long)diy;
		}
	}
	else{
		for( *year = 1969 ; ; (*year)-- ){
			diy = ISLEAP(*year) ? 366:365;
			*doy += (long)diy;
			if( *doy >= 0L ) break;
		}
	}
	*doy = *doy + 1L;
	month_day(*year, *doy, month, day);
}

void month_day(int year, int doy, int *month, int *day)
{
	int i,dim,leap;

	leap = ISLEAP(year);
	*day = doy;
	for( i = 0 ; i < 12 ; i ++ ){
		dim = days_in_month[i];
		if( leap && i == 1 ) dim++;
		if( *day <= dim ) break;
		*day -= dim;
	}
	*month = i + 1;
}


/*
 * convert julian date to epoch time
 */
double dtoepoch(int year, int doy)
{
	register long	cnt;
	long	days;

	cnt = (long)year;
	days = 0L;
	if (cnt > 1970L)
		while (--cnt >= 1970L)
			days += ISLEAP(cnt) ? 366 : 365;
	else if (cnt < 1970L)
		while (cnt < 1970L) {
			days -= ISLEAP(cnt) ? 366 : 365;
			cnt++;
		}
	return( (days + ((doy - 1) )) * 86400.);
}

void setepoch(int year, int doy, int hour, int minute, int second, int msec, double *epoch)
{
	/* first test */
	*epoch = 0;
	if(year == -12345 || doy == -12345 ||
		hour == -12345 || minute == -12345 ||
		second == -12345 || msec == -12345) return;
	/* now convert to epoch */
	htoe1(year, doy, hour, minute, second, msec, epoch);
}



void printkzdatestr(int year, int doy, char *str)
{
int month, day;
month_day(year, doy, &month, &day);

if(year != -12345 && doy != -12345)
sprintf(str,"KZDATE   %3s %2.2d (%3.3d), %4.4d",month_name[month -1], day, doy, year);
else
	str[0]='\0';
}

void printkztimestr(int hour, int minute, int second, int msec, char *str)
{
if(hour != -12345 && minute != -12345 && second != -12345 && msec != -12345)
sprintf(str,"KZTIME         %2.2d:%2.2d:%2.2d.%3.3d",hour,minute,second,msec);
else
	str[0]='\0';
}
void printkdatestr(int year, int doy, char *str)
{
int month, day;
month_day(year, doy, &month, &day);

if(year != -12345 && doy != -12345)
sprintf(str,"%3s %2.2d (%3.3d), %4.4d",month_name[month -1], day, doy, year);
else
	str[0]='\0';
}

void printktimestr(int hour, int minute, int second, int msec, char *str)
{
if(hour != -12345 && minute != -12345 && second != -12345 && msec != -12345)
sprintf(str,"%2.2d:%2.2d:%2.2d.%3.3d",hour,minute,second,msec);
else
	str[0]='\0';
}

void printtimestr(double epoch, char *str)
{
int year, doy, month, day, hour, minute,second,msec;
etoh( epoch, &year, &doy, &month, &day, &hour, &minute, &second, &msec);
sprintf(str,"%3s %2.2d (%3.3d), %4.4d %2.2d:%2.2d:%2.2d.%3.3d",month_name[month -1], day, doy, year,hour,minute,second,msec);
}
