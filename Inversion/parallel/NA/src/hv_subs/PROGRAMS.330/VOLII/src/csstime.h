#ifndef _CSSTIME_H
#define _CSSTIME_H

#define ISLEAP(yr)	((!(yr % 4) && yr % 100) || !(yr % 400))

struct date_time{
	double epoch;
	long date;
	int year;
	int month;
	char mname[4];
	int day;
	long doy;
	int hour;
	int minute;
	float second;
};

void htoe(struct date_time *dt) ;	/* convert from human to epoch	*/
double dtoepoch(long date) ;		/* convert julian date to epoch	*/
void month_day(struct date_time *dt) ;	/* from epoch fill in monty/day	*/
void etoh(struct date_time *dt) ;	/* epoch to human		*/
void mdtodate(struct date_time *dt);	/* from epoch to YYYY DOY */
void timestr(struct date_time *dt, char *str) ;
					/* 1999 12 31 23:59:59.999	*/
void timeprintstr(struct date_time *dt,char *str) ;
					/* epoch jday mon 12,1999 23:59:59.999*/

#endif
