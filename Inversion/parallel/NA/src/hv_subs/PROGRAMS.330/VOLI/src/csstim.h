#ifndef _CSSTIM_H
#define _CSSTIM_H

#define ISLEAP(yr)	((!((yr) % 4) && (yr) % 100) || !((yr) % 400))

struct date_time {
	double epoch;
	int date;
	int year;
	int month;
	char mname[4];
	int day;
	long doy;
	int hour;
	int minute;
	int second;
	int msec;
};



void htoe1(int nzyear, int nzjday, int nzhour, int nzmin, int nzsec, int nzmsec, double *epoch);	/* convert from human to epoch YYYYDDDHHMMSS.sss */
void htoe2(int nzyear, int nzmonth, int nzday, int nzhour, int nzmin, int nzsec, int nzmsec, double *epoch);	/* convert from human to epoch YYYYMMDDHHMMSS.sss  */
void etoh(double epoch,int *year,int *doy,int *month,int *day,int *hour,int *minute,int *second,int *msec);
int mdtodate(int nzyear, int nzmonth, int nzday);
double dtoepoch(int year, int doy);
void month_day(int year, int doy, int *month, int *day);
void printkzdatestr(int year, int doy, char *str);
void printkztimestr(int hour, int minute, int second, int msec, char *str);
void printkdatestr(int year, int doy, char *str);
void printktimestr(int hour, int minute, int second, int msec, char *str);
void setepoch(int year, int doy, int hour, int minute, int second, int msec, double *epoch);
void printtimestr(double epoch, char *str);


#endif
