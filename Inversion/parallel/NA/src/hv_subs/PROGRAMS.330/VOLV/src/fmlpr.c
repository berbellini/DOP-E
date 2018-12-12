
/* originally developed by c j ammon pennsylvania state university
   modified by RBHerrmann of Saint Louis University to make more compatible with
   Berkeley announcement 

   02 JAN 2008 - modified to permit more comments, e.g., filtering and 
                 stations used.  Also modified to permit inputting 
                 the moment tensor directly
   05 JAN 2009 - the Ammon lprmech is now called fmlpr, and creates two files:
		 a .msg file which is the printer plot beach ball, and the .ndk
		 which is the source information in the Global CMT ndk format:

http://www.seismology.harvard.edu/projects/CMT/catalog/allorder.ndk_explained

*/

/* includes */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "msg.h"
#include "ndk.h"

/* function prototypes */
void PrintPradiation(FILE *outStream, float Mij[3][3], float plnt,float stkt,float plnp,float stkp);
void phase(float x,float y,float z,float *del,float *bet);
void trans(float dip,float stk,float slip,float *stkt,float *plnt,float *stkp,float *plnp, float *plnn, float *stkn);
void tpdss(float stkt,float plnt,float stkp,float plnp,
	float *stk0, float *dip0, float *rak0,
	float *stk1, float *dip1, float *rak1) ;
void outmsg( int year, int month, int day, int hour, int minute, int second, struct msg_ *m, char *file1, char *file3);
void outndk( int year, int month, int day, int hour, int minute, int second, struct ndk_ *n);


main(int argc, char **argv)
{
	char title[132];
	char datetimeloc[132];
	char geog[132];
	char file1[132];
	char file2[132];
	char file3[132];
	float strike, dip, rake, Mw;
	float Mo;
	int i,j;

	int exponent;

	int year, month, day, hour, minute, second, millisecond;
	float evla, evlaerr, evlo, evloerr, evdp, evdperr;
	float depth;
	float mag1;

	struct msg_ m;
	struct ndk_ n;
	time_t secs;
	struct tm *gmt;

	FILE  *file2id;
	

	/* get input */
	if(argc != 12)
	{
		fprintf(stderr,
			"usage: fmlpr TITLE DATETIME GEOG EVDP MW strike dip rake FILE1 FILE2 FILE3 \n\n");
		fprintf(stderr,
			"  TITLE DATETIME and GEOG are enclosed in quotes\n");
		fprintf(stderr,
			"  FILE1 - name of file with stations used\n");
		fprintf(stderr,
			"  FILE2 - file with data summary: NSTA NCMP MINPER\n");
		fprintf(stderr,
			"  FILE3 - file with filtering commands\n");
		fprintf(stderr,
			"  If FILE1 FILE2 FILE3 do not exist indicate by pair of double quotes\n");
		exit(1);
	}
	strcpy(title,argv[1]);
	strcpy(datetimeloc,argv[2]);
	strcpy(geog,argv[3]);
	sscanf(argv[4],"%f",&evdp);
	sscanf(argv[5],"%f",&Mw);
	sscanf(argv[6],"%f",&strike);
	sscanf(argv[7],"%f",&dip);
	sscanf(argv[8],"%f",&rake);
	strcpy(file1,argv[9]);
	strcpy(file2,argv[10]);
	strcpy(file3,argv[11]);

	/* parse the datetimeloc */
	sscanf(datetimeloc,"%d %d %d %d %d %d %d %f %f %f %f",
		&year,&month,&day, &hour,&minute,&second,&millisecond,
		&evla, &evlo, &depth, &mag1);
	/* end get input */

	 Mo = 1.5*Mw + 16.10;
	 Mo = pow (10.0, Mo);


	/* create the MSG formatted output */
	m.evla = evla;
	m.evlo = evlo;
	m.depth = depth;
	m.mag1 = mag1;
	m.strike = strike;
	m.dip = dip;
	m.rake = rake;
	m.evdp = evdp;
	m.Mw = Mw;
	m.Mo = Mo;
	strcpy(m.title, title);
	strcpy(m.auth,"ENS ");
	strcpy(m.geog,geog);
	sprintf(m.evdate,"%4.4d/%2.2d/%2.2d",year,month,day);
	sprintf(m.evtime,"%2.2d:%2.2d:%2.2d:%1.1d",hour,minute,second,millisecond/100);
	outmsg(  year,  month,  day,  hour,  minute,  second, &m, file1, file3 );
	 
	 

	/* create the NDK formatted output, using some of the MSG computations */
	/* First adjust the everything according to the exponent */
	exponent = 0;
	n.m0 = Mo;
	while (n.m0 > 10){
		exponent++;
		n.m0/=10.0 ;
		m.mhrv[0][0]/=10.0;
		m.mhrv[0][1]/=10.0;
		m.mhrv[0][2]/=10.0;
		m.mhrv[1][0]/=10.0;
		m.mhrv[1][1]/=10.0;
		m.mhrv[1][2]/=10.0;
		m.mhrv[2][0]/=10.0;
		m.mhrv[2][1]/=10.0;
		m.mhrv[2][2]/=10.0;
		
	}

	/* fill up the ndk_ structure */
	/* Line 1 */
	strcpy(n.auth,"ENS ");
	sprintf(n.evdate,"%4.4d/%2.2d/%2.2d",year,month,day);
	sprintf(n.evtime,"%2.2d:%2.2d:%2.2d.%1.1d",hour,minute,second,millisecond/100);
	n.evla = evla;
	n.evlo = evlo;
	n.depth = depth;
	n.mag1 = mag1;
	n.mag2 = 0.0;
	strcpy(n.geog,geog);
	/* Line 2 */
	sprintf(n.cmtname,"C%4.4d%2.2d%2.2d%2.2d%2.2d%2.2dA",year,month,day,hour,minute,second);
	n.bsta = 0;
	n.bcomp = 0;
	n.bper = 0.;
	n.ssta = 0 ;
	n.scomp =  0;
	n.sper =  0.0;
	if(strlen(file2) > 0){
		if((file2id = fopen(file2, "r")) != NULL){
			if(fscanf(file2id,"%d %d %f",
				&n.ssta,&n.scomp,&n.sper)!=3){
				n.ssta = 0 ;
				n.scomp =  0;
				n.sper =  0.0;
			}
			fclose(file2id);
		}
	}
	n.msta = 0;
	n.mcomp = 0;
	n.mper = 0.;
	n.cmpsrc=2;
	strcpy(n.momratefunc,"TRIHD:");
	n.halfdur= 1.05*pow(Mo/1.0e24,0.33333);
	/* Line 3 */
        n.cmttime = 0.0;
        n.cmttimeerr = 0.0;
        n.clat = evla;
        n.claterr = 0.0;
        n.clon = evlo;
        n.clonerr = 0.0;
        n.cdepth = evdp;
        n.cdeptherr = 0.0 ;
	strcpy(n.depthtype,"FREE");
	time(&secs);
	gmt = gmtime(&secs);
	sprintf(n.timestamp,"S-%4.4d%2.2d%2.2d%2.2d%2.2d%2.2d",
		1900+gmt->tm_year,1+gmt->tm_mon,gmt->tm_mday,gmt->tm_hour,gmt->tm_min,gmt->tm_sec);
	/* Line 4 */
	n.exponent = exponent;
	n.Mrr = m.mhrv[0][0];
	n.EMrr = 0.0 ;
	n.Mtt = m.mhrv[1][1];
	n.EMtt = 0.0 ;
	n.Mpp = m.mhrv[2][2];
	n.EMpp = 0.0 ;
	n.Mrt = m.mhrv[0][1];
	n.EMrt = 0.0 ;
	n.Mrp = m.mhrv[0][2];
	n.EMrp = 0.0 ;
	n.Mtp = m.mhrv[1][2];
	n.EMtp = 0.0 ;
	/* Line 5 */
	strcpy(n.version,"V10");
	n.ev1 = n.m0;
	n.pl1 = m.plnt ;
	n.az1 = m.stkt;
	n.ev2 = 0.0;
	n.pl2 = m.plnn ;
	n.az2 = m.stkn;
	n.ev3 = - n.m0;
	n.pl3 = m.plnp ;
	n.az3 = m.stkp;
	n.strike1 = m.stk0;
	n.dip1 = m.dip0;
	n.rake1 =  m.rak0;
	n.strike2 = m.stk1;
	n.dip2 = m.dip1;
	n.rake2 = m.rak1;
	
	

	/* create the NDK formatted output */
	outndk(  year,  month,  day,  hour,  minute,  second, &n);


	exit(0);

}
void outmsg( int year, int month, int day, int hour, int minute, int second,  struct msg_ *m, char *file1, char *file3)
{
FILE *msgid;
char msgname[20];
char buffer[132];
float stkt, stkp, stkn, plnt, plnp, plnn;
float stk0, dip0, rak0;
float stk1, dip1, rak1;
int i, j;
FILE *fid;
	/* create event based file namefor MSG output */
	sprintf(msgname,"%4.4d%2.2d%2.2d%2.2d%2.2d%2.2d.msg",
		year,month,day,hour,minute,second);
	msgid = fopen(msgname,"w+");

	fprintf(msgid," %s\n",m->title);
	fprintf(msgid," %4s %10s %10s %6.2f %7.2f %5.1f %3.1f %s\n",
		m->auth,m->evdate,m->evtime,m->evla,m->evlo,m->depth,m->mag1,m->geog);
	if(strlen(file1) > 0){
		if((fid = fopen(file1, "r")) != NULL){
			fprintf(msgid," \n");
			fprintf(msgid," Stations used:\n");
			while(fgets(buffer, 132, fid) !=NULL)
				fprintf(msgid,"   %s",buffer);
			fclose(fid);
		}
	}
	if(strlen(file3) > 0){
		if((fid = fopen(file3, "r")) != NULL){
			fprintf(msgid," \n");
			fprintf(msgid," Filtering commands used:\n");
			while(fgets(buffer, 132, fid) !=NULL)
				fprintf(msgid,"   %s",buffer);
			fclose(fid);
		}
	}
	fprintf(msgid," \n");
	fprintf(msgid," Best Fitting Double Couple\n");
	fprintf(msgid,"  Mo = %8.2e dyne-cm\n",m->Mo);
	fprintf(msgid,"  Mw = %4.2f \n",m->Mw);
	fprintf(msgid,"  Z  = %d km\n",(int)m->evdp);

	trans( m->dip, m->strike, m->rake, &stkt, &plnt, &stkp, &plnp,  &plnn,  &stkn);
	m->stkt = stkt;
	m->plnt = plnt;
	m->stkn = stkn;
	m->plnn = plnn;
	m->stkp = stkp;
	m->plnp = plnp;
	tpdss( stkt, plnt, stkp, plnp,
	 	&stk0,  &dip0,  &rak0,
	 	&stk1,  &dip1,  &rak1) ;
	m->stk0 = stk0;
	m->dip0 = dip0;
	m->rak0 = rak0;
	m->stk1 = stk1;
	m->dip1 = dip1;
	m->rak1 = rak1;
	fprintf(msgid,"  Plane   Strike  Dip  Rake\n");
	fprintf(msgid,"   NP1      %3.0f   %3.0f   %3.0f\n",m->stk0,m->dip0,m->rak0);
	fprintf(msgid,"   NP2      %3.0f   %3.0f   %3.0f\n",m->stk1,m->dip1,m->rak1);
	fprintf(msgid,"  Principal Axes:\n");
	fprintf(msgid,"   Axis    Value   Plunge  Azimuth\n");
	fprintf(msgid,"    T  %9.2e    %3.0f     %3.0f\n",m->Mo,m->plnt,m->stkt);
	fprintf(msgid,"    N  %9.2e    %3.0f     %3.0f\n",0.0,m->plnn,m->stkn);
	fprintf(msgid,"    P  %9.2e    %3.0f     %3.0f\n",- m->Mo,m->plnp,m->stkp);

	sdr_to_mij(m->strike,m->dip,m->rake,m->mar);
	for(i=0 ; i < 3 ; i++){
		 for(j=0;j<3;j++){
			 m->mar[i][j] *= m->Mo;
			/* build up a copy for mhrv */
			m->mhrv[i][j] = m->mar[i][j];
		}
	}

	fprintf(msgid,"\n");
	fprintf(msgid," Moment Tensor: (dyne-cm)\n");
	fprintf(msgid,"    Component   Value\n");
	fprintf(msgid,"       Mxx    %9.2e\n",m->mar[0][0]);
	fprintf(msgid,"       Mxy    %9.2e\n",m->mar[0][1]);
	fprintf(msgid,"       Mxz    %9.2e\n",m->mar[0][2]);
	fprintf(msgid,"       Myy    %9.2e\n",m->mar[1][1]);
	fprintf(msgid,"       Myz    %9.2e\n",m->mar[1][2]);
	fprintf(msgid,"       Mzz    %9.2e\n",m->mar[2][2]);

	 
	PrintPradiation(msgid, m->mhrv,plnt,stkt,plnp,stkp);
	ar_to_hrv_mij(m->mhrv);
	 
	fprintf(msgid," Global CMT Convention Moment Tensor:\n");
	fprintf(msgid,"      R          T          P\n");
	PrintMTensor(msgid, m->mhrv);

	fclose(msgid);
}

void outndk( int year, int month, int day, int hour, int minute, int second, struct ndk_ *n)
{
FILE *ndkid;
char ndkname[20];

	/* create event based file namefor NDK output */

	sprintf(ndkname,"%4.4d%2.2d%2.2d%2.2d%2.2d%2.2d.ndk",
		year,month,day,hour,minute,second);
	ndkid = fopen(ndkname,"w+");

	/* LINE 1 */
	fprintf(ndkid,"%-4s %10s %10s %6.2f %7.2f %5.1f %3.1f %3.1f %s\n",
		n->auth,n->evdate,n->evtime,n->evla,n->evlo,n->depth,n->mag1,n->mag2,n->geog);

	/* LINE 2 */
	fprintf(ndkid,"%-16s B:%3d %4d %3.0f S:%3d %4d %3.0f M:%3d %4d %3.0f CMT: %1d %6s%5.1f\n",
		n->cmtname, n->bsta, n->bcomp, n->bper, n->ssta, n->scomp, n->sper, 
		n->msta, n->mcomp, n->mper,n->cmpsrc,n->momratefunc,n->halfdur);



	/* LINE 3 */
	fprintf(ndkid,"CENTROID: %8.1f %3.1f %6.2f %4.2f %7.2f %4.2f %5.1f %4.1f %-4s %-16s\n",

		n->cmttime,n->cmttimeerr,n->clat,n->claterr,n->clon,n->clonerr,
		n->cdepth,n->cdeptherr,n->depthtype,n->timestamp);

	/* LINE 4 */
	fprintf(ndkid,"%2d %6.3f %5.3f %6.3f %5.3f %6.3f %5.3f %6.3f %5.3f %6.3f %5.3f %6.3f %5.3f\n",

		n->exponent,n->Mrr,n->EMrr, n->Mtt,n->EMtt, n->Mpp,n->EMpp, 
		n->Mrt,n->EMrt, n->Mrp,n->EMrp, n->Mtp,n->EMtp);

	/* LINE 5 */
	fprintf(ndkid,"%3s %7.3f %2.0f %3.0f %7.3f %2.0f %3.0f %7.3f %2.0f %3.0f %7.3f %3.0f %2.0f %4.0f %3.0f %2.0f %4.0f\n",

		n->version,n->ev1,n->pl1,n->az1,n->ev2,n->pl2,n->az2,n->ev3,n->pl3,n->az3,
		n->m0,n->strike1,n->dip1,n->rake1,n->strike2,n->dip2,n->rake2);

	fclose(ndkid);
}

#define ABS(a)   ( (a) > (0) ? (a): -(a) )
#define SGN(a,b) ( (b) > (0) ? (ABS(a)): -(ABS(a)))

void trans(float dip,float stk,float slip,float *stkt,float *plnt,float *stkp,float *plnp, float *plnn, float *stkn)
{
/*
c-----
c     This takes strike, dip, and slip angles and defines a coordinate
c     transformation in which the slip vector, fault normal vector
c     and the normal vector of the auxiliary plane in terms of the
c     local NS, EW, and vertical cartesian coordinate system.
c-----
c     dip   angle of dip measured from horizontal plane.
c           0 < dip < 90
c     stk   direction of nodal plane strike. When looking
c           along stirke, fault plane dips downward to the
c           right
c           0 < stk < 360
c     slip  direction of movement along nodal plane.
c           If -180 < slip < 180, then the P-wave
c           vertically downward will a first motion polarity
c           equal to the sign of slip
c-----
c
c     REFERENCE:
c
c     Herrmann, R. B. (1975). A student's guide to the use of
c     P and S wave data, Earthquake Notes 46, 29-39
c
c-----
c     The X, Y and Z axes form a right handed coordinate
c     system, such that the compressional quadrant is
c     occurs whenever the xy > 0. The Z axis is the
c     null axis. The Y-axis is normal to the fault
c     plane and the X-axis is normal to the auxiliary plane
c     Note that for the input convention on dip, the
c     Y-axis will always point in the negative z-direction
c
c-----
*/
	float degrad, sins, coss, sind, cosd, sinf, cosf;
	float a11, a12, a13, a21, a22, a23, a31, a32, a33;
	float t1, t2, t3, p1, p2, p3;
	float xy, z3;
	float tmp1, tmp2, tmp3;
	float delx, dely, delz, betx, bety, betz;
	float delt, delp, bett, betp;
	float pnorm, tnorm, xnorm, ynorm, znorm;
	degrad=0.017452329;
	sins=sin(stk*degrad);
	coss=cos(stk*degrad);
	sind=sin(dip*degrad);
	cosd=cos(dip*degrad);
	sinf=sin(slip*degrad);
	cosf=cos(slip*degrad);
	/*
		X-axis
	*/
	a11=cosf*coss+sinf*cosd*sins;
	a12=cosf*sins-sinf*cosd*coss;
	a13=-sinf*sind;
	/*
		Y-axis
	*/
	a21=-sins*sind;
	a22=coss*sind;
	a23=-cosd;
	/*
		Z-axis
	*/
	a31=coss*sinf-cosd*cosf*sins;
	a32=cosd*cosf*coss+sinf*sins;
	a33=cosf*sind;
	t1=a11+a21;
	t2=a12+a22;
	t3=a13+a23;
	tnorm=sqrt(2.);
	if(t3 <= 0.0) tnorm=-tnorm;
	t1=t1/tnorm;
	t2=t2/tnorm;
	t3=t3/tnorm;
	p1=a11-a21;
	p2=a12-a22;
	p3=a13-a23;
	pnorm=sqrt(2.);
	if(p3 <= 0.0) pnorm=-pnorm;
	p1=p1/pnorm;
	p2=p2/pnorm;
	p3=p3/pnorm;
	xnorm=1.0;
	if(a13 < 0.0) xnorm=-1.0;
	ynorm=1.0;
	if(a23 < 0.0) ynorm=-1.0;
	znorm=1.0;
	if(a33 < 0.0) znorm=-1.0;
	a11=a11/xnorm;
	a12=a12/xnorm;
	a13=a13/xnorm;
	a21=a21/ynorm;
	a22=a22/ynorm;
	a23=a23/ynorm;
	a31=a31/znorm;
	a32=a32/znorm;
	a33=a33/znorm;
/*
c------
c     now all vectors point into the lower hemisphere.
c     To get the correct orientations, we require
c     that if the center of the focal sphere is a compression
c     that the X,Y, Z axes form a right handed coordinate system in the
c     lower hemisphere, otherwise it will form a left handed
c     coordinate system
c-----
c     obtain P-wave polarity at the center of the focal
c     sphere
c-----
*/
	xy = a13*a23;
	/*
		determine if right handed or left handed coordinate system
	*/
	z3 = a11*a22 - a12*a21;
	/*
		make right handed coordinate system
	*/
	if(z3 < 0.0){
		tmp1=a11 ;
		tmp2=a12 ;
		tmp3=a13 ;
		a11=a21 ;
		a12=a22 ;
		a13=a23 ;
		a21=tmp1 ;
		a22=tmp2 ;
		a23=tmp3 ;
	}
	if(SGN(1.0,xy) != SGN(1.0,slip)){
		tmp1=a11 ;
		tmp2=a12 ;
		tmp3=a13 ;
		a11=a21 ;
		a12=a22 ;
		a13=a23 ;
		a21=tmp1 ;
		a22=tmp2 ;
		a23=tmp3 ;
	}
	phase(a11,a12,a13,&delx,&betx);
	phase(a21,a22,a23,&dely,&bety);
	phase(a31,a32,a33,&delz,&betz);
	phase(t1,t2,t3,&delt,&bett);
	phase(p1,p2,p3,&delp,&betp);
/*	plx = 90-betx;
	ply = 90-bety;
	plz = 90-betz;
	plt=90-bett;
	plp=90-betp;
*/
	*stkt=delt;
	*plnt=90.0-bett;
	*stkp=delp;
	*plnp=90.0-betp;
	*stkn=delz;
	*plnn=90.0-betz;
}

void phase(float x,float y,float z,float *del,float *bet)
{
/*
c-----
c     This expresses a three component vector in terms
c     of two angles in a spherical coordinate system
c-----
c     x     x-coordinate of vector
c     y     y-coordinate of vector
c     z     z-coordinate of vector
c-----
c     del   direction of projection of the vector on a
c           horizontal plane, measured positively, clockwise
c           from north - the trend of the vector.
c     bet   angle that the 3-D vector makes with the
c           downward vertical. This is just 90 - plunge
c-----
*/
	float degrad;
	degrad=0.017452329;
	*bet=acos(z);
	*bet/=degrad;
	if(x < 0.0){
		*del=atan((y/x))/degrad+180. ;
	} else if(x == 0.0){
		if(y < 0.0)
			*del = 270;
		else if(y == 0.0)
			*del = 0;
		else
			*del = 90.;
	} else {
		if(y < 0.0)
			*del = atan((y/x))/degrad + 360.0;
		else if(y >= 0.0)
			*del=atan((y/x))/degrad;
		
	}
}

void tpdss(float stkt,float plnt,float stkp,float plnp,
	float *stk0, float *dip0, float *rak0,
	float *stk1, float *dip1, float *rak1) 
{
/*
	This takes the orientations of the tension and pressure 
	axes, and determines the two nodal planes 
	convert (P,T) to (dip,slip,stk) and (X,Y,Z). 
*/
	float t1,t2,t3,p1,p2,p3,f1,f2,f3,xn1,xn2,xn3 ;
	float deg;
	float stk,slip,dip,stkk,slipp,dipp,xx,yy ;
	deg=3.141592653/180.0 ;
/*
	strike angle is measured from the north clockwise. 
	plunge angle is measured from the horizon downward. 
*/
	stkp=stkp*deg ;
	plnp=plnp*deg ;
	stkt=stkt*deg ;
	plnt=plnt*deg ;
	t1=cos(plnt)*cos(stkt) ;
	t2=cos(plnt)*sin(stkt) ;
	t3=sin(plnt) ;
	p1=cos(plnp)*cos(stkp) ;
	p2=cos(plnp)*sin(stkp) ;
	p3=sin(plnp) ;
	f1=t1+p1 ;
	f2=t2+p2 ;
	f3=t3+p3 ;
	yy=sqrt(f1*f1+f2*f2+f3*f3) ;
	f1=f1/yy ;
	f2=f2/yy ;
	f3=f3/yy ;
	xn1=t1-p1 ;
	xn2=t2-p2 ;
	xn3=t3-p3 ;
	yy=sqrt(xn1*xn1+xn2*xn2+xn3*xn3) ;
	xn1=xn1/yy ;
	xn2=xn2/yy ;
	xn3=xn3/yy ;
	if(xn3 > 0.0){ 
		xn1=-xn1 ;
		xn2=-xn2 ;
		xn3=-xn3 ;
		f1=-f1 ;
		f2=-f2 ;
		f3=-f3 ;
	}
	dip=acos(-xn3) ;
	stk=atan2(-xn1,xn2) ;
	xx=f1*xn2-f2*xn1 ;
	slip=atan2(-f3,xx) ;
	*dip0=dip/deg ;
	*rak0=slip/deg ;
	if(*rak0 < -180.)*rak0+=180. ;
	if(*rak0 > 180.)*rak0-=360. ;
	*stk0=stk/deg ;
	if(*stk0 < 0.0) *stk0+=360. ;

	if(f3 > 0.0){
		xn1=-xn1 ;
		xn2=-xn2 ;
		xn3=-xn3 ;
		f1=-f1 ;
		f2=-f2 ;
		f3=-f3 ;
	}
	dip=acos(-f3) ;
	stk=atan2(-f1,f2) ;
	xx=f2*xn1-f1*xn2 ;
	slip=atan2(-xn3,xx) ;
	*dip1=dip/deg ;
	*rak1=slip/deg ;
	if(*rak1 < -180.)*rak1+=180. ;
	if(*rak1 > 180.)*rak1-=360. ;
	*stk1=stk/deg ;
	if(*stk1 < 0.0) *stk1+=360. ;
}
