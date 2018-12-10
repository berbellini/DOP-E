#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#define ABS(a  ) ( (a) >  0  ? (a):-(a))
#define MIN(a,b) ( (b) > (a) ? (a):(b) )
#define MAX(a,b) ( (b) < (a) ? (a):(b) )
#define PI  3.141592653589793238512808959406186204433
#define PI2 1.570796326794896619256404479703093102216

/* CHANGES
	03 SEP 2008 - permit use of GCARC for DELdEG and DIST for DELKM
*/


void delaz(float evla,float evlo,float stla,float stlo,float *del,float *az,float *baz,float *delkm,float *saz,float *caz);
void dircos(double lat,double lon,double *stla,double *clat,double *stlo,double *clon,double *aa,double *bb,double *cc);
void gcmdln(int argc, char **argv, float *evla,float *evlo,float *stla,float *stlo,int *dodelkm,int *dodeldg,int *doaz,int *dobaz, int *donl);
void usage();
double ellip(double k2, double t1, double t2);

main(int argc, char **argv)
{
	float evla, evlo, stla, stlo;
	float del, az, baz, delkm, saz, caz;
	int dodelkm, dodeldg, doaz, dobaz, donl;
	gcmdln(argc, argv, &evla,&evlo,&stla,&stlo,&dodelkm,&dodeldg,
		&doaz,&dobaz,&donl);
	delaz(evla,evlo,stla,stlo, &del, &az, &baz, &delkm, &saz, &caz);
	if(doaz == 0 && dobaz == 0 && dodeldg == 0 && dodelkm ==0){
		printf("DELDEG=%f DELKM=%f AZ=%f BAZ=%f\n",del,delkm,az,baz);
	} else {
		if(dodelkm)
			printf("%f ",delkm);
		if(dodeldg)
			printf("%f ",del);
		if(doaz)
			printf("%f ",az);
		if(dobaz)
			printf("%f ",baz);
		if(donl)
			printf("\n");

	}
}

double degrad = 0.0174532925199433;	/* degrees to radians */
double a = 6378.137 ;			/* equatorial radius */
double b = 6356.752 ;			/* polar radius = a(1-f) */
double f = 0.003352810664747481 ;	/* f = 1 - b/a = 
						 1/298.257223563 */
double e = 0.081819191 ;			/* eccentricity (A.38)
						e^2 = (a^2 - b^2)/a^2
						e^2 = 2f - f^2 */
double e2 = 0.006694379990141317 ;	/* e^2 */
double eps = 0.006739496742276435 ;	/* second ellipticity (A.42)
						eps = (a^2 - b^2)/b^2
						= e^2/(1-e^2) */
static float RE=6371.003;

void delaz(float evla,float evlo,float stla,float stlo,float *del,float *az,float *baz,float *delkm,float *saz,float *caz)
{

/*
	evla	float	Epicenter latitude (degrees)
	evlo	float	Epicenter longitude (degrees)
	stla	float	Station latitude (degrees)
	stlo	float	Station longitude (degrees)
	del	float *	Epicentral Distance (degrees)
	az	float *	Epicenter to Station Azimuth
	baz	float *	Station to Epicenter Backazimuth
	delkm	float *	Epicentral Distance (km)
	saz	float *	Sine of Azimuth
	caz	float *	Cosine of Azimuth
*/
	double devla, devlo, dstla, dstlo;
	double sevla, cevla, sevlo, cevlo, ea, eb, ec;
	double sstla, cstla, sstlo, cstlo, sa, sb, sc;
	double cdel, fac, tmp;
	double cbz, sbz;

	int i, ni, ang;
	double c,s,t1,t2,t3,xp,zp;
	double x,y,z;

	double xp1, xp2, xp3;
	double yp1, yp2, yp3;
	double zp1, zp2, zp3;
	double xpe, zpe;
	double xps, zps;
	double re, rs;		/* radial distance from center of sphere to surface */
	double bp, ap;
	double te, ts;
	double e2p;

	devla = (double)evla;
	devlo = (double)evlo;
	dstla = (double)stla;
	dstlo = (double)stlo;

	/* major stability hack */
	if(devla ==  90.0)devla =  89.99999 ;
	if(devla == -90.0)devla = -89.99999 ;
	if(dstla ==  90.0)dstla =  89.99999 ;
	if(dstla == -90.0)dstla = -89.99999 ;
	if(devla == 0.0 && dstla == 0.0){
		dstla =  0.000001 ;
		devla = -0.000001 ;
	}
	if( ( dstla == - devla ) && ( ABS( devla - dstlo) == 180.0)){
		dstla += 0.000001 ;
		dstlo += 0.000001 ;
	}

	dircos(devla,devlo,&sevla,&cevla,&sevlo,&cevlo,&ea,&eb,&ec);
	dircos(dstla,dstlo,&sstla,&cstla,&sstlo,&cstlo,&sa,&sb,&sc);

	/*
	compute distance
	Choose correct formula for short and large distances
	*/
	cdel = ea*sa + eb*sb + ec*sc;
	/*
		if DEL == 0
	*/
	if(cdel ==  1.0){
		*del = 0.0;
	/*
		if DEL = [0,20)
	*/
	} else if(cdel >  0.9396){
		fac = (ea-sa)*(ea-sa) + (eb-sb)*(eb-sb) + (ec-sc)*(ec-sc);
		fac = sqrt(fac)/2.0;
		*del = 2.0*asin(fac);
	/*
		if DEL = [20,160]
	*/
	} else if(cdel <= 0.9396 && cdel >= -0.9396){
		*del = acos(cdel);
	/*
		if DEL = (160,180]
	*/
	} else {
		fac = (ea+sa)*(ea+sa) + (eb+sb)*(eb+sb) + (ec+sc)*(ec+sc);
		fac = sqrt(fac)/2.0;
		*del = 2.0*acos(fac);
	}
	
	/*
	check for station or epicenter at pole
	*/
	if( *del == 0.0 ){
		*az = 0.0 ;
		*baz = 0.0 ;
		*delkm = 0.0;
		return;
	}
	if(devla == 90.0 ){
		*az = 360.0 - dstlo ;
		*baz = 0.0;
		*saz = sin(degrad*(360 - dstlo));
		*caz = cos(degrad*(360 - dstlo));
	} else if(dstla == 90.0 ){
		*az = 0.0;
		*baz = 360.0 - devlo;
		*saz =  0.0;
		*caz =  1.0;
	} else if(devla == -90.0 ){
		*az =  stlo;
		*baz = 180.0;
		*saz = sin(degrad*(dstlo)); 
		*caz = cos(degrad*(dstlo)); 
	} else if(dstla == -90.0 ){
		*az = 180.0 ;
		*baz = devlo ;
		*saz =  0.0 ;
		*caz = -1.0 ;
	} else {
		*saz = cevla*(cstla * sin(degrad*(dstlo - devlo)));
		*caz = (sstla - cdel*sevla);
		fac = sqrt((*saz)*(*saz) + (*caz)*(*caz));
		if(fac > 0.0){
			*saz = *saz / fac;
			*caz = *caz / fac;
			*az = atan2(*saz,*caz);
		
			sbz = - cstla*(cevla * sin(degrad*(dstlo - devlo)));
			cbz = (sevla - cdel*sstla);
			*baz = atan2(sbz,cbz);
		} else {
			*az = 0.0;
			*caz = 1.0;
			*saz = 0.0;
			*baz = 180.0;
		}
 		*az = *az / degrad;
		*baz = *baz / degrad;
	}

	*delkm = (*del) * RE;
	*del = (*del) / degrad;
	/* put az and baz in the range [0,360) */
	if( *az < 0.0)*az = *az + 360.0 ;
	if( *az >= 360.0) *az =  *az - 360.0 ;
	if( *baz < 0.0)*baz = *baz + 360.0 ;
	if( *baz >= 360.0) *baz = *baz - 360.0 ;

	/* to get this to work for very small distances, we will only do
		the elliptic integral is *del > 0.001 degrees
	*/
	if( *delkm > 1.0 ) {

		/* to compute the distance in kilometers we
			follow the idea of Rudoe (Bomford, 1980) but do it
			algebraically as much as possible.
			1. E and S define a plane which goes through the origin
			2. The intersection of this plane with the spheroid is
			  an ellipse whose major axis is 'a' but whose minor axis is
			  b <=  minor_axis <= a
			3. Now we just define  latitudes on the ellipse 
			4. Get the distance by an elliptic integral 
				or approximate formula
		*/
		/* define the plane:  alpha x + beta y + gamma z = 0 since
			plane goes through the origin. The vector normal to the
			plane is ( alpha, beta, gamma) = E x S
		*/ 
		re = a*sqrt(1.0 -e2*sevla*sevla);
		rs = a*sqrt(1.0 -e2*sstla*sstla);
		yp1 =  eb * sc - sb * ec ;
		yp2 =  ec * sa - ea * sc ;
		yp3 =  ea * sb - eb * sa ;
		/* safety */
		fac = sqrt ( yp1*yp1 + yp2*yp2 + yp3*yp3);
		if(ABS(yp1) <  1.0e-6 * fac)
			yp1 =  1.0e-6 ;
		if(ABS(yp2) <  1.0e-6 * fac)
			yp2 = -1.0e-6 ;
		fac = sqrt ( yp1*yp1 + yp2*yp2 + yp3*yp3);
		yp1 /= fac ;
		yp2 /= fac ;
		yp3 /= fac;
		/* since this defines the normal, now define other unit vectors
		*/
		xp1 =  yp2;
		xp2 = -yp1 ;
		xp3 = 0.0 ;
		fac = sqrt ( xp1*xp1 + xp2*xp2 + xp3*xp3);
		xp1 /= fac ;
		xp2 /= fac ;
		xp3 /= fac;

		zp1 = - yp1*yp3 ;
		zp2 = - yp2*yp3 ;
		zp3 = yp1*yp1 + yp2*yp2 ;
		fac = sqrt ( zp1*zp1 + zp2*zp2 + zp3*zp3);
		zp1 /= fac ;
		zp2 /= fac ;
		zp3 /= fac;

		/*  by construction the xp axis is horizontal - this has
			an ellipse width of 2*ap = 2*a
		    the xp axis is the director of the minor axis and
			has length of 2*bp
			we solve the ellipse equation
		*/
		ap = a ;
		bp = 1./sqrt( (zp1/a)*(zp1/a) + (zp2/a)*(zp2/a) + (zp3/b)*(zp3/b));
		e2p = (1.0 - bp/ap)*(1.0 + bp/ap);

		/* now get the points on the new ellipse E(3d) -> E'(2d)
			S(3d) -> S'(2d) */

		xpe = ea*xp1 + eb*xp2 + ec*xp3 ;
		zpe = ea*zp1 + eb*zp2 + ec*zp3 ;
		te = atan2(zpe/bp , xpe/ap);
		xps = sa*xp1 + sb*xp2 + sc*xp3 ;
		zps = sa*zp1 + sb*zp2 + sc*zp3 ;

		ts = atan2(zps/bp , xps/ap);
/*
		printf("te %f ",te/degrad);
		printf("ts %f ",ts/degrad);
		printf("e2p: %f ",e2p);
*/

		/* compute the linear distance by integrating along the path.
			However be careful to use the minor arc and to use
			increasing angle */
		/* first order */
		if ( ts > te ){
			tmp = ts ;
			ts = te ;
			te = tmp ;
		}
/*
printf("ts %f te %f\n",ts/degrad,te/degrad);
*/
		/* now check minor arc */
		if( (te -ts ) <= PI ){
			/* minor arc */
			*delkm = a * ellip(e2p, (ts)/degrad,(te)/degrad) ;
		} else {
			*delkm = a * ellip(e2p, (te)/degrad,180.0) 
				+ a * ellip(e2p, -180.0,(ts)/degrad) ;
		}
	}
}

#define N 500
double ellip(double k2, double t1, double t2)
{
	double degrad = 0.0174532925199433;
	/* compute incomplete elliptic integral 
	   t2          2    2    1/2
        INT   ( 1.0 - k cos t  )  dt
           t1
	NOTE I USE cos AND NOT sin BECAUSE OF MY COORDINATES
	*/
	double dt, t, ct, st, dct, dst;
	int i;
	double sum;
	double st2;
	double ct2;
	double tmp;
	double k;

/*
printf("Integrate k2=%f (%f,%f)\n",k2,t1,t2);
*/

/*
	special care for a small interval
*/
	if(ABS(t2 - t1) < 0.01){
		k = sqrt(k2);
		ct = cos(0.5*(t1+t2)*degrad);
		sum = (t2-t1)*degrad * sqrt(1.0-ct*k)*sqrt(1.0+ct*k);
	} else {
		sum = 0.00;
		dt = (t2-t1)*degrad/N ;
		st = sin(t1*degrad);
		ct = cos(t1*degrad);
		dct = cos(dt);
		dst = sin(dt);
		t = t1 ;
		for(i=0;i<=N;i++){
			/*
			Use rectangular rule
			*/
			st2 = st*st ;
			ct2 = ct*ct ;
			if(i==0)
				sum += 0.5*dt*sqrt(1.0-ct2*k2);
			else if(i==N)
				sum += 0.5*dt*sqrt(1.0-ct2*k2);
			else
				sum += dt*sqrt(1.0-ct2*k2) ;
			tmp = ct*dct - st*dst;
			st  = st*dct + ct*dst;
			ct = tmp;
			t+=dt;
		}
	}
	return (sum);
	
}

	
void dircos(double lat,double lon,double *slat,double *clat,double *slon,double *clon,double *aa,double *bb,double *cc)
{
/*
	convert geographic latitude to geocentric
	
	The realtion between geocentric and geographic latitude is
	tan phi c = ( 1 - f)^2 tan phi g

	To avoid problems at poles, define sin phi c and cos phi c
	so that the relation holds and also that s^2 + c^2 = 1

*/
	double c,s,e4,fac;
	c = cos(degrad*lat);
	s = sin(degrad*lat);
	e4 = (1.-e2)*(1.-e2);
	fac = sqrt(e4 + (1.0-e4)*c*c);
	*slat = (1.-e2) * s /fac;
	*clat =      c /fac;
	*slon = sin(degrad*lon);
	*clon = cos(degrad*lon);

	*aa = (*clat) * (*clon);
	*bb = (*clat) * (*slon);
	*cc = (*slat);
}

void gcmdln(int argc, char **argv, float *evla,float *evlo,float *stla,float *stlo,int *dodelkm,int *dodeldg,int *doaz,int *dobaz, int *donl)
{
	char *cp;

	*evla = 0.0;
	*evlo = 0.0;
	*stla = 0.0;
	*stlo = 0.0;
	*dodelkm = 0;
	*dodeldg = 0;
	*doaz = 0;
	*dobaz = 0;
	*donl = 0;

	if(argc == 1)usage();

	while(argc-- > 1 ) {
		if(*argv[1] == '-'){
			cp = argv[1];
			cp++;
			if(strcmp("ELAT",cp) == 0 || strcmp("evla",cp)==0 || 
				strcmp("elat",cp) == 0 || strcmp("EVLA",cp)==0){
				argv++;
				argc--;
				cp = argv[1];
				*evla = atof(cp);
			} else if(strcmp("ELON",cp) == 0 || strcmp("evlo",cp)==0||
				strcmp("elon",cp) == 0 || strcmp("EVLO",cp)==0){
				argv++;
				argc--;
				cp = argv[1];
				*evlo = atof(cp);
			} else if(strcmp("SLAT",cp) == 0 || strcmp("stla",cp)==0 ||
				strcmp("slat",cp) == 0 || strcmp("STLA",cp)==0){
				argv++;
				argc--;
				cp = argv[1];
				*stla = atof(cp);
			} else if(strcmp("SLON",cp) == 0 || strcmp("stlo",cp)==0||
				strcmp("slon",cp) == 0 || strcmp("STLO",cp)==0){
				argv++;
				argc--;
				cp = argv[1];
				*stlo = atof(cp);
			} else if(strcmp("DELKM",cp) == 0 || strcmp("delkm",cp)==0){
				*dodelkm = 1;
			} else if(strcmp("DIST",cp) == 0 || strcmp("dist",cp)==0){
				*dodelkm = 1;
			} else if(strcmp("DELDEG",cp) == 0 || strcmp("deldeg",cp)==0){
				*dodeldg = 1;
			} else if(strcmp("GCARC",cp) == 0 || strcmp("gcarc",cp)==0){
				*dodeldg = 1;
			} else if(strcmp("AZ",cp) == 0 || strcmp("az",cp)==0){
				*doaz = 1;
			} else if(strcmp("BAZ",cp) == 0 || strcmp("baz",cp)==0){
				*dobaz = 1;
			} else if(strcmp("NL",cp) == 0 || strcmp("nl",cp)==0){
				*donl = 1;
			} else if(strcmp("h",cp) == 0 || strcmp("?",cp)==0){
				usage();
			}
			argv++;
	}
}
}

void usage()
{
fprintf(stderr,"Usage: udelaz -ELAT evla -ELON -evlo -SLAT stla -SLON stlo\n");
fprintf(stderr,"              [-DELKM -DELDEG -AZ -BAZ] [-NL]\n");
fprintf(stderr," -ELAT elat or -EVLA evla  (default 0.0) Event    latitude N=+\n");
fprintf(stderr," -SLAT slat or -STLA stla  (default 0.0) Station  latitude N=+\n");
fprintf(stderr," -ELON elon or -EVLO evlo  (default 0.0) Event   longitude E=+\n");
fprintf(stderr," -SLON slon or -STLO stlo  (default 0.0) Station longitude E=+\n");
fprintf(stderr,"  Output will be of form\n");
fprintf(stderr,"  DELDEG=0.000000 DELKM=0.000000 AZ=0.000000 BAZ=0\n");
fprintf(stderr,"  unless one of the following flags is used, in which case\n");
fprintf(stderr,"  number is returned but no text, for use in SHELL script:\n");
fprintf(stderr,"  AZ=`udelaz -ELAT 15 -ELON -20 -SLAT 40 -SLON -30 -AZ`\n");
fprintf(stderr," -AZ             Output only event -> station azimuth\n");
fprintf(stderr," -BAZ            Output only station -> event azimuth\n");
fprintf(stderr," -DELDEG or -GCARC Output only great circle distance in degrees\n");
fprintf(stderr," -DELKM or -DIST  Output only great circle distance in km\n");
fprintf(stderr," -NL    (default false) Output newline (not required for shells cripts\n");
exit(0);
}


