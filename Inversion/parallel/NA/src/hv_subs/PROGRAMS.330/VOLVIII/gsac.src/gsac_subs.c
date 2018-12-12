/* changes
 * ensured that BAZ and AZ are correcly computed for EVLA=STLA EVLO=STLO
 * this is done by the *rad in lines 93
 * */
#include "gsac.h"
#include	<math.h>
static void dircos(float lat,float lon,double *slat,double *clat,double *slon, double *clon, double *aa,double *bb,double *cc);
void delaz(float elat, float elon, float slat, float slon, float *deldeg, float *az, float *baz, float *delkm);
double ellip(double k2, double t1, double t2);

/* The Earth Ellipsoid parameters are from
        WGS84 http://en.wikipedia.org/wiki/Figure_of_the_Earth 
        a = 6378.137;
        f = 1./298.257223563;

   the results agree with Equation 1b of
        Thomas, P. D. (1965). Geodesic arc length on the reference
        ellipsoid to second-order terms in the flattening,
        J. Geophys. Res. 70, 3331-3340.
    when I use the WGS84 values of a and f.

    Compared to udelaz (which is also inside of gsac_subs.c) I find the
        following differences:
        STLA    STLO   EVLA   EVLO   udelaz  Thomas
        0       0       90      90      10001.965820 10001.958678 NS
        0       0       0       90      10018.753906 10018.754171 Equator
        -90     90      90      90      20003.929688 20003.915653 NS
        0       0       0       179.99  20036.392578 20036.395759 Equator
        0       0       0.01    0       1.105761 1.098228      NS
        0       0.01    0       0       1.113213 1.113195      Equator

        The difference for the Equator 0 0.01 0 0 may be due to the
        integration not being fine enough

        The Thomas formula is OK and a lot simpler but requires care
        when DELTA = 0 or 180 degrees because of division by
        (1+cos d) and (1-cos d). A taylor series expansion might

*/


/* general support routines 
 * distance computation
 * */
#define PI  3.141592653589793238512808959406186204433
#define PI2 1.570796326794896619256404479703093102216
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
void delaz(float evla,float evlo,float stla,float stlo,float *del,float *az,float *baz,float *delkm)
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
	double sevla, cevla, sevlo, cevlo, ea, eb, ec;
	double sstla, cstla, sstlo, cstlo, sa, sb, sc;
	double cdel, fac, tmp;
	double cbz, sbz;
	float saz, caz;

	double xp1, xp2, xp3;
	double yp1, yp2, yp3;
	double zp1, zp2, zp3;
	double xpe, zpe;
	double xps, zps;
	double re, rs;		/* radial distance from center of sphere to surface */
	double bp, ap;
	double te, ts;
	double e2p;

	/* major stability hack */
	if(evla ==  90.0)evla =  89.9999 ;
	if(evla == -90.0)evla = -89.9999 ;
	if(stla ==  90.0)stla =  89.9999 ;
	if(stla == -90.0)stla = -89.9999 ;
	if(evla == 0.0 && stla == 0.0){
		stla =  0.0001 ;
		evla = -0.0001 ;
	}
	if( ( stla == - evla ) && ( ABS( evla - stlo) == 180.0)){
		stla += 0.001 ;
		stlo += 0.001 ;
	}

	dircos(evla,evlo,&sevla,&cevla,&sevlo,&cevlo,&ea,&eb,&ec);
	dircos(stla,stlo,&sstla,&cstla,&sstlo,&cstlo,&sa,&sb,&sc);

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
	if(evla == 90.0 ){
		*az = 360.0 - stlo ;
		*baz = 0.0;
		saz = sin(degrad*(360 - stlo));
		caz = cos(degrad*(360 - stlo));
	} else if(stla == 90.0 ){
		*az = 0.0;
		*baz = 360.0 - evlo;
		saz =  0.0;
		caz =  1.0;
	} else if(evla == -90.0 ){
		*az =  stlo;
		*baz = 180.0;
		saz = sin(degrad*(stlo)); 
		caz = cos(degrad*(stlo)); 
	} else if(stla == -90.0 ){
		*az = 180.0 ;
		*baz = evlo ;
		saz =  0.0 ;
		caz = -1.0 ;
	} else {
		saz = cevla*(cstla * sin(degrad*(stlo - evlo)));
		caz = (sstla - cdel*sevla);
		fac = sqrt((saz)*(saz) + (caz)*(caz));
		if(fac > 0.0){
			saz = saz / fac;
			caz = caz / fac;
			*az = atan2(saz,caz);
		
			sbz = - cstla*(cevla * sin(degrad*(stlo - evlo)));
			cbz = (sevla - cdel*sstla);
			*baz = atan2(sbz,cbz);
		} else {
			*az = 0.0;
			caz = 1.0;
			saz = 0.0;
			*baz = 180.0;
		}
 		*az = *az / degrad;
		*baz = *baz / degrad;
	}

	*delkm = *del * re;
	*del = *del / degrad;
	/* put az and baz in the range [0,360) */
	if( *az < 0.0)*az = *az + 360.0 ;
	if( *az >= 360.0) *az =  *az - 360.0 ;
	if( *baz < 0.0)*baz = *baz + 360.0 ;
	if( *baz >= 360.0) *baz = *baz - 360.0 ;

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

/*
printf("Integrate (%f,%f)\n",t1,t2);
*/

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
	return (sum);
	
}

	
void dircos(float lat,float lon,double *slat,double *clat,double *slon,double *clon,double *aa,double *bb,double *cc)
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
