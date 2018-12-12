/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: RLINER                                                c
c                                                                     c
c      COPYRIGHT (C)  1986, 1989 R. B. Herrmann                       c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
*/
/* 
	line plotting algorithm Integer Digital Differential Analyzer
	Berger, M. (1986). Computer Graphics with Pascal, The
	Benjamin/Cummings Publishing Company, Inc.,
	pp 39-41, converted from Pascal to C by RBHerrmann
								*/

#ifndef INT
#ifdef MSDOS
#define INT long
#else
#define INT int
#endif
#endif

extern int dc_curpen;	/* current pen number for gray area plot */
extern INT dc_left, dc_bottom, dc_right, dc_top;
extern INT dc_rotate;
extern INT dc_shadeoff;

extern INT dv_lineclip(INT *x0,INT *z0,INT *x1,INT *z1,
	INT MinX,INT MinY,INT MaxX,INT MaxY);
extern void dv_zzpoint(INT x,INT y);
extern void dv_zpoint(INT x,INT y);

void dv_rliner(INT x0,INT z0,INT x1,INT z1)
{
	register x,y;
	register count,error;
	register ypx, ymx;
	INT deltax,deltay;
	extern INT dv_lineclip();
	INT tmp;
	/* invoke rotation if required */
	if(dc_rotate){
		tmp = z0;
		z0 = x0;
		x0 = dc_right - tmp;
		tmp = z1;
		z1 = x1 ;
		x1 = dc_right - tmp;
	}
	if(dv_lineclip(&x0,&z0,&x1,&z1,dc_left,dc_bottom,dc_right,dc_top))
		return;
	error = 0;
	deltax = x1 - x0;
	deltay = z1 - z0;
	if(deltay < 0)
	{
		tmp = x0;
		x0 = x1;
		x1 = tmp;
		tmp = z0;
		z0 = z1;
		z1 = tmp;
		deltax = -deltax;
		deltay = -deltay;
	}
	ypx = deltax+deltay;
	ymx = deltay-deltax;
	x = x0;
	y = z0;
	x = x0;
	y = z0;
	dv_zzpoint (x,y);
	if(deltax >= 0){			/* positive slope */
		if(deltax > deltay){		/* 0 < slope < 1 */
			for(count =1;count < deltax;count++){
				if(error <=  0){
					x++;
					dv_zzpoint(x,y);
					error+= deltay;
				}
				else{
					x++;
					y++;
					dv_zzpoint(x,y);
					error+=ymx;
				}
			}
		}
		else if(deltax < deltay) {	/* slope > 1 */
			for(count=1;count < deltay;count++){
				if(error <  0){
					x++;
					y++;
					dv_zzpoint(x,y);
					error+= ymx;
				}
				else {
					y++;
					dv_zzpoint(x,y);
					error-=deltax;
				}
			}
		}
		else {				/* slope = 1 */
			x++;
			y++;
			while(x<x1){
				dv_zzpoint(x++,y++);
			}
		}
	}
	else {					/* negative slope */
		if(-deltax > deltay){		/* -1 < slope < 0 */
			for(count=1;count < -deltax;count++){
				if(error <=  0){
					x--;
					dv_zzpoint(x,y);
					error+=deltay;
				}
				else {
					x--;
					y++;
					dv_zzpoint(x,y);
					error+= ypx;
				}
			}
		}
		else if (-deltax < deltay) {				/* slope < -1 */
			for(count=1;count<deltay;count++){
				if(error <  0){
					x--;
					y++;
					dv_zzpoint(x,y);
					error+= ypx;
				}
				else {
					y++;
					dv_zzpoint(x,y);
					error+=deltax;
				}
			}
		}
		else {
			x--;
			y++;
			while(x>x1){
				dv_zzpoint(x--,y++);
			}
		}
	}
	dv_zzpoint(x1,z1);
}
