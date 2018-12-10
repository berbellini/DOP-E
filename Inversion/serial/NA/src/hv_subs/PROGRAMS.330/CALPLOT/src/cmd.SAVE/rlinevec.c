/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: RLINER                                                c
c                                                                     c
c      COPYRIGHT (C)  1997       R. B. Herrmann                       c
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
#ifdef MSDOS
#define INT long
#else
#define INT int
#endif

/* prototype definitions */
void dv_rliner(INT x0,INT z0,INT x1,INT z1);
extern INT dv_lineclip(INT *x0,INT *z0,INT *x1,
		INT *z1,INT dc_left,INT dc_bottom,INT dc_right,INT dc_top);
extern void dv_zzpoint(INT x, INT y);
extern void dv_concur(INT isx,INT isy);
extern void dv_movcur(INT isx,INT isy);

#include	<math.h>
#define		ABS(x)	((x) > 0 ? (x) : (-x))

extern int dc_curpen;	/* current pen number for gray area plot */
extern INT dc_left, dc_bottom, dc_right, dc_top;
extern INT dc_oldx1, dc_oldy1;
extern INT dc_rotate;
extern INT dc_hardwarefill;
extern INT dc_shdon;

extern INT dc_ClipRight, dc_ClipTop, dc_ClipBottom, dc_ClipLeft;

void dv_rliner(INT x0,INT z0,INT x1,INT z1)
{
	INT tmp;
	INT i,maxv,minv;
	/* invoke rotation if required */
	if(dc_rotate){
		tmp = x0;
		x0 = z0 ;
		z0 = dc_top - tmp;
		tmp = x1;
		x1 = z1 ;
		z1 = dc_top - tmp;
	}
	if(dv_lineclip(&x0,&z0,&x1,&z1,dc_ClipLeft,dc_ClipBottom,dc_ClipRight,dc_ClipTop))
		return;
	if(dc_oldx1 != x0 || dc_oldy1 != z0)	/* cut down on move eg	*/
		dv_movcur(x0,z0);		/* extra pen up/down for*/
					/* mechanical plotters  */
	/* if dc_curpen > 1000 impose a gray level shading */
	if(dc_curpen < 1000){
		dv_concur(x1,z1);
	}
	else {
		if(x1 == x0){
			if(z1 > z0){
				maxv = z1;
				minv = z0;
			}
			else {
				maxv = z0;
				minv = z1;
			}
			if(dc_hardwarefill || dc_shdon  )
				dv_concur(x1,z1);
			else {
				for(i=minv;i<= maxv;i++){
					dv_zzpoint(x1,i);
				}
			}
		}
		else if(z1 == z0){
			if(x1 > x0){
				minv = x0;
				maxv = x1;
			}
			else {
				minv = x1;
				maxv = x0;
			}
			if(dc_hardwarefill || dc_shdon)
				dv_concur(x1,z1);
			else {
				for(i=minv;i <= maxv;i++){
					dv_zzpoint(i,z1);
				}

			}
			
		}
		else {
			dv_concur(x1,z1);
		}
	}
	dc_oldx1 = x1;
	dc_oldy1 = z1;
}
