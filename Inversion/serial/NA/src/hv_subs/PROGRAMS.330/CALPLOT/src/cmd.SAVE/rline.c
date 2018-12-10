/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: RLINE                                                 c
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

#define		ABS(x)	((x) > 0 ? (x) : (-x))

/* prototypes */

static INT do_d(INT x,INT y,INT wx,INT wy);
void dv_rline(INT X0,INT Y0,INT X1,INT Y1);
extern void dv_rliner(INT x0,INT z0,INT x1,INT z1);

/* external variable definitions */

extern INT	dc_rotate;
extern INT dc_shadeoff;
extern INT dc_xlinewidth;	/* current line width  in x device units */
extern INT dc_ylinewidth;	/* current line width  in y device units */
extern INT dc_linewidth;	/* current line width  in plot units */
extern INT dc_herelinewidth; /* set linewidth here */

void dv_rline(INT X0,INT Y0,INT X1,INT Y1)
{
	INT x, y, wx, wy;
	INT d1, d2, di_d();
	INT bearing = 0;
	if(dc_linewidth == 0 || dc_shadeoff ==1 || dc_herelinewidth == 0)
		dv_rliner(X0,Y0,X1,Y1);
	else {
		if(X1 > X0)
			bearing |= 01;
		if(Y1 > Y0)
			bearing |= 010;
		if(X1 == X0 )
			bearing |= 0100;
		if(Y1 == Y0)
			bearing |= 01000;
		wx = dc_xlinewidth/2;
		wy = dc_ylinewidth/2;
		y =   0;
		x = -wx;
		dv_rliner(X0+x,Y0+y,X1+x,Y1+y);
		dv_rliner(X0-x,Y0-y,X1-x,Y1-y);
		while(x < 0){
			d1 = do_d(x+1, y,wx,wy);
			d2 = do_d(x, y+1,wx,wy);
			if( d1 < d2)
				x++;
			else
				y++;
			/* recall x is negative and y positive */
			/* 11 to NE, 10 to NW, 00 to SW, 01 to SE */
			if(bearing == 00){ 
				dv_rliner(X0+x,Y0+y,X1+x,Y1+y);
				dv_rliner(X0-x,Y0-y,X1-x,Y1-y);
				dv_rliner(X0+wx+x,Y0+wy-y,X1+x,Y1-y);
				dv_rliner(X0-x,Y0+y,X1-wx-x,Y1-wy+y);
			} else if(bearing == 01){
				dv_rliner(X0+x,Y0-y,X1+x,Y1-y);
				dv_rliner(X0-x,Y0+y,X1-x,Y1+y);
				dv_rliner(X0+x,Y0+y,X1+wx+x,Y1-wy+y);
				dv_rliner(X0+wx+x,Y0+wy-y,X1-x,Y1-y);
			} else if(bearing == 010){
				dv_rliner(X0+x,Y0-y,X1+x,Y1-y);
				dv_rliner(X0-x,Y0+y,X1-x,Y1+y);
				dv_rliner(X0-x,Y0-y,X1+wx+x,Y1+wy-y);
				dv_rliner(X0+wx+x,Y0-wy+y,X1+x,Y1+y);
			} else if(bearing == 011){ 
				dv_rliner(X0+x,Y0+y,X1+x,Y1+y);
				dv_rliner(X0-x,Y0-y,X1-x,Y1-y);
				dv_rliner(X0+x,Y0-y,X1+wx+x,Y1+wy-y);
				dv_rliner(X0-wx-x,Y0-wy+y,X1-x,Y1+y);
			/* order is important */
			} else if(bearing == 01100){ /* point */
				dv_rliner(X0+x,Y0+y,X1-x,Y0+y);
				dv_rliner(X0+x,Y0-y,X1-x,Y0-y);
			} else if(bearing & 0100){
				if(bearing & 010){
					dv_rliner(X0+x,Y0-y,X1+x,Y1+y);
					dv_rliner(X0-x,Y0-y,X1-x,Y1+y);
				} else {
					dv_rliner(X0+x,Y0+y,X1+x,Y1-y);
					dv_rliner(X0-x,Y0+y,X1-x,Y1-y);
				}
			} else if(bearing & 01000){
				if(bearing &01){
					dv_rliner(X0+x,Y0+y,X1-x,Y0+y);
					dv_rliner(X0+x,Y0-y,X1-x,Y0-y);
				} else {
					dv_rliner(X0-x,Y0+y,X1+x,Y1+y);
					dv_rliner(X0-x,Y0-y,X1+x,Y1-y);
				}
			}
		}
		
	}
}

static INT do_d(INT x,INT y,INT wx,INT wy)
{
	double tmp;
	double dx, dy, dwx,dwy;
	dx = x; dy = y ; dwx = wx ; dwy = wy;
	
	tmp = dwy*dwy*dx*dx + dwx*dwx*dy*dy - dwx*dwx*dwy*dwy;
	return( (INT) ABS(tmp) );
}
