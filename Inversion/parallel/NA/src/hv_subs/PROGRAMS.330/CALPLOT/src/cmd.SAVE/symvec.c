/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: SYMVEC                                                c
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
#include	<math.h>


#include <stdio.h>
#ifdef MSDOS
#define INT long
#else
#define INT int
#endif

extern INT dc_curfont;

void dv_symvec(INT xloc,INT yloc,INT ht,char *s,INT ang,INT nochar);
static void scont(INT x,INT y);
static void smove(INT x,INT y);
static void symloc(INT jndex,double *scal,double *shifx,double *shify);
static void symscl(double x1,double y1,double *x,double *y,
	int *ipen,double angle,double sinth,double costh,
	double height,double scal,double shx,double shy);
static void symrot(double *x,double *y, double s, double c);

/* mapping of character into table entry */
/* SYMBOLS and ASCII */
static int atab[] = {
   		  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,
  		 10,  11,  12,  13,  14,  15,  16,  17,  18,  19,
  		 20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
  		 30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
  		 40,  41,  42,  43,  44,  45,  46,  47,  48,  49,
  		 50,  51,  52,  53,  54,  55,  56,  57,  58,  59,
  		 60,  61,  62,  63,  64,  65,  66,  67,  68,  69,
  		 70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
  		 80,  81,  82,  83,  84,  85,  86,  87,  88,  89,
  		 90,  91,  92,  93,  94,  95,  96,  97,  98,  99,
 		100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
 		110,  16 
} ;

/* SYMBOLS and GREEK */
static int gtab[] = {
   		  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,
  		 10,  11,  12,  13,  14,  15,  16,  17, 111,  19,
  		112,  21,  22, 113,  24,  25,  26,  27,  28,  29,
  		 30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
  		 40,  41,  42,  43,  44,  45,  46,  47, 114,  49,
  		 50,  72, 115,  53, 116, 117,  56,  57, 118,  59,
  		119,  61,  62,  63, 120, 121,  64, 122,  68,  73,
  		123, 124, 125, 126,  74,  75, 127,  77, 128,  79,
  		 80, 129, 130, 131, 132, 133, 134, 135, 136, 137,
  		138, 139, 140, 141, 142,  95, 143, 144, 145, 146,
 		147, 148, 149, 150, 151, 152, 153, 107, 108, 109,
 		110,  16 
} ;
/* symbol from the CALPLOT symbol routine, but in C */
static int chartb[][9] = {
	 822,4240,   4,4442,2200,   0,   0,   0,   0,		/* 0 square */
	1222,4241,3010, 103,1434,4342,2200,   0,   0,		/* octagon */
	 622,4210,1442,2200,   0,   0,   0,   0,   0,		/* triangle */
	 722,2420,2202,4222,   0,   0,   0,   0,   0,		/* plus */
	 722,4004,2200,4422,   0,   0,   0,   0,   0,		/* X */
	 722,4220, 224,4222,   0,   0,   0,   0,   0,		/* diamond */
	 722, 242,2024,4222,   0,   0,   0,   0,   0,		/* hat */
	 622,  44,4004,2200,   0,   0,   0,   0,   0,		/* table */
	1122,   4,  44,4044,2221,2322,   0,   0,   0,		/* Z */
	 722,4422,4022, 222,   0,   0,   0,   0,   0,		/* Y */
	1422,3331,4031,1100,1113, 413,3344,2200,   0,		/* sq-X */
	1422,4202,2220,2422,4004,2200,4422,   0,   0,		/* asterisk */
	 722,  44,4004,  22,   0,   0,   0,   0,   0,		/* hr-glass */
	 422,4202,2200,   0,   0,   0,   0,   0,   0,		/* vertlin */
	1122,4210,1442,3230, 234,3222,   0,   0,   0,		/* star */
	 322,2422,   0,   0,   0,   0,   0,   0,   0,		/* dash */
	   0,   0,   0,   0,   0,   0,   0,   0,   0,		/* space */
	 462,3217, 200,   0,   0,   0,   0,   0,   0,		/* ! */
	 461,4148,6300,   0,   0,   0,   0,   0,   0,		/* " */
	 861,   7,6345,4425,2400,   0,   0,   0,   0,		/* # */
	1454,6361,5040,3133,2414, 301,1007,6200,   0,		/* $ */
	1260,5051,6160,1814, 403,1305,6400,   0,   0,		/* % */
	1004,4050,6152,2010, 102,2400,   0,   0,   0,		/* & */
	 463,6242,6300,   0,   0,   0,   0,   0,   0,		/* ' */
	 664,6352,1203, 400,   0,   0,   0,   0,   0,		/* ( */
	 660,6152,1201,   0,   0,   0,   0,   0,   0,		/* ) */
	 812,5255,1439,3015,5400,   0,   0,   0,   0,		/* * */
	 412,5235,3400,   0,   0,   0,   0,   0,   0,		/* + */
	 422,2302,2200,   0,   0,   0,   0,   0,   0,		/* , */
	 221,2400,   0,   0,   0,   0,   0,   0,   0,		/* - */
	 510, 112,2110,   0,   0,   0,   0,   0,   0,		/* . */
	 200,6400,   0,   0,   0,   0,   0,   0,   0,		/* / */
	 963,6150,1001, 314,5463,   0,   0,   0,   0,		/* 0 */
	 501, 302,6251,   0,   0,   0,   0,   0,   0,		/* 1 */
	 904,  21,3344,5463,6150,   0,   0,   0,   0,		/* 2 */
	 910, 103,1434,4337,6460,   0,   0,   0,   0,		/* 3 */
	 561,3034,5803,   0,   0,   0,   0,   0,   0,		/* 4 */
	 910, 103,1434,4340,6064,   0,   0,   0,   0,		/* 5 */
	 962,2010, 103,1424,3331,   0,   0,   0,   0,		/* 6 */
	 750,6064,5432,1101,   0,   0,   0,   0,   0,		/* 7 */
	1631,4050,6163,5444,3331,2010, 103,1424,3300,		/* 8 */
	 901,4454,6361,5040,3133,   0,   0,   0,   0,		/* 9 */
	1043,4232,3343,2712,1323,2200,   0,   0,   0,		/* : */
	 943,4232,3343,2702,2322,   0,   0,   0,   0,		/* ; */
	 354,3014,   0,   0,   0,   0,   0,   0,   0,		/* < */
	 441,4429,2100,   0,   0,   0,   0,   0,   0,		/* = */
	 350,3410,   0,   0,   0,   0,   0,   0,   0,		/* > */
	1202, 313,1202,2744,5463,6150,4000,   0,   0,		/* ? */
	1623,4352,4131,2223,3454,6361,5020, 103,1400,		/* @ */
	 500,6204,2623,   0,   0,   0,   0,   0,   0,		/* A */
	1300, 314,2433,3133,4454,6360,6101,   0,   0,		/* B */
	 854,6361,5010, 103,1400,   0,   0,   0,   0,		/* C */
	 800, 314,5463,6061, 100,   0,   0,   0,   0,		/* D */
	 764,6030,3330,   4,   0,   0,   0,   0,   0,		/* E */
	 664,6030,3330,   0,   0,   0,   0,   0,   0,		/* F */
	1054,6361,5010, 103,1434,3300,   0,   0,   0,		/* G */
	 600,6030,3464, 400,   0,   0,   0,   0,   0,		/* H */
	 601, 302,6261,6300,   0,   0,   0,   0,   0,		/* I */
	 510, 103,1464,   0,   0,   0,   0,   0,   0,		/* J */
	 600,6020,6442, 400,   0,   0,   0,   0,   0,		/* K */
	 360,   4,   0,   0,   0,   0,   0,   0,   0,		/* L */
	 500,6032,6404,   0,   0,   0,   0,   0,   0,		/* M */
	 400,6004,6400,   0,   0,   0,   0,   0,   0,		/* N */
	 910,5061,6354,1403, 110,   0,   0,   0,   0,		/* O */
	 801,6160,6354,4433,3100,   0,   0,   0,   0,		/* P */
	1003, 110,5061,6354,1427, 400,   0,   0,   0,		/* Q */
	1001,6160,6354,4433,3133, 400,   0,   0,   0,		/* R */
	1254,6361,5040,3133,2414, 301,1000,   0,   0,		/* S */
	 460,6462, 200,   0,   0,   0,   0,   0,   0,		/* T */
	 664,1403, 110,6000,   0,   0,   0,   0,   0,		/* U */
	 364, 260,   0,   0,   0,   0,   0,   0,   0,		/* V */
	 564, 332, 160,   0,   0,   0,   0,   0,   0,		/* W */
	 460, 405,6400,   0,   0,   0,   0,   0,   0,		/* X */
	 564,3202,3260,   0,   0,   0,   0,   0,   0,		/* Y */
	 404,  64,6000,   0,   0,   0,   0,   0,   0,		/* Z */
	 462,6000, 200,   0,   0,   0,   0,   0,   0,		/* [ */
	 260, 400,   0,   0,   0,   0,   0,   0,   0,		/* \ */
	 462,6404, 200,   0,   0,   0,   0,   0,   0,		/* ] */
	 320,5224,   0,   0,   0,   0,   0,   0,   0,		/* ^ */
	 200, 400,   0,   0,   0,   0,   0,   0,   0,		/* _ */
	 461,6242,6100,   0,   0,   0,   0,   0,   0,		/* ` */
	1041,4334, 419, 301,1021,2400,   0,   0,   0,		/* a */
	1061, 126,3233,2414, 302,1100,   0,   0,   0,		/* b */
	 644,4130,1001, 400,   0,   0,   0,   0,   0,		/* c */
	1064, 424,3332,2111, 203,1400,   0,   0,   0,		/* d */
	 920,2434,4341,3010, 103,   0,   0,   0,   0,		/* e */
	 702,5263,5444,3633,   0,   0,   0,   0,   0,		/* f */
	1362,6354,4433,3241,5162,4914, 301,   0,   0,		/* g */
	 761, 126,3233,2404,   0,   0,   0,   0,   0,		/* h */
	 771,6145,4101, 502,   0,   0,   0,   0,   0,		/* i */
	 672,6215, 112,5200,   0,   0,   0,   0,   0,		/* j */
	 661, 136, 426,4300,   0,   0,   0,   0,   0,		/* k */
	 561,6202, 603,   0,   0,   0,   0,   0,   0,		/* l */
	1030,  25,3122, 227,3324, 400,   0,   0,   0,		/* m */
	 731, 126,3233,2404,   0,   0,   0,   0,   0,		/* n */
	 932,3324,1403, 211,2132,   0,   0,   0,   0,		/* o */
	 861, 166,6354,4433,3100,   0,   0,   0,   0,		/* p */
	 963, 304,3831,4050,6163,   0,   0,   0,   0,		/* q */
	 631, 126,3233,2400,   0,   0,   0,   0,   0,		/* r */
	 843,4130,2123,1403, 100,   0,   0,   0,   0,		/* s */
	 662,1203,1446,4300,   0,   0,   0,   0,   0,		/* t */
	 734, 436,1102, 314,   0,   0,   0,   0,   0,		/* u */
	 340, 244,   0,   0,   0,   0,   0,   0,   0,		/* v */
	 930,1001,1232,1703,1434,   0,   0,   0,   0,		/* w */
	 431, 406,3400,   0,   0,   0,   0,   0,   0,		/* x */
	 944,3332,4161,6914, 302,   0,   0,   0,   0,		/* y */
	 431,3401, 400,   0,   0,   0,   0,   0,   0,		/* z */
	 964,6352,4231,2212, 304,   0,   0,   0,   0,		/* { */
	 462,4227, 200,   0,   0,   0,   0,   0,   0,		/* | */
	 960,6152,4233,2212, 100,   0,   0,   0,   0,		/* } */
	 420,3123,3400,   0,   0,   0,   0,   0,   0,		/* 110 ~ */
	/* special greek symbols */
	 560, 264,3633,   0,   0,   0,   0,   0,   0,		/* 111-18 for all */
	 760,6434,3034, 400,   0,   0,   0,   0,   0,		/* 112-20 existance */
	1160,6253,4434,3134,2413, 200,   0,   0,   0,		/* 113-23 member of */ 
	 800, 325,2345,5142,5300,   0,   0,   0,   0,		/* 114-48 appr equal */
	 400,6204,0000,   0,   0,   0,   0,   0,   0,		/* 115-52 DELTA */
	1421,3041,4334,2321,3066,6362, 201, 300,   0,		/* 116-54 PHI */
	 400,6064,5400,   0,   0,   0,   0,   0,   0,		/* 117-55 GAMMA */
	1030,4111, 203,1454,6352,4400,   0,   0,   0,		/* 118-58 theta */ 
	 300,6204,   0,   0,   0,   0,   0,   0,   0,		/* 119-60 LAMBDA */
	 400,6064, 400,   0,   0,   0,   0,   0,   0,		/* 120-64 PI */
	1110,5061,6354,1403, 110,3633,   0,   0,   0,		/* 121-65 THETA */
	 714, 400,3260,6454,   0,   0,   0,   0,   0,		/* 122-67 SIGMA */
	 963,6150,3021,2314, 302,   0,   0,   0,   0,		/* 123-70 sigma1 */
	1200, 111,2050,6163,5424,1303, 400,   0,   0,		/* 124-71 OMEGA */
	1010,0004,1455,6064,5436,3300,   0,   0,   0,		/* 125-72 XI */
	1350,3021,2334,5467,6163,6202, 103,   0,   0,		/* 126-73 PSI */
	1500,1011, 100, 813,1404, 337,4243,3332,   0,		/* 127-76 therefore */
	 500, 252, 204,   0,   0,   0,   0,   0,   0,		/* 128-78 perpendicular to */
	1144,2304,2302, 110,3041,4223,   0,   0,   0,		/* 129-81 alpha */
	1501,2112,1324,3443,4243,5464,7372,6101,   0,		/* 130-82 beta */
	 650,6103,1405,6400,   0,   0,   0,   0,   0,		/* 131-83 chi */
	1353,6261,5041,3324,1403, 110,2032,   0,   0,		/* 132-84 delta */
	1043,4130,2123,2110, 103,1400,   0,   0,   0,		/* 133-85 epsilon */
	1221,3040,5153,4434,2321,2202,7200,   0,   0,		/* 134-86 phi 1 */
	 750,6162,3303,3364,   0,   0,   0,   0,   0,		/* 135-87 gamma */
	 760,6121,6163,5404,   0,   0,   0,   0,   0,		/* 136-88 eta */
	 431,4202,1300,   0,   0,   0,   0,   0,   0,		/* 137-89 iota */
	1002,6263,5434,2321,3050,6100,   0,   0,   0,		/* 138-90 phi 2 */
	 900,5030,3153,5409, 331,   0,   0,   0,   0,		/* 139-91 kappa */
	 760,6103, 403,3200,   0,   0,   0,   0,   0,		/* 140-92 lambda */
	1000,6030,2122,3323,2423,6300,   0,   0,   0,		/* 141-93 mu */
	 540,4102,3444,   0,   0,   0,   0,   0,   0,		/* 142-94 nu */
	 714, 343,4440,4101,   0,   0,   0,   0,   0,		/* 143-96 pi */
	1130,3353,6261,5010, 102,1333,   0,   0,   0,		/* 144-97 theta */
	 900,5061,6253,3322,2130,   0,   0,   0,   0,		/* 145-98 rho */
	 944,4130,1001, 213,3342,   0,   0,   0,   0,		/* 146-99 sigma */
	 630,4144,4202,1300,   0,   0,   0,   0,   0,		/* 147-100 tau */
	 830,4111, 203,1434,4300,   0,   0,   0,   0,		/* 148-101 upsilon */
	1341,3010, 112,3212, 314,3443,4440,   0,   0,		/* 149-102 omega 1 */
	1141,3010, 112,3212, 314,3443,   0,   0,   0,		/* 150-103 omega 2 */
	1471,6263,6251,4243,4130,2011,1303, 200,   0,		/* 151-104 xi */
	 850,3021,2334,5467, 200,   0,   0,   0,   0,		/* 152-105 psi */
	1072,6364,6331,2112,1404, 300,   0,   0,   0 		/* 153-106 zeta */
} ;

extern INT dc_sdx;

void dv_symvec(INT xloc,INT yloc,INT ht,char *s,INT ang,INT nochar)
{

	double x0, y0;
	double rangle, scal, shifx, shify, shx, shy;
	double angle, height;
	double sinth, costh;
	double x1, y1, x, y, x2, y2;
	int ipen, index, jndex, uplim, j, l, m, loops;



	x0 = xloc;
	y0 = yloc;
	height = (double)ht / 6.0;
	angle = (double)ang;
	rangle = -(angle +0.01)/57.29578;
	sinth = sin(rangle);
	costh = cos(rangle);

	ipen = 3;
	if(nochar < 0){
		if(nochar < -1)
			di_cont((INT)xloc,(INT)yloc);
		uplim = 1;
	} else {
		uplim = nochar;
	}
	for(j=0;j<uplim;j++){
		jndex = (char)s[j];
		scal = 1.0;
		shifx = 0.0;
		shify = 0.0;
		/* obtain index into character generation table
			as well as symbol positioning and scaling 
			The character is assumed to be ASCII, e.g.,
			that 32 = space, 127 = not printable  */	
		if(nochar < 0){
			if(dc_curfont == 4){
				index = gtab[jndex];
			} else {
				index = atab[jndex];
			}
			if(jndex < 16){	
				symloc(jndex   ,&scal,&shifx,&shify);
			} else {
				symloc(jndex+16,&scal,&shifx,&shify);
			}
		} else {
			if(jndex < 040 || jndex > 0176)
				index = 16;	/* force space */
			else {
				if(dc_curfont == 4)
					index = gtab[jndex - 16];
				else
					index = atab[jndex - 16];
			}
			symloc(jndex   ,&scal,&shifx,&shify);
		}

		/* plot the symbol */
		if(index == 16){	/* force a space */
			if(nochar > -1){
				y0=y0-6.0*height*sinth;
				x0=x0+6.0*height*costh;
				smove((INT)x0,(INT)y0);
			}
		} else {	/* parse table and plot symbol */
		/*
       		the first two digits  of the current character's first 
       		four digit integer are now isolated and used to determine 
       		the number of loops which must be done to plot all of the 
       		vectors which make up the character.  this is accomplished 
       		by allowing an integer division by 100 to truncate off the 
       		following two digits. 
  								*/
			loops=chartb[index][0]/100 ;
		/*
		all other digits are yx values, if x >=5 pen up
		with x coordinate given by x mod 5 (or x-5 here )
       		the coordinates of the first point are now obtained. 
       		they are then enlarged by the scaling factor, "height", 
       		rotated if the angle from horizontal is noticeably 
       		large, and moved to via dark vector since it is the 
       		starting point of the character. 
 		*/
		/* 	Initialize shift */
			shy = shify*height ;
			shx = shifx*height ;
			x1=(double)(chartb[index][0]-chartb[index][0]/10*10) ;
			y1=(double)((chartb[index][0]-chartb[index][0]/100*100)/10) ;
			symscl(x1,y1,&x,&y,&ipen,angle,sinth,costh,
				height,scal,shx,shy);
			smove((INT)(x+x0),(INT)(y+y0));
			
		/*
       		looping is initiated to treat the subsequent vectors of 
       		the current character plot. 
 
 
       		the particular four digit integer of the current character 
       		to be treated is found here by "m" and treated in a way 
       		basically like the first one was:,,the digits are isolated 
       		and tested to determine whether dark or not, then scaled to 
       		"height", rotated if necessary, and plotted in the desig- 
       		nated mode. 
 
 		*/
			for(m=1;m<loops;m+=2){
				l=m/2+1 ;
				y1=(chartb[index][l]/1000) ;
				x1=(chartb[index][l]/100-y1*10);
				symscl(x1,y1,&x,&y,&ipen,angle,sinth,costh,
					height,scal,shx,shy);
				if(ipen==2)
					scont((INT)(x+x0),(INT)(y+y0));
				else
					smove((INT)(x+x0),(INT)(y+y0));
				if(m != (loops-1) ){
					y2=(chartb[index][l]/10-y1*100-x1*10) ;
					x2=(chartb[index][l]-y1*1000-x1*100-y2*10) ;
				symscl(x2,y2,&x,&y,&ipen,angle,sinth,costh,
					height,scal,shx,shy);
				if(ipen==2)
					scont((INT)(x+x0),(INT)(y+y0));
				else
					smove((INT)(x+x0),(INT)(y+y0));
				}	
			} /* end plot character */
			if(nochar > -1){
				y0=y0-6.0*height*sinth;
				x0=x0+6.0*height*costh;
				smove((INT)x0,(INT)y0);
			}
		} /* end plot non-blank */
	} /* end plot string */
}

static void symrot(double *x,double *y, double s, double c)
{
	double x1, y1;
	x1 = *x;
	y1 = *y;
	*x =  c * x1 + s * y1;
	*y = -s * x1 + c * y1;
}

struct pos {
	INT shift;
	INT scale;
	float shify;
};

static struct pos ascpos[] = {
	{ 0147,	0141,	2.33 },
	{ 0152,	0143,	2.33 },
	{ 0160, 0145,	3.00 },
	{ 0161, 0151,	3.00 },
	{ 0171, 0163,	3.00 },
	{ 0000,	0166,	0.00 },
	{   -1,   -1,	0.00 }
};

static struct pos grkpos[] = {
	{ 86,	107,	2.00 },	/* sigma1 */
	{ 98,	0,	1.00 },	/* beta */
	{ 99,	0,	2.00 },	/* chi */
	{102,	0,	2.00 },	/* phi1 */
	{103,	0,	2.00 },	/* gamma */
	{104,	0,	2.00 },	/* eta */
	{106,	0,	2.00 },	/* phi2 */
	{109,	0,	2.00 },	/* mu */
	{114,	0,	2.00 },	/* rho */
	{121,	0,	2.00 },	/* psi */
	{120,	0,	1.00 }, /* xi */
	{122,	0,	1.00 }, /* zeta */
	{ -1,  -1,      0.0   }
};
/*	this subroutine determines if the particular character
	letter needs to be shifted or scaled down		*/
static void symloc(INT jndex,double *scal,double *shifx,double *shify)
{

	/*
		shift down g,j,p,q,y
		scale down a,c,e,i,s,v
		Note that we use ASCII character equivalent as the test 
 	*/

struct pos *poshere;
	INT i;

	if(dc_curfont == 4){
		poshere = &grkpos[0];
	} else {
		poshere = &ascpos[0];
	}

	if(jndex < 16){
		*scal = 1.5;
		*shifx = 3.0;
		*shify = 3.0;
	} else {
		for(i=0;poshere[i].shift >= 0;i++){
	            if (jndex ==  poshere[i].shift){
			*shify = poshere[i].shify;
			*scal = 1.0;
	                *shifx = 1./7. ;
	            }
	            if (jndex ==  poshere[i].scale) {
	                *scal  =  0.8 ;
	                *shifx = -1.0 ;
	            }
		}
	}
}

static void symscl(double x1,double y1,double *x,double *y,
	int *ipen,double angle,double sinth,double costh,
	double height,double scal,double shx,double shy)
{
	double fabs();
	double fx, fy;
	if(x1 >= 5){
		fx=(x1-5)*height*scal-shx;
		*ipen = 3;
	} else {
		fx=x1*height*scal - shx;
		*ipen = 2;
	}
	fy = y1*height*scal - shy;
	if(dc_curfont == 2){	/* italics */
		fx = 0.3 * fy + fx;
	}
	*x = fx;
	*y = fy;
	if(fabs(angle) > 2.0)
		symrot(x,y,sinth,costh);
}

/* routines to permit a fake bold */
static INT smx, smy;
static void smove(INT x,INT y)
{
	smx = x;
	smy = y;
}

static void scont(INT x,INT y)
{
	if(dc_curfont == 3){
		di_move(smx+dc_sdx,smy);
		di_cont(  x+dc_sdx,  y);
	}
	di_move(smx,smy);
	di_cont(  x,  y);
	smx = x;
	smy = y;
}

