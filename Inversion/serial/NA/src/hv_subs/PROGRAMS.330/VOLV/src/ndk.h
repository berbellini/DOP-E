struct ndk_ {
	/* line 1 */
	char auth[5];
	char evdate[11];
	char evtime[11];
	float evla; 
	float evlo; 
	float depth; 
	float mag1;
	float mag2;
	char geog[25];
	/* line 2 */
	char cmtname[17];
	int bsta; 
	int  bcomp;
	float bper;
	int ssta; 
	int  scomp;
	float sper;
	int msta; 
	int  mcomp;
	float mper;
	int cmpsrc;
	char momratefunc[7];
	float halfdur;
	/* line 3*/
	float cmttime;
	float  cmttimeerr;
	float clat;
	float claterr;
	float clon;
	float clonerr;
	float cdepth;
	float cdeptherr; 
	char depthtype[5];
	char timestamp[17];
	/* line 4*/
	int exponent; 
	float Mrr; 
	float EMrr; 
	float Mtt; 
	float EMtt; 
	float Mpp; 
	float EMpp; 
	float Mrt; 
	float EMrt; 
	float Mrp; 
	float EMrp; 
	float Mtp; 
	float EMtp;
	/* line 5*/
	char version[4];
	float ev1; 
	float pl1; 
	float az1; 
	float ev2; 
	float pl2; 
	float az2; 
	float ev3; 
	float pl3; 
	float az3; 
	float m0;
	float strike1; 
	float dip1; 
	float rake1; 
	float strike2; 
	float dip2; 
	float rake2; 
};
