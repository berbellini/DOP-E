struct msg_ {
	char title[132];
	char geog[132];
	char auth[5];		/* event location information */
	char evdate[11]; 	/* event location information */
	char evtime[11];	/* event location information */
	float evla;		/* event location information */
	float evlo;		/* event location information */
	float depth;		/* event location information */
	float mag1;		/* event location information */
	float Mw;		/* moment tensor Mw */
	float Mo;		/* moment tensor Mo */
	float evdp;		/* moment tensor depth */
	float strike;		/* moment tensor strike */
	float dip;
	float rake;
	float stkt;
	float plnt;
	float stkp;
	float plnp;
	float stkn;
	float plnn;
	float stk0;
	float dip0;
	float rak0;
	float stk1;
	float dip1;
	float rak1;
	float mar[3][3];	/* moment tensor in Aki-Richards format */
	float mhrv[3][3];	/* moment tensor in CMT format */
	
};
