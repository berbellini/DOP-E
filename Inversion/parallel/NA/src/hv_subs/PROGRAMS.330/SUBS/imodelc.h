/* include file for igetmod and iputmod */

#define NL 200

struct isomod {
	float d[NL];
	float a[NL];
	float b[NL];
	float rho[NL];
	float qa[NL];
	float qb[NL];
	float etap[NL];
	float etas[NL];
	float frefp[NL];
	float frefs[NL];
	float refdep;
} ;

#define MODEL_ISO  0
#define MODEL_TI   1
#define MODEL_ANI  2
#define MODEL_KGS  0
#define MODEL_FLAT 0
#define MODEL_SPH  1
#define MODEL_1D   1
#define MODEL_2D   2
#define MODEL_3D   3
#define MODEL_CVEL 0
#define MODEL_VVEL 1


void getmod(char *mname, int *mmax,char *title,int *iunit,int *iiso,int *iflsph,
     int *idimen,int *icnvel,int *ierr, int listmd, struct isomod *ismod);

void putmod(char *mname, int mmax,char *title,int iunit,int iiso,int iflsph,
     int idimen,int icnvel, int listmd, struct isomod ismod);

void getdval(float depth, int mmax, float *a, float *b, float *rho, float *qa, float *qb,
	float *etap, float *etas, float *frefp, float *frefs, int *ret, struct isomod *ismod);
