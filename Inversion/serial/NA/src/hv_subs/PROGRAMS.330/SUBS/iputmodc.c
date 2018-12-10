/* CHANGES:
	27 MAY 2010 - created
	04 OCT 2010 - corrected logic in lgstr
		introduced use of isblank()
*/
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "imodelc.h"
#define ABS(a  ) ( (a) >  0  ? (a):-(a))

FILE *fmodel;
char *cmt="      H(KM)   VP(KM/S)   VS(KM/S) RHO(GM/CC)     QP         QS       ETAP       ETAS      FREFP      FREFS\n";

static int lgstr(char *str);

void putmod(char *mname, int mmax,char *title,int iunit,int iiso,int iflsph,
     int idimen,int icnvel, int listmd, struct isomod ismod)
{
/*
c-----
c       CHANGES
c       03 MAY 2002 permit write to standard output
c-----
c       General purpose model input
c       This model specification is designed to be as 
c           general as possible
c
c       Input lines
c       Line 01: MODEL
c       Line 02: Model Name
c       Line 03: ISOTROPIC or ANISOTROPIC or 
c           TRANSVERSELY ANISOTROPIC
c       Line 04: Model Units, First character is length (k for kilometer
c           second is mass (g for gm/cc), third is time (s for time)
c       Line 05: FLAT EARTH or SPHERICAL EARTH
c       Line 06: 1-D, 2-D or 3-D
c       Line 07: CONSTANT VELOCITY
c       Line 08: open for future use
c       Line 09: open for future use
c       Line 10: open for future use
c       Line 11: open for future use
c       Lines 12-end:   These are specific to the model
c           For ISOTROPIC the entries are
c           Layer Thickness, P-velocity, S-velocity, Density, Qp, Qs,
c           Eta-P, Eta S (Eta is frequency dependence), 
c           FreqRefP, FreqRefP
c-----
cMODEL
cTEST MODEL.01
cISOTROPIC
cKGS
cFLAT EARTH
c1-D
cCONSTANT VELOCITY
cLINE08
cLINE09
cLINE10
cLINE11
c H  VP  VS   RHO   QP  QS   ETAP   ETAS REFP  REFS
c1.0    5.0 3.0 2.5 0.0 0.0 0.0 0.0 1.0 1.0
c2.0    5.1 3.1 2.6 0.0 0.0 0.0 0.0 1.0 1.0
c7.0    6.0 3.5 2.8 0.0 0.0 0.0 0.0 1.0 1.0
c10.0   6.5 3.8 2.9 0.0 0.0 0.0 0.0 1.0 1.0
c20.0   7.0 4.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0
c40.0   8.0 4.7 3.3 0.0 0.0 0.0 0.0 1.0 1.0
c-----
c-----
c       wlun    I*4 - logical unit for writing model file. This
c                 unit is released after the use of this routine
c       mname   C*(*)   - model name
c       mmax    I*4 - number of layers in the model, last layer is
c                    halfspace
c       title   C*(*)   - title of the model file
c       iunit   I*4 - 0 Kilometer, Gram, Sec
c       iiso    I*4 - 0 isotropic 
c                 1 transversely anisotropic 
c                 2 general anisotropic 
c       iflsph  I*4 - 0 flat earth model
c                 1 spherical earth model
c       idimen  I*4 - 1 1-D
c               - 2 2-D
c               - 3 3-D
c       icnvel  I*4 - 0 constant velocity
c                 1 variable velocity
c       lverby  L   - .false. quiet output
c------
*/
	int j, ls;
	float curdep,dout;

	if(strncmp(mname,"stdout",6)==0){
		fmodel = stdout;
        } else {
		fmodel = fopen(mname,"w+");
	}

	/* LINE 01 */
	fprintf(fmodel,"MODEL.01\n");

	/* LINE 02 */
	/* do not output excessive blanks */
	ls = lgstr(title);
	if(ls < strlen(title))
		title[ls] = '\0';
	fprintf(fmodel,"%s\n",title);

	/* LINE 03 */
        if(iiso == MODEL_ISO)
            fprintf(fmodel,"ISOTROPIC\n");
        else if(iiso == MODEL_TI )
            fprintf(fmodel,"TRANSVERSE ANISOTROPIC\n");
        else if( iiso == MODEL_ANI )
            fprintf(fmodel,"ANISOTROPIC\n");

	/* LINE 04 */
        fprintf(fmodel,"KGS\n");

	/* LINE 05 */
        if(iflsph == MODEL_FLAT)
            fprintf(fmodel,"FLAT EARTH\n");
        else if(iflsph == MODEL_SPH)
            fprintf(fmodel,"SPHERICAL EARTH\n");

	/* LINE 06 */
	if( idimen == MODEL_1D)
		fprintf(fmodel,"1-D\n");
	else if( idimen == MODEL_2D)
		fprintf(fmodel,"2-D\n");
	else if( idimen == MODEL_3D)
		fprintf(fmodel,"3-D\n");

	/* LINE 07 */
	if(icnvel == MODEL_CVEL)
            fprintf(fmodel,"CONSTANT VELOCITY\n");
	else if(icnvel == MODEL_VVEL)
            fprintf(fmodel,"VARIABLE VELOCITY\n");

	/* put lines 8 through 11 */
        fprintf(fmodel,"LINE08\n");
        fprintf(fmodel,"LINE09\n");
        fprintf(fmodel,"LINE10\n");
        fprintf(fmodel,"LINE11\n");

	/* put model specifically for 1-D flat isotropic */
	/* put out comment line */
	fprintf(fmodel,"%s\n",cmt);

	/* output the actual model */
        
	curdep = 0.0 ;
	for(j=0;j<mmax;j++){
		curdep +=  ABS(ismod.d[j]) ;
		if(curdep <= ismod.refdep)
			dout = - ismod.d[j]  ;
		else
			dout = ismod.d[j];
fprintf(fmodel,"%10.4f %10.4f %10.4f %10.4f %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g\n",
		dout,
		ismod.a[j],
		ismod.b[j],
		ismod.rho[j],
		ismod.qa[j],
		ismod.qb[j],
		ismod.etap[j],
		ismod.etas[j],
		ismod.frefp[j],
		ismod.frefs[j]);
	}
	if(fmodel != stdout)fclose(fmodel);
}

static int lgstr(char *str)
{
	int n, i;
	int lg;
        n = strlen(str) ;
	lg = n;
	i = n-1 ;

        while(isblank(str[i]) && i>=0){
		lg = i;
		i--;
	}
        return (lg);
}
