/* CHANGES:
	26 MAY 2010 - igetmod.c created

        --------------------------------------------
	To call this, do this in the main program:

	#include "imodel.h"
	char mname[100];
	char title[100];
	struct isomod isomod;
	main()
	{

        int mmax, iunit, iiso, iflsph, idimen, icnvel, ierr, listmd;

        strcat(mname,"CUS.mod");
        listmd = 1;
        getmod(mname,&mmax,title,&iunit,&iiso,&iflsph,
           &idimen,&icnvel,&ierr,listmd,&isomod);
	for(i=0 ; i < mmax; i++){
		printf("%9.4f\n",isomod.a[i]);
        --------------------------------------------
*/




#include <stdio.h>
#include <string.h>

#include "imodelc.h"


#define LSTR 132
static char instr[LSTR];


void getmod(char *mname, int* mmax,char *title,int *iunit,int *iiso,int *iflsph,
     int *idimen,int *icnvel,int *ierr, int listmd, struct isomod  *isomod)
{
	
/*
c-----
c       HISTORY
c
c       09 08 2000  gave ierr an initial default value for g77
c       01 13 2001  put in close(lun) if file is not model file
c       03 MAY 2002     Modify to permit read from standard input
c       06 JUL 2005 moved inquire to permit use of STDIN
c
c	25 MAY 2010 Fortran version translated to this C version
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
c       Line 04: Model Units, First character 
c           is length (k for kilometer
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
c       rlun    I*4 - logical unit for reading model file. This
c                 unit is released after the use of this routine
c       mname   C*(*)   - model name - if this is stdin or 
c           STDIN just read
c                 from standard input
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
c       ierr    I*4 - 0 model file correctly read in
c               - -1 file does not exist
c               - -2 file is not a model file
c                 -3 error in the model file
c       listmd  L   - .true. list the model
c------
*/

	FILE *fmodel;

	float x0, x1, x2, x3, x4, x5,x6,x7,x8,x9;

        int  j, i, irefdp ;

	/* test to see of the filename is stdin, or if not, if
		the file exists for reading */
        *ierr = 0 ;
	if(strncmp(mname,"stdin",5) == 0 ){
		fmodel = stdin;
	} else {
		if((fmodel = fopen(mname, "r") ) == NULL){
			*ierr = -1;
			return ;
		};
	};
	/* read in the model */
	/* read LINE 01 */
	fgets(instr, LSTR, fmodel);
	if(strncmp(instr,"MODEL",5) != 0 && strncmp(instr,"model",5) != 5){
		*ierr = -2 ;
		if(fmodel != stdin)fclose(fmodel);
		return;
	}

	/* read LINE 02 */
	fgets(instr, LSTR, fmodel);
	strcpy(title,instr);

	/* read LINE 03 */
	fgets(instr, LSTR, fmodel);
	*iiso = -1;
	if(strncmp(instr,"ISO",3) == 0 || strncmp(instr,"iso",3) == 0){
		*iiso = MODEL_ISO;
	} else if(strncmp(instr,"TRA",3) == 0 || strncmp(instr,"tra",3) == 0){
		*iiso = MODEL_TI ;
	} else if(strncmp(instr,"ANI",3) == 0 || strncmp(instr,"ani",3) == 0){
		*iiso = MODEL_ANI;
	}

	/* read LINE 04 */
	fgets(instr, LSTR, fmodel);
	if(strncmp(instr,"KGS",3) == 0 || strncmp(instr,"kgs",3) == 0){
		*iunit = MODEL_KGS;
	}

	/* read LINE 05 */
	fgets(instr, LSTR, fmodel);
	if(strncmp(instr,"FLA",3) == 0 || strncmp(instr,"fla",3) == 0){
		*iflsph = MODEL_FLAT;
	} else if(strncmp(instr,"SPH",3) == 0 || strncmp(instr,"sph",3) == 0){
		*iflsph = MODEL_SPH;
	}

	/* read LINE 06 */
	fgets(instr, LSTR, fmodel);
	if(strncmp(instr,"1-D",3) == 0 || strncmp(instr,"1-d",3) == 0){
		*idimen = MODEL_1D;
	} else if(strncmp(instr,"2-D",3) == 0 || strncmp(instr,"2-d",3) == 0){
		*idimen = MODEL_2D;
	} else if(strncmp(instr,"3-D",3) == 0 || strncmp(instr,"3-d",3) == 0){
		*idimen = MODEL_3D;
	}

	/* read LINE 07 */
	fgets(instr, LSTR, fmodel);
	if(strncmp(instr,"CON",3) == 0 || strncmp(instr,"con",3) == 0){
		*icnvel = MODEL_CVEL;
	} else if(strncmp(instr,"VAR",3) == 0 || strncmp(instr,"var",3) == 0){
		*icnvel = MODEL_VVEL;
	}

	/* skip ines 8 through 11 */
	for(i=8 ; i <= 11 ; i++)
		fgets(instr, LSTR, fmodel);

	/* now the 1-D model */
	/* read the comment line */
	fgets(instr, LSTR, fmodel);

	/* read in the model */

	j = 0;
	isomod->refdep = 0.0 ;
	irefdp = -1;

        *mmax = 0 ;
        if(*iiso == 0){
		while(fgets(instr,LSTR,fmodel) != NULL){
            		j = *mmax ;
		sscanf(instr,"%g %g %g %g %g %g %g %g %g %g",&x0,&x1,&x2,&x3,&x4,&x5,&x6,&x7,&x8,&x9);
			isomod->d[j] = x0;
			isomod->a[j] = x1;
			isomod->b[j] = x2;
			isomod->rho[j] = x3;
			isomod->qa[j] = x4;
			isomod->qb[j] = x5;
			isomod->etap[j] = x6;
			isomod->etas[j] = x7;
			isomod->frefp[j] = x8;
			isomod->frefs[j] = x9;
		if(isomod->d[j] < 0 ){
			isomod->d[j] = -isomod->d[j];
                   	isomod->refdep +=  isomod->d[j];
                    	irefdp = j;
		}
			j++;
            		*mmax = j;
		}
	if( *mmax > 0 ){
            if(listmd != 0){
            	*ierr = 0 ;
		printf(" LAYER             H      P-VEL     S-VEL   DENSITY  \n");
		for(i= 0 ; i < *mmax ; i++){
		printf(" %5d     %10.3f%10.3f%10.3f%10.3f\n",i+1,isomod->d[i],isomod->a[i],isomod->b[i],isomod->rho[i]);
		if(i == irefdp )printf(" -SURFACE - - - - - - - - - - - - - - - - - - - - \n");
		}
            }
        	} else { 
            		*ierr = -3 ;
			fprintf(stderr,"Error in model file\n");
		}
	}
	/* beware this closes stdin */
	
	if(fmodel != stdin)fclose(fmodel);

}

void getdval(float depth, int mmax, float *a, float *b, float *rho, float *qa, float *qb,
        float *etap, float *etas, float *frefp, float *frefs, int *ret, struct isomod *ismod)
{
	/* search through the layers to find the material properties for a
		given depth - here the depth is absolute, e.g.,
		a negative value indicates a point above
	*/
	int i, k;
	float dmin, dmax;
	dmin = - ismod->refdep;
	if(mmax == 1){
		k = 0;
	} else {
		if( depth < dmin){
			k = 0;
		} else {
			for(i = 0 ; i < mmax ; i ++){
				dmax = dmin + ismod->d[i];
				if(depth >= dmin && depth < dmax)
					k = i;
				dmin = dmax;
			}
		}
	}
	if(depth >= dmax){
		k = mmax -1;	/* assume bottom layer is halfspace */
	}
	*a = ismod->a[k];
	*b = ismod->b[k];
	*rho = ismod->rho[k];
	*qa = ismod->qa[k];
	*qb = ismod->qb[k];
	*frefp = ismod->frefp[k];
	*frefs = ismod->frefs[k];
	*etap = ismod->etap[k];
	*etas = ismod->etas[k];
}
