/* This requires more work so that
 * 1. We test to see if it is a valid SAC file
 * 2. We test to see if it is a BYTE swapped valid SAC file
 */
/* CHANGES
 09 JAN 2005 - added  nerr = -1; in bwsach
 03 MAR 2009 - added logic in brsac and brsach to avoid reading non-sac files
		note this logic is safely duplicated since I usually do a
		gsac_valid_sacfile which calls brsach before one is 
		actually does a brsac 
		Eventually get the file size from the system and check
		file size versus what is estimated from header
 10 SEP 2009 - corrected valid sac file test to include a spectra file, e.g., IXY
		Change
		if(has12345 == 0 || sachdr.ihdr[H_NVHDR] != 6 || sachdr.ihdr[H_IFTYPE] !=ENUM_ITIME){
to
		if(has12345 == 0 || sachdr.ihdr[H_NVHDR] != 6 || (sachdr.ihdr[H_IFTYPE] !=ENUM_ITIME && sachdr.ihdr[H_IFTYPE] !=ENUM_IXY)){
		
*/

#include	<stdio.h>
#include	 "gsac.h"
#include	 "gsac_sac.h"
#include	"csstim.h"

static float *x = (float *)NULL;

/* these are absolutely necessary because the addition to cut in
	absoltue time */
#define CUT_GMT		-20
#define CUT_CAL		-21


int brsac(char *fname,struct  sachdr_ *sachdrret, float **data)
{
	int maxpts, nread;
	int newnpts;		/* for cut new number of points */
	float newb, newe;	/* new values for b and e if cut */
	int newbeg, newend;		/* index of first point in new series */
	float oldb, olde;
	float delta;
	float t0, t1, tv;
	double tzref;
	int i;
	int i1, i2;
	int j1, j2, j3;
	int nerr;
	int has12345;
	int rhas12345, ihas12345, chas12345;
	FILE *fptr;
	struct sachdr_ sachdr;
	/*
		determine if the file exists for reading
	*/
	if((fptr=fopen(fname,"rb")) == NULL){
		sachdr.ihdr[H_NPTS] = 0;
		nerr = -1;
		perror("fopen error in brsac:");
		return(nerr);
	}
#ifdef WIN32
	setmode(fileno(fptr), O_BINARY);
#endif
	/*
		file exists get header information only
	*/
	
	nerr = 0;
	if(fread(&sachdr,sizeof( struct sachdr_),1,fptr)!=1){
		nerr=-3;
		fclose(fptr);
		return(nerr);
	}
	maxpts = sachdr.ihdr[H_NPTS] ;
			nerr = 0 ;

	/* there must be at least one -12345 for this to be a sac file!! */
	has12345 = 0;
	rhas12345 = 0;
	ihas12345 = 0;
	chas12345 = 0;
	for (i=0;i < 70 ; i++){
		if(sachdr.rhdr[i] == -12345.0) {
			has12345++ ; 
			rhas12345++;
		}
	}
	for (i=0;i < 40 ; i++){
		if(sachdr.ihdr[i] == -12345){
			has12345++; 
			ihas12345++;
		}
	}
	for (i=0;i < 24 ; i++){
		if(strncmp(sachdr.chdr[i],  "-12345  ",8)==0){
			has12345++;
			chas12345++;
		}
	}
	if(has12345 == 0 || sachdr.ihdr[H_NVHDR] != 6 || 
		(sachdr.ihdr[H_IFTYPE] !=ENUM_ITIME && sachdr.ihdr[H_IFTYPE] !=ENUM_IXY)){
	if(has12345 == 0)
		fprintf(stderr,"Problem: %s is not a Sac file \n",fname);
	else
		fprintf(stderr,"Problem: %s may be a byte-swapped Sac file. Separately run 'saccvt -I'\n",fname);
		nerr = -3;
		fclose(fptr);
		return(nerr);
	}
	/* get tzref - needed to cut CAL or cut GMT */
	
	if(sachdr.ihdr[H_NZYEAR] == -12345 ||
		sachdr.ihdr[H_NZJDAY] == -12345 ||
		sachdr.ihdr[H_NZHOUR] == -12345 ||
		sachdr.ihdr[H_NZMIN] == -12345 ||
		sachdr.ihdr[H_NZSEC] == -12345 ||
		sachdr.ihdr[H_NZMSEC] == -12345){
		/*
		printf("Using 1970 000 00 00 00 000\n");
		*/
		tzref = 0.0;
	} else {
	htoe1(sachdr.ihdr[H_NZYEAR], 
		sachdr.ihdr[H_NZJDAY], 
		sachdr.ihdr[H_NZHOUR],
		sachdr.ihdr[H_NZMIN],
		sachdr.ihdr[H_NZSEC],
		sachdr.ihdr[H_NZMSEC],
		&tzref);
	}

	if(gsac_control.docut && sachdr.rhdr[gsac_control.cutint[0]]!= -12345.
			&&sachdr.rhdr[gsac_control.cutint[1]]!= -12345.){

		oldb = sachdr.rhdr[H_B];
		olde = sachdr.rhdr[H_E];
		for(i=0;i<2;i++){
			if(gsac_control.cutint[i] == CUT_GMT){
				tv = (float)(gsac_control.cutepoch[i]-tzref);
			} else if(gsac_control.cutint[i] == CUT_CAL){
				tv = (float)(gsac_control.cutepoch[i]-tzref);
			} else {
				tv = sachdr.rhdr[gsac_control.cutint[i]]
					+gsac_control.cutoff[i] - oldb;
			}
		if(i==0)
			t0 = tv;
		else
			t1 = tv;
		}
		delta = sachdr.rhdr[H_DELTA];
		newb = MIN(t0,t1);
		newe = MAX(t0,t1);
		newnpts = 1 + (int)((newe - newb)/delta + 0.49);
		if(newb > 0)
			newbeg =  (int)((newb)/delta +0.49);
		else
			newbeg =  (int)((newb)/delta -0.49);
/*
printf("  CUT oldb %f olde %f   newb %f newe %f\n maxpts %d	newnpts %d newbeg %d\n",oldb,olde,newb,newe,maxpts,newnpts,newbeg);
*/
		newb = oldb + newbeg*delta;
		newe = newb + (newnpts -1)*delta;
		newend = newbeg + newnpts -1 ;
/*
printf("   newb %f newe %f newend %f\n",newb,newe,newend);
*/
	} else {
		oldb = sachdr.rhdr[H_B];
		olde = sachdr.rhdr[H_E];
		newb = oldb;
		newe = olde;
		newbeg = 0;
		newnpts = sachdr.ihdr[H_NPTS] ;
		newend = newbeg + newnpts  ;
/*
printf("NOCUT oldb %f olde %f   newb %f newe %f\n maxpts %d	newnpts %d newbeg %d\n",oldb,olde,newb,newe,maxpts,newnpts,newbeg);
*/
	}	
	/* now update the header values */
	sachdr.rhdr[H_B] = newb;
	sachdr.rhdr[H_E] = newe;
	sachdr.ihdr[H_NPTS] = newnpts;
	/* now read into the data array with card, but we need an
	 * offset in one case to skip a few */
	

	/* there are several possibilities here 
	 * if the data are within the original window we cant to go from
	 *   OFF1 to OFF2 where 0 <=OFF1 and OFF2 <= maxpts
	 * if we read off the end then we find that we can only read to 
	 * 	maxpts
	 * if OFF1 < 0 then we fill the data array starting at maxbeg, or
	 * like
	 * 	fread(&(*data)[maxbeg-,sizeof(float),MIN(maxpts,newnpts)
	 * 	et cetera et cetera et cetera */



	i1 = newbeg;
	i2 = newend;
	j1 = MAX(0,-i1);
	j2 = MAX(0, i1);
	j3 = MIN(i2-j2+1,maxpts);
/*
printf("i1 %d i2 %d j1 %d j2 %d j3 %d\n (*data)[%d + i] = x[%d + i ], i=0;i< %d;i++;\n",i1,i2,j1,j2,j3,j1,j2,j3);
*/
	/* allocate the data values setting them to zero */
	if(*data == NULL){
		if((*data = (float *)calloc(newnpts,sizeof(float))) == NULL){
			nerr = -2;
			fclose(fptr);
			return(nerr);
		}
	} else {
		/* set to zero */
		if((*data = (float *)realloc(*data,newnpts*sizeof(float))) == NULL){
			nerr = -2;
			fclose(fptr);
			return(nerr);
		}
	}
	/* initialize to zero since realloc does not do this */
	for(i=0;i<newnpts;i++)
		(*data)[i] = 0.0;
	/* allocate the temporary array to be read in */
	if(x == (float *)NULL)
		x = (float *)calloc(maxpts,sizeof(float));
	else
		x = (float *)realloc(x,maxpts*sizeof(float));
	nread = fread(x,sizeof(float),maxpts,fptr);
	i = 0;
	while(i < j3){
		if((j2+i) >= 0 && (j2+i)< maxpts){
			/* KLUDGE TO MAKE THIS WORK FOR e -1 e 0.05 */
		(*data)[j1 + i] = x[j2 + i ];
/*
		printf(" i %d (*data)[%d+%d] = %f x[%d+%d]=%f\n",i,j1,i,(*data)[j1 + i],j2,i, x[j2 + i ]);
*/

		}
		i++;
	}
	/* now we will fill the data array */
	/*
	nread = fread(*data,sizeof(float),maxpts,fptr);
	*/
	fclose (fptr) ;
	*sachdrret = sachdr;
	return(nerr);
}

int brsach(char *fname,struct  sachdr_ *sachdrret)
{
	int maxpts;
	int i;
	int nerr;
	int has12345;
	int rhas12345, ihas12345, chas12345;
	FILE *fptr;
	struct sachdr_ sachdr;
	/*
		determine if the file exists for reading
	*/
	if((fptr=fopen(fname,"rb")) == NULL){
		sachdr.ihdr[H_NPTS] = 0;
		nerr = -1;
		perror("fopen error in brsac:");
		return(nerr);
	}
#ifdef WIN32
	setmode(fileno(fptr), O_BINARY);
#endif
	/*
		file exists get header information only
	*/
	
	nerr = 0;
	if(fread(&sachdr,sizeof( struct sachdr_),1,fptr)!=1){
		nerr=-3;
		fclose(fptr);
		return(nerr);
	}
	maxpts = sachdr.ihdr[H_NPTS] ;
			nerr = 0 ;

	/* there must be at least one -12345 for this to be a sac file!! */
	has12345 = 0;
	rhas12345 = 0;
	ihas12345 = 0;
	chas12345 = 0;
	for (i=0;i < 70 ; i++){
		if(sachdr.rhdr[i] == -12345.0) {
			has12345++ ; 
			rhas12345++;
		}
	}
	for (i=0;i < 40 ; i++){
		if(sachdr.ihdr[i] == -12345){
			has12345++; 
			ihas12345++;
		}
	}
	for (i=0;i < 24 ; i++){
		if(strncmp(sachdr.chdr[i],  "-12345  ",8)==0){
			has12345++;
			chas12345++;
		}
	}
	/* File 
		must have at least one -12345 sequence, 
		must have the correct * header version, and 
		must be equally spaced time series 
	*/
	if(has12345 == 0 || sachdr.ihdr[H_NVHDR] != 6 || 
		(sachdr.ihdr[H_IFTYPE] !=ENUM_ITIME && sachdr.ihdr[H_IFTYPE] !=ENUM_IXY)){
	if(has12345 == 0)
		fprintf(stderr,"Problem: %s is not a Sac file \n",fname);
	else
		fprintf(stderr,"Problem: %s may be a byte-swapped Sac file. Separately run 'saccvt -I'\n",fname);
		nerr = -3;
		fclose(fptr);
		return(nerr);
	}
	
	fclose (fptr) ;
	*sachdrret = sachdr;
	return(nerr);
}

int bwsac(char *fname,struct  sachdr_ sachdr, float *data){
	FILE *fptr;
	int maxpts, nwrite, nerr;
	if((fptr=fopen(fname,"wb")) == NULL){
		nerr = -1;
		perror("fopen error in bwsac:");
		return(nerr);
	}
#ifdef WIN32
	setmode(fileno(fptr), O_BINARY);
#endif
	fseek(fptr, 0L, SEEK_SET);
	if(fwrite(&sachdr,sizeof( struct sachdr_),1,fptr)!=1){
		nerr=-3;
		fclose(fptr);
		return(nerr);
	}
	maxpts = sachdr.ihdr[H_NPTS] ;
	nwrite = fwrite(data,sizeof(float),maxpts,fptr);
	fclose(fptr);
	return(0);
}

int bwsach(char *fname,struct  sachdr_ sachdr){
	int nerr;
	FILE *fptr;
	if((fptr=fopen(fname,"r+b")) == NULL){
		perror("fopen error in bwsach:");
		nerr = -1;
		return(nerr);
	}
#ifdef WIN32
	setmode(fileno(fptr), O_BINARY);
#endif
	fseek(fptr, 0L, SEEK_SET);

	if(fwrite(&sachdr,sizeof( struct sachdr_),1,fptr)!=1){
		nerr=-3;
		fclose(fptr);
		return(nerr);
	}
	fclose(fptr);
	return(0);
}


int gsac_valid_sacfile(char *name)
{
        struct sacfile_ tsacdata ;
        int iret;
        /* examine the trace header to determine if the
         * file is a valid SAC file
         * return YES or NO
         * The check is performed by opening the file, reading the
         * 3 headers, and determining the following:
         * Is the Version correct, is there one -12345. or -12345
         * in the real/integer header. Also determine if the
         * byte order is reversed
         * */
        strcpy(tsacdata.sac_ifile_name , name);
        strcpy(tsacdata.sac_ofile_name , name);
        iret = brsach(tsacdata.sac_ifile_name,&tsacdata.sachdr);
        if(iret < 0)
                return NO;
        else
                return YES;
}

struct sacfile_ *sacdata ;
int *sortptr;
float *sortfloat;


void gsac_alloc_trace(int oldmax ){
/* allocate a data structure */
/* we need counter here for the number of current successes */
int k;
	if(oldmax == 0 && gsac_control.max_number_traces == 0 ){
		sacdata = (struct sacfile_ *)calloc(1, sizeof(struct sacfile_));
		sortptr = (int *)calloc(1, sizeof(int));
		sortfloat = (float *)calloc(1, sizeof(float));
		sortptr[0] = 0;
		sacdata[0].sac_data = (float *)NULL;
		sacdata[0].sac_spectra = (float *)NULL;
		gsac_control.max_number_traces++;
		sortptr[0] = 0;
	} else {
		if(gsac_control.max_number_traces >(oldmax-1))
		sacdata = realloc(sacdata, sizeof(struct sacfile_)*(gsac_control.max_number_traces +1));
		sortptr = realloc(sortptr, sizeof(int)*(gsac_control.max_number_traces +1));
		sortfloat = realloc(sortfloat, sizeof(float)*(gsac_control.max_number_traces +1));
		sacdata[gsac_control.max_number_traces].sac_data = (float *)NULL;
		sacdata[gsac_control.max_number_traces].sac_spectra = (float *)NULL;
		gsac_control.max_number_traces++;
		for(k=0;k<gsac_control.max_number_traces;k++){
			sortptr[k] = k;
		}
	}
}

