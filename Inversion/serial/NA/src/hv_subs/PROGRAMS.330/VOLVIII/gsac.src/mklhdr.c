#include <stdio.h>

/*
struct

key rik position display user_change
*/
/* string patterns */

static int  NRSTR=70;
char *rstr[] = {
"DELTA   ", "DEPMIN  ", "DEPMAX  ", "SCALE   ", "ODELTA  ", 
"B       ", "E       ", "O       ", "A       ", "FMT     ", 
"T0      ", "T1      ", "T2      ", "T3      ", "T4      ", 
"T5      ", "T6      ", "T7      ", "T8      ", "T9      ", 
"F       ", "RESP0   ", "RESP1   ", "RESP2   ", "RESP3   ", 
"RESP4   ", "RESP5   ", "RESP6   ", "RESP7   ", "RESP8   ", 
"RESP9   ", "STLA    ", "STLO    ", "STEL    ", "STDP    ", 
"EVLA    ", "EVLO    ", "EVEL    ", "EVDP    ", "FHDR40  ", 
"USER0   ", "USER1   ", "USER2   ", "USER3   ", "USER4   ",
"USER5   ", "USER6   ", "USER7   ", "USER8   ", "USER9   ", 
"DIST    ", "AZ      ", "BAZ     ", "GCARC   ", "SB      ", 
"SDELTA  ", "DEPMEN  ", "CMPAZ   ", "CMPINC  ", "XMINIMUM", 
"XMAXIMUM", "YMINIMUM", "YMAXIMUM", "ADJTM   ", "TIMMAX  ", 
"TIMMIN  ", "FHDR67  ", "FHDR68  ", "FHDR69  ", "FHDR70  " 
	};
static int NISTR=40;
char *istr[] = {
"NZYEAR  ", "NZJDAY  ", "NZHOUR  ", "NZMIN   ", "NZSEC   ", 
"NZMSEC  ", "NVHDR   ", "NINF    ", "NHST    ", "NPTS    ", 
"NSNPTS  ", "NSN     ", "NXSIZE  ", "NYSIZE  ", "NHDR15  ", 
"IFTYPE  ", "IDEP    ", "IZTYPE  ", "IHDR4   ", "IINST   ", 
"ISTREG  ", "IEVREG  ", "IEVTYP  ", "IQUAL   ", "ISYNTH  ", 
"IDHR11  ", "IDHR12  ", "IDHR13  ", "IDHR14  ", "IDHR15  ", 
"IDHR16  ", "IDHR17  ", "IDHR18  ", "IDHR19  ", "IDHR20  ", 
"LEVEN   ", "LPSPOL  ", "LOVROK  ", "LCALDA  ", "LHDR5   " 
	};
static int NCSTR=24;
char *cstr[] = {
"KSTNM   ", "KEVNM   ", "KEVNMC  ", "KHOLE   ", 
"KO      ", "KA      ", "KT0     ", "KT1     ", 
"KT2     ", "KT3     ", "KT4     ", "KT5     ", 
"KT6     ", "KT7     ", "KT8     ", "KT9     ", 
"KF      ", "KUSER0  ", "KUSER1  ", "KUSER2  ", 
"KCMPNM  ", "KNETWK  ", "KDATRD  ", "KINST   "
	};

static int NESTR=50;
char *estr[] = {
	"ITIME   ", "IRLIM   ", "IAMPH   ", "IXY     ", "IUNKN   ", 
	"IDISP   ", "IVEL    ", "IACC    ", "IB      ", "IDAY    ", 
	"IO      ", "IA      ", "IT0     ", "IT1     ", "IT2     ", 
	"IT3     ", "IT4     ", "IT5     ", "IT6     ", "IT7     ", 
	"IT8     ", "IT9     ", "IRADNV  ", "ITANNV  ", "IRADEV  ", 
	"ITANEV  ", "INORTH  ", "IEAST   ", "IHORZA  ", "IDOWN   ", 
	"IUP     ", "ILLLBB  ", "IWWSN1  ", "IWWSN2  ", "IHGLP   ", 
	"ISRO    ", "INUCL   ", "IPREN   ", "IPOSTN  ", "IQUAKE  ", 
	"IPREQ   ", "IPOSTQ  ", "ICHEM   ", "IOTHER  ", "IGOOD   ", 
	"IGLCH   ", "IDROP   ", "ILOWSN  ", "IRLDTA  ", "IVOLTS  "
	} ;
static int NIVAL=50;
char *Istr[] = {
	"IFTYPE  ", "IFTYPE  ", "IFTYPE  ", "IFTYPE  ", "IDEP    ",
	"IDEP    ", "IDEP    ", "IDEP    ", "IZTYPE  ", "IZTYPE  ",
	"IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ",
	"IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ",
	"IZTYPE  ", "IZTYPE  ", "IRADNV  ", "ITANNV  ", "IRADEV  ",
	"ITANEV  ", "INORTH  ", "IEAST   ", "IHORZA  ", "IDOWN   ",
	"IUP     ", "ILLLBB  ", "IWWSN1  ", "IWWSN2  ", "IHGLP   ",
	"ISRO    ", "IEVTYP  ", "IEVTYP  ", "IEVTYP  ", "IEVTYP  ",
	"IEVTYP  ", "IEVTYP  ", "IEVTYP  ", "IQUAL   ", "IQUAL   ",
	"IQUAL   ", "IQUAL   ", "IQUAL   ", "ISYNTH  ", "IDEP    "    
};

main()
{
	int i;
	for (i=0 ; i < NRSTR ; i++)
		printf("{\"%s\", RHDR, %2d, YES, YES }, \n",rstr[i],i);
	for (i=0 ; i < 15 ; i++)
		printf("{\"%s\", IHDR, %2d, YES, YES }, \n",istr[i],i);
	for (i=16 ; i < 25 ; i++)
		printf("{\"%s\", EHDR, %2d, YES, YES }, \n",istr[i],i);
	for (i=26 ; i < 35 ; i++)
		printf("{\"%s\", IHDR, %2d, YES, YES }, \n",istr[i],i);
	for (i=36 ; i < 40 ; i++)
		printf("{\"%s\", LHDR, %2d, YES, YES }, \n",istr[i],i);
	for (i=0 ; i < NCSTR ; i++)
		printf("{\"%s\", CHDR, %2d, YES, YES }, \n",cstr[i],i);
	/*
		printf("%d\n",i);
		printf("%s %d\n",rstr[i],i);
		*/
}

