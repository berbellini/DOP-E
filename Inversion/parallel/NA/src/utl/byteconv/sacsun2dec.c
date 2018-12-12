/*** sacsun2dec - 
 *   writes a byte swaped little endian Seismic Analysis Code (SAC) file
 *   from a big endian SAC file. 
 *   SAC (c) LLNL
 *   Run from little endian machine, DEC alpha. 
 *   saclinux2sun is gcc complied G. A. Ichinose (1996)
 *   sacsun2dec is revised by T. Shibutani, 1999-08-08.
 ***/
#include <stdio.h>
#include <string.h>
#include "sac.h"
#define NUMHEADBYTES 632 /* floats and ints are 440 rest are characters */

int main(int ac, char **av)
{
        struct sac_header sp;
	float *data, fbuf;
	int lbuf;
	int i, j;
	FILE *fp;
	char cbuf[NUMHEADBYTES];
	float float_swap(char cbuf[]);
	int  int_swap(char cbuf[]);
	
	if (ac == 1) { 
		fprintf(stderr, "Usage: %s [sacfile(s)]\n", av[0]);
		exit(-1);
	}
	
	for( i = 1;  i < ac; i++) {
		if ( (fp = fopen(av[i], "rb")) == NULL) {
			fprintf(stderr, "%s Error Opening file %s\n", av[0], av[i]);
		}

	/* set some sac header defaults */
		sp = sac_null;

	/** read in header **/
		fread(cbuf, 440*sizeof(char), 1, fp);
		sp.delta  = float_swap(cbuf);
		sp.depmin = float_swap(cbuf+4);
		sp.depmax = float_swap(cbuf+8);
		sp.scale  = float_swap(cbuf+12);
		sp.odelta = float_swap(cbuf+16);
		sp.b      = float_swap(cbuf+20);
		sp.e      = float_swap(cbuf+24);
		sp.o      = float_swap(cbuf+28);
		sp.a      = float_swap(cbuf+32);
		sp.internal1 = float_swap(cbuf+36);
		sp.t0     = float_swap(cbuf+40);
		sp.t1     = float_swap(cbuf+44);
		sp.t2     = float_swap(cbuf+48);
		sp.t3     = float_swap(cbuf+52);
		sp.t4     = float_swap(cbuf+56);
		sp.t5     = float_swap(cbuf+60);
		sp.t6     = float_swap(cbuf+64);
		sp.t7     = float_swap(cbuf+68);
		sp.t8     = float_swap(cbuf+72);
		sp.t9     = float_swap(cbuf+76);
		sp.f      = float_swap(cbuf+80);
		sp.resp0  = float_swap(cbuf+84);
		sp.resp1  = float_swap(cbuf+88);
		sp.resp2  = float_swap(cbuf+92);
		sp.resp3  = float_swap(cbuf+96);
		sp.resp4  = float_swap(cbuf+100);
		sp.resp5  = float_swap(cbuf+104);
		sp.resp6  = float_swap(cbuf+108);
		sp.resp7  = float_swap(cbuf+112);
		sp.resp8  = float_swap(cbuf+116);
		sp.resp9  = float_swap(cbuf+120);
		sp.stla   = float_swap(cbuf+124);
		sp.stlo   = float_swap(cbuf+128);
		sp.stel   = float_swap(cbuf+132);
		sp.stdp   = float_swap(cbuf+136);
		sp.evla   = float_swap(cbuf+140);
		sp.evlo   = float_swap(cbuf+144);
		sp.evel   = float_swap(cbuf+148);
		sp.evdp   = float_swap(cbuf+152);
		sp.unused1 = float_swap(cbuf+156);
		sp.user0  = float_swap(cbuf+160);
		sp.user1  = float_swap(cbuf+164);
		sp.user2  = float_swap(cbuf+168);
		sp.user3  = float_swap(cbuf+172);
		sp.user4  = float_swap(cbuf+176);
		sp.user5  = float_swap(cbuf+180);
		sp.user6  = float_swap(cbuf+184);
		sp.user7  = float_swap(cbuf+188);
		sp.user8  = float_swap(cbuf+192);
		sp.user9  = float_swap(cbuf+196);
		sp.dist   = float_swap(cbuf+200);
		sp.az     = float_swap(cbuf+204);
		sp.baz    = float_swap(cbuf+208);
		sp.gcarc  = float_swap(cbuf+212);
		sp.internal2 = float_swap(cbuf+216);
		sp.internal3 = float_swap(cbuf+220);
		sp.depmen = float_swap(cbuf+224);
		sp.cmpaz  = float_swap(cbuf+228);
		sp.cmpinc = float_swap(cbuf+232);
		sp.unused2 = float_swap(cbuf+236);
		sp.unused3 = float_swap(cbuf+240);
		sp.unused4 = float_swap(cbuf+244);
		sp.unused5 = float_swap(cbuf+248);
		sp.unused6 = float_swap(cbuf+252);
		sp.unused7 = float_swap(cbuf+256);
		sp.unused8 = float_swap(cbuf+260);
		sp.unused9 = float_swap(cbuf+264);
		sp.unused10 = float_swap(cbuf+268);
		sp.unused11 = float_swap(cbuf+272);
		sp.unused12 = float_swap(cbuf+276);
		sp.nzyear = int_swap(cbuf+280);
		sp.nzjday = int_swap(cbuf+284);
		sp.nzhour = int_swap(cbuf+288);
		sp.nzmin  = int_swap(cbuf+292);
		sp.nzsec  = int_swap(cbuf+296);
		sp.nzmsec = int_swap(cbuf+300);
		sp.internal4 = int_swap(cbuf+304);
		sp.internal5 = int_swap(cbuf+308);
		sp.internal6 = int_swap(cbuf+312);
		sp.npts   = int_swap(cbuf+316);
		sp.internal7 = int_swap(cbuf+320);
		sp.internal8 = int_swap(cbuf+324);
		sp.unused13 = int_swap(cbuf+328);
		sp.unused14 = int_swap(cbuf+332);
		sp.unused15 = int_swap(cbuf+336);
		sp.iftype = int_swap(cbuf+340);
		sp.idep   = int_swap(cbuf+344);
		sp.iztype = int_swap(cbuf+348);
		sp.unused16 = int_swap(cbuf+352);
		sp.iinst  = int_swap(cbuf+356);
		sp.istreg = int_swap(cbuf+360);
		sp.ievreg = int_swap(cbuf+364);
		sp.ievtyp = int_swap(cbuf+368);
		sp.iqual  = int_swap(cbuf+372);
		sp.isynth = int_swap(cbuf+376);
		sp.unused17 = int_swap(cbuf+380);
		sp.unused18 = int_swap(cbuf+384);
		sp.unused19 = int_swap(cbuf+388);
		sp.unused20 = int_swap(cbuf+392);
		sp.unused21 = int_swap(cbuf+396);
		sp.unused22 = int_swap(cbuf+400);
		sp.unused23 = int_swap(cbuf+404);
		sp.unused24 = int_swap(cbuf+408);
		sp.unused25 = int_swap(cbuf+412);
		sp.unused26 = int_swap(cbuf+416);
		sp.leven  = int_swap(cbuf+420);
		sp.lpspol = int_swap(cbuf+424);
		sp.lovrok = int_swap(cbuf+428);
		sp.lcalda = int_swap(cbuf+432);
		sp.unused27 = int_swap(cbuf+436);
		fread(cbuf, (632-440)*sizeof(char), 1, fp);
		strncpy( sp.kstnm, cbuf, 7);
		strncpy( sp.kevnm, cbuf+ 8, 15);
		strncpy( sp.khole, cbuf+24, 7);
		strncpy( sp.ko, cbuf+32, 7);
		strncpy( sp.ka, cbuf+40, 7);
		strncpy( sp.kt0, cbuf+48, 7);
		strncpy( sp.kt1, cbuf+56, 7);
		strncpy( sp.kt2, cbuf+64, 7);
		strncpy( sp.kt3, cbuf+72, 7);
		strncpy( sp.kt4, cbuf+80, 7);
		strncpy( sp.kt5, cbuf+88, 7);
		strncpy( sp.kt6, cbuf+96, 7);
		strncpy( sp.kt7, cbuf+104, 7);
		strncpy( sp.kt8, cbuf+112, 7);
		strncpy( sp.kt9, cbuf+120, 7);
		strncpy( sp.kf, cbuf+128, 7);
		strncpy( sp.kuser0, cbuf+136, 7);
		strncpy( sp.kuser1, cbuf+144, 7);
		strncpy( sp.kuser2, cbuf+152, 7);
		strncpy( sp.kcmpnm, cbuf+160, 7);
		strncpy( sp.knetwk, cbuf+168, 7);
		strncpy( sp.kdatrd, cbuf+176, 7);
		strncpy( sp.kinst, cbuf+184, 7);

		fprintf(stderr, "sta=%-4.4s dt=%g b=%g e=%g\n", sp.kstnm, sp.delta, sp.b, sp.e);
		fprintf(stderr, "%d %dj %dh %dm %ds %dms   n=%d\n",
			sp.nzyear, sp.nzjday, sp.nzhour, sp.nzmin, sp.nzsec, sp.nzmsec, sp.npts);

	/** read the data **/
		data = (float *)malloc(sp.npts*sizeof(float));
		for ( j=0; j<sp.npts; j++) {
			fread(cbuf, sizeof(char), 4, fp);
			data[j] = (float) float_swap(cbuf);
		}
		fclose(fp);

        /** write the header **/
		if ( (fp = fopen(av[i], "w")) == NULL ) {
			fprintf(stdout, "Error in opening output file: %s\n",av[i] );
		}
		fwrite(&sp, sizeof(struct sac_header), 1, fp);
		fwrite(&data[0], sp.npts*sizeof(float), 1, fp);

		fclose(fp);
		free(data);
	}
	return 0;
}

int int_swap(char cbuf[])
{
        union {
                char cval[4];
                int lval;
        } l_union;

	l_union.cval[3] = cbuf[0];
	l_union.cval[2] = cbuf[1];
	l_union.cval[1] = cbuf[2];
	l_union.cval[0] = cbuf[3];
        return(l_union.lval);
}

float float_swap(char cbuf[])
{
        union {
                char cval[4];
                float fval;
        } f_union;

        f_union.cval[3] = cbuf[0];
	f_union.cval[2] = cbuf[1];
	f_union.cval[1] = cbuf[2];
	f_union.cval[0] = cbuf[3];
        return(f_union.fval);
}

