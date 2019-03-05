/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: ZZPOINT                                               c
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
#ifdef MSDOS
#define INT long
#else
#define INT int
#endif

extern INT dc_curpen;
extern INT dc_mapcurpen;

static INT dim = 16;
INT dc_dim2 = 255;
static INT bit[16][16] = {
/* 4x4 magic square dithered INTo 16x16 */
{  1,225, 49,209, 15,239, 63,223,  4,228, 52,212, 14,238, 62,222},
{177, 81,129, 97,191, 95,143,111,180, 84,132,100,190, 94,142,110},
{193, 33,241, 17,207, 47,255, 31,196, 36,244, 20,206, 46,254, 30},
{113,145, 65,161,127,159, 79,175,116,148, 68,164,126,158, 78,174},
{ 12,236, 60,220,  6,230, 54,214,  9,233, 57,217,  7,231, 55,215},
{188, 92,140,108,182, 86,134,102,185, 89,137,105,183, 87,135,103},
{204, 44,252, 28,198, 38,246, 22,201, 41,249, 25,199, 39,247, 23},
{124,156, 76,172,118,150, 70,166,121,153, 73,169,119,151, 71,167},
{ 13,237, 61,221,  3,227, 51,211, 16,240, 64,224,  2,226, 50,210},
{189, 93,141,109,179, 83,131, 99,192, 96,144,112,178, 82,130, 98},
{205, 45,253, 29,195, 35,243, 19,208, 48,256, 32,194, 34,242, 18},
{125,157, 77,173,115,147, 67,163,128,160, 80,176,114,146, 66,162},
{  8,232, 56,216, 10,234, 58,218,  5,229, 53,213, 11,235, 59,219},
{184, 88,136,104,186, 90,138,106,181, 85,133,101,187, 91,139,107},
{200, 40,248, 24,202, 42,250, 26,197, 37,245, 21,203, 43,251, 27},
{120,152, 72,168,122,154, 74,170,117,149, 69,165,123,155, 75,171}
};



void dv_zzpoint(INT x,INT y)
{
	static INT ipen, jpen;
	INT i, ix, iy, ixm, iym;
	if(dc_curpen < 1000)
		dv_zpoint(x,y);
	else {
		if(bit[x%dim][y%dim] >= dc_mapcurpen)
			dv_zpoint(x,y);
	}
}