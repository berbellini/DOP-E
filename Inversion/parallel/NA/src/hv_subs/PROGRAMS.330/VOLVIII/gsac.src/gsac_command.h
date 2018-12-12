#ifndef GSAC_COMMAND_H_
#define GSAC_COMMAND_H_
#include	"gsac_docommand.h"
#include	"gsac_help.h"
#define NOP 	0
#define EXIT 	1
#define	LP 	2
#define	HP 	3
#define	BP 	4
#define	IN 	5
#define WH 	6
#define CH 	7
#define READ 	8
#define WRITE 	9
#define PLOT 	10
#define PLOTPK 	11
#define LH 	12
#define HELP 	13
#define RTR 	14
#define RMEAN 	15
#define BG 	16
#define ADD 	17
#define SUB 	18
#define MUL 	19
#define DIV 	20
#define INT 	21
#define DIF 	22
#define BR 	23
#define FFT 	24
#define PLOTSP 	25
#define SORT 	26
#define DELETE 	27
#define SYNCHRONIZE 28
#define COLOR 	29
#define FG 	30
#define CUTERR 31
#define ROTATE 32
#define CUT 33
#define TRANSFER 34
#define CD 35
#define FILTER 36
#define HILBERT 37
#define ENVELOPE 38
#define ABSF 39
#define SQR 40
#define SQRT 41
#define QDP 42
#define PRS 43
#define MERGE 45
#define YLIM 46
#define AGC 47
#define PCTL 48
#define HOLD 49
#define PAUSE 50
#define ECHO 51
#define XLIM 52
#define EXP 53
#define LOG 54
#define LOG10 55
#define EXP10 56
#define XLIN 57
#define XLOG 58
#define YLIN 59
#define YLOG 60
#define LINLIN 61
#define LINLOG 62
#define LOGLIN 63
#define LOGLOG 64
#define ROTATE3 65
#define SGN 66
#define CONVOLVE 67
#define CORRELATE 68
#define STACK 69
#define REVERSE 70
#define MARKTIMES 71
#define TAPER 72
#define HISTORY 73
#define WHITEN 74
#define FILEID 75
#define SMOOTH 76
#define REFRACTION 77
#define MULF 78
#define DIVF 79
#define ADDF 80
#define SUBF 81
#define WRITESPEC 82
#define DECIMATE 83
#define MAP 84
#define TRIANGLE 85
#define BOXCAR 86
#define TRAPEZOID 87
#define VERSION 88
#define XGRID 89
#define YGRID 90
#define GRID 91
#define BACKGROUND 92
#define OUTCSV 93
#define READHDR 94
#define TITLE 95
#define MOMENTTENSOR 96
#define SHIFT 97
#define RICKER 98
int gsac_docommand(int nstr, char  **cmdstr, char *input_lineptr);
typedef void (*gsac_cmd)(void);
typedef void (*gsac_set_cmd)(int ncmd, char **cmdstr);

struct gsac_command {
	int cmd;
	char *str;
	gsac_set_cmd gsac_set_param;
	gsac_cmd gsac_exec;
	char **helper;
	} gsac_command_list[]= {
		{NOP, "", &gsac_set_param_nop, &gsac_exec_nop, help_default},
		{EXIT, "EXIT", &gsac_set_param_nop, &gsac_exec_nop, help_default},
		{EXIT, "QUIT", &gsac_set_param_nop, &gsac_exec_nop, help_default},
		{EXIT, "Q", &gsac_set_param_nop, &gsac_exec_nop, help_default},
		{LP, "LP", &gsac_set_param_lp, &gsac_exec_lp, help_lp},
		{LP, "LOWPASS", &gsac_set_param_lp, &gsac_exec_lp, help_lp},
		{HP, "HP", &gsac_set_param_hp, &gsac_exec_hp, help_hp},
		{HP, "HIGHPASS", &gsac_set_param_hp, &gsac_exec_hp, help_hp},
		{BP, "BP", &gsac_set_param_bp, &gsac_exec_bp, help_bp},
		{BP, "BANDPASS", &gsac_set_param_bp, &gsac_exec_bp, help_bp},
		{BR, "BR", &gsac_set_param_br, &gsac_exec_br, help_br},
		{BR, "BANDREJECT", &gsac_set_param_br, &gsac_exec_br, help_br},
		{IN, "INTERP", &gsac_set_param_in, &gsac_exec_in, help_in},
		{IN, "INTERPOLATE", &gsac_set_param_in, &gsac_exec_in, help_in},
		{READ, "R", &gsac_set_param_read, &gsac_exec_read, help_read},
		{READ, "READ", &gsac_set_param_read, &gsac_exec_read, help_read},
		{WRITE, "W", &gsac_set_param_write, &gsac_exec_write, help_write},
		{WRITE, "WRITE", &gsac_set_param_write, &gsac_exec_write, help_write},
		{PLOT, "PLOT2", &gsac_set_param_plot, &gsac_exec_plot, help_plot},
		{PLOT, "PLOT1", &gsac_set_param_plot, &gsac_exec_plot, help_plot},
		{PLOT, "PLOT ", &gsac_set_param_plot, &gsac_exec_plot, help_plot},
		{PLOT, "P2", &gsac_set_param_plot, &gsac_exec_plot, help_plot},
		{PLOT, "P1", &gsac_set_param_plot, &gsac_exec_plot, help_plot},
		{PLOT, "P", &gsac_set_param_plot, &gsac_exec_plot, help_plot},
		{PLOTPK, "PLOTPK", &gsac_set_param_plotpk, &gsac_exec_plotpk, help_plotpk},
		{PLOTPK, "PPK", &gsac_set_param_plotpk, &gsac_exec_plotpk, help_plotpk},
		{PLOTSP, "PLOTSP", &gsac_set_param_plotsp, &gsac_exec_plotsp, help_plotsp},
		{PLOTSP, "PSP", &gsac_set_param_plotsp, &gsac_exec_plotsp, help_plotsp},
		{HELP, "HELP", &gsac_set_param_nop, &gsac_exec_nop, help_default},
		{HELP, "H", &gsac_set_param_nop, &gsac_exec_nop, help_default},
		{LH, "LH", &gsac_set_param_lh, &gsac_exec_lh, help_lh},
		{LH, "LISTHDR", &gsac_set_param_lh, &gsac_exec_lh, help_lh},
		{LH, "LISTHEADER", &gsac_set_param_lh, &gsac_exec_lh, help_lh},
		{RTR, "RTREND", &gsac_set_param_rtr, &gsac_exec_rtr, help_rtr},
		{RTR, "RTR", &gsac_set_param_rtr, &gsac_exec_rtr, help_rtr},
		{RMEAN, "RMEAN", &gsac_set_param_rmean, &gsac_exec_rmean, help_rmean},
		{BG, "BG", &gsac_set_param_bg, &gsac_exec_bg, help_bg},
		{BG, "BEGINGRAPHICS", &gsac_set_param_bg, &gsac_exec_bg, help_bg},
		{BG, "BD", &gsac_set_param_bg, &gsac_exec_bg, help_bg},
		{BG, "BEGINDEVICES", &gsac_set_param_bg, &gsac_exec_bg, help_bg},
		{ADD, "ADD", &gsac_set_param_add, &gsac_exec_add, help_add},
		{SUB, "SUB", &gsac_set_param_sub, &gsac_exec_sub, help_sub},
		{MUL, "MUL", &gsac_set_param_mul, &gsac_exec_mul, help_mul},
		{DIV, "DIV", &gsac_set_param_div, &gsac_exec_div, help_div},
		{INT, "INT", &gsac_set_param_int, &gsac_exec_int, help_int},
		{DIF, "DIF", &gsac_set_param_dif, &gsac_exec_dif, help_dif},
		{FFT, "FFT", &gsac_set_param_fft, &gsac_exec_fft, help_fft},
		{FFT, "DFT", &gsac_set_param_fft, &gsac_exec_fft, help_fft},
		{CH , "CHNHDR", &gsac_set_param_ch, &gsac_exec_ch, help_ch},
		{CH , "CHANGEHEADER", &gsac_set_param_ch, &gsac_exec_ch, help_ch},
		{CH , "CH" , &gsac_set_param_ch, &gsac_exec_ch, help_ch},
		{SORT, "SORT", &gsac_set_param_sort, &gsac_exec_sort, help_sort},
		{DELETE, "DELETE", &gsac_set_param_del, &gsac_exec_del, help_del},
		{DELETE, "DEL", &gsac_set_param_del, &gsac_exec_del, help_del},
		{SYNCHRONIZE, "SYNCHRONIZE", &gsac_set_param_sync, &gsac_exec_sync, help_sync},
		{SYNCHRONIZE, "SYNC", &gsac_set_param_sync, &gsac_exec_sync, help_sync},
		{COLOR, "COLOR", &gsac_set_param_color, &gsac_exec_color, help_color},
		{FG, "FG", &gsac_set_param_fg, &gsac_exec_fg, help_fg},
		{FG, "FUNCGEN", &gsac_set_param_fg, &gsac_exec_fg, help_fg},
		{CUTERR, "CUTERR", &gsac_set_param_cuterr, &gsac_exec_cuterr, help_cuterr},
		{ROTATE, "ROT", &gsac_set_param_rot, &gsac_exec_rot, help_rot},
		{ROTATE, "ROTATE", &gsac_set_param_rot, &gsac_exec_rot, help_rot},
		{ROTATE3, "ROTATE3", &gsac_set_param_rot3, &gsac_exec_rot3, help_rot3},
		{ROTATE3, "ROT3", &gsac_set_param_rot3, &gsac_exec_rot3, help_rot3},
		{CUT, "CUT", &gsac_set_param_cut, &gsac_exec_cut, help_cut},
		{TRANSFER, "TRANSFER", &gsac_set_param_trans, &gsac_exec_trans, help_trans},
		{TRANSFER, "TRANS", &gsac_set_param_trans, &gsac_exec_trans, help_trans},
		{CD, "CD", &gsac_set_param_cd, &gsac_exec_cd, help_cd},
		{FILTER, "FILTER", &gsac_set_param_filter, &gsac_exec_filter, help_filter},
		{FILTER, "FILT", &gsac_set_param_filter, &gsac_exec_filter, help_filter},
		{HILBERT, "HILBERT", &gsac_set_param_hilb, &gsac_exec_hilb, help_hilb},
		{ENVELOPE, "ENVELOPE", &gsac_set_param_env, &gsac_exec_env, help_env},
		{ENVELOPE, "ENV", &gsac_set_param_env, &gsac_exec_env, help_env},
		{ABSF, "ABS", &gsac_set_param_math, &gsac_exec_math, help_abs},
		{SQR, "SQR", &gsac_set_param_math, &gsac_exec_math, help_sqr},
		{SQRT, "SQRT", &gsac_set_param_math, &gsac_exec_math, help_sqrt},
		{EXP, "EXP", &gsac_set_param_math, &gsac_exec_math, help_exp},
		{EXP10, "EXP10", &gsac_set_param_math, &gsac_exec_math, help_exp10},
		{LOG, "LOG", &gsac_set_param_math, &gsac_exec_math, help_log},
		{LOG10, "LOG10", &gsac_set_param_math, &gsac_exec_math, help_log10},
		{XLOG, "XLOG", &gsac_set_param_plotctl, &gsac_exec_plotctl, help_xlog},
		{XLIN, "XLIN", &gsac_set_param_plotctl, &gsac_exec_plotctl, help_xlin},
		{YLOG, "YLOG", &gsac_set_param_plotctl, &gsac_exec_plotctl, help_ylog},
		{YLIN, "YLIN", &gsac_set_param_plotctl, &gsac_exec_plotctl, help_ylin},
		{LINLIN, "LINLIN", &gsac_set_param_plotctl, &gsac_exec_plotctl, help_linlin},
		{LINLOG, "LINLOG", &gsac_set_param_plotctl, &gsac_exec_plotctl, help_linlog},
		{QDP, "QDP", &gsac_set_param_qdp, &gsac_exec_qdp, help_qdp},
		{PRS, "PRS", &gsac_set_param_prs, &gsac_exec_prs, help_prs},
		{PRS, "PLOTRECORDSECTION", &gsac_set_param_prs, &gsac_exec_prs, help_prs},
		{MERGE, "MERGE", &gsac_set_param_merge, &gsac_exec_merge, help_merge},
		{YLIM, "YLIM", &gsac_set_param_ylim, &gsac_exec_ylim, help_ylim},
		{AGC, "AGC", &gsac_set_param_dagc, &gsac_exec_dagc, help_dagc},
		{PCTL, "PCTL", &gsac_set_param_pctl, &gsac_exec_pctl, help_pctl},
		{HOLD, "HOLD", &gsac_set_param_hold, &gsac_exec_hold, help_hold},
		{PAUSE, "PAUSE", &gsac_set_param_pause, &gsac_exec_pause, help_pause},
		{ECHO, "ECHO", &gsac_set_param_echo, &gsac_exec_echo, help_echo},
		{XLIM, "XLIM", &gsac_set_param_xlim, &gsac_exec_xlim, help_xlim},
		{ROTATE3, "ROTATE3", &gsac_set_param_rot3, &gsac_exec_rot3, help_rot3},
		{SGN, "SGN", &gsac_set_param_sgn, &gsac_exec_sgn, help_sgn},
		{SGN, "SIGN", &gsac_set_param_sgn, &gsac_exec_sgn, help_sgn},
		{CONVOLVE, "CONVOLVE", &gsac_set_param_conv, &gsac_exec_conv, help_conv},
		{CONVOLVE, "CON", &gsac_set_param_conv, &gsac_exec_conv, help_conv},
		{CONVOLVE, "CONV", &gsac_set_param_conv, &gsac_exec_conv, help_conv},
		{CORRELATE, "CORRELATE", &gsac_set_param_corr, &gsac_exec_corr, help_corr},
		{CORRELATE, "COR", &gsac_set_param_corr, &gsac_exec_corr, help_corr},
		{STACK, "STACK", &gsac_set_param_stack, &gsac_exec_stack, help_stack},
		{REVERSE, "REVERSE", &gsac_set_param_rev, &gsac_exec_rev, help_rev},
		{REVERSE, "REV", &gsac_set_param_rev, &gsac_exec_rev, help_rev},
		{MARKTIMES, "MARKTIMES", &gsac_set_param_markt, &gsac_exec_markt, help_markt},
		{MARKTIMES, "MARKT", &gsac_set_param_markt, &gsac_exec_markt, help_markt},
		{TAPER, "TAPER", &gsac_set_param_taper, &gsac_exec_taper, help_taper},
		{HISTORY, "HISTORY", &gsac_set_param_hist, &gsac_exec_hist, help_hist},
		{FILEID, "FILEID", &gsac_set_param_fileid, &gsac_exec_fileid, help_fileid},
		{SMOOTH, "SMOOTH", &gsac_set_param_smth, &gsac_exec_smth, help_smth},
		{REFRACTION, "REFRACTION", &gsac_set_param_refr, &gsac_exec_refr, help_refr},
		{REFRACTION, "REFR", &gsac_set_param_refr, &gsac_exec_refr, help_refr},
		{MULF, "MULF", &gsac_set_param_mulf, &gsac_exec_mulf, help_mulf},
		{DIVF, "DIVF", &gsac_set_param_divf, &gsac_exec_divf, help_divf},
		{ADDF, "ADDF", &gsac_set_param_addf, &gsac_exec_addf, help_addf},
		{SUBF, "SUBF", &gsac_set_param_subf, &gsac_exec_subf, help_subf},
		{WRITESPEC, "WRITESPEC", &gsac_set_param_writesp, &gsac_exec_writesp, help_writesp},
		{WRITESPEC, "WRITESP", &gsac_set_param_writesp, &gsac_exec_writesp, help_writesp},
		{WRITESPEC, "WSP", &gsac_set_param_writesp, &gsac_exec_writesp, help_writesp},
		{DECIMATE, "DECIMATE", &gsac_set_param_dec, &gsac_exec_dec, help_dec},
		{DECIMATE, "DEC", &gsac_set_param_dec, &gsac_exec_dec, help_dec},
		{MAP, "MAP", &gsac_set_param_map, &gsac_exec_map, help_map},
		{TRIANGLE, "TRIANGLE", &gsac_set_param_triangle, &gsac_exec_triangle, help_triangle},
		{TRAPEZOID, "TRAPEZOID", &gsac_set_param_trapezoid, &gsac_exec_trapezoid, help_trapezoid},
		{BOXCAR, "BOXCAR", &gsac_set_param_boxcar, &gsac_exec_boxcar, help_boxcar},
		{VERSION, "VERSION", &gsac_set_param_v, &gsac_exec_v, help_v},
		{XGRID, "XGRID", &gsac_set_param_xgrid, &gsac_exec_xgrid, help_xgrid},
		{YGRID, "YGRID", &gsac_set_param_ygrid, &gsac_exec_ygrid, help_ygrid},
		{GRID, "GRID", &gsac_set_param_grid, &gsac_exec_grid, help_grid},
		{BACKGROUND, "BACKGROUND", &gsac_set_param_back, &gsac_exec_back, help_back},
		{OUTCSV, "OUTCSV", &gsac_set_param_outcsv, &gsac_exec_outcsv, help_outcsv},
		{READHDR, "READHDR", &gsac_set_param_rh, &gsac_exec_rh, help_rh},
		{READHDR, "RH", &gsac_set_param_rh, &gsac_exec_rh, help_rh},
		{RICKER, "RICKER", &gsac_set_param_ricker, &gsac_exec_ricker, help_ricker},
		{RICKER, "RI", &gsac_set_param_ricker, &gsac_exec_ricker, help_ricker},
		{TITLE, "TITLE", &gsac_set_param_title, &gsac_exec_title, help_title},
		{TITLE, "T", &gsac_set_param_title, &gsac_exec_title, help_title},
		{MOMENTTENSOR, "MOMENTTENSOR", &gsac_set_param_mt, &gsac_exec_mt, help_mt},
		{MOMENTTENSOR, "MT", &gsac_set_param_mt, &gsac_exec_mt, help_mt},
		{SHIFT, "SHIFT", &gsac_set_param_shift, &gsac_exec_shift, help_shift},
		{WHITEN, "WHITEN", &gsac_set_param_whit, &gsac_exec_whit, help_whit},
		{WH, "WH", &gsac_set_param_nop, &gsac_exec_wh, help_wh},
		{WH, "WRITEHDR", &gsac_set_param_nop, &gsac_exec_wh, help_wh},
		{WH, "WRITEHEADER", &gsac_set_param_nop, &gsac_exec_wh, help_wh}
		};

#endif



