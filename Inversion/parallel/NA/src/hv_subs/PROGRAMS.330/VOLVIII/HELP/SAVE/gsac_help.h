char *help_default [] = {
""
};

char *help_read[] = {
""
};

char *help_write[] = {
""
};

char *help_fileid[] = {
""
};

char *help_whit[] = {
""
};

char *help_hist[] = {
""
};

char *help_taper[] = {
""
};

char *help_markt[] = {
""
};

char *help_rev[] = {
""
};

char *help_stack[] = {
""
};

char *help_smth[] = {
"\n",
"Apply a smoothing operator\n",
"\n",
"SMOOTH [MEAn|MEDian] [Halfwidth n ] [Pass p] [Default}\n",
"\n",
"MEAN (default): Apply an averaging operator\n",
"MEDIAN: Apply a median filter\n",
"Halfwidth  n:    Smoothing  operator consists of 2n+1 points (de-\n",
"fault 1)\n",
"Pass p: Apply the operator p times\n",
"Default:        Reset to MEAN Halfwidth 1 Pass 1\n",
"\n",
"\n",
"\n",
"The Pass option is new.  SMOOTH MEAN HALFWIDTH 1 PASS 1\n",
"\n",
"ABS, ENV\n",
""
};

char *help_xlim[] = {
""
};

char *help_echo[] = {
""
};

char *help_pause[] = {
""
};

char *help_hold[] = {
""
};

char *help_pctl[] = {
""
};

char *help_dagc[] = {
""
};

char *help_ylim[] = {
""
};

char *help_merge[] = {
""
};

char *help_lp[] = {
"LowPass filter traces\n",
"\n",
"LowPass [options]\n",
"\n",
"where options is one or more of the following: [Butter ] [ Corner\n",
"fc ] [ NPoles npoles ] [ Passes npass] Butter  : Butterworth fil-\n",
"ter\n",
"Corner  : Corner frequency (R) range 0 - Nyquist\n",
"NPoles  : Number of poles  (I) range 1 - 10\n",
"Passes  : Number of passes (I) range 1 - 2 Lowpass filter using a\n",
"BI-LINEAR Z-transformation implementation of a lowpass filter.  A\n",
"bi-linear method is chosen since this is easily implemented alge-\n",
"braically. Passes = 1 gives a causal  filter  while  Passes  =  2\n",
"gives a zero-phase filter with a 6db point at the corner frequen-\n",
"cy.  USER1 = permin, USER2=permax, where permin=1.0/filt_fh,  and\n",
"permax=  0.01/(npts  *  dt)  for  use  by  sacmft96  adn sacpom96\n",
"        HIGHPASS, BANDPASS, BANDREJECT\n",
""
};

char *help_hp[] = {
""
};

char *help_bp[] = {
""
};

char *help_br[] = {
""
};

char *help_lh[] = {
""
};

char *help_bg[] = {
""
};

char *help_rtr[] = {
""
};

char *help_rmean[] = {
""
};

char *help_add[] = {
""
};

char *help_sub[] = {
""
};

char *help_div[] = {
""
};

char *help_mul[] = {
""
};

char *help_int[] = {
""
};

char *help_dif[] = {
""
};

char *help_fft[] = {
""
};

char *help_plot[] = {
""
};

char *help_plotpk[] = {
""
};

char *help_plotsp[] = {
""
};

char *help_ch [] = {
""
};

char *help_sort [] = {
""
};

char *help_del [] = {
""
};

char *help_sync [] = {
""
};

char *help_wh [] = {
""
};

char *help_color [] = {
""
};

char *help_fg [] = {
""
};

char *help_cuterr [] = {
""
};

char *help_rot [] = {
""
};

char *help_rot3[] = {
""
};

char *help_cut [] = {
""
};

char *help_in [] = {
""
};

char *help_trans [] = {
""
};

char *help_cd [] = {
""
};

char *help_filter [] = {
""
};

char *help_hilb [] = {
""
};

char *help_env [] = {
""
};

char *help_abs [] = {
"Take the absolute value of each data point ABS\n",
""
};

char *help_sqr [] = {
""
};

char *help_sqrt [] = {
""
};

char *help_exp [] = {
""
};

char *help_log [] = {
""
};

char *help_exp10 [] = {
""
};

char *help_log10 [] = {
""
};

char *help_ylin [] = {
""
};

char *help_ylog [] = {
""
};

char *help_linlin [] = {
""
};

char *help_linlog [] = {
""
};

char *help_qdp [] = {
""
};

char *help_prs [] = {
""
};

char *help_conv[] = {
""
};

char *help_corr[] = {
""
};

char *help_sgn[] = {
""
};

