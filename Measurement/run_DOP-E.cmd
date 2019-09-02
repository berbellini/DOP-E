#!/bin/sh
# schimmel@ictja.csic.es
#==============================================#
#
#################################
#   PROVIDE 3-C data (SAC files): 
#################################
Z=Z.sac
N=N.sac
E=E.sac


##############################
#   SET PARAMETERS FOR POLFRE:
##############################
##### DOP power:
pow=3     # 3
##### freq. dependent stability window for DOP:
wlen=17     # 17
##### minimum DOP:
dopm=0.75   # 0.75
##### Gauss window for ST: cycle*T=2*std
cycle=1  #1
##### frequency range:
f1=0.033   
f2=0.5    	

##### neighbouring frequencies to average: 2+nflen+1:
nflen=2     #2
##### numb. of frequencies in band f1-f2:
nfr=100  
##### max number samples to process:
nsp=8192
##### average spectral matrix rather than spectra (keep as it is):
ave=ave
##### frequencies are spaced on a log scale (default is linear)
flog=flog

 
par1=" wlenf="$wlen" pow="$pow" "$ave" dopm="$dopm" nsp="$nsp" "
par2=" f1="$f1" f2="$f2" nflen="$nflen" cycle="$cycle" nfr="$nfr""

echo $Z $N $E
#############
# EXECUTE PG
#############
\rm -f azi_dopm.asc
bin/DOP-E_v1.2 $Z $N $E $par1 $par2 hv wdeg=10 zdeg=10
cat azi_dopm.asc >> output_file.asc


exit






