import os
import sys
import numpy as np


def gm_periodfile(working_folder, periodfile):
	os.chdir(working_folder)
	os.system("cp " + periodfile + " ./periods.txt")
	os.system('sprep96  -M ' + 'model.d' + ' -NMOD 1  -PARR periods.txt' + '  -R   > /dev/null')
	os.system('sdisp96  > sdisp96.out')
	os.system('sregn96  > sregn96.out')
	command = "sdpegn96 -R -U -ASC > ciccio"
	os.system(command)
	lines = open("SREGN.ASC","r").readlines()
	TT, EE = [], []
	for i in range(1,len(lines)):
        	T = float(lines[i].split()[2])
        	E = np.log10(abs(float(lines[i].split()[8])))
        	TT.append(T)
        	EE.append(E)

    	return TT, EE



working_folder = sys.argv[1]
periodfile = sys.argv[2]

T0, E = gm_periodfile(working_folder, periodfile)

out = open("./predicted_ellipticity", "w")
for i in range(0,len(E)):
	out.write(str(T0[i])+'\t'+str(round(E[i],5))+'\n')
out.close()




