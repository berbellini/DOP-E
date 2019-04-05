import os
import sys
import matplotlib.pylab as plt
sys.path.insert(0, './plot_routines')
from plot_results import *



inversion_name = "DEMO1.1"   # choose a name for your inversion


# preparing folders
os.system("mkdir Results")
out_folder = "Results/" + inversion_name
os.system("mkdir " + out_folder)
os.system("mkdir models")
os.system("mkdir out_models")

exists = os.path.isfile(out_folder + "/full_list.txt") # check if the result file exists, to avoid overwriting it



if exists == True:
	print "\n\n*** Warning: result file already exists! ***"
	print "*** Change inversion name or remove full_list.txt file from " + out_folder + " folder ***"
else:
	os.chdir("NA/src")
	os.system("make clean")    # recompile code
	os.system("make all")
	os.chdir("../data")


	#************************
	# RUN INVERSION
	os.system("../bin/hv_na")
	sys.exit()
	#************************
	os.chdir("../../")


	# Move results to result folder
	os.system("mv full_list.txt " + out_folder)

	# Save inversion parameters to keep track
	os.system("cp NA/data/na.in " + out_folder)
	os.system("cp NA/data/hv_files/hv.in " + out_folder)
	os.system("cp NA/data/hv_files/NA_MDL/hv_param " + out_folder)

	lines = open("NA/data/hv_files/hv.in").readlines()

	# Save observed curve to results folder
	observed_data = lines[6].split()[0]
	os.system("cp NA/data/hv_files/OBS/" + observed_data + " " + out_folder)

	# Save model used for litho to results folder
	ref_model = lines[8].split()[0]
	os.system("cp NA/data/hv_files/REF_MDL/" + ref_model + " " + out_folder)


	# Plot results
	plot_results(out_folder, observed_data, ref_model)









