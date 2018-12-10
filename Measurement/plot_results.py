import os
import sys
import matplotlib.pylab as plt
import numpy as np 
from collections import defaultdict
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['xtick.minor.width'] = 1.5


results_file = "CCD_data/total_CCD_1.2009.asc"



#load measurements from earthquakes
lines = open("CCD_ellipticity_earthquakes.txt").readlines()
Teq, Eeq, Sdeq = [],[],[]
for i in range(0,len(lines)):
	Teq.append(float(lines[i].split()[0]))
	Eeq.append(float(lines[i].split()[1]))
	Sdeq.append(float(lines[i].split()[2]))

#load theo from CRUST1.0
lines = open("CCD.ell").readlines()
Ttheo, Etheo= [],[]
for i in range(0,len(lines)):
	Ttheo.append(float(lines[i].split()[0]))
	Etheo.append(np.log10(float(lines[i].split()[1])))


dop_threshold = 0.9  #0.9

hv_dic 	= defaultdict(list)
baz_dic = defaultdict(list)
baz_all = []
freq_list = []

# plot month 
lines = open(results_file).readlines()
for i in range(0,len(lines)):
	if len(lines[i].split()) != 6:
		freq = float(lines[i].split()[1])
		T = 1/freq
		dop = float(lines[i].split()[2])
		baz = float(lines[i].split()[0])
		if freq not in freq_list :
			freq_list.append(freq)
		hv = np.log10(float(lines[i].split()[5])) 
		if dop >= dop_threshold:
			hv_dic[freq].append(hv) 
			baz_dic[freq].append(baz)
			baz_all.append(baz)


f_out = open("periods_list.txt","w")
f_out.write("Periods\n")
for f in sorted(freq_list, reverse=True):
	f_out.write(str(round(1/f,3))+"\n")
f_out.close()



# calc mean and std and plot
fig = plt.figure(1, figsize=(11.69, 8.27))
T_list, hv_mean, hv_std, hv_median, hv_err_perc = [],[],[], [],[]
for freq in sorted(freq_list):
	T = 1/freq
	if T <= 15.0:
		hv_mean_freq = np.mean(hv_dic[freq])
		hv_sd_freq = np.std(hv_dic[freq])


		print T, freq,  hv_mean_freq

		TT = [T] * len(hv_dic[freq])
		#plt.scatter(TT, hv_dic[freq], s=1, color="0.75", zorder=0)

		T_list.append(T)
		hv_mean.append(np.mean(hv_dic[freq]))
		hv_std.append(np.std(hv_dic[freq]))

		median = np.median(hv_dic[freq])
		percentile_min = np.percentile(hv_dic[freq], 15.9)
		percentile_max = np.percentile(hv_dic[freq], 84.1)
		err_inf = abs(percentile_min - median)
		err_sup = abs(percentile_max - median)
		error = (err_inf + err_sup) / 2.0
		hv_median.append(median)
		hv_err_perc.append(error)

#plt.scatter(T_list, hv_mean, color="black", zorder=2, label="Meas. from ambient noise")
#plt.errorbar(T_list, hv_mean, yerr=hv_std, fmt=" ", color="black", zorder=1)

plt.scatter(T_list, hv_median, color="black", zorder=2, label="Ellipticity from DOP-E")
plt.errorbar(T_list, hv_median, yerr=hv_err_perc, fmt=" ", color="black", zorder=1)


plt.errorbar(Teq, Eeq, yerr=Sdeq, fmt = " ", color="red")
plt.scatter(Teq, Eeq, color="red", label="Ellipticity from earthquakes")

plt.plot(Ttheo, Etheo, color="green", linewidth=2, linestyle="--", label="Theo. ell. from CRUST1.0 + LITHO1.0")

plt.legend(fontsize=15)
plt.xlabel("Period (s)", size=20)
plt.ylabel("Log(H/V)", size=20)
plt.xscale("log")
plt.xlim(1.8,60)
plt.ylim(-1,1.5)

plt.savefig("Results_CCD_1.2009.png")

