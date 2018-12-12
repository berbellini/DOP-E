import os
import sys
import matplotlib.pylab as plt
from operator import itemgetter
from M_function import *
from matplotlib.pylab import cm
from matplotlib.ticker import FormatStrFormatter
import math
import matplotlib.ticker as ticker

def normalize_misfit(m, m1, m2, g1, g2):
	d = m1 - g1
	m1s = 0
	m2s = m2 - m1
	ms = m - m1
	g1s = 0
	g2s = g2 - g1 
	x = g1 + ms * g2s / m2s
	
	return x

def prepare4plot(tthick, vvs, vvp, rrho):
	Z, VS, VP, RHO = [],[],[],[]

	Z.append(0.0)
	VP.append(vvp[0])
	VS.append(vvs[0])
	RHO.append(rrho[0])
	
	z = 0.0
	for i in range(0,len(vvs)-1):
		h = tthick[i]
		Vp = vvp[i]
		Vs = vvs[i]
		rho = rrho[i]
		
		z += h 

		Z.append(z)
		Z.append(z)

		VP.append(Vp)
		VS.append(Vs)
		RHO.append(rho)

		VP.append(vvp[i+1])
		VS.append(vvs[i+1])
		RHO.append(rrho[i+1])

	Z.append(z)
	VP.append(vvp[i+1])
	VS.append(vvs[i+1])
	RHO.append(rrho[i+1])

	return Z, VS, VP, RHO


def plot_results(result_folder, observed_data, input_model):

	plt.rcParams['xtick.labelsize'] = 15
	plt.rcParams['ytick.labelsize'] = 15
	plt.rcParams['xtick.major.size'] = 6
	plt.rcParams['xtick.minor.size'] = 6
	plt.rcParams['xtick.major.width'] = 2
	plt.rcParams['xtick.minor.width'] = 2

	result_file = result_folder + "/full_list.txt"
	observed_data = result_folder + "/" + observed_data
	input_model = result_folder + "/" + input_model

	lines = open(observed_data).readlines()
	TTobs, EEobs, ssdobs = [],[],[]
	for i in range(1,len(lines)):
		Tobs = float(lines[i].split()[0])
		Eobs = float(lines[i].split()[1])
		sdobs = float(lines[i].split()[2])
		TTobs.append(Tobs)
		EEobs.append(Eobs)
		ssdobs.append(sdobs)


	Zi, VPi, VSi, RHOi = plot_array_from_model(input_model)


	model_list = []
	lines = open(result_file).readlines()
	n_model = 1
	iindex = []
	nmod = 0
	ccost, vvs1, vvs2, vvs3, vvs4, vvs5, vvs6 = [], [],[],[],[],[],[]
	hh1, hh2, hh3 = [],[],[]
	for i in range(0,len(lines)):
		if lines[i].strip() and lines[i].split()[0] == "mft:":
			mft = float(lines[i].split()[1])
			rough = float(lines[i].split()[3])
			cost = float(lines[i].split()[5])
			cost = mft

			EE, TT = [], []
			for n in range(3,1000):
				if lines[i+n].split()[0] != "!":
					T = float(lines[i+n].split()[0])
					E = float(lines[i+n].split()[1])
					EE.append(E)
					TT.append(T)
				else:
					break


			tthick, vvs, vvp, rrho = [],[],[],[]
			for m in range(n+2, 1000):
				if lines[i+m][1] !="#":
					#print lines[i+m]
					thick = float(lines[i+m].split()[0])
					vs = float(lines[i+m].split()[2]) 
					vp = float(lines[i+m].split()[1])
					rho = float(lines[i+m].split()[3])
					tthick.append(thick)
					vvs.append(vs)
					vvp.append(vp)
					rrho.append(rho)
				else:
					break

			iindex.append(nmod)
			nmod+=1
			model_list.append([n_model, cost, TT, EE, tthick, vvs, vvp, rrho, i])
			n_model += 1
			vs1 = vvs[0]
			vs2 = vvs[1]
			vs3 = vvs[2]
			vs4 = vvs[3]
			vs5 = vvs[4]
			vs6 = vvs[5]

			h1 = tthick[0]
			h2 = tthick[1]
			h3 = tthick[2]

			ccost.append(cost)
			vvs1.append(vs1)
			vvs2.append(vs2)
			vvs3.append(vs3)
			vvs4.append(vs4)
			vvs5.append(vs5)
			vvs6.append(vs6)

			hh1.append(h1)
			hh2.append(h2)
			hh3.append(h3)


	model_list_sorted = sorted(model_list, key = itemgetter(1), reverse=True)
	min_cost = model_list_sorted[-1][1]
	max_cost = model_list_sorted[0][1]




	perc = 20	#  20%

	threshold = min_cost * (perc + 100)/100.

	TTbest = model_list_sorted[-1][2]
	EEbest = model_list_sorted[-1][3]
	tthick = model_list_sorted[-1][4]
	vvs = model_list_sorted[-1][5]
	vvp = model_list_sorted[-1][6]
	rrho = model_list_sorted[-1][7]
	line_best = model_list_sorted[-1][8]

	out = open(result_folder+"/best_model.d","w")
	for i in range(line_best, line_best+1000):
	#	print lines[i]
		if lines[i].split()[0] == "model:":
			for j in range(1,11):
				out.write(lines[i+j])
			break
	out.close()

	mean = (vvs[0]*tthick[0] + vvs[1]*tthick[1])/(tthick[0]+tthick[1])

	Zbest, VSbest, VPbest, RHObest = prepare4plot(tthick, vvs, vvp, rrho)


	cmap = cm.hot
	fig = plt.figure(1, figsize=(8.27, 11.69))
	fig.subplots_adjust( wspace=1.)
	ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=4)
	ax2 = plt.subplot2grid((2, 4), (1, 0), colspan=2)
	ax3 = plt.subplot2grid((2, 4), (1, 2), colspan=2)




	for model in  sorted(model_list_sorted, key = itemgetter(1), reverse=True):
		cost = model[1]
		TT = model[2]
		EE = model[3]
		tthick = model[4]
		vvs = model[5]
		vvp = model[6]
		rrho = model[7]
		Z, VS, VP, RHO = prepare4plot(tthick, vvs, vvp, rrho)

		if cost <= threshold:
			colorVal = normalize_misfit(cost, min_cost, threshold, 0.5,1)

			ax2.plot(VS,Z, color=str(colorVal), linewidth=2, zorder=0)
			ax3.plot(VS,Z, color=str(colorVal), linewidth=2, zorder=0)
			ax1.plot(TT,EE, color=str(colorVal), linewidth=2, zorder=9)

	#plt.subplot(223)
	ax2.plot(VSbest,Zbest, color="red", zorder=2, label="Best model", linewidth=2)
	ax2.plot(VSi, Zi, color="black", label= "Litho1.0", linewidth=2, zorder=2, linestyle=":")
	ax2.set_ylim(45,0)
	ax2.set_xlim(0.2,5.0)
	ax2.set_xlabel("Vs (km/s)",size=15)
	ax2.set_ylabel("Depth (km)",size=15)
	ax2.legend(loc=3, fontsize=13)
	ax2.xaxis.set_major_locator(ticker.FixedLocator([1,2,3,4,5]))
	#plt.axhline(6,color="0.5",linestyle="--",linewidth=0.5)
	#ax2.tick_params(labelsize=13)


	#plt.subplot(224)
	ax3.plot(VSbest,Zbest, color="red", zorder=2, label="Best model", linewidth=2)
	ax3.plot(VSi, Zi, color="black", label= "Litho1.0", linewidth=2, zorder=2, linestyle=":")
	ax3.set_ylim(6,0)
	ax3.set_xlim(0., 4)
	ax3.set_xlabel("Vs (km/s)",size=15)
	ax3.set_ylabel("Depth (km)",size=15)
	ax3.legend(loc=3, fontsize=13)

	ax1.errorbar(TTobs, EEobs, yerr=ssdobs, color="black", fmt=" ",zorder=0, alpha=0.3)
	ax1.scatter(TTobs, EEobs, color="black", s=30, label="Observed data", zorder=0)
	ax1.plot(TTbest,EEbest, color="red", label="Theo. ellipticity from best model",zorder=10, linewidth=2)
	ax1.set_xscale("log")
	ax1.set_xlim(0.9*min(TT),max(TT)*1.1)
	ax1.set_ylim(-1,1)
	ax1.set_xlabel("Period (s)",size=15)
	ax1.set_ylabel("Log(H/V)",size=15)
	ax1.legend(loc=4, fontsize=15)
	#ax1.tick_params(labelsize=13)
	ax1.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
	ax1.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
	#ax1.xaxis.set_major_locator(ticker.FixedLocator([2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90]))
	#ax1.set_xticks([9])

	ax1.set_xticklabels([1,1,1,1,1,1,1,1,2,3,4,5,6,7,8,"",20,"",40,"",60,"",80,""],minor=True)


	plt.suptitle("Real data inversion\nStation: CCD\nMisfit threshold = " +str(perc)+"%",size=18)
	plt.savefig(result_folder + "/results_CCD", dpi=200)
	plt.close()


	#----------------------------------------------------------

	fig = plt.figure(1, figsize=(8.27, 11.69))
	fig.subplots_adjust(wspace=0.5, hspace=0.4, top=0.9)
	plt.subplot(9,1,1)
	plt.scatter(iindex, zip(*model_list)[1], color="black", s=1)
	plt.ylabel("Misfit", size=13)
	plt.yscale("log")
	plt.ylim(min(zip(*model_list)[1])*0.5, max(zip(*model_list)[1]))

	plt.subplot(9,1,2)
	plt.scatter(iindex, vvs1, color="black", s=1)
	plt.ylabel("Vs 1\n(km/s)", size=13)
	
	plt.subplot(9,1,3)
	plt.scatter(iindex, hh1, color="black", s=1)
	plt.ylabel("H 1\n(km)", size=13)

	plt.subplot(9,1,4)
	plt.scatter(iindex, vvs2, color="black", s=1)
	plt.ylabel("Vs 2\n(km/s)", size=13)
	
	plt.subplot(9,1,5)
	plt.scatter(iindex, hh2, color="black", s=1)
	plt.ylabel("H 2\n(km)", size=13)
	
	plt.subplot(9,1,6)
	plt.scatter(iindex, hh3, color="black", s=1)
	plt.ylabel("H 3\n(km)", size=13)

	plt.subplot(9,1,7)
	plt.scatter(iindex, vvs4, color="black", s=1)
	plt.ylabel("Vs 4\n(km/s)", size=13)

	plt.subplot(9,1,8)
	plt.scatter(iindex, vvs5, color="black", s=1)
	plt.ylabel("Vs 5\n(km/s)", size=13)

	plt.subplot(9,1,9)
	plt.scatter(iindex, vvs6, color="black", s=1)
	plt.ylabel("Vs 6\n(km/s)", size=13)
	
	plt.xlabel("# Model", size=13)


	plt.suptitle("Inversion evolution", size=18)
	plt.savefig(result_folder + "/convergence.png")
	plt.close()


	#======================================================================




	plt.figure(figsize=(15,15))
	plt.subplots_adjust(left=0.1, right = 0.9, top=0.9, bottom=0.1, hspace=0.2, wspace=0.2)

	ccost, vvs1, vvs2, vvs3, vvs4, vvs5, vvs6 = [], [],[],[],[],[],[]
	for model in  sorted(model_list_sorted, key = itemgetter(1), reverse=True):
		cost = model[1]
		TT = model[2]
		EE = model[3]
		tthick = model[4]
		vvs = model[5]
		vvp = model[6]
		rrho = model[7]

		vs1 = vvs[0]
		vs2 = vvs[1]
		vs3 = vvs[2]
		vs4 = vvs[3]
		vs5 = vvs[4]
		vs6 = vvs[5]

		ccost.append(cost)
		vvs1.append(vs1)
		vvs2.append(vs2)
		vvs3.append(vs3)
		vvs4.append(vs4)
		vvs5.append(vs5)
		vvs6.append(vs6)








	#====================================================


	vmin = np.log10(min_cost)
	vmax = np.log10(max_cost)
	cmap = cm.jet
	m="*"
	fc="red"
	eg="black"
	s=400
	plt.subplot(5,5,1)
	cp = plt.scatter(vvs1, vvs2, c=np.log10(ccost),s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	plt.ylabel("Vs2", fontsize=15)
	plt.xlim(min(vvs1), max(vvs1))
	plt.ylim(min(vvs2), max(vvs2))


	plt.subplot(5,5,6)
	cp = plt.scatter(vvs1, vvs3, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	plt.ylabel("Vs3", fontsize=15)
	#plt.xlim(min(vvs1), max(vvs1))
	#plt.ylim(min(vvs3), max(vvs3))

	plt.subplot(5,5,7)
	cp = plt.scatter(vvs2, vvs3, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	#plt.xlim(min(vvs2), max(vvs2))
	#plt.ylim(min(vvs3), max(vvs3))

	plt.subplot(5,5,11)
	cp = plt.scatter(vvs1, vvs4, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	plt.xlim(min(vvs1), max(vvs1))
	plt.ylim(min(vvs4), max(vvs4))
	plt.ylabel("Vs4", fontsize=15)

	plt.subplot(5,5,12)
	cp = plt.scatter(vvs2, vvs4, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	plt.xlim(min(vvs2), max(vvs2))
	plt.ylim(min(vvs4), max(vvs4))

	plt.subplot(5,5,13)
	cp = plt.scatter(vvs3, vvs4, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	#plt.xlim(min(vvs3), max(vvs3))
	#plt.ylim(min(vvs4), max(vvs4))


	plt.subplot(5,5,16)
	cp = plt.scatter(vvs1, vvs5, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	plt.xlim(min(vvs1), max(vvs1))
	plt.ylim(min(vvs5), max(vvs5))
	plt.ylabel("Vs5", fontsize=15)
	

	plt.subplot(5,5,17)
	cp = plt.scatter(vvs2, vvs5, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	plt.xlim(min(vvs2), max(vvs2))
	plt.ylim(min(vvs5), max(vvs5))
	plt.xlabel("Vs2", fontsize=15)

	plt.subplot(5,5,18)
	cp = plt.scatter(vvs3, vvs5, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	#plt.xlim(min(vvs3), max(vvs3))
	#plt.ylim(min(vvs5), max(vvs5))
	plt.xlabel("Vs3", fontsize=15)

	plt.subplot(5,5,19)
	cp = plt.scatter(vvs4, vvs5, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	plt.xlim(min(vvs4), max(vvs4))
	plt.ylim(min(vvs5), max(vvs5))
	plt.xlabel("Vs4", fontsize=15)

	plt.subplot(5,5,21)
	cp = plt.scatter(vvs1, vvs6, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	plt.xlim(min(vvs1), max(vvs1))
	plt.ylim(min(vvs6), max(vvs6))
	plt.xlabel("Vs1", fontsize=15)
	plt.ylabel("Vs6", fontsize=15)

	plt.subplot(5,5,22)
	cp = plt.scatter(vvs2, vvs6, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	plt.xlim(min(vvs2), max(vvs2))
	plt.ylim(min(vvs6), max(vvs6))
	plt.xlabel("Vs2", fontsize=15)

	plt.subplot(5,5,23)
	cp = plt.scatter(vvs3, vvs6, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	#plt.xlim(min(vvs3), max(vvs3))
	#plt.ylim(min(vvs6), max(vvs6))
	plt.xlabel("Vs3", fontsize=15)

	plt.subplot(5,5,24)
	cp = plt.scatter(vvs4, vvs6, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	plt.xlim(min(vvs4), max(vvs4))
	plt.ylim(min(vvs6), max(vvs6))
	plt.xlabel("Vs4", fontsize=15)

	plt.subplot(5,5,25)
	cp = plt.scatter(vvs5, vvs6, c=np.log10(ccost), s=20, linewidth=0, vmin=vmin, vmax=vmax, cmap=cmap)
	plt.xlim(min(vvs5), max(vvs5))
	plt.ylim(min(vvs6), max(vvs6))
	plt.xlabel("Vs5", fontsize=15)


	plt.savefig(result_folder + "/correlation.png")
	plt.close()

	return
