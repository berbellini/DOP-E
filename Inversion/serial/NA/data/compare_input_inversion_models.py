import os
import sys
import matplotlib.pylab as plt
from M_function import *
import time



model1 = "../../synthetic_test/input_models/PRMA_MAMBo-E.d"
model2 = "model.d"

obs  = "hv_files/OBS/PRMA.6.2010.0.2_ell.txt"
pred = "hv_files/SYNT/predicted_ellipticity"


for n in range(0,100000):

	Z1, VP1, VS1, RHO1 = plot_array_from_model(model1)
	Z2, VP2, VS2, RHO2 = plot_array_from_model(model2)

	Tobs, Eobs, sdobs  = [],[],[]
	lines = open(obs).readlines()
	for i in range(1,len(lines)):
		Tobs.append(float(lines[i].split()[0]))
		Eobs.append(float(lines[i].split()[1]))
		sdobs.append(float(lines[i].split()[2]))

	Ts, Es = [],[]
	lines = open(pred).readlines()
	for i in range(1,len(lines)):
		Ts.append(float(lines[i].split()[0]))
		Es.append(float(lines[i].split()[1]))



	plt.figure(1, figsize=(8.27, 11.69))
	plt.subplot(234)
	plt.plot(VS1, Z1, color="black", linewidth=2)
	plt.plot(VS2, Z2, color="red", linewidth=2)
	plt.ylim(5,0)
	plt.xlim(0.0,4)
	plt.ylabel("Depth (km)")
	plt.xlabel("Vs (km/s)")

	plt.subplot(235)
	plt.plot(VP1, Z1, color="black", linewidth=2)
	plt.plot(VP2, Z2, color="red", linewidth=2)
	plt.ylim(5,0)
	plt.xlabel("Vp (km/s)")

	plt.subplot(236)
	plt.plot(RHO1, Z1, color="black", linewidth=2)
	plt.plot(RHO2, Z2, color="red", linewidth=2)
	plt.ylim(5,0)
	plt.xlabel("Rho (g/cm$^3$)")

	plt.subplot(211)
	plt.scatter(Tobs, Eobs, color="black",  label="observed")
	plt.errorbar(Tobs, Eobs, yerr=sdobs, color="0.5", fmt=" ", zorder=0)
	plt.plot(Ts, Es, color="red", label="predicted")
	plt.xscale("log")
	plt.legend()
	plt.savefig("tmp.png")
	plt.close()

	time.sleep(2)


