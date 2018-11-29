import os
import sys
import glob
from calendar import monthrange


# Created by Andrea Berbellini, November 2018. 
# This program performs the DOP-E analysis on 1 year of data in January 2009 from station CCD (Concordia, Antarctica).
# Data are distributed by ORFEUS data center and belong to the GEOSCOPE network:
#
# Data reference:
# Institut De Physique Du Globe De Paris (IPGP), ; Ecole Et Observatoire Des Sciences De La Terre De Strasbourg (EOST). 
# (1982). GEOSCOPE, French Global Network of broad band seismic stations. Institut de Physique du Globe de Paris (IPGP).
#  https://doi.org/10.18715/geoscope.g  



# This code is free to download from the GJI website. 
# For more details, tests, examples plese refer to:
# Berbellini, A., Schimmel, M., Ferreira, A.M.G., Morelli, A.,
#    Constraining S-wave velocity using Rayleigh wave ellipticity
#    from polarization analysis of seismic noise, Geophysical Journal
#    International, 2018, in press.
# 
# DOP-E method is an update of the previous code by Martin Schimmel, for more details please refer to:
# Schimmel M., and J. Gallart, The use of instantaneous 
#    polarization attributes for seismic signal detection 
#    and image enhancement , Geophys.J.Int.,, 155, 653-668, 
#    doi:10.1046/j.1365-246X.2003.02077.x, 2003.
# Schimmel & Gallart, Degree of polarization filter for
#    frequency-dependent signal enhancement through noise
#    suppression, Bull.Seism.Soc.Am., 94, 1016-1035, 
#    doi: 10.1785/0120030178, 2004.
# Schimmel & Gallart, The inverse S Transform in filters 
#    with time-frequency localization , IEEE Transactions on 
#    Signal Processing, 53 (11), 4417 - 4422, 
#    doi:10.1109/TSP.2005.857065, 2005.
# Stutzmann, E., Schimmel, M., Patau, G., Maggi, A., Global 
#    climate imprint on seismic noise , Geochem. Geophys. Geosyst., 
#    10, Q11004, doi:10.1029/2009GC002619, 2009
# Schimmel, M., Stutzmann, E., Ardhuin, F., Gallart, J., 
#    Polarized Earth's Ambient Microseismic Noise , Geochem. Geophys. 
#    Geosyst., doi:10.1029/2011GC003661, 2011.
# Obrebski, M.J., Ardhuin, F., Stutzmann, E., Schimmel, M., 
#    How moderate sea states can generate loud seismic noise in 
#    the deep ocean, Geophys. Res. Lett., 39, L11601, 
#    doi: 10.1029/2012GL051896, 2012.
# Sergeant A., Stutzmann E., Maggi A., Schimmel M., Ardhuin F., 
#    Obrebski M., Frequency-dependent noise sources in the North 
#    Atlantic Ocean, Geochem. Geophys. Geosyst., 14, 
#    doi:10.1002/2013GC004905, 2013.



data_folder = "./CCD_data/"
year = "2009"
measurement_code = "demo"



for m in range(1,2):   #12):
	month = str(m)
	ndays = monthrange(int(year), int(month))[1]
	for d in range(1,ndays+1):
		day = str(d)
		day_folder = data_folder + month + "." + year + "/" + day
		output_file = day_folder + "/" + day + "." + month + "." + year + "."+ measurement_code +".asc"
		if not os.path.isfile(output_file):

			filelist = glob.glob(day_folder+"/*BHZ*chunk*")

			for i in range(0,len(filelist)):
				f = filelist[i]
				Zfile = f
				Nfile = f.replace("BHZ","BHN")
				Efile = f.replace("BHZ","BHE")
			
				if os.path.isfile(Zfile) and os.path.isfile(Zfile) and os.path.isfile(Zfile):
					print Zfile
					print Nfile
					print Efile
					print "---------------"

					os.system("cp " + Zfile + " Z.sac")
					os.system("cp " + Nfile + " N.sac")
					os.system("cp " + Efile + " E.sac")

					os.system("sh run_DOP-E.cmd")
					os.system("rm Z.sac")
					os.system("rm N.sac")
					os.system("rm E.sac")
					sys.exit()
			os.system("cp output_file.asc " + output_file)
			os.system("rm output_file.asc")
			os.system("cp run_DOP-E.cmd " + day_folder)
			
			
		else:
			pass
		print output_file
		os.system("cat " + output_file+ " >> " + data_folder + "total_CCD_1.2009.asc")













