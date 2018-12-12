
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import sys

def brocher_Vp_Vs(Vs):
	Vp = 0.9409 + 2.0947*Vs - 0.8206*(Vs**2) + 0.2683*(Vs**3) - 0.0251*(Vs**4)
	return Vp

def brocher_rho_Vp(Vp):
	rho = 1.6612*Vp - 0.4721*(Vp**2) + 0.0671*(Vp**3) - 0.0043*(Vp**4) + 0.000106*(Vp**5)
	return rho

def brocher_rho_Vs(Vs):
	Vp = brocher_Vp_Vs(Vs)
	rho = brocher_rho_Vp(Vp)
	return rho

def prem(Z):
	radius=6371000.0  # average radius of the Earth
	#prem=('/localstorage/alessandro/Phase_velocity/test/prem_old')
	#prem = '/home/andrea/Copy/Ellipticity_UCL/lib/PREM_ANISOTROPIC_IDV.csv'
	prem = '/home/andrea/Dropbox/Work/Ellipticity_UCL/lib/PREM_ISOTROPIC_NOOCEAN.csv'
	# read reference prem

	ref_in=open(prem,'r')
	linesref=ref_in.readlines()
	radius_ref=[]
	depth=[]
	density=[]
	Vp = []
	Vs = []
	Vps=[]
	Vsv=[]
	Qk=[]
	Qsh=[]
	Vph=[]
	Vsh=[]
	etha=[]
	for i in reversed(range(3,len(linesref))):
		args0= (linesref[i].split())
		args=[]
		for j in range(len(args0)):        
		    if len(args0[j])>0 :
		       args.append(args0[j])

		radius_ref.append(float(args[0]))
		depth.append((radius - float(args[0])) / 1000.0)
		density.append(float(args[1])/1000.0) 
		Vp.append(float(args[2])/1000.0) 
		Vs.append(float(args[3])/1000.0)
	#	Vps.append(float(args[2])/1000.0) 
	#	Vsv.append(float(args[3])/1000.0)
	#	Vph.append(float(args[6])/1000.0)
	#	Vsh.append(float(args[7])/1000.0 )
		etha.append(float(args[8].strip('\n')))
		Qk.append(float(args[4]))
		Qsh.append(float(args[5]))

	density_func= interp1d(depth, density, kind='linear')
	Vs_func= interp1d(depth, Vs, kind='linear')
	Vp_func= interp1d(depth, Vp, kind='linear')
	#Vps_func= interp1d(depth, Vps, kind='linear')
	#Vsv_func= interp1d(depth, Vsv, kind='linear')
	#Vph_func= interp1d(depth, Vph, kind='linear')
	#Vsh_func= interp1d(depth, Vsh, kind='linear')
	etha_func= interp1d(depth, etha, kind='linear')
	Qk_func= interp1d(depth, Qk, kind='linear')
	Qsh_func= interp1d(depth, Qsh, kind='linear')


	return  density_func(Z),Vp_func(Z),Vs_func(Z), etha_func(Z),Qk_func(Z),Qsh_func(Z) #,Vph_func(Z),Vsh_func(Z),etha_func(Z),Qk_func(Z),Qsh_func(Z)	


def modify_model(model_pre, model_post, parameter, layer, dx):
	lines = open(model_pre, "r").readlines()
	out = open(model_post, "w")
	header_lines = 12
	k = layer + header_lines - 1
	l = lines[k].split()
	

	H = float(l[0])
	Vp = float(l[1])
	Vs = float(l[2])
	rho = float(l[3])
	qp = float(l[4])
	qs = float(l[5])
	etap = float(l[6])
	etas = float(l[7])
	frefp = float(l[8])
	frefs = float(l[9])
	#VpVs = Vp / Vs
	
	
	if parameter == 'H':
		H = float(l[0]) + (dx/100.0) * float(l[0])
		delta_par = H - float(l[0])
	if parameter == 'Vp':
		Vp = float(l[1]) + (dx/100.0) * float(l[1])
		delta_par = Vp - float(l[1])
	if parameter == 'Vs':
		Vs = float(l[2]) + (dx/100.0) * float(l[2])
		delta_par = Vs - float(l[2])
	if parameter == 'rho':
		rho = float(l[3]) + (dx/100.0) * float(l[3])
		delta_par = rho - float(l[3])
	if parameter == 'VpVs':
		VpVs_new  = VpVs + (dx/100.0) * VpVs
		delta_par = VpVs_new - VpVs
		Vp = Vs * VpVs_new

	for i in range(0,len(lines)):
		if i != k:
			out.write(lines[i])
		else:
			out.write(str(H)+'\t'+str(Vp)+'\t'+str(Vs)+'\t'+str(rho)\
			+'\t'+str(qp)+'\t'+str(qs)+'\t'+str(etap)+'\t'+str(etas)\
			+'\t'+str(frefp)+'\t'+str(frefs)+'\n')
			
	return delta_par


def modify_model_log(model_pre, model_post, parameter, layer, dx):
	lines = open(model_pre, "r").readlines()
	out = open(model_post, "w")
	header_lines = 12
	k = layer + header_lines - 1
	l = lines[k].split()

	H = float(l[0])
	Vp = float(l[1])
	Vs = np.log10(float(l[2]))
	rho = float(l[3])
	qp = float(l[4])
	qs = float(l[5])
	etap = float(l[6])
	etas = float(l[7])
	frefp =+ float(l[8])
	frefs = float(l[9])
	
	
	if parameter == 'H':
		H = float(l[0]) + (dx/100.0) * float(l[0])
	if parameter == 'Vp':
		Vp = float(l[1]) + (dx/100.0) * float(l[1])
	if parameter == 'Vs':
		Vs = float(l[2]) + (dx/100.0) * float(l[2])
		
		delta_Vs = np.log10(Vs) - np.log10(float(l[2]))

	if parameter == 'rho':
		rho = float(l[3]) + (dx/100.0) * float(l[3])

	for i in range(0,len(lines)):
		if i != k:
			out.write(lines[i])
		else:
			out.write(str(H)+'\t'+str(Vp)+'\t'+str(Vs)+'\t'+str(rho)\
			+'\t'+str(qp)+'\t'+str(qs)+'\t'+str(etap)+'\t'+str(etas)\
			+'\t'+str(frefp)+'\t'+str(frefs)+'\n')
			
	return delta_Vs



def write_model(old_model_path, new_model_path, m_new, parameter_list):
	num_layer = len(m_new)/len(parameter_list)
	new = open(new_model_path, "w") 
	old_lines = open(old_model_path, "r").readlines()

	if 'Vs' in parameter_list:
		index_Vs = parameter_list.index('Vs')
	if 'H' in parameter_list:
		index_H = parameter_list.index('H')
	if 'H' in parameter_list and 'Vs' in parameter_list:
		index_Vs = parameter_list.index('Vs')
		index_H = parameter_list.index('H')
	if 'Vp' in parameter_list:
		index_Vp = parameter_list.index('Vp')
	if 'rho' in parameter_list:
		index_rho = parameter_list.index('rho')
	

	for i in range(0,12):
		new.write(old_lines[i])

	for i in range(12,12+num_layer):
		j = i-12
		H = float(old_lines[i].split()[0])
		Vp = float(old_lines[i].split()[1])
		Vs = float(old_lines[i].split()[2])
		rho = float(old_lines[i].split()[3])
		A = float(old_lines[i].split()[4])
		B = float(old_lines[i].split()[5])
		C = float(old_lines[i].split()[6])
		D = float(old_lines[i].split()[7])
		E = float(old_lines[i].split()[8])
		F = float(old_lines[i].split()[9])

		if "Vs" in parameter_list: 
			Vs = m_new[j + num_layer * index_Vs]
		if "Vp" in parameter_list: 
			Vp = m_new[j + num_layer * index_Vp]
		if "H" in parameter_list: 
			H = m_new[j + num_layer * index_H]
		if "Vs" in parameter_list and "H" in parameter_list:
			Vs = m_new[j + num_layer * index_Vs]
			H = m_new[j + num_layer * index_H]

		A = float(old_lines[i].split()[4])
		B = float(old_lines[i].split()[5])
		C = float(old_lines[i].split()[6])
		D = float(old_lines[i].split()[7])
		E = float(old_lines[i].split()[8])
		F = float(old_lines[i].split()[9])


		new.write(str("%.3f" %H)+'\t'+str("%.3f" %Vp)+'\t'+str("%.3f" %Vs)+'\t'+str(rho)+'\t'+\
			str(A)+'\t'+str(B)+'\t'+str(C)+'\t'+str(D)+'\t'+str(E)+'\t'+str(F)+'\n')

	for i in range(12+len(m_new),len(old_lines)):
		new.write(old_lines[i])
	
	return


def write_model_log(old_model_path, new_model_path, m_new, parameter_list):
	num_layer = len(m_new)/len(parameter_list)
	new = open(new_model_path, "w") 
	old_lines = open(old_model_path, "r").readlines()

	if 'Vs' in parameter_list:
		index_Vs = parameter_list.index('Vs')
	if 'H' in parameter_list:
		index_H = parameter_list.index('H')
	if 'H' in parameter_list and 'Vs' in parameter_list:
		index_Vs = parameter_list.index('Vs')
		index_H = parameter_list.index('H')
	if 'Vp' in parameter_list:
		index_Vp = parameter_list.index('Vp')
	if 'rho' in parameter_list:
		index_rho = parameter_list.index('rho')
	

	for i in range(0,12):
		new.write(old_lines[i])

	for i in range(12,12+num_layer):
		j = i-12
		H = float(old_lines[i].split()[0])
		Vp = float(old_lines[i].split()[1])
		Vs = float(old_lines[i].split()[2])
		rho = float(old_lines[i].split()[3])
		A = float(old_lines[i].split()[4])
		B = float(old_lines[i].split()[5])
		C = float(old_lines[i].split()[6])
		D = float(old_lines[i].split()[7])
		E = float(old_lines[i].split()[8])
		F = float(old_lines[i].split()[9])

		if "Vs" in parameter_list: 
			Vs = m_new[j + num_layer * index_Vs]
		if "Vp" in parameter_list: 
			Vp = m_new[j + num_layer * index_Vp]
		if "H" in parameter_list: 
			H = m_new[j + num_layer * index_H]
		if "Vs" in parameter_list and "H" in parameter_list:
			Vs = m_new[j + num_layer * index_Vs]
			H = m_new[j + num_layer * index_H]

		A = float(old_lines[i].split()[4])
		B = float(old_lines[i].split()[5])
		C = float(old_lines[i].split()[6])
		D = float(old_lines[i].split()[7])
		E = float(old_lines[i].split()[8])
		F = float(old_lines[i].split()[9])


		new.write(str("%.3f" %H)+'\t'+str("%.3f" %Vp)+'\t'+str("%.3f" %Vs)+'\t'+str(rho)+'\t'+\
			str(A)+'\t'+str(B)+'\t'+str(C)+'\t'+str(D)+'\t'+str(E)+'\t'+str(F)+'\n')

	for i in range(12+len(m_new),len(old_lines)):
		new.write(old_lines[i])
	
	return


def plot_array_from_model(model):
	Z = []
	VP = []
	VS = []
	RHO = []
	z= 0.0
	model_file = open(model, 'r')
	model_lines = model_file.readlines()


	Z.append(0.0)
	VP.append(float(model_lines[12].split()[1]))
	VS.append(float(model_lines[12].split()[2]))
	RHO.append(float(model_lines[12].split()[3]))
	
	
	for i in range(12,len(model_lines)-1):
		h = float(model_lines[i].split()[0])
		Vp = float(model_lines[i].split()[1])
		Vs = float(model_lines[i].split()[2])
		rho = float(model_lines[i].split()[3])
		
		z += h 

		Z.append(z)
		Z.append(z)

		VP.append(Vp)
		VS.append(Vs)
		RHO.append(rho)

		VP.append(float(model_lines[i+1].split()[1]))
		VS.append(float(model_lines[i+1].split()[2]))
		RHO.append(float(model_lines[i+1].split()[3]))

	Z.append(z)
	VP.append(float(model_lines[i+1].split()[1]))
	VS.append(float(model_lines[i+1].split()[2]))
	RHO.append(float(model_lines[i+1].split()[3]))

	
	return Z, VP, VS, RHO


def m_2_model(m, output_file):

	#spessore degli strati del modello interpolato
	thick = [0.1,0.1,0.1,0.1,0.1,0.5,0.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,\
		2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,\
		2.0,2.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,10.0,10.0,10.0,\
		20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,\
		20.0,20.0,20.0,20.0,20.0,20.0]

	header = 'MODEL.01\nMODELLO DEMO\nISOTROPIC\nKGS\nSPHERICAL EARTH\n1-D\nCONSTANT VELOCITY\nLINE08\nLINE09\nLINE10\nLINE11\nH\tVP\tVS\tRHO\tQP\tQS\tETAP\tETAS\tFREFP\tFREFS\n'
	residual = "\t1450.0\t600.0\t1.0\t1.0\t1.0\t1.0\n"

	x = []
	y = []
	depth = 0.0
	y.append(depth)
	
	j = 0
	for i in range(0,len(m)/3):
		H = m[j]
		Vs = 10**(m[j+1])
		Vs_grad = m[j+2]
		j += 3
		Vs2 = H * Vs_grad + Vs
		depth += H

		x.append(Vs)
		x.append(Vs2)
		y.append(depth)
		y.append(depth)

	del y[-1]
	y[-1] = 500.0
		
	f = interp1d(y, x, kind='linear')
	ynew = np.arange(0,500,0.001)

	H_sum = 0.0
	Vs_array = []
	d = 0.0
	for t in thick:
		half = t/2.0
		d = H_sum + half
		H_sum += t
		
		Vs_int = f(d)
		Vs_array.append(Vs_int)
		
	output_model = open(output_file, "w")
	output_model.write(header)


	for i in range(0,len(thick)):
		Vs = Vs_array[i]
		
		l = str(thick[i]) + '\t'+ str("%.3f"%brocher_Vp_Vs(Vs)) + '\t' + str("%.3f" %Vs) + '\t' + str("%.3f"%brocher_rho_Vs(Vs)) + residual
		output_model.write(l)
	
	l = str(500.0) + '\t'+ str("%.3f"%brocher_Vp_Vs(Vs)) + '\t' + str("%.3f" %Vs) + '\t' + str("%.3f"%brocher_rho_Vs(Vs)) + residual
	output_model.write(l)	
	l = str(0.1) + '\t'+ str("%.3f"%brocher_Vp_Vs(Vs)) + '\t' + str("%.3f" %Vs) + '\t' + str("%.3f"%brocher_rho_Vs(Vs)) + residual
	output_model.write(l)
	output_model.close()

	return


def m_2_mod_generic(m, which_list, output_file):

	#thick = [1.0,10.806, 8.606, 16.322, 43.99, 500.0]
	thick_4_gradient = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,\
			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,\
			5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,\
			5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,\
			5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,\
			10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,\
			10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,\
			10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0]

	thick = [10.0, 800.0]
	header = 'MODEL.01\nMODELLO DEMO\nISOTROPIC\nKGS\nSPHERICAL EARTH\n1-D\nCONSTANT VELOCITY\nLINE08\nLINE09\nLINE10\nLINE11\nH\tVP\tVS\tRHO\tQP\tQS\tETAP\tETAS\tFREFP\tFREFS\n'
	residual = "\t1450.0\t600.0\t1.0\t1.0\t1.0\t1.0\n"

	vvs = []
	hh = []
	vvp = []
	rrho = []
	gradVs = []
	deltaVs = []
	j = 0
	rhotemp = 0.0
	for i in range(0,len(m)):
		if which_list[i] == 'H':
			hh.append(m[i])	
		if which_list[i] == 'Vs':
			vvs.append(m[i])
		if which_list[i] == 'logVs':
			vvs.append(10**m[i])
		if which_list[i] == 'Vp':
			vvp.append(m[i])
		if which_list[i] == 'logVp':
			vvp.append(10**m[i])
		if which_list[i] == 'rho':
			rrho.append(m[i])
			rhotemp = m[i]
		if which_list[i] == 'deltarho':
			rhotemp += m[i]
			rrho.append(rhotemp)
		if which_list[i] == 'Vp/Vs':
			if which_list[i-1] == 'logVs':
				vp = 10**m[i-1] * m[i]
			if which_list[i-1] == 'Vs':
				vp = m[i-1] * m[i]
			vvp.append(vp)
		if which_list[i] == 'gradVs':
			gradVs.append(m[i])
		if which_list[i] == 'deltaVs':
			deltaVs.append(m[i])	
		if which_list[i] == 'X':
			rrho.append(brocher_rho_Vs(vvs[-1]))
		j += 1


	

	out = open(output_file, "w")
	out.write(header)

	if hh:
		pass
	else:
		hh = list(thick)

	if vvp:
		pass
	else:	
		for i in vvs:
			vvp.append(brocher_Vp_Vs(i))

	if rrho:
		pass
	else:	
		for i in vvs:
			rrho.append(brocher_rho_Vs(i))

	
	for i in range(0,len(hh)):
		line = str(hh[i])+'\t'+str('%.4f'%vvp[i])+'\t'+str('%.4f'%vvs[i])+'\t'+str('%.4f'%rrho[i])+residual	
	
		out.write(line)

	#	last_line = str(800.0)+'\t'+str('%.4f'%vvp[i])+'\t'+str('%.4f'%vvs[i])+'\t'+str('%.4f'%rrho[i])+residual	
	#	out.write(last_line)
		
	


	out.close()



		
		































	



























