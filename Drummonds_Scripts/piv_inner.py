import pandas as pd
from pandas import DataFrame
import numpy as np
import PIV
import h5py
import matplotlib.pyplot as plt
import hotwire as hw

### PIV INNER PLOTS AND ANALYSIS ########
#########################################

def piv_inner(date, num_tests, utau, legend1):
	################################################
	#   PURPOSE
	#   1. Inner Normalize
	#   2. PLOT
	#       uplus vs yplus
	#       urmsplus vs yplus
	#       uvprimeplus vs yplus
	##################################################
	##note- vel and axis are flipped to properlly calc delta

	## Initalize variables
	umean = dict()
	vmean = dict()
	urms = dict()
	vrms = dict()
	uvprime = dict()
	x = dict()
	y = dict()
	conf = dict()
	uplus = dict()
	yplus = dict()
	urmsplus = dict()
	vrmsplus = dict()
	uvprimeplus = dict()
	#read in variables
	for j in range(0, num_tests):
		name = 'data/PIV_' + date + '_' +str(j) +'.h5'
		umean[j] = np.array(pd.read_hdf(name, 'umean_profile_avg'))[:-3]
		vmean[j] = np.array(pd.read_hdf(name, 'vmean_profile_avg'))[:-3]
		urms[j] = np.array(pd.read_hdf(name, 'urms_profile_avg'))[:-3]
		vrms[j] = np.array(pd.read_hdf(name, 'vrms_profile_avg'))[:-3]
		uvprime[j] = np.array(pd.read_hdf(name, 'uvprime_profile_avg'))[:-3]
		x[j] = np.array(pd.read_hdf(name, 'xaxis'))[:-3]
		y[j] = np.array(pd.read_hdf(name, 'yaxis'))[:-3]
		#conf[j] = pd.read_hdf(name, 'confidence')[:-3]
	#air_prop = hw.air_prop(23)
	####Vincenti data###
	name = '../outside_data/FPF_Vincenti_Data.h5'
	unorm = np.array(pd.read_hdf(name, 'unorm'))
	ynorm = np.array(pd.read_hdf(name, 'ynorm'))
	Re_tau = ['1450', '2180', '2280', '3270', '5680', '6430', '10770', '15480', '15740']
	urms_vincenti = pd.read_csv('../outside_data/urmsplus_fpf_primary.csv', delimiter=',')
	yplus_vincenti = pd.read_csv('../outside_data/yplus_fpf_primary.csv', delimiter=',')
	# ### WUMOIN DATA #######
	wumoin = dict()
	wumoin = pd.read_hdf('../outside_data/WuMoin2010.h5', 'zpg_re1840')
	# ## Reza Data ##
	# # #Reza PIV Velocity data
	name = '../outside_data/PIV_ZPG_071016.h5'
	urms_reza =np.array(pd.read_hdf(name, 'rms'))
	ynorm_reza = np.array(pd.read_hdf(name, 'yplus'))

	### INNER NORMALIZE #
	#####################
	#approximate utau, calculated and then adjusted
	#utau = [.133]
	for j in range(0, num_tests):
		uplus[j] = umean[j]/utau
		yplus[j] = (y[j]*utau)/(1.538*10**(-5))
		urmsplus[j] = urms[j]/utau
		vrmsplus[j] = vrms[j]/utau
		uvprimeplus[j] = uvprime[j]/utau**2
		#conf[j, 'uplus'] = conf[j, 'u']/utau
		#conf[j, 'urmsplus'] = conf[j, 'urms']/utau
		#conf[j, 'vrmsplus'] = conf[j, 'vrms']/utau
		#conf[j, 'uvprimeplus'] = conf[j, 'uvprime']/utau**2

	### Plot figure of heated and unheated PIV data
	###############################################
	marker_u = ['-xr', '-or','-sr']
	marker_v = ['-xb', '-ob','-sb']

	legend2 = [r'$Re_\tau=$1450', r'$Re_\tau=$2180', r'$Re_\tau=$2280', r'$Re_\tau=$3270', r'$Re_\tau=$5680', r'$Re_\tau=$6430', r'$Re_\tau=$10770', r'$Re_\tau=$15480', r'$Re_\tau=$15740'] + legend1
	marker_V = ['-xr', '-xg', '-xb', '-xm', '-xk', '-xc', '-sr', '-sg', '-sb']
	plt.figure()
	##Uplus Yplus
	#plot vincenti data
	for j in range(0, 9):
		plt.semilogx(ynorm[j], unorm[j], marker_V[j])
	#plot PIV data
	for j in range(0, num_tests):
		plt.semilogx(yplus[j], uplus[j], marker_u[j])
	plt.legend(legend2, loc=2, fontsize=10)
	plt.ylabel('$U^+$', fontsize=20)
	plt.xlabel('$y^+$', fontsize=20)
	plt.grid(True)
	#plt.axis([.0001, 30000, -30, 35])
	plt.show()

	##Urmsplus Yplus
	plt.figure()
	legend2 =  [r'$Re_\tau=$6430, Hot-Wire', r'$Re_{\theta}=$1840, WuMoin'] + legend1
	#plt.semilogx(ynorm_reza[3][:-1], urms_reza, 'k')
	plt.semilogx(yplus_vincenti['yplus_6430'], urms_vincenti['urmsplus_6430'], 'k')
	plt.semilogx(wumoin['yplus'], wumoin['urms_plus'], '--k')
	for j in range(0, num_tests):	
		plt.semilogx(yplus[j], urmsplus[j], marker_u[j])
	#plt.fill_between(yplus, urmsplus[:,0], urmsplus[:,0] + np.array(conf['urmsplus']), color='r', alpha=.5)
	#plt.fill_between(yplus, urmsplus[:,0], urmsplus[:,0] - np.array(conf['urmsplus']), color='r', alpha=.5)
	plt.legend(legend2, loc=2, fontsize=10)
	plt.ylabel('$u_{rms}^+$', fontsize=20)
	plt.xlabel('$y^+$', fontsize=20)
	plt.grid(True)
	#plt.axis([.001, 30000, 0, 35])
	plt.show()

	##Vrmsplus Yplus
	plt.figure()
	legend2 =  [r'$Re_\tau=$6430, Hot-Wire', r'$Re_{\theta}=$1840, Wumoin'] + legend1
	#plt.semilogx(ynorm_reza[3][:-1], urms_reza, 'k')
	plt.semilogx(yplus_vincenti['yplus_6430'], urms_vincenti['urmsplus_6430'], 'k')
	plt.semilogx(wumoin['yplus'], wumoin['vrms_plus'], '--k')
	for j in range(0, num_tests):
		plt.semilogx(yplus[j], vrmsplus[j], marker_v[j])
	#plt.fill_between(yplus, urmsplus[:,0], urmsplus[:,0] + np.array(conf['urmsplus']), color='r', alpha=.5)
	#plt.fill_between(yplus, urmsplus[:,0], urmsplus[:,0] - np.array(conf['urmsplus']), color='r', alpha=.5)
	plt.legend(legend2, loc=2, fontsize=10)
	plt.ylabel('$v_{rms}^+$', fontsize=20)
	plt.xlabel('$y^+$', fontsize=20)
	plt.grid(True)
	#plt.axis([.001, 30000, 0, 35])
	plt.show()

	##UVprimeplus Yplus
	legend2 =  [r'$Re_{\theta}=$1840, WuMoin'] + legend1
	plt.figure()
	#WuMoin
	plt.semilogx(wumoin['yplus'], wumoin['-uv_plus'], '--k')
	#PIV data
	for j in range(0, num_tests):
		plt.semilogx(yplus[j], -1*uvprimeplus[j], marker_u[j])
	#plt.fill_between(yplus, uvprimeplus[:,0], uvprimeplus[:,0] + np.array(conf['uvprimeplus']), color='r', alpha=.5)
	#plt.fill_between(yplus, uvprimeplus[:,0], uvprimeplus[:,0] - np.array(conf['uvprimeplus']), color='r', alpha=.5)
	plt.legend(legend2, loc=2, fontsize=10)
	plt.ylabel('$(u^,v^,)^+$', fontsize=16)
	plt.xlabel('$y^+$', fontsize=16)
	plt.axis([1, 10000, 0, 5])
	plt.grid(True)
	plt.show()

	return
