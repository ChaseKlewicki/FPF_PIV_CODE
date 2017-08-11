import pandas as pd
from pandas import DataFrame
import numpy as np
import PIV
import h5py
import matplotlib.pyplot as plt
import hotwire as hw

################################################
#   PURPOSE
#   1. Inner Normalize
#   2. PLOT
#       uplus vs yplus
#       urmsplus vs yplus
#       uvprimeplus vs yplus
##################################################
##note- vel and axis are flipped to properlly calc delta

#Set Parameters
def inner(date,num_tests,sizex,sizey,folder,utau):

	## Initalize variables
	conf = dict()
	# yplus = dict()
	# urmsplus = dict()
	# uvprimeplus = dict()
	#read in variables
	name = folder + 'data/PIV_' + date + '.h5'
	umean = np.array(pd.read_hdf(name, 'umean_profile_avg'))[:-3]
	vmean = np.array(pd.read_hdf(name, 'vmean_profile_avg'))[:-3]
	urms = np.array(pd.read_hdf(name, 'urms_profile_avg'))[:-3]
	vrms = np.array(pd.read_hdf(name, 'vrms_profile_avg'))[:-3]
	uvprime = np.array(pd.read_hdf(name, 'uvprime_profile_avg'))[:-3]
	x = np.array(pd.read_hdf(name, 'xaxis'))[:-3]
	y = np.array(pd.read_hdf(name, 'yaxis'))[:-3]
	conf = pd.read_hdf(name, 'confidence')[:-3]
	#air_prop = hw.air_prop(23)
	####Vincenti data###
	name = folder + 'data/outside_data/FPF_Vincenti_Data.h5'
	unorm = np.array(pd.read_hdf(name, 'unorm'))
	ynorm = np.array(pd.read_hdf(name, 'ynorm'))
	Re_tau = ['1450', '2180', '2280', '3270', '5680', '6430', '10770', '15480', '15740']
	name = folder + 'data/outside_data/urmsplus_fpf_primary.csv'
	urms_vincenti = pd.read_csv(name , delimiter=',')
	name = folder + 'data/outside_data/urmsplus_fpf_primary.csv'
	yplus_vincenti = pd.read_csv(name, delimiter=',')
	# ### WUMOIN DATA #######
	wumoin = dict()
	name = folder + 'data/outside_data/yplus_versus_-uv_plus_Re_1840.dat'
	temp = pd.read_csv(name , delimiter=' ')
	temp = temp.shift(1)
	temp.columns = ['0', 'yplus','uvplus']
	temp['yplus'][0] = 0.2558821659990199
	temp['uvplus'][0] = 0.00009276450027256462
	wumoin['yplus_1840'] = temp['yplus']
	wumoin['uvplus_1840'] = temp['uvplus']
	# ## Reza Data ##
	# # #Reza PIV Velocity data
	name = folder + 'data/outside_data/PIV_ZPG_071016.h5'
	urms_reza =np.array(pd.read_hdf(name, 'rms'))
	ynorm_reza = np.array(pd.read_hdf(name, 'yplus'))

	### INNER NORMALIZE #
	#####################
	#approximate utau, calculated and then adjusted
	#utau = [.11]
	uplus = umean/utau
	yplus = (y*utau)/(1.538*10**(-5))
	urmsplus = urms/utau
	vrmsplus = vrms/utau
	uvprimeplus = uvprime/utau[0]**2
	conf['uplus'] = conf['u']/utau
	conf['urmsplus'] = conf['urms']/utau
	conf['vrmsplus'] = conf['vrms']/utau
	conf['uvprimeplus'] = conf['uvprime']/utau[0]**2

	### Plot figure of heated and unheated PIV data
	###############################################
	legend1 = [r'$Re_\tau=$1450', r'$Re_\tau=$2180', r'$Re_\tau=$2280', r'$Re_\tau=$3270', r'$Re_\tau=$5680', r'$Re_\tau=$6430', r'$Re_\tau=$10770', r'$Re_\tau=$15480', r'$Re_\tau=$15740', r'$Re_\tau=$7510, PIV', 'Conf Int.']
	marker_V = ['-xr', '-xg', '-xb', '-xm', '-xk', '-xc', '-sr', '-sg', '-sb']
	# marker = ['-or', '-ob']
	plt.figure()
	##Uplus Yplus
	#plot vincenti data
	for j in range(0, 9):
		plt.semilogx(ynorm[j], unorm[j], marker_V[j])
	#plot PIV data
	plt.semilogx(yplus, uplus, '-or')
	plt.fill_between(yplus, uplus[:,0], uplus[:,0] + np.array(conf['uplus']), color='r', alpha=.5)
	plt.fill_between(yplus, uplus[:,0], uplus[:,0] - np.array(conf['uplus']), color='r', alpha=.5)
	plt.legend(legend1, loc=0, fontsize=10)
	plt.ylabel('$U^+$', fontsize=20)
	plt.xlabel('$y^+$', fontsize=20)
	plt.grid(True)
	plt.axis([.001, 30000, 0, 35])
	plt.show()

	##Urmsplus Yplus
	plt.figure()
	legend1 =  [r'$Re_\tau=$6430, Hot-Wire', r'$Re_\tau=$7510, PIV', r'Conf Int.']
	#plt.semilogx(ynorm_reza[3][:-1], urms_reza, 'k')
	plt.semilogx(yplus_vincenti, urms_vincenti, 'k')
	plt.semilogx(yplus, urmsplus, '-xb')
	plt.fill_between(yplus, urmsplus[:,0], urmsplus[:,0] + np.array(conf['urmsplus']), color='r', alpha=.5)
	plt.fill_between(yplus, urmsplus[:,0], urmsplus[:,0] - np.array(conf['urmsplus']), color='r', alpha=.5)
	plt.legend(legend1, loc=0, fontsize=10)
	plt.ylabel('$u_{rms}^+$', fontsize=20)
	plt.xlabel('$y^+$', fontsize=20)
	plt.grid(True)
	plt.tight_layout()
	#plt.axis([.001, 30000, 0, 35])
	plt.show()

	##Vrmsplus Yplus
	plt.figure()
	legend1 =  [r'$Re_\tau=$7510, PIV', r'Conf Int.']
	plt.semilogx(yplus, vrmsplus, '-xb')
	plt.fill_between(yplus, vrmsplus[:,0], vrmsplus[:,0] + np.array(conf['vrmsplus']), color='r', alpha=.5)
	plt.fill_between(yplus, vrmsplus[:,0], vrmsplus[:,0] - np.array(conf['vrmsplus']), color='r', alpha=.5)
	plt.legend(legend1, loc=0, fontsize=10)
	plt.ylabel('$v_{rms}^+$', fontsize=20)
	plt.xlabel('$y^+$', fontsize=20)
	plt.grid(True)
	plt.tight_layout()
	#plt.axis([.001, 30000, 0, 35])
	plt.show()

	##UVprimeplus Yplus
	legend1 =  [r'$Re_{\theta}=$1840, WuMoin', r'$Re_{\theta}= 30288$, PIV', r'Conf Int']
	plt.figure()
	#WuMoin
	plt.semilogx(wumoin['yplus_1840'], wumoin['uvplus_1840'], 'k')
	#PIV data
	plt.semilogx(yplus, uvprimeplus, '-xr')
	plt.fill_between(yplus, uvprimeplus[:,0], uvprimeplus[:,0] + np.array(conf['uvprimeplus']), color='r', alpha=.5)
	plt.fill_between(yplus, uvprimeplus[:,0], uvprimeplus[:,0] - np.array(conf['uvprimeplus']), color='r', alpha=.5)
	plt.legend(legend1, loc=0, fontsize=10)
	plt.ylabel('$(u^,v^,)^+$', fontsize=16)
	plt.xlabel('$y^+$', fontsize=16)
	plt.axis([1, 10000, 0, 5])
	plt.grid(True)
	plt.show()

	# data = pd.read_csv('data/urms_side_error.csv')
	# legend1 = ['Reza, PIV','$96\%$FOV (standard)', '$37\%$FOV', '$45\%$FOV', '$53\%$FOV', '$61\%$FOV', '$68\%$FOV', '$76\%$FOV', '$84\%$FOV', '$92\%$FOV']
	# plt.figure()
	# plt.semilogx(ynorm_reza[3][:-1], urms_reza, 'k')
	# for j in range(0, 9):
	#     pos1 = str(j)
	#     plt.semilogx(data['y'], data[pos1], '-x')
	# plt.ylabel('$u_{rms}^+$', fontsize=16)
	# plt.xlabel('$y^+$', fontsize=16)
	# plt.legend(legend1, loc=0, fontsize=12)
	# plt.show()


	### caluser CHART
	# air_prop = hw.air_prop(23)
	# k=.387
	# B = 4.32
	# u_clauser = dict()
	# utau_range = np.arange(.101, .16, .005)
	# #create vector for legend which displays all utau values
	# utau_text = list()
	# utau_text = np.str(utau_range[0])
	# for j in range(1, len(utau_range)):
	#     utau_text = np.append(utau_text, np.str(utau_range[j]))
	# utau_text = np.append(utau_text, 'data')
	# #count through all utau values and seperate vector for which data will then
	# # be plotted against
	# count = 0
	# for j in utau_range:
	#     u_clauser[count] = (1/k)*(j/4.5)*np.log(y*4.5/air_prop['nu']) + (1/k)*(j/4.5)*np.log(j/4.5) + B*(j/4.5)
	#     count+=1
	# plt.figure()
	# for j in range(0, count):
	#      plt.semilogx(y*4.5/air_prop['nu'], u_clauser[j])
	# plt.semilogx(y*4.5/air_prop['nu'], umean/4.5, '-kx')
	# plt.legend(utau_text, loc=0, fontsize=12)
	# plt.ylabel('$u/u_{\infty}$', fontsize=16)
	# plt.xlabel(r'$\frac{y u_{\infty}}{\nu}$', fontsize=16)
	# #plt.axis([4, 10, .4, 1.6])
	# plt.grid(True)
	# plt.show()
	print('Done!')
	return
