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
date = '072017'
num_tests = 10
sizex = 129  #Num vertical points
sizey = 129  #num horizontal points
## Initalize variables
conf = dict()
# yplus = dict()
# urmsplus = dict()
# uvprimeplus = dict()
#read in variables
name = 'D:/test_072017/PIV_' + date + '.h5'
umean = np.array(pd.read_hdf(name, 'umean_profile_avg'))*-1
vmean = np.array(pd.read_hdf(name, 'vmean_profile_avg'))
urms = np.array(pd.read_hdf(name, 'urms_profile_avg'))
vrms = np.array(pd.read_hdf(name, 'vrms_profile_avg'))
uvprime = np.array(pd.read_hdf(name, 'uvprime_profile_avg'))
x = np.array(pd.read_hdf(name, 'xaxis'))
y = np.array(pd.read_hdf(name, 'yaxis'))
conf = pd.read_hdf(name, 'confidence')
#air_prop = hw.air_prop(23)
####Vincenti data###
name = 'Drummonds_Scripts/data/outside_data/FPF_Vincenti_Data.h5'
unorm = np.array(pd.read_hdf(name, 'unorm'))
ynorm = np.array(pd.read_hdf(name, 'ynorm'))
Re_tau = ['1450', '2180', '2280', '3270', '5680', '6430', '10770', '15480', '15740']
# ### WUMOIN DATA #######
wumoin = dict()
temp = pd.read_csv('Drummonds_Scripts/data/outside_data/yplus_versus_-uv_plus_Re_1840.dat', delimiter=' ')
temp = temp.shift(1)
temp.columns = ['0', 'yplus','uvplus']
temp['yplus'][0] = 0.2558821659990199
temp['uvplus'][0] = 0.00009276450027256462
wumoin['yplus_1840'] = temp['yplus']
wumoin['uvplus_1840'] = temp['uvplus']
# ## Reza Data ##
# # #Reza PIV Velocity data
name = 'Drummonds_Scripts/data/outside_data/PIV_ZPG_071016.h5'
urms_reza =np.array(pd.read_hdf(name, 'rms'))
ynorm_reza = np.array(pd.read_hdf(name, 'yplus'))

### INNER NORMALIZE #
#####################
#approximate utau, calculated and then adjusted
utau = [.11]
uplus = umean/utau
yplus = (y*utau)/(1.538*10**(-5))
urmsplus = urms/utau
uvprimeplus = uvprime/utau
conf['uplus'] = conf['u']/utau
conf['urmsplus'] = conf['urms']/utau
conf['uvprimeplus'] = conf['uvprime']/utau

### Plot figure of heated and unheated PIV data
###############################################
legend1 = [r'$Re_\tau=$1450', r'$Re_\tau=$2180', r'$Re_\tau=$2280', r'$Re_\tau=$3270', r'$Re_\tau=$5680', r'$Re_\tau=$6430', r'$Re_\tau=$10770', r'$Re_\tau=$15480', r'$Re_\tau=$15740', r'$Re_\tau=$7510, PIV', 'Conf Int.']
marker_V = ['-xr', '-xg', '-xb', '-xm', '-xk', '-xc', '-sr', '-sg', '-sb']
# marker = ['-or', '-ob']
plt.figure(figsize = (10,10))
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
plt.figure(figsize = (10,10))
legend1 =  ['Reza, PIV', r'$Re_\tau=$7510, PIV', r'Conf Int.']
plt.semilogx(ynorm_reza[3][:-1], urms_reza, 'k')
plt.semilogx(yplus, urmsplus, '-xb')
plt.fill_between(yplus, urmsplus[:,0], urmsplus[:,0] + np.array(conf['urmsplus']), color='r', alpha=.5)
plt.fill_between(yplus, urmsplus[:,0], urmsplus[:,0] - np.array(conf['urmsplus']), color='r', alpha=.5)
plt.legend(legend1, loc=0, fontsize=10)
plt.ylabel('$u_{rms}^+$', fontsize=20)
plt.xlabel('$y^+$', fontsize=20)
plt.grid(True)
#plt.axis([.001, 30000, 0, 35])
plt.show()

##UVprimeplus Yplus
legend1 =  [r'$Re_{\theta}=$1840, WuMoin', r'$Re_{\theta}= 30288$, PIV', r'Conf Int']
plt.figure(figsize = (10,10))
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
