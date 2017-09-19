import pandas as pd
from pandas import DataFrame
import numpy as np
import PIV
import h5py
import matplotlib.pyplot as plt
import hotwire as hw

################################################
#   PURPOSE
#   1. Compute Integral Parameters
#   2. Outer Normalize
#   3. Plot
##################################################
#note- vel and axis are flipped to properlly calc delta


def piv_outer(date, num_tests):
	#read in variables
	name = 'data/PIV_' + date + '.h5'
	umean_fov = np.array(pd.read_hdf(name, 'umean'))
	vmean_fov = np.array(pd.read_hdf(name, 'vmean'))
	umean = np.array(pd.read_hdf(name, 'umean_profile_avg'))
	vmean = np.array(pd.read_hdf(name, 'vmean_profile_avg'))
	urms = np.array(pd.read_hdf(name, 'urms_profile_avg'))
	vrms = np.array(pd.read_hdf(name, 'vrms_profile_avg'))
	uvprime = np.array(pd.read_hdf(name, 'uvprime_profile_avg'))
	x = np.array(pd.read_hdf(name, 'xaxis'))
	y = np.array(pd.read_hdf(name, 'yaxis'))
	num_tests = len(umean)

	###1.  Outer Normalize #############
	###################################

	###2.  Outer Normalize #############
	###################################

	###3.  PLOTS ######################
	###################################

	#mean profiles
	#U vs y
	plt.figure()
	plt.plot(umean, y, '-xr')
	plt.xlabel('U (m/sec)', fontsize=14)
	plt.ylabel('Wall Normal Position (m)', fontsize=14)
	plt.legend([r'$Re_{\theta}=$30288, PIV Data'], loc=0)
	plt.show()

	#V vs y
	plt.figure()
	plt.plot(vmean, y, '-xb')
	plt.xlabel('V (m/sec)', fontsize=14)
	plt.ylabel('Wall Normal Position (m)', fontsize=14)
	plt.legend([r'$Re_{\theta}=$30288, PIV Data'], loc=0)
	plt.show()


	#urms vs y
	plt.figure()
	plt.plot(urms, y, '-xr')
	plt.xlabel('$U_{rms}$ (m/sec)', fontsize=20)
	plt.ylabel('Wall Normal Position (m)', fontsize=14)
	plt.legend([r'$Re_{\theta}=$30288, PIV Data'], loc=0)
	plt.show()

	#vrms vs y
	plt.figure()
	plt.plot(vrms, y, '-xb')
	plt.xlabel('$V_{rms}$ (m/sec)', fontsize=20)
	plt.ylabel('Wall Normal Position (m)', fontsize=14)
	plt.legend([r'$Re_{\theta}=$30288, PIV Data'], loc=0)
	plt.show()

	#uprime vs y
	plt.figure()
	plt.plot(uvprime, y, '-xr')
	plt.xlabel('$u^,v^,$', fontsize=20)
	plt.ylabel('Wall Normal Position (m)', fontsize=14)
	plt.legend([r'$Re_{\theta}=$30288, PIV Data'], loc=0)
	plt.show()

	### Mean Vecotr plot
	skip_num = 5

	umean_fov2 = umean_fov[:, 0:-1:skip_num]
	vmean_fov2 = vmean_fov[:, 0:-1:skip_num]
	x2 = x[0:-1:skip_num]
	y2 = y

	Y = np.tile(y2, (len(x2), 1))
	Y = np.transpose(Y)
	X = np.tile(x2-.0543, (len(y2), 1))
	mean_fov2 = (umean_fov2**2 + vmean_fov2**2)**(1/2)

	contour_levels = np.arange(0, 5, .05)
	plt.figure()
	c = plt.contourf(X, Y, mean_fov2, levels = contour_levels, linewidth=40, alpha=.6)
	cbar = plt.colorbar(c)
	cbar.ax.set_ylabel('Velocity (m/sec)')
	plt.hold(True)
	q = plt.quiver(X, Y, umean_fov2, vmean_fov2, angles='xy', scale=50, width=.0025)
	p = plt.quiverkey(q, .11, -.025, 4,"4 m/s",coordinates='data',color='r')
	plt.axis([0, .1, 0, .2])
	plt.ylabel('Wall Normal Position, $y/\delta$', fontsize=18)
	plt.xlabel('Streamwise Position, x (m)', fontsize=14)
	plt.title('Mean PIV Vector Field', fontsize=14)
	plt.show()
	print('Done!')
	return
