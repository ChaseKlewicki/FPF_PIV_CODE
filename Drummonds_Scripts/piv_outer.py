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


def piv_outer(date, num_tests, legend1):
	#initalize variables
	umean_fov = dict()
	vmean_fov = dict()
	umean = dict()
	vmean = dict()
	urms = dict()
	vrms = dict()
	uvprime = dict()
	x = dict()
	y = dict()
	for j in range(0, num_tests):
		#read in variables
		name = 'data/PIV_' + date + '_' +str(j) + '.h5'
		umean_fov[j] = np.array(pd.read_hdf(name, 'umean'))
		vmean_fov[j] = np.array(pd.read_hdf(name, 'vmean'))
		umean[j] = np.array(pd.read_hdf(name, 'umean_profile_avg'))
		vmean[j] = np.array(pd.read_hdf(name, 'vmean_profile_avg'))
		urms[j] = np.array(pd.read_hdf(name, 'urms_profile_avg'))
		vrms[j] = np.array(pd.read_hdf(name, 'vrms_profile_avg'))
		uvprime[j] = np.array(pd.read_hdf(name, 'uvprime_profile_avg'))
		x[j] = np.array(pd.read_hdf(name, 'xaxis'))
		y[j] = np.array(pd.read_hdf(name, 'yaxis'))

	###2.  Outer Normalize #############
	###################################

	###3.  PLOTS ######################
	###################################
	marker_u = ['-xr', '-or','-sr']
	marker_v = ['-xb', '-ob','-sb']
	#mean profiles
	#U vs y
	plt.figure()
	for j in range(0, num_tests):
		plt.plot(y[j], umean[j], marker_u[j])
	plt.ylabel('U (m/sec)', fontsize=14)
	plt.xlabel('Wall Normal Position (m)', fontsize=14)
	plt.legend(legend1, loc=0)
	plt.show()

	#V vs y
	plt.figure()
	for j in range(0, num_tests):
		plt.plot(y[j], vmean[j], marker_v[j])
	plt.ylabel('V (m/sec)', fontsize=14)
	plt.xlabel('Wall Normal Position (m)', fontsize=14)
	plt.legend(legend1, loc=0)
	plt.show()


	#urms vs y
	plt.figure()
	for j in range(0, num_tests):
		plt.plot(y[j], urms[j], marker_u[j])
	plt.ylabel('$U_{rms}$ (m/sec)', fontsize=20)
	plt.xlabel('Wall Normal Position (m)', fontsize=14)
	plt.legend(legend1, loc=0)
	plt.show()

	#vrms vs y
	plt.figure()
	for j in range(0, num_tests):
		plt.plot(y[j], vrms[j], marker_v[j])
	plt.ylabel('$V_{rms}$ (m/sec)', fontsize=20)
	plt.xlabel('Wall Normal Position (m)', fontsize=14)
	plt.legend(legend1, loc=0)
	plt.show()

	#uprime vs y
	plt.figure()
	for j in range(0, num_tests):
		plt.plot(y[j], uvprime[j], marker_u[j])
	plt.ylabel('$u^,v^,$', fontsize=20)
	plt.xlabel('Wall Normal Position (m)', fontsize=14)
	plt.legend(legend1, loc=0)
	plt.show()

	### Mean Vecotr plot
	skip_num = 5
	umean_fov2 = umean_fov[0]
	vmean_fov2 = vmean_fov[0]
	x2 = x[0]
	umean_fov2 = umean_fov2[:, 0:-1:skip_num]
	vmean_fov2 = vmean_fov2[:, 0:-1:skip_num]
	x2 = x2[0:-1:skip_num]
	y2 = y[0]

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
