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
def outer(date,num_tests,sizex,sizey,folder):
	#Parameter set
	#date = '072117'

	#read in variables
	name = folder + 'data/PIV_' + date + '.h5'
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
	plt.plot(umean, y, '-xb')
	plt.xlabel('U (m/sec)')
	plt.ylabel('Wall Normal Position (m)')
	plt.legend(['400rpm'], loc=0)
	plt.show()

	#V vs y
	plt.figure()
	plt.plot(vmean, y, '-xr')
	plt.xlabel('V (m/sec)')
	plt.ylabel('Wall Normal Position (m)')
	plt.legend(['400rpm'], loc=0)
	plt.show()


	#urms vs y
	plt.figure()
	plt.plot(urms, y, '-xb')
	plt.xlabel('$U_{rms}$ (m/sec)')
	plt.ylabel('Wall Normal Position (m)')
	plt.legend(['400rpm'], loc=0)
	plt.show()

	#vrms vs y
	plt.figure()
	plt.plot(vrms, y, '-xr')
	plt.xlabel('$V_{rms}$ (m/sec)')
	plt.ylabel('Wall Normal Position (m)')
	plt.legend(['400rpm'], loc=0)
	plt.show()

	#uprime vs y
	plt.figure()
	plt.plot(uvprime, y, '-xr')
	plt.xlabel('$u^,v^,$')
	plt.ylabel('Wall Normal Position (m)')
	plt.legend(['400rpm'], loc=0)
	plt.show()

	 ### Mean Vecotr plot
	skip_num = 3

	umean_fov2 = umean_fov[:, 0:-1:skip_num]
	vmean_fov2 = vmean_fov[:, 0:-1:skip_num]
	x2 = x[0:-1:skip_num]
	y2 = y

	Y = np.tile(y2/.82, (len(x2), 1))
	Y = np.transpose(Y)
	X = np.tile(x2-.0543, (len(y2), 1))
	mean_fov2 = (umean_fov2**2 + vmean_fov2**2)**(1/2)

	contour_levels = np.arange(0, 3.3, .05)
	plt.figure()
	c = plt.contourf(X, Y, mean_fov2, levels = contour_levels, linewidth=40, alpha=.6)
	cbar = plt.colorbar(c)
	cbar.ax.set_ylabel('Velocity (m/sec)')
	plt.hold(True)
	q = plt.quiver(X, Y, umean_fov2, vmean_fov2, angles='xy', scale=50, width=.0025)
	p = plt.quiverkey(q, .11, -.025, 4,"4 m/s",coordinates='data',color='r')
	plt.axis([0, .1, 0, .246])
	plt.ylabel('Wall Normal Position, $y/\delta$')
	plt.xlabel('Streamwise Position, x (m)')
	plt.title('Mean PIV Vector Field')
	plt.show()

	##need to run readin program before running below plot
	#im_num = 768
	#Umask2 = Umask[0, im_num, :, 0:-1:skip_num]
	#Vmask2 = Vmask[0, im_num, :, 0:-1:skip_num]
	#mean_fov1 = (Umask2**2 + Vmask2**2)**(1/2)

	#contour_levels = np.arange(0, 4.6, .05)
	#plt.figure()
	#c = plt.contourf(X, Y, mean_fov1, levels = contour_levels, linewidth=40, alpha=.6)
	#cbar = plt.colorbar(c)
	#cbar.ax.set_ylabel('Velocity (m/sec)')
	#plt.hold(True)
	#q = plt.quiver(X, Y, -1*Umask2, Vmask2, angles='xy', scale=50, width=.0025)
	#p = plt.quiverkey(q, .11, -.025, 4,"4 m/s",coordinates='data',color='r')
	#plt.axis([0, .1, 0, .246])
	#plt.ylabel('Wall Normal Position, $y/\delta$')
	#plt.xlabel('Streamwise Position (m)')
	#plt.title('Instantaneuos PIV Vector Field')
	#plt.show()
	print('Done!')
	return
