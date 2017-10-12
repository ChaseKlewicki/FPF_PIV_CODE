import pandas as pd
import numpy as np
import PIV as piv
import time
import sys
import h5py
from scipy.signal import medfilt
import matplotlib.pyplot as plt
import hotwire as hw

### PIV READ IN CODE ##############
####################################
def piv_readin(date, file_name, base_name, num_images, data_delimiter, sizex, sizey, walloffset, side_error):
	#Initalize variables
	num_tests = len(base_name)
	u = np.ndarray([num_tests, num_images, sizey, sizex])
	v = np.ndarray([num_tests, num_images, sizey, sizex])
	umean = np.ndarray([num_tests, sizey, sizex])
	#vmean1 = np.ndarray([num_tests, sizey, sizex])
	#mask = np.zeros([num_tests, 3])
	umean_profile = dict()
	vmean_profile = dict()
	urms_profile = dict()
	vrms_profile = dict()
	uvprime_profile = dict()
	#determine file name
	#file_name = dict()
	#for j in range(1, num_images+1):
	#    file_name[j] = '/B' + str('{0:05}'.format(i)) + '.txt'
	for j in base_name:
		#Read in
		[x, y, u[j], v[j]] = piv.piv_readin_mod(base_name[j], file_name, data_delimiter, num_images+1, sizey, sizex)
		#Obtain mean vector field
		umean[j] = np.nanmean(u[j, :], axis=0)
	#determine mask position
	tempmask = piv.mask_loc(umean[j])
	mask = list(tempmask)
	#use this to find the mean vel in each image, and look for bad images
	## Resize vecotor field to crop out masked areas and
	# create new vectors which take out the masked areas and any side errors
	sizex_mask = mask[3] - mask[2] - side_error*2
	sizey_mask = mask[1] - mask[0]
	Umask = np.ndarray([num_tests, num_images, sizey_mask, sizex_mask])
	Vmask = np.ndarray([num_tests, num_images, sizey_mask, sizex_mask])
	Ufilt = np.ndarray([num_tests, num_images, sizey_mask, sizex_mask])
	Vfilt = np.ndarray([num_tests, num_images, sizey_mask, sizex_mask])
	umean = np.ndarray([num_tests, sizey_mask, sizex_mask])
	vmean = np.ndarray([num_tests, sizey_mask, sizex_mask])
	for j in base_name:
		Umask[j] = u[j][:, mask[0]:mask[1], int(mask[2]+side_error):int(mask[3]-side_error)]
		Vmask[j] = v[j][:, mask[0]:mask[1], int(mask[2]+side_error):int(mask[3]-side_error)]
		## FILTER IMAGES
		#calc u infinity
		umean[j] = np.nanmean(Umask[j], axis=0)
		uinfinity = np.mean(umean[j, 0:20, :])
		#Filter Images
		[Ufilt[j], Vfilt[j], bad_im_count] = piv.filt_images(Umask[j], Vmask[j], uinfinity, sizey)
		umean[j] = np.nanmean(Ufilt[j], axis=0)
		vmean[j] = np.nanmean(Vfilt[j], axis=0)
	## Determine RMS quantities ##
	uprime = np.ndarray([num_tests, num_images, sizey_mask, sizex_mask])
	vprime = np.ndarray([num_tests, num_images, sizey_mask, sizex_mask])
	uvprime  = np.ndarray([num_tests, num_images, sizey_mask, sizex_mask])
	uvprime_mean = np.ndarray([num_tests, sizey_mask, sizex_mask])
	urms = np.ndarray([num_tests, sizey_mask, sizex_mask])
	vrms = np.ndarray([num_tests, sizey_mask, sizex_mask])
	for j in range(0, num_tests):
		for jj in range(0, num_images):
			uprime[j, jj] = ((Ufilt[j][jj]-umean[j]))
			vprime[j, jj] = ((Vfilt[j][jj]-vmean[j]))
			uvprime[j, jj] = uprime[j, jj]*vprime[j, jj]
		uvprime_mean[j] = np.nanmean(uvprime[j], axis=0)
		urms[j] = np.nanmean(uprime[j]**2, axis=0)**(1/2)
		vrms[j] = np.nanmean(vprime[j]**2, axis=0)**(1/2)

	## wall position adjustment ###########
	#convert to m and take off wall position as seen in images
	x = (x)/1000
	y = (y-walloffset)/1000
	xmask = x[ (mask[2]+side_error):(mask[3]-side_error) ]
	ymask = y[ mask[0]:mask[1] ]
	## Create Mean Profiles for each data set#######
	for j in range(0, num_tests):
	    umean_profile[j] = np.mean(umean[j], axis=1)
	    vmean_profile[j] = np.mean(vmean[j], axis=1)
	    urms_profile[j] = np.mean(urms[j], axis=1)
	    vrms_profile[j] = np.mean(vrms[j], axis=1)
	    uvprime_profile[j] = np.mean(uvprime_mean[j], axis=1)

	## Average multiple profiles together
	#use this if multiple tests are performed at the same condition
	umean_profile_avg = np.zeros(len(umean_profile[0]))
	vmean_profile_avg = np.zeros(len(umean_profile[0]))
	urms_profile_avg = np.zeros(len(umean_profile[0]))
	vrms_profile_avg = np.zeros(len(umean_profile[0]))
	uvprime_profile_avg = np.zeros(len(umean_profile[0]))
	#average datasets together
	for j in range(0, num_tests):
	    umean_profile_avg = umean_profile_avg + umean_profile[j]
	    vmean_profile_avg = vmean_profile_avg + vmean_profile[j]
	    urms_profile_avg = urms_profile_avg + urms_profile[j]
	    vrms_profile_avg = vrms_profile_avg + vrms_profile[j]
	    uvprime_profile_avg = uvprime_profile_avg + uvprime_profile[j]
	#divide profiles by number of tests which were combined
	umean_profile_avg = umean_profile_avg / num_tests
	vmean_profile_avg = vmean_profile_avg / num_tests
	urms_profile_avg = urms_profile_avg / num_tests
	vrms_profile_avg = vrms_profile_avg / num_tests
	uvprime_profile_avg = uvprime_profile_avg / num_tests

	##calculate conf interval
	conf = dict()
	Neff = 75
	conf['u'] =  (np.nanmean(np.nanmean(np.nanvar(Umask, axis=1), axis=0), axis=1))**(1/2) * (1/Neff)**(1/2)
	conf['v'] =  (np.nanmean(np.nanmean(np.nanvar(Vmask, axis=1), axis=0), axis=1))**(1/2) * (1/Neff)**(1/2)
	conf['urms'] =  (np.nanmean(np.nanvar(urms, axis=0), axis=1))**(1/2) * (1/(2*Neff-1))**(1/2)
	conf['vrms'] =  (np.nanmean(np.nanvar(vrms, axis=0), axis=1))**(1/2) * (1/(2*Neff-1))**(1/2)
	sigma_u = (np.nanmean(np.nanvar(Umask, axis=1), axis=0))**(1/2)
	sigma_v = (np.nanmean(np.nanvar(Vmask, axis=1), axis=0))**(1/2)
	conf['uvprime'] = np.nanmean(sigma_u * sigma_v * (1+ (np.nanmean(uvprime_mean, axis=0)/(sigma_u * sigma_v))**2 / (Neff - 1))**(1/2), axis=1)

	###  WRITE OUT DATA
	####################
	#open hdf5 file
	hdf = pd.HDFStore('data/PIV_' + date + '.h5')
	hdf.put('umean', pd.DataFrame(umean[0]))
	hdf.put('vmean', pd.DataFrame(vmean[0]))
	hdf.put('umean_profile_avg', pd.DataFrame(umean_profile_avg))
	hdf.put('vmean_profile_avg', pd.DataFrame(vmean_profile_avg))
	hdf.put('urms_profile_avg', pd.DataFrame(urms_profile_avg))
	hdf.put('vrms_profile_avg', pd.DataFrame(vrms_profile_avg))
	hdf.put('uvprime_profile_avg', pd.DataFrame(uvprime_profile_avg))
	hdf.put('confidence', pd.DataFrame(conf))
	hdf.put('xaxis', pd.Series(xmask))
	hdf.put('yaxis', pd.Series(ymask))
	hdf.put('mask', pd.DataFrame(mask))
	hdf.close()

	print('Data Saved!')
	return(Ufilt, Vfilt, xmask, ymask, bad_im_count)
