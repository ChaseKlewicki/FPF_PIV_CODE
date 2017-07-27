import pandas as pd
import numpy as np
import PIV
import time
import sys
import h5py
from scipy.signal import medfilt
import matplotlib.pyplot as plt
import hotwire as hw

################################################################
#   PURPOSE
#   1. Readin in PIV data_sets
#   2. Find Mask
#   3. Resize
#   4. Determine mean profiles and HO quantities
################################################################
#Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total:
        print()

def piv_readin_mod(base_name_input, data_sets, sizex, sizey):
    #initalize data
    temp_u = np.ndarray([data_sets-1, sizex, sizey])
    temp_v = np.ndarray([data_sets-1, sizex, sizey])
    count = 0
    x_range = np.arange(1, data_sets)
    #setup progressbar
    printProgressBar(0, len(x_range), prefix = 'Reading In:', suffix = 'Complete', length = 50)
    for i in x_range:
        #create file name for each txt file
        loc = base_name_input + '/B' + str('{0:05}'.format(i)) + '.txt'
        #read in txt file but skip first row
        temp = pd.read_csv(loc, sep='\t', skiprows=1, header=None)
        #rename columns to designated davis output
        temp.columns = ['Xlocation (mm)', 'Ylocation (mm)', 'U (m/sec)', 'V (m/sec)']
        #for j in range(0, len(temp['U (m/sec)'])):
            #temp['U (m/sec)'][j] = float(temp['U (m/sec)'][j].replace(',','.'))
            #temp['V (m/sec)'][j] = float(temp['V (m/sec)'][j].replace(',','.'))
            #temp['Xlocation (mm)'][j] = float(temp['Xlocation (mm)'][j].replace(',','.'))
            #temp['Ylocation (mm)'][j] = float(temp['Ylocation (mm)'][j].replace(',','.'))
        #reorganize into seperate arrays
        temp_x = np.array(np.reshape(temp['Xlocation (mm)'], (sizex, sizey)))
        temp_y = np.array(np.reshape(temp['Ylocation (mm)'], (sizex, sizey)))
        temp_u[count] = np.array(np.reshape(temp['U (m/sec)'], (sizex, sizey)))
        temp_v[count] = np.array(np.reshape(temp['V (m/sec)'], (sizex, sizey)))
        count+=1
        printProgressBar(i, len(x_range), prefix = 'Reading In:', suffix = 'Complete', length = 50)
    x_axis = temp_x[0]
    y_axis = temp_y[:,0]
    print('Done Read in!')
    #sys.stdout.write("\n")
    return(x_axis, y_axis, temp_u, temp_v)

def filt_images(tempU):
    Umean = np.mean(np.mean(tempU))
    STDmean = np.mean(np.var(tempU, axis=0))**(1/2)
    len1 = len(tempU)
    for j in range(0, len1):
        temp_mean = np.mean(np.mean(tempU[j]))
        if temp_mean > Umean + 2*STDmean:
            tempU[j] = tempU[j]*np.nan
        if temp_mean < Umean - 2*STDmean:
            tempU[j] = tempU[j]*np.nan
    return(tempU)


## to work on ##
################
# -work to read in top line of PIV file such that size and units and be predetermined
# determine better way to store PIV datasets

## INITIAL CODE USED FOR READING IN

#Parameter set
date = '061617'
filter_width = 21
num_images = 1000
sizex = 129
sizey = 129
walloffset = 2.22
side_error = 2

#list name of data set folders
base_name = dict()
#List the base name for each test to be read in and analyzed, names taken directly from folder
base_name[0] = 'D:/test_061617/Cam_Date=170616_Time=144917_TR_SeqPIV_MP(2x16x16_50ov_ImgCorr)=unknown'
base_name[1] = 'D:/test_061617/Cam_Date=170616_Time=145207_TR_SeqPIV_MP(2x16x16_50ov_ImgCorr)=unknown'
base_name[2] = 'D:/test_061617/Cam_Date=170616_Time=145516_TR_SeqPIV_MP(2x16x16_50ov_ImgCorr)=unknown'
base_name[3] = 'D:/test_061617/Cam_Date=170616_Time=145802_TR_SeqPIV_MP(2x16x16_50ov_ImgCorr)=unknown'
base_name[4] = 'D:/test_061617/Cam_Date=170616_Time=150030_TR_SeqPIV_MP(2x16x16_50ov_ImgCorr)=unknown'
base_name[5] = 'D:/test_061617/Cam_Date=170616_Time=150311_TR_SeqPIV_MP(2x16x16_50ov_ImgCorr)=unknown'
base_name[6] = 'D:/test_061617/Cam_Date=170616_Time=150551_TR_SeqPIV_MP(2x16x16_50ov_ImgCorr)=unknown'
base_name[7] = 'D:/test_061617/Cam_Date=170616_Time=150828_TR_SeqPIV_MP(2x16x16_50ov_ImgCorr)=unknown'
base_name[8] = 'D:/test_061617/Cam_Date=170616_Time=151133_TR_SeqPIV_MP(2x16x16_50ov_ImgCorr)=unknown'
base_name[9] = 'D:/test_061617/Cam_Date=170616_Time=151411_TR_SeqPIV_MP(2x16x16_50ov_ImgCorr)=unknown'

#Initalize variables
num_tests = len(base_name)
U = np.ndarray([num_tests, num_images-1, sizey, sizex])
V = np.ndarray([num_tests, num_images-1, sizey, sizex])
u_filt = np.ndarray([num_tests, num_images-1, sizey, sizex])
v_filt = np.ndarray([num_tests, num_images-1, sizey, sizex])
umean = np.ndarray([num_tests, sizey, sizex])
vmean = np.ndarray([num_tests, sizey, sizex])
mask = np.zeros([num_tests, 3])
umean_profile = dict()
vmean_profile = dict()
urms_profile = dict()
vrms_profile = dict()
uvprime_profile = dict()
for j in base_name:
    #Read in
    [x, y, U[j], V[j]] = PIV.piv_readin_mod(base_name[j], num_images, sizey, sizex)
    #filter out images that contain spurious vectors
    U[j] = filt_images(U[j])
    V[j] = filt_images(V[j])
    #Filter data sets (median filter along axis=0 (through images))
    #u_filt[j] = medfilt(U[j, :], kernel_size=[filter_width, 1, 1])
    #v_filt[j] = medfilt(V[j, :], kernel_size=[filter_width, 1, 1])
    u_filt[j] = U[j, :]
    v_filt[j] = V[j, :]
    #Obtain mean vecotr field
    umean[j] = np.nanmean(u_filt[j], axis=0)
    vmean[j] = np.nanmean(v_filt[j], axis=0)
    #print('Done Reading')
    #bar.finish()
    #determine mask position
tempmask = PIV.mask_loc(umean[j])
mask = list(tempmask)
## Resize vecotr field to crop out masked areas and
# create new vectors which take out the masked areas and any side errors
sizex_mask = mask[3] - mask[2] - side_error*2
sizey_mask = mask[1] - mask[0]
Umask = np.ndarray([num_tests, num_images-1, sizey_mask, sizex_mask])
Vmask = np.ndarray([num_tests, num_images-1, sizey_mask, sizex_mask])
umean = np.ndarray([num_tests, sizey_mask, sizex_mask])
vmean = np.ndarray([num_tests, sizey_mask, sizex_mask])
for j in base_name:
    Umask[j] = -1*U[j][:, mask[0]:mask[1], int(mask[2]+side_error):int(mask[3]-side_error)]
    Vmask[j] = V[j][:, mask[0]:mask[1], int(mask[2]+side_error):int(mask[3]-side_error)]
    umean[j] = np.nanmean(Umask[j], axis=0)
    vmean[j] = np.nanmean(Vmask[j], axis=0)

## Determine RMS quantities ##
uprime = np.ndarray([num_tests, num_images-1, sizey_mask, sizex_mask])
vprime = np.ndarray([num_tests, num_images-1, sizey_mask, sizex_mask])
uvprime  = np.ndarray([num_tests, num_images-1, sizey_mask, sizex_mask])
uvprime_mean = np.ndarray([num_tests, sizey_mask, sizex_mask])
urms = np.ndarray([num_tests, sizey_mask, sizex_mask])
vrms = np.ndarray([num_tests, sizey_mask, sizex_mask])
for j in range(0, num_tests):
    for jj in range(0, 499):
        uprime[j, jj] = ((Umask[j][jj]-umean[j]))
        vprime[j, jj] = ((Vmask[j][jj]-vmean[j]))
        uvprime[j, jj] = uprime[j, jj]*vprime[j, jj]
    uvprime_mean[j] = np.nanmean(uvprime[j], axis=0)
    urms[j] = np.nanmean(uprime[j]**2, axis=0)**(1/2)
    vrms[j] = np.nanmean(vprime[j]**2, axis=0)**(1/2)

## wall position adjustment ###########
#convert to m and take off wall position as seen in images
x = (x)/1000
y = (y-walloffset)/1000
xmask = x[ mask[2]:mask[3] ]
ymask = y[ mask[0]:mask[1] ]
## Create Mean Profiles #######
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
Neff = 10
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
hdf = pd.HDFStore('D:/test_061617/PIV_' + date + '.h5')
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
