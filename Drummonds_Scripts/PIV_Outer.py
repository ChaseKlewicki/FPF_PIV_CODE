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
def outer(date,num_tests,sizex,sizey,name):
    ## Initalize variables

    #read in variables
    umean = np.array(pd.read_hdf(name, 'umean_profile_avg'))*-1
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
    plt.figure(figsize = (10,10))
    plt.plot(umean, y, '-xb')
    plt.xlabel('U (m/sec)')
    plt.ylabel('Wall Normal Position (m)')
    plt.legend(['400rpm'], loc=0)
    plt.show()

    #V vs y
    plt.figure(figsize = (10,10))
    plt.plot(vmean, y, '-xr')
    plt.xlabel('V (m/sec)')
    plt.ylabel('Wall Normal Position (m)')
    plt.legend(['400rpm'], loc=0)
    plt.show()


    #urms vs y
    plt.figure(figsize = (10,10))
    plt.plot(urms, y, '-xb')
    plt.xlabel('$U_{rms}$ (m/sec)')
    plt.ylabel('Wall Normal Position (m)')
    plt.legend(['400rpm'], loc=0)
    plt.show()

    #vrms vs y
    plt.figure(figsize = (10,10))
    plt.plot(vrms, y, '-xr')
    plt.xlabel('$V_{rms}$ (m/sec)')
    plt.ylabel('Wall Normal Position (m)')
    plt.legend(['400rpm'], loc=0)
    plt.show()

    #uprime vs y
    plt.figure(figsize = (10,10))
    plt.plot(uvprime, y, '-xr')
    plt.xlabel('$u^,v^,$')
    plt.ylabel('Wall Normal Position (m)')
    plt.legend(['400rpm'], loc=0)
    plt.show()



    print('Done!')
