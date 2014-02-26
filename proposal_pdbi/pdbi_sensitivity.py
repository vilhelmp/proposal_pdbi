#!/bin/python
"""
Script to calculate various useful things when preparing 
proposals for the Plateau de Bure Interferometer (PdBI)

Magnus Vilhelm Persson (magnusp@vilhelm.nu)

"""
from scipy.constants import c as sol
from scipy import array
from bisect import bisect_left, bisect_right
########################################################################
#                   PARAMETERS
#
# all you need to change below!
#
winter = False       # 0 for summer, 1 for winter
nu = 260.3E9        # frequency of line to observe in Hz
Na = 5              # number of antennas (6 winter, 5 summer atm)
Nc = 2              # number of configurations/tracks
Ton = 5             # on-source integration time
                    # observation time is ~1.6*Ton
dv = 0.1            # The requested velocity resolution
Np = 2              # Number of polarizations with the same setup
########################################################################
########################################################################
#                           CONSTANTS
# double check against each new announcment, this is from 2014 Summer Call
band_limits = array([[80, 116],
                     [129, 174],
                     [201, 267],
                     [277, 371]])*1E9
# conversion factor from Kelvin to Jansky (Jy/K)
JypK_summer = array([22, 29, 35, 45])
JypK_winter = array([22, 29, 35, 45])
# array defining the T_sys for typical season conditions
# and sources where decliation > 20 degrees
# each position in the array corresponds to a certain frequency
# interval:
#           <100 GHz
#           <116 GHz (115 GHz)
#           <150 GHz
#           <174 GHz
#           <267 GHz
#           <371 GHz
# whatever limit that comes first, top to bottom
T_sys_summer = array([100, 180, 150, 200, 250, 500])
T_sys_winter = array([100, 170, 130, 170, 200, 370])
T_sys_frequencies = array([110, 116, 150, 174, 267, 371])*1E9
# efficiency for the different bands
eta_summer = array([0.9, 0.8, 0.6, 0.5])
eta_winter = array([0.9, 0.85, 0.8, 0.7])
########################################################################
########################################################################
########################################################################
#                           PARSING
# get the Tsys for the observations
index = bisect_left(T_sys_frequencies, nu)
if winter:
    T_sys = T_sys_winter[index]
else:
    T_sys = T_sys_summer[index]
# get the efficiency
index  = (nu >= band_limits[:, 0]) * (nu <= band_limits[:, 1])
if winter:
    eta = eta_winter[index][0]
    JypK = JypK_winter[index][0]
    JypK = JypK_winter[index][0]
else:
    eta = eta_summer[index][0]
    JypK = JypK_summer[index][0]
########################################################################
#                           PRINTING
# pretty printing function
CSI = "\x1B["
end = CSI+'m'
blue = lambda s: '{0}{1}{2}{3}'.format(CSI, ';34m', s, end) 
# The speed of light in km/s
sol *= 1e-3 
# calculate the rms from the detection equation
# taken from the call for proposals doc. 
rms = JypK*T_sys/(eta*(Na*(Na-1)*Nc*Ton*3600*dv/sol*nu*Np)**0.5)*1e3
# calculate the rms in 1 km/s bins
sig = rms*(1./dv)**0.5*dv
# pretty print the information
print blue('\t  Sensitivity - PdBI')
print blue('#'*40)
print 'Frequency \t: {0:2.5f} Ghz'.format(nu*1e-9)
print 'Velocity res.\t: {0:2.2f} km/s ({1:2.3f} MHz)'.format(dv, dv/sol*nu*1e-6)
print 'RMS noise \t: {0:2.3f} mJy/channel ({1:2.3f} mJy/km/s)'.format(rms, sig)
print 'T_on  \t\t: {0:} hrs'.format(Ton)
print 'Tracks \t\t: {0}'.format(Nc)
print 'Season \t\t: {0}'.format(['Summer', 'Winter'][winter])
print blue('#'*40)
