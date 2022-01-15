#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 18:15:12 2021

@author: carlos
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os
# Constants definition
h = 6.626e-34 # J/s
c = 2.998e8 # m/s
# Problem data
D = 1.5 # m
r = D/2 # m
lamb = 550e-9 # m
pixel = 0.33 #arcsec/pixel

#----------------------FLUX OF A  0 MAGNITUDE OBJECT---------------------------
# In-band flux
F_ib0 = 3.565e-8 # W/m**2/\mu m**-1 = J/s/m**2/\mu m**-1  
# Average energy of a photon in the detector:
E_0 = h*c/lamb # J
# Number of photons per second
F_0 = F_ib0/E_0 # 1/s/m**2/\mu m**-1  
B = 0.1 # \mu m. Bandwidth
Q = 0.2 # photon/electron convertion
# Total flux at the focus of the telescope 
P_0 = F_0*B*Q*np.pi*r**2 # electrons/s

print('Total flux at the focus of the telescope of a 0 magnitude star: ',P_0)
mag_15 = True
if mag_15==True: # Case election (a or b)
    m15,m14 = 15.000, 14.999
elif mag_15==False:
    m15,m14 = 4.000, 3.999
FWHM = 1 

# -----------------------FLUX OF OUR STARS AND SKY-----------------------------
P = lambda m : P_0*10**(-0.4*m) # Function to compute the photons per second 

# For our magnitudes:
P_15,P_14 = P(m15),P(m14) # electrons/s
P_sky = P(22.5) # Sky background flux in electrons/arcsec^2/s

# electrons of sky per pixel:
P_s_p = P_sky*pixel*pixel
# diffraction limit
dif_lim=1.22*lamb/2/r #rad
dif_lim = dif_lim*180*3600/np.pi # arcsec

# sky contribution to each star (sky contribution to the considered apperture):
aperture_diameter = 1.2 # arcsec
aperture = np.pi*(aperture_diameter/2)**2 # arcsec**2
aperture_pixel = aperture*pixel**2 # pixel**2
sky_background = P_sky*aperture



#%% 


sigma = FWHM/2.355 # arcsec
sig_pix = sigma/pixel

func = lambda r,sig: 1-np.exp(-r**2/2/sig**2) # Function of the fraction of the 
# total area contained in a circle of radious r for a Gaussian distribution.

# Integral over 1.2":
rad = 0.6*FWHM
fraction = func(rad,sigma) # Fraction of the total counts contained
integrated_counts_14, integrated_counts_15 = P_14*fraction, P_15*fraction # Total counts. 



#-------------Total exposure time to reach the required limit.-----------------
t = (integrated_counts_14+integrated_counts_15+2*sky_background)/\
    (integrated_counts_14**2+integrated_counts_15**2-2*integrated_counts_14*integrated_counts_15)
t_ = (integrated_counts_14+integrated_counts_15+2*sky_background)/\
    (integrated_counts_14-integrated_counts_15)**2
dist = lambda total,r,sig: total/2/np.pi/sig**2*np.exp(-r**2/2/sig**2)

# Counts in the central pixel considering a gaussian distribution:
C_cent_14, C_cent_15 = dist(P_14,0,sig_pix), dist(P_15,0,sig_pix)
# We use sigma in pixels in order to get the distribution in this unit.
# consideing the sky
C_cent_14 += P_s_p
C_cent_15 += P_s_p

A, B, C = C_cent_14**2, -C_cent_14*(1+2e4), 1e8

t_max = (-B-np.sqrt(B**2-4*A*C))/2/A
print(C_cent_14*abs(t_max)+np.sqrt(C_cent_14*abs(t_max)))
# Number of exposures:

n = t/t_max

