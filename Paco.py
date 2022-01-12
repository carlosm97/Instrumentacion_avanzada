#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 20 20:57:50 2021

@author: carlos
"""
#%% ------------------------------DATA-----------------------------------------
import matplotlib.pyplot as plt
import numpy as np
# Telescope:
D = 0.5 # m
f_D = 15 
F = f_D*D
# Lents and collimator:
f_1 = 220e-3 # m
d_1 = 38e-3 # m 
f_d_1 = f_1/d_1
f_2 = 220e-3 # m
d_2 = 38e-3 # m 
f_d_2 = f_2/d_2
#print(f_D_l)
# Grid
lenth = 40 # mm
sigma_1 = 300 # lines/mm
Blaze = 5 # º
Blaze_SI = Blaze*np.pi/180
alpha = 0 
seeing = 1 # "
lam = 5000 # \AA
# Detector: 
L = 150 # pixels 
l = 15 # microns/pixel
l *= 1e-6 # m/pixel
# In[i] Ancho del filtro

sigma = 1/sigma_1

beta = 2*Blaze_SI

# Orden en que lam tiene la máxima eficiencia 

m = sigma*1e7/lam*(np.sin(alpha)+np.sin(2*Blaze_SI+alpha))
lam_max = sigma*1e7/int(m)*(np.sin(alpha)+np.sin(2*Blaze_SI+alpha))

# Pureza espectral libre, de acuerdo con la f: 
    
D_lam = 2*lam/(2*int(m)+3)

# Considerando un criterio de ancho de filtro de +-1/3 o +-1/4:
print('Ancho del filtro ',D_lam,'\AA')


# In[ii] Pureza espectral dada por la red
Beta = np.arcsin(int(m)*lam*1e-7/sigma)

dl = lambda r_1, w_1, A_1: r_1*w_1/A_1/f_1

# Ángulo subtendido por la rendija sobre el cielo:
Phi = seeing/3600 # deg
Phi *= np.pi/180 # rad

# CONSIDERANDO QUE LA RED NO ESTÁ TOTALMENTE ILUMINADA:
    
W = f_1/f_D # m. Diámetro del cono de luz que ilumina el colimador. 
Beta = np.arcsin(int(m)*lam*1e-7/sigma) # Ángulo Beta, sacado de la ecuación de red. 

R = W*(np.sin(Beta)+np.sin(alpha))/Phi/D # Ec. 13.2.5 del Astronomical Optics
d_lam = lam/R
print('Pureza espectral dada por la red: ',d_lam,'\AA')



# In[iii]
# Dispersión angular:
#A1 = (np.sin(alpha)+np.sin(beta))/np.cos(beta)/lam # rad/\AA
A = int(m)/np.cos(beta)/sigma #rad/mm
A *= 1e-7 # rad/\AA
# Dispersión lineal:
P_1 = A*f_2 # m/\AA

# Tamaño del detector CCD:
T = L*l # m

# Rango espectral cubierto por el detector:
D_lam = T/P_1
print('Rango espectral cubierto por el detector CCD: ', D_lam,'\AA')

D_lam_pix = D_lam/L
# ¿Cumple el teorema del muestreo?
print('Si cumple el teorema del muestreo: ',2*D_lam_pix<d_lam)