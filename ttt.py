# -*- coding: utf-8 -*-
"""
Created on Wed May 31 17:03:43 2017

@author: Jonas
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, RealData
import matplotlib as mpl
mpl.rc("figure", figsize=(12, 6))

T=np.loadtxt('Tundp.txt')[:26]
chi=np.loadtxt('Tundp.txt')[26:]

plt.errorbar(T,chi,yerr=0.01,xerr=0.04,fmt='x',label='gemessene Werte')
plt.xlabel('T in K')
plt.ylabel('$\chi$ ~ U in V')
plt.grid(True)

def Curie(beta,T):
    return beta[0]+beta[1]/T
    
def CurieWeiss(beta,T):
    return beta[0]+beta[1]/(T-beta[2])
    
data = RealData(T,chi,0.04,0.01)
model = Model(Curie)
odr = ODR(data, model, beta0=[0,6])
#odr.set_job(fit_type=2)
outputC = odr.run()
paramsC=outputC.beta
errorC=outputC.sd_beta

plt.plot(T,Curie(paramsC,T),label='Fit Curie-Gesetz')

model = Model(CurieWeiss)
odr = ODR(data, model, beta0=[0,6,0])
#odr.set_job(fit_type=2)
outputCW= odr.run()
paramsCW=outputCW.beta
errorCW=outputCW.sd_beta

plt.plot(T,CurieWeiss(paramsCW,T),label='Fit Curie-Weiss-Gesetz')
plt.legend(loc='best')
#plt.savefig('chifit.pdf')
plt.show()