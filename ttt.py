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
chi=np.loadtxt('Tundp.txt')[26:52]

plt.errorbar(T,chi,yerr=0.01,xerr=0.04,fmt='x',label='Messwerte')
plt.xlabel('T in K')
plt.ylabel('$\chi$ ~ U in V')
plt.grid(True)

def Curie(beta,T):
    return beta[0]+beta[1]/T
    
def CurieWeiss(beta,T):
    return beta[0]+beta[1]/(T-beta[2])
    
data = RealData(T,chi,0.04,0.01)
model = Model(Curie)
odr = ODR(data, model, beta0=[0.686,5.53])
#odr.set_job(fit_type=2)
outputC = odr.run()
paramsC=outputC.beta
errorC=outputC.sd_beta

plt.plot(T,Curie(paramsC,T),label='Fit Curie-Gesetz')

model = Model(CurieWeiss)
odr = ODR(data, model, beta0=[0.184,4.861,0.175])
#odr.set_job(fit_type=2)
outputCW= odr.run()
paramsCW=outputCW.beta
errorCW=outputCW.sd_beta

plt.plot(T,CurieWeiss(paramsCW,T),label='Fit Curie-Weiss-Gesetz')
plt.legend(loc='best')
#plt.savefig('chifit.pdf')
plt.show()

'''
aufwärmen
'''

p=np.loadtxt('Tundp.txt')[52:68]
chi2=np.loadtxt('Tundp.txt')[68:]

def CW(chi):
    return paramsCW[1]/(chi-paramsCW[0])+paramsCW[2]

T2=CW(chi2)
plt.errorbar(T2,p,yerr=1,xerr=0.04,fmt='x',label='Messwerte')
plt.xlabel('T in K')
plt.ylabel('p in Torr')
plt.grid(True)

def lin(beta,T):
    return beta[0]*T+beta[1]
    
d=10 #denn Punkt 11 liegt in beiden fits über dem lambda-Punkt
    
data = RealData(T2[:d],p[:d],0.04,1)
model = Model(lin)
odr = ODR(data, model, beta0=[67.8,-110])
#odr.set_job(fit_type=2)
outputL = odr.run()
paramsL=outputL.beta
errorL=outputL.sd_beta

l=np.linspace(1.71,2.25)

plt.plot(l,lin(paramsL,l),label='Fit p-T unterhalb $\lambda$-Punkt',zorder=10,color='r')

data = RealData(T2[d:],p[d:],0.04,1)
model = Model(lin)
odr = ODR(data, model, beta0=[8183,-17652])
#odr.set_job(fit_type=2)
outputL2 = odr.run()
paramsL2=outputL2.beta
errorL2=outputL2.sd_beta

plt.plot(l,lin(paramsL2,l),label='Fit p-T oberhalb $\lambda$-Punkt',zorder=10,color='g')
plt.xlim([1.7,2.3])
plt.ylim([0,450])
plt.legend(loc='best')
#plt.savefig('pT.pdf')
plt.show()

print('T_lambda = ',(paramsL2[1]-paramsL[1])/(paramsL[0]-paramsL2[0]), 'Messpunkt 11: ',T2[10])

def deltaT(chi):
    return np.sqrt(errorCW[2]**2+(errorCW[1]/(chi-paramsCW[0]))**2+(paramsCW[1]*0.01/(chi-paramsCW[0])**2)**2+(paramsCW[1]*errorCW[0]/(chi-paramsCW[0])**2)**2)