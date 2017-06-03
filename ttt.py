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
plt.xlabel('T in K', size=15)
plt.ylabel('$\chi$ ~ U in V', size=15)
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
#plt.savefig('graphics/chifit.pdf')
plt.show()

'''
aufwärmen
'''

p=np.loadtxt('Tundp.txt')[52:69]
chi2=np.loadtxt('Tundp.txt')[69:]

def CW(chi):
    return paramsCW[1]/(chi-paramsCW[0])+paramsCW[2]

def deltaT(chi):
    return np.sqrt(errorCW[2]**2+(errorCW[1]/(chi-paramsCW[0]))**2+(paramsCW[1]*0.01/(chi-paramsCW[0])**2)**2+(paramsCW[1]*errorCW[0]/(chi-paramsCW[0])**2)**2)

T2=CW(chi2)
dT=deltaT(chi2)
plt.errorbar(T2,p,yerr=1,xerr=dT,fmt='x',label='Messwerte')
plt.xlabel('T in K', size=15)
plt.ylabel('p in Torr', size=15)
plt.grid(True)

def lin(beta,T):
    return beta[0]*T+beta[1]
    
d=11 #denn Punkt 12 liegt in beiden fits über dem lambda-Punkt
    
data = RealData(T2[:d],p[:d],dT[:d],1)
model = Model(lin)
odr = ODR(data, model, beta0=[65.8,-110])
#odr.set_job(fit_type=2)
outputL = odr.run()
paramsL=outputL.beta
errorL=outputL.sd_beta

l=np.linspace(1.71,2.25)

plt.plot(l,lin(paramsL,l),label='Fit p-T unterhalb $\lambda$-Punkt',zorder=10,color='r')

data = RealData(T2[d:],p[d:],dT[d:],1)
model = Model(lin)
odr = ODR(data, model, beta0=[8269,-17842])
#odr.set_job(fit_type=2)
outputL2 = odr.run()
paramsL2=outputL2.beta
errorL2=outputL2.sd_beta

plt.plot(l,lin(paramsL2,l),label='Fit p-T oberhalb $\lambda$-Punkt',zorder=10,color='g')
plt.xlim([1.7,2.3])
plt.ylim([0,450])
plt.legend(loc='best')
#plt.savefig('graphics/pT.pdf')
plt.show()

print('T_lambda = ',(paramsL2[1]-paramsL[1])/(paramsL[0]-paramsL2[0]), 'Messpunkt 12: ',T2[11])

DTl=np.sqrt((errorL2[1]/(paramsL[0]-paramsL2[0]))**2+(errorL[1]/(paramsL[0]-paramsL2[0]))**2+((paramsL2[1]-paramsL[1])*errorL[0]/(paramsL[0]-paramsL2[0])**2)**2+((paramsL2[1]-paramsL[1])*errorL2[0]/(paramsL[0]-paramsL2[0])**2)**2)

print('Delta T_lambda = ',DTl)

Tc=([1.7,1.77,1.81,1.85,1.90,1.95,2  ,2.03,2.06,2.1,2.18,2.20,2.21])
c =([2.4,2.9 ,3.2 ,3.5 ,4.1 ,4.6 ,5.4,6   ,6.7 ,8.4,3.2 ,3.1 ,3])

m=0.145*25.3*np.pi*4*4
x=0
for i in range(len(Tc)-1):
    x+=(c[i+1]+c[i])/2*(Tc[i+1]-Tc[i])