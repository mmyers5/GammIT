import numpy as np
import scipy.integrate as integrate

def integrand_tau(x,SFR,tau):
  integrand = np.exp(-x/tau)
  return integrand

def integrate_tau(age1,age2,SFR,tau):
  I = integrate.quad(integrand_tau,age1,age2,args=(SFR,tau))
  return I[0]

def calc_age(tau,M,SFR):
  age = tau*np.log(M/(tau*SFR) + 1)
  return age

def calc_age2(M,SFR):
  age = M/SFR
  return age

def calc_constant(age2,tau,M):
  B = tau*(np.exp(age2/tau)-1)
  C = M/B
  return C

def calc_SFR(age2,tau,M):
  SFR = calc_constant(age2,tau,M)*np.exp(0/tau)
  return SFR
  '''
  now = 13.8e9
  A = calc_constant(age2,tau,M)
  age = now-age2
  SFR = A*np.exp(-(now-age)/tau)
  return SFR
  '''
