import numpy as np
from scipy import signal
from scipy import fft
from scipy import special

from tools import *

def R_22(omega, a1, Uc, sigma1c, ke): #same as R_33
    kc = Uc / omega
    res = a1 * 6*special.gamma(17/6.)*sigma1c**2*ke**(2/3)*(3*ke**2 + 8*(a1*kc)**2) / (Uc*np.sqrt(np.pi)*55*special.gamma(1/3.) * (ke**2 + (a1*kc)**2)**(11/6))
    return(kc, res)

def R_11(omega, a1, Uc, sigma1c, ke):
    kc = Uc / omega
    res = a1 * 36*special.gamma(17/6.)*sigma1c**2*ke**(2/3) / (Uc*np.sqrt(np.pi)*55*special.gamma(1/3.) * (ke**2 + (a1*kc)**2)**(5/6))
    return(kc, res)

def Lij(const, kt, eps_ij):
    res = const*kt**(2/3)/eps_ij
    return(res)


def sigma1c(sigma1, w1, sigma2, w2, sigma3, w3):
    res = w1*sigma1**2 + w2*sigma2**2 + w3*sigma3**2
    return(res)


def sigma(datas, axis = "streamwise"):
    if axis == "streamwise" or axis == "wallnormal":
        res = fft.fft(datas[:,:,:], axis=0)
        res = np.mean(res, axis=2)
        
    if axis == "spanwise":
        res = fft.fft(datas[:,:,:], axis=0)
        res = np.mean(res, axis=1)
    
    return(res)
        

    
    


def Lambda2_22(a1, a2, ke, omega, Uc):
    kc = Uc / omega
    res = a2 * 55 * np.sqrt(np.pi) * special.gamma(1/3.) * (3*ke**2 + 11*(a1*kc)**2) * (ke**2 + (a1*kc)**2)**(-1/2.) /(108 * special.gamma(17/6.) * (3*ke**2 + 8*(a1*kc)**2))
    return(kc, res)
    
    
def Lambda2_33(a1, a2, ke, omega, Uc):
    kc = Uc / omega
    res = 55 * np.sqrt(np.pi) * special.gamma(1/3.) * (a1**3*a2*kc**2) * (ke**2 + (a1*kc)**2)**(-1/2.) /(a1*6 * special.gamma(17/6.) * (3*ke**2 + 8*(a1*kc)**2))
    return(kc, res)


def Lambda2_11(a1, a2, ke, omega, Uc):
    kc = Uc / omega
    res = 55 * np.sqrt(np.pi) * special.gamma(1/3.) * (a1*a2) * (ke**2 + (a1*kc)**2)**(-1/2.) /(216 * a1 * special.gamma(17/6.))
    return(kc, res)
    