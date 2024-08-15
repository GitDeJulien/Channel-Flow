import numpy as np
from scipy import fft
from scipy import special

from tools import *

def phi_22(omega, a1, Uc, sigma1c, Le): #same as R_33
    nw = omega.shape[0]
    res = np.zeros((nw))
    kc = omega / Uc
    #ke = 1./Le
    for w in range(nw):
        res[w] = a1 * 6*special.gamma(17/6.)*sigma1c**2*Le[w]*(3 + 8*(a1*kc[w]*Le[w])**2) / (np.sqrt(np.pi)*55*special.gamma(1/3.) * (1 + (a1*kc[w] * Le[w])**2)**(11/6.))
    return(kc, res)

def phi_11(omega, a1, Uc, sigma1c, Le):
    nw = omega.shape[0]
    res = np.zeros((nw))
    kc = omega / Uc
    ke = 1./Le
    for w in range(nw):
        res[w] = a1 * 36*special.gamma(17/6.)*sigma1c**2*Le[w] / (np.sqrt(np.pi)*55*special.gamma(1/3.) * (1 + (a1*kc[w]*Le[w])**2)**(5/6.))
    return(kc, res)

def L(C, kt, omega):
    nw = omega.shape[0]
    res = np.zeros((nw))
    eps = omega * 0.09 * kt
    for w in range(nw):
        res[w] = C*kt[w]**(2./3)/eps[w]
    return(res)


def sigma_1c(sigma1, w1, sigma2, w2, sigma3, w3):
    res = w1*sigma1 + w2*sigma2 + w3*sigma3
    # nw = sigma1.shape[0]
    # res = np.zeros_like(sigma1)
    # for w in range(nw):
    #     res[w] = w1*sigma1[w] + w2*sigma2[w] + w3*sigma3[w]
    return(res)


def sigma_power2(data_fluct, axis = "streamwise"):
    if axis == "streamwise" or axis == "wallnormal":
        res = np.mean(np.mean(np.mean(data_fluct[:,:,:]**2, axis=0), axis=1))
        
    if axis == "spanwise":
        res = np.mean(np.mean(np.mean(data_fluct[:,:,:]**2, axis=0), axis=0))
        
    return(res)
        
        
######################################
######################################


def Lambda2_22(a1, a2, ke, omega, Uc):
    kc = Uc * omega
    res = a2 * 55 * np.sqrt(np.pi) * special.gamma(1/3.) * (3*ke**2 + 11*(a1*kc)**2) * (ke**2 + (a1*kc)**2)**(-1/2.) /(108 * special.gamma(17/6.) * (3*ke**2 + 8*(a1*kc)**2))
    return(kc, res)
    
    
def Lambda2_33(a1, a2, ke, omega, Uc):
    kc = Uc * omega
    res = 55 * np.sqrt(np.pi) * special.gamma(1/3.) * (a1**3*a2*kc**2) * (ke**2 + (a1*kc)**2)**(-1/2.) /(a1*6 * special.gamma(17/6.) * (3*ke**2 + 8*(a1*kc)**2))
    return(kc, res)


def Lambda2_11(a1, a2, ke, omega, Uc):
    kc = Uc * omega
    res = 55 * np.sqrt(np.pi) * special.gamma(1/3.) * (a1*a2) * (ke**2 + (a1*kc)**2)**(-1/2.) /(216 * a1 * special.gamma(17/6.))
    return(kc, res)
    