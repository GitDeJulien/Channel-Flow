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
        res[w] = a1 * sigma1c**2 * Le[w]*special.gamma(5./6) / (np.sqrt(np.pi)*special.gamma(1./3) * (1 + (a1*kc[w]*Le[w])**2)**(5/6.)) * (4./3 - 5/(6*(1+(a1*kc[w]*Le[w])**2)))
    return(kc, res)

def phi_11(omega, a1, Uc, sigma1c, Le):
    nw = omega.shape[0]
    res = np.zeros((nw))
    kc = omega / Uc
    #ke = 1./Le
    for w in range(nw):
        res[w] = a1 * sigma1c**2 * Le[w]*special.gamma(5./6) / (np.sqrt(np.pi)*special.gamma(1./3) * (1 + (a1*kc[w]*Le[w])**2)**(5/6.))
    return(kc, res)

# def phi_22(omega, a1, Uc, sigma1c, Le): #same as R_33
#     nw = omega.shape[0]
#     res = np.zeros((nw))
#     kc = omega / Uc
#     #ke = 1./Le
#     for w in range(nw):
#         res[w] = a1*sigma1c**2*Le[w]*(1 + 8./3*(special.gamma(1./3)*a1*kc[w]*Le[w])**2) / (np.pi* (1 + (special.gamma(1./3)*a1*kc[w]*Le[w])**2)**(11/6.))
#     return(kc, res)

# def phi_11(omega, a1, Uc, sigma1c, Le):
#     nw = omega.shape[0]
#     res = np.zeros((nw))
#     kc = omega / Uc
#     #ke = 1./Le
#     for w in range(nw):
#         res[w] = a1 * sigma1c**2*Le[w] / (np.pi * (1 + (1.339*a1*kc[w]*Le[w])**2)**(5/6.))
#     return(kc, res)

# def phi_22(omega, a1, Uc, sigma1c, Le): #same as R_33
#     nw = omega.shape[0]
#     res = np.zeros((nw))
#     kc = omega / Uc
#     #ke = 1./Le
#     for w in range(nw):
#         res[w] = a1 * 6*special.gamma(17/6.)*sigma1c**2*Le[w]*(3 + 8*(a1*kc[w]*Le[w])**2) / (np.sqrt(np.pi)*55*special.gamma(1/3.) * (1 + (a1*kc[w] * Le[w])**2)**(11/6.))
#     return(kc, res)

# def phi_11(omega, a1, Uc, sigma1c, Le):
#     nw = omega.shape[0]
#     res = np.zeros((nw))
#     kc = omega / Uc
#     #ke = 1./Le
#     for w in range(nw):
#         res[w] = a1 * 36*special.gamma(17/6.)*sigma1c**2*Le[w] / (np.sqrt(np.pi)*55*special.gamma(1/3.) * (1 + (a1*kc[w]*Le[w])**2)**(5/6.))
#     return(kc, res)

def L(C, kt, omega):
    nw = omega.shape[0]
    res = np.zeros((nw))
    eps = omega * 0.09 * kt
    for w in range(nw):
        res[w] = C*kt[w]**(3./2)/eps[w]
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
    