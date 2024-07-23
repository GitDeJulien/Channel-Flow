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
    