import numpy as np
from tqdm import tqdm
from scipy import integrate

from tools import *
from plot import *



def frozen_turbulence(datas, zplan, z, nt, split_time, dt, n1, dx = None, tEnd = None, tStart = None, ch = "spectra"):
    
    if split_time == 'Y':
        split_t = 1000
        num_split_t = nt // split_t
    
        if ch == "spectra":

            var1 = 0.
            inte_space = 0.
            inte_time = 0.
            
            for n in tqdm(range(1,num_split_t), desc=f'PSD (z={z[zplan]:.2f})', colour= 'GREEN'):

                var1 += np.mean(np.mean(np.var(datas[(n-1)*split_t:n*split_t,:,:], axis=0), axis=-1))

                freq, psd, time_spectra = WelshPowerSpectra(datas[(n-1)*split_t:n*split_t,:,:], dt, dx, geom = 'plan', axis='time')
                omega = 2*np.pi*freq
                inte_time += integrate.simpson(time_spectra, freq)
                freq, psd, space_spectra = WelshPowerSpectra(datas[(n-1)*split_t:n*split_t,:,:], dt, dx, geom = 'plan', axis='streamwise')
                k = 2*np.pi*freq
                inte_space += integrate.simpson(space_spectra, freq)
            
            
            var1 /= (num_split_t-1)
            inte_time /= (num_split_t-1)
            inte_space /= (num_split_t-1)
            print('<uiui>=',var1)
            print('Inte_time:', inte_time)
            print('Inte_space:', inte_space)
            
            return(omega, k, time_spectra, space_spectra)
            
        if ch == 'corr':
            
            Dt = np.linspace(0,(tEnd-tStart)*dt, split_t)
            Dx = np.linspace(0,np.pi,n1)
            
            R_time = np.zeros((split_t))
            R_space = np.zeros((n1))
                
            for n in tqdm(range(1,num_split_t), desc=f'Corr 1D (z={z[zplan]:.2f})', colour= 'GREEN'):
                
                R_time += Autocorrelation_1d(datas[(n-1)*split_t:n*split_t,:,:], mode_corr='half', geom='plan', axis='time')
                R_space += Autocorrelation_1d(datas[(n-1)*split_t:n*split_t,:,:], mode_corr='half', geom='plan', axis='streamwise')
        
            R_space /= (num_split_t-1)
            R_time /= (num_split_t-1)
            
            return(Dt, Dx, R_time, R_space)
        
        
    if split_time == 'n':
        
        if ch == "spectra":
            var1 = np.mean(np.mean(np.var(datas[:,:,:], axis=0), axis=-1))

            freq, psd, time_spectra = WelshPowerSpectra(datas[:,:,:], dt, dx, geom = 'plan', axis='time')
            omega = 2*np.pi*freq
            inte_time = integrate.simpson(time_spectra, freq)
            freq, psd, space_spectra = WelshPowerSpectra(datas[:,:,:], dt, dx, geom = 'plan', axis='streamwise')
            k = 2*np.pi*freq
            inte_space = integrate.simpson(space_spectra, freq)
            
            print('<uiui>=',var1)
            print('Inte_time:', inte_time)
            print('Inte_space:', inte_space)
            
            return(omega, k, time_spectra, space_spectra)
            
        if ch == 'corr':
            
            R_time = np.zeros((split_t))
            R_space = np.zeros((n1))
                
            R_time = Autocorrelation_1d(datas[:,:,:], mode_corr='half', geom='plan', axis='time')
            R_space = Autocorrelation_1d(datas[:,:,:], mode_corr='half', geom='plan', axis='streamwise')

            return(Dt, Dx, R_time, R_space)
        
    
    