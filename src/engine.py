import numpy as np
from tqdm import tqdm
from scipy import integrate

from tools import *
from plot_figures import *
from parameters import *



def frozen_turbulence(datas, zplan, z, nt, split_time, dt, n1, dx = None, x1 = None, Uc = None, delta_x=None, ch = "spectra", scaling='density'):
    
    if split_time == 'Y':
        num_split_t = nt // split_t
    
        if ch == "spectra":

            var1 = 0.
            inte_space = 0.
            inte_time = 0.
            sum_space = 0.
            
            for n in tqdm(range(1,num_split_t), desc=f'PSD (z={z[zplan]:.2f})', colour= 'GREEN'):

                var1 += np.mean(np.mean(np.var(datas[(n-1)*split_t:n*split_t,:,:], axis=0), axis=-1))

                freq, psd, time_spectra = WelshPowerSpectra(datas[(n-1)*split_t:n*split_t,:,:], dt, dx, geom = 'plan', axis='time', scaling=scaling)
                omega = 2*np.pi*freq
                inte_time += integrate.simpson(time_spectra, freq)
                freq, psd, space_spectra = WelshPowerSpectra(datas[(n-1)*split_t:n*split_t,:,:], dt, dx, geom = 'plan', axis='streamwise', scaling=scaling)
                k = 2*np.pi*freq
                inte_space += integrate.simpson(space_spectra, freq)
            
            
            var1 /= (num_split_t-1)
            inte_time /= (num_split_t-1)
            inte_space /= (num_split_t-1)
            sum_space /= (num_split_t-1)
            print('<uiui>=',var1)
            print('Inte_time:', inte_time)
            print('Inte_space:', inte_space)
            
            return(omega, k, time_spectra, space_spectra)
            
        if ch == 'corr':
            
            #Dt = np.linspace(0,(tEnd-tStart)*dt, split_t)
            Dt = np.linspace(0, split_t*dt, split_t)
            Dx = np.linspace(0,xlen//2,n1)
            
            R_time = np.zeros((split_t))
            R_space = np.zeros((n1))
                
            for n in tqdm(range(1,num_split_t), desc=f'Corr 1D (z={z[zplan]:.2f})', colour= 'GREEN'):
                
                R_time += Autocorrelation_1d(datas[(n-1)*split_t:n*split_t,:,:], mode_corr='half', geom='plan', axis='time')
                R_space += Autocorrelation_1d(datas[(n-1)*split_t:n*split_t,:,:], mode_corr='half', geom='plan', axis='streamwise')
        
            R_space /= (num_split_t-1)
            R_time /= (num_split_t-1)
            
            ind1 = 0
            ind2 = 0
            for i in range(Dt.shape[0]):
                if R_time[i]<0.0:
                    ind1 = i
                    break
            
            for i in range(Dx.shape[0]):
                if R_space[i]<0.0:
                    ind2 = i
                    break
                
            if ind1 == 0:
                ind1 = Dt.shape[0]
            
            if ind2 == 0:
                ind2 = Dx.shape[0]
                
            
            return(ind1, ind2, Dt, Dx, R_time, R_space)
        
        if ch == 'corr2d':
            
            R = np.zeros((split_t, n1))  
            R_full = np.zeros((2*split_t, 2*n1)) 
            
            Dt = np.linspace(-split_t*dt,split_t*dt,2*split_t)

            Dx = np.linspace(-xlen, xlen, 2*n1)
            
            # Dt = np.linspace((tEnd-tStart)*dt//2, (tEnd-tStart)*dt//2, split_t)

            # Dx = np.linspace(-xlen//2, xlen//2, n1)
            
            for n in tqdm(range(1,num_split_t), desc=f'Autocorr 2D', colour= 'GREEN'):
                
                R += Correlation_2d(datas[(n-1)*split_t:n*split_t,:,:], geom='plan', axis='streamwise')[1]
                #R_full += Correlation_2d(datas[(n-1)*split_t:n*split_t,:,:], geom='plan', axis='streamwise')[1]
                
            R /= (num_split_t-1)
            
            R_full[0:split_t, 0:n1] = R
            R_full[split_t:2*split_t, 0:n1] = R
            R_full[0:split_t, n1:2*n1] = R
            R_full[split_t:2*split_t, n1:2*n1] = R
            
            levels = [0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]
            coef = get_ellispses_slop(R, levels , 0.01, Dt, Dx)
            
            
            return(Dt, Dx, R_full, coef)
            #return(Dt, Dx, R, coef)
        
        if ch == 'gamma':
            frequency,_,_,_,_,_ = Crosscorrelation_2d(datas[0:split_t,:,:], delta_x, x1, dt, Uc, geom = "plan", axis = "streamwise")
            length = frequency.shape[0]
            funct = np.zeros((length))
            kc = np.zeros((length))
            
            
            for n in tqdm(range(1,num_split_t), desc=f'Crosscorr', colour= 'GREEN'):
                funct += Crosscorrelation_2d(datas[(n-1)*split_t:n*split_t,:,:], delta_x, x1, dt, Uc, geom = "plan", axis = "streamwise")[3]
                kc += Crosscorrelation_2d(datas[(n-1)*split_t:n*split_t,:,:], delta_x, x1, dt, Uc, geom = "plan", axis = "streamwise")[4]
            
            funct /= (num_split_t-1)
            kc /= (num_split_t-1)
            omega = 2*np.pi*frequency
            
            #omega_lim = 2*np.pi*Uc/dx
            omega_lim = 80 ##have to be computed (mesh size?)
            for i in range(omega.shape[0]):
                if omega[i]>omega_lim:
                    ind1 = i
                    break
                
            omega_lim = 40
            for i in range(omega.shape[0]):
                if omega[i]>omega_lim:
                    ind2 = i
                    break
            
            
            return(ind1, ind2, kc, omega, funct)

        
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
            
            Dt = np.linspace(0, nt*dt, nt)
            Dx = np.linspace(0,xlen//2,n1)
            
            R_time = np.zeros((nt))
            R_space = np.zeros((n1))
                
            R_time = Autocorrelation_1d(datas[:,:,:], mode_corr='half', geom='plan', axis='time')
            R_space = Autocorrelation_1d(datas[:,:,:], mode_corr='half', geom='plan', axis='streamwise')

            ind1 = 0
            ind2 = 0
            for i in range(Dt.shape[0]):
                if R_time[i]<0.0:
                    ind1 = i
                    break
            
            for i in range(Dx.shape[0]):
                if R_space[i]<0.0:
                    ind2 = i
                    break
                
            if ind1 == 0:
                ind1 = Dt.shape[0]
            
            if ind2 == 0:
                ind2 = Dx.shape[0]
                
            
            return(ind1, ind2, Dt, Dx, R_time, R_space)
        
        if ch == 'corr2d':
            
            R = np.zeros((nt, n1))  
            R_full = np.zeros((2*nt, 2*n1)) 
            
            Dt = np.linspace(-nt*dt,nt*dt,2*nt)

            Dx = np.linspace(-xlen, xlen, 2*n1)
            
                
            R += Correlation_2d(datas[:,:,:], geom='plan', axis='streamwise')[1]
                
            R_full[0:nt, 0:n1] = R
            R_full[nt:2*nt, 0:n1] = R
            R_full[0:nt, n1:2*n1] = R
            R_full[nt:2*nt, n1:2*n1] = R
            
            levels = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95]
            coef = get_ellispses_slop(R, levels, 0.001, Dt, Dx)
            
            return(Dt, Dx, R_full, coef)
        
        
        if ch == 'gamma':
            
            frequency = Crosscorrelation_2d(datas[:,:,:], dx, x1, dt, Uc, geom = "plan", axis = "streamwise")[1]
            funct = Crosscorrelation_2d(datas[:,:,:], dx, x1, dt, Uc, geom = "plan", axis = "streamwise")[3]
            kc = Crosscorrelation_2d(datas[:,:,:], dx, x1, dt, Uc, geom = "plan", axis = "streamwise")[4]
            

            omega = 2*np.pi*frequency
            
            #omega_lim = 2*np.pi*Uc/dx
            omega_lim = 100
            for i in range(omega.shape[0]):
                if omega[i]>omega_lim:
                    ind = i
                    break
            
            return(ind, kc, omega, funct)
        

