import numpy as np
from scipy import signal
from scipy import fft
from tqdm import tqdm


from read_fpar import *

font = dict(family="serif",
        size=16,
        color="black")

#====================================================================
#========= Get number of intervals for welch method =================
#====================================================================

def get_nperseg(len_sig):
    """ Get the appropriate length of the welch segment
    INPUT: - len_sig (int) : length of the signal

    OUTPUT: - nperseg (int) : signal length to be used
    """
    nperseg_min = 512 #minimum signal length
    # recommended signal length if the signal is longer than the minimum length
    nperseg_des = 2**(np.ceil(np.log2(len_sig))-3)
    if(nperseg_des>=nperseg_min):
        nperseg = nperseg_des
    else:
        nperseg=len_sig

    return int(nperseg)

#====================================================================
#===================== Get ellipses slop ============================
#====================================================================

def get_ellispses_slop(R, levels, eps, Dt, Dx, n1, split_t = None):
    """ Get the appropriate slop of ellispses for the 2D Correlation contour plot 
    INPUT: - R (np.array) : 2D correlation
        - levels (list) : list of relevent ellipse values
        - eps (float) : precision of ellipse points
        - Dt, Dx (np.array) : time and space series
        - n1 (int) : number of points
        - split_t (None or int) : number of time iteration if the time is splited
        
    OUTPUT: - coef (tuple (float, float)) : a and b coeffitient of the fited line y = a*x+b
    """
    
    x_slop = []
    y_slop = []
    
    for level in levels:
        imax = 0
        jmax = 0
        for i in range(R.shape[0]//2):
            for j in range(R.shape[1]//2):
                if R[i,j] > level - eps and R[i,j] < level + eps:
                    if i>imax:
                        imax = i
                    if j>jmax:
                        jmax = j
        y_slop.append(Dt[imax])
        x_slop.append(Dx[jmax])
        
    
    # for level in levels:
    #     imin = R.shape[0]
    #     jmin = R.shape[1]
    #     for i in range(R.shape[0]//2 , R.shape[0]):
    #         for j in range(R.shape[1]//2 , R.shape[1]):
    #             if R[i,j] > level - eps and R[i,j] < level + eps:
    #                 if i<imin:
    #                     imin = i
    #                 if j<jmin:
    #                     jmin = j
    #     y_slop.append(Dt[split_t-1 + imin])
    #     x_slop.append(Dx[n1-1+jmin])

        
    coef = np.polyfit(x_slop, y_slop, 1)
    
    return(coef)

#====================================================================
#======================= Welsh Power Spectra ========================
#====================================================================

def WelshPowerSpectra(data, dt, dx, geom = 'plan', axis ='time'):
    """ Compute the frequency spectra psd(w, x) and when averaged psd_avg(w) 
    or the spatial spectra psd(t, k), when averaged psd_avg(k) """ 
    
    if geom == 'line':
        if axis == 'time':
            niters, npts  = data.shape
            nperseg_time = get_nperseg(niters)
            freq, Pxx = signal.welch(data, fs = 1/dt, window='hann', nperseg=nperseg_time, scaling='density', axis=0)
            psd = Pxx
            psd_avg = np.mean(psd, axis=1)
            
        elif axis == 'streamwise' or axis == 'spanwise' or axis == 'wallnormal':
            niters, npts  = data.shape
            nperseg_space = get_nperseg(npts)
            freq, Pxx = signal.welch(data, fs = 1/dx, window='hann', nperseg=nperseg_space, scaling='density', axis=1)
            psd = Pxx
            psd_avg = np.mean(psd, axis=0)
            
        else:
            print('Error: The argument axis is incorrect. Please choose between \'time\', \'streamwise\', \'spanwise\' and \'wallnormal\'.')
            exit(-1)
            
    
    elif geom == 'plan':
        if axis == 'time':
            niters, npts, nlines = data.shape
            nperseg_time = get_nperseg(niters)
            freq, Pxx = signal.welch(data, fs = 1/dt, window='hann', nperseg=nperseg_time, scaling='density', axis=0)
            psd = np.mean(Pxx, axis=2)
            psd_avg = np.mean(psd, axis=1)  
            
        elif axis == 'streamwise' or axis == 'wallnormal':
            niters, npts, nlines = data.shape
            nperseg_space = get_nperseg(npts)
            freq, Pxx = signal.welch(data, fs = 1/dx, window='hann', nperseg=nperseg_space, scaling='density', axis=1)
            psd = np.mean(Pxx, axis=2)
            psd_avg = np.mean(psd, axis=0) 
            
        elif axis == 'spanwise':
            niters, nlines, npts = data.shape
            nperseg_space = get_nperseg(npts)
            freq, Pxx = signal.welch(data, fs = 1/dx, window='hann', nperseg=nperseg_space, scaling='density', axis=2)
            psd = np.mean(Pxx, axis=1)
            psd_avg = np.mean(psd, axis=0)
            
        else:
            print('Error: The argument axis is incorrect. Please choose between \'time\', \'streamwise\', \'spanwise\' and \'wallnormal\'.')
            exit(-1)
            
            
    else:
        print('Error: Argument \'geom\' is not valid. Please select \'line\' or \'plan\'.')
        exit(0)
        
    return(freq, psd, psd_avg)


#====================================================================
#======================= Autocorrelation 1D =========================
#====================================================================

def Autocorrelation_1d(data, geom = 'plan', mode_corr = 'half', axis = 'time'):
    
    if geom == 'line':
        if axis == 'time':
            niters, npts = data.shape
            full_corr = np.zeros((2*niters-1))
            R = np.zeros((niters))
            for pts in range(npts):
                full_corr[:] += signal.correlate(data[:,pts], data[:,pts], mode='full', method='auto')
            full_corr /= (npts)
            full_corr /= max(full_corr)
            R[:] = full_corr[full_corr.size //2 :]
            
        elif axis == 'streamwise' or axis == 'spanwise' or axis == 'wallnormal':
            niters, npts = data.shape
            full_corr = np.zeros((2*npts-1))
            R = np.zeros((npts))
            for iters in range(niters):
                full_corr[:] += signal.correlate(data[iters,:], data[iters,:], mode='full', method='auto')
            full_corr /= (niters)
            full_corr /= max(full_corr)
            R[:] = full_corr[full_corr.size //2 :]
            
        else:
            print('Error: The argument axis is incorrect. Please choose between \'time\', \'streamwise\', \'spanwise\' and \'wallnormal\'.')
            exit(-1)
        
        if mode_corr == 'full':
            return(full_corr)
        elif mode_corr == 'half':
            return(R)
        else:
            print('Error: The argument mode_corr is incorrect. Please choose between \'full\' and \'half\'.')
            exit(-1)
        
    elif geom == 'plan':
        if axis == 'time':
            niters, npts, nlines = data.shape
            full_corr = np.zeros((2*niters-1))
            R = np.zeros((niters))
            for lines in range(nlines):
                for pts in range(npts):
                    full_corr[:] += signal.correlate(data[:,pts,lines], data[:,pts,lines], mode='full', method='auto')
            full_corr /= (npts*nlines)
            full_corr /= max(full_corr)
            R[:] = full_corr[full_corr.size //2 :]
            
        elif axis == 'streamwise' or axis == 'wallnormal':
            niters, npts, nlines = data.shape
            full_corr = np.zeros((2*npts-1))
            R = np.zeros((npts))
            for lines in range(nlines):
                for iters in range(niters):
                    full_corr[:] += signal.correlate(data[iters,:,lines], data[iters,:,lines], mode='full', method='auto')
            full_corr /= (niters*nlines)
            full_corr /= max(full_corr)
            R[:] = full_corr[full_corr.size //2 : ]
            
        elif axis == 'spanwise':
            niters, nlines, npts = data.shape
            full_corr = np.zeros((2*npts-1))
            R = np.zeros((npts))
            for lines in range(nlines):
                for iters in range(niters):
                    full_corr[:] += signal.correlate(data[iters,lines,:], data[iters,lines,:], mode='full', method='auto')
            full_corr /= (niters*nlines)
            full_corr /= max(full_corr)
            R[:] = full_corr[full_corr.size //2 :]
            
        else:
            print('Error: The argument axis is incorrect. Please choose between \'time\', \'streamwise\', \'spanwise\' and \'wallnormal\'.')
            exit(-1)
        
        if mode_corr == 'full':
            return(full_corr)
        elif mode_corr == 'half':
            return(R)
        else:
            print('Error: The argument geom is incorrect. Please choose between \'line\' and \'plan\'.')
            exit(-1)
            
    else:
        print('Error: The argument mode_corr is incorrect. Please choose between \'full\' and \'half\'.')
        exit(-1)
        


#====================================================================
#======================= Autocorrelation 2D =========================
#====================================================================        
def Correlation_2d(data, geom = 'plan', axis = 'streamwise'):
    
    if geom == 'line':
        niters, npts = data.shape
        phi = np.zeros((niters, npts))
        R = np.zeros((niters, npts))
        fourier = fft.fft2(data[:,:], workers = 3)
        phi += np.real(fourier * np.conj(fourier))
        R = np.real(fft.ifft2(phi, workers=3))
        R /= np.max(R, axis=(0,1))
        
    if geom == 'plan':
        if axis == 'streamwise' or axis == 'wallnormal':
            niters, npts, nlines = data.shape
            phi = np.zeros((niters, npts), dtype=complex)
            R = np.zeros((niters, npts))
            for lines in range(nlines):
                fourier = fft.fft2(data[:,:,lines], workers=3)
                phi += fourier * np.conj(fourier)
            phi /= nlines
            R = np.real(fft.ifft2(phi, workers=3))
            R /= np.max(R, axis=(0,1))
            
        if axis == 'spanwise':
            niters, nlines, npts = data.shape
            phi = np.zeros((niters, npts))
            R = np.zeros((niters, npts))
            for lines in range(nlines):
                fourier = fft.fft2(data[:,lines,:], workers=3)
                phi += np.real(fourier * np.conj(fourier))
            phi /= nlines
            R = np.real(fft.ifft2(phi, workers=3))
            R /= np.max(R, axis=(0,1))
        
    return(phi, R)


#====================================================================
#======================= Crosscorelation 2D =========================
#==================================================================== 

def Crosscorrelation_2d(datas, delta_x, coords, dt, Uc, geom = "plan", axis = "streamwise"):
    
    if geom == 'line':
        niters, npts = datas.shape
        nperseg_time = get_nperseg(niters)
        
        r = np.abs(coords[delta_x] - coords[0])

        frequency, psd = signal.welch(datas[:, :], fs=1/dt, window='hann', nperseg=nperseg_time//8, scaling='density', axis=0)
        psd = np.mean(psd, axis=1)

        avg_corr = np.zeros((niters))
        cnt = 0
        for x in range(0, npts-delta_x-1, 2):
            corr = signal.correlate(datas[:, x], datas[:, x+delta_x], mode='full')
            avg_corr += corr[corr.size//2:]
            cnt += 1
        avg_corr /= cnt
        n = frequency.shape[0]
        phi = np.real(fft.fft(avg_corr, n=n, workers=3))
        
        # phi, R = Correlation_2d(datas, geom = geom, axis = axis)
        # phi = np.mean(phi, axis=0)
        
        kc = 2*np.pi*frequency/Uc

        funct = np.abs(phi*np.exp(-1j * kc * r) / psd)
        funct = np.log(funct)
            
    if geom == 'plan':
        if axis == 'streamwise':
            niters, npts, nlines = datas.shape
            nperseg_time = get_nperseg(niters)
            
            r = np.abs(coords[delta_x] - coords[0])

            frequency, psd = signal.welch(datas[:, :, :], fs=1/dt, window='hann', nperseg=nperseg_time//8, scaling='density', axis=0)
            psd = np.mean(np.mean(psd, axis=2), axis=1)

            avg_corr = np.zeros((niters))
            cnt = 0
            for line in range(nlines):
                for x in range(0, npts-delta_x-1, 2):
                    corr = signal.correlate(datas[:, x, line], datas[:, x+delta_x, line], mode='full')
                    avg_corr += corr[corr.size//2:]
                    cnt += 1
            avg_corr /= cnt
            n = frequency.shape[0]
            phi = np.real(fft.fft(avg_corr, n=n, workers=3))
            
            # phi, R = Correlation_2d(datas, geom = geom, axis = axis)
            # phi = np.mean(phi, axis=0)
            
            kc = 2*np.pi*frequency/Uc

            funct = np.abs(phi[:n]*np.exp(-1j * kc * r) / psd)
            funct = np.log(funct)
            
        if axis == 'spanwise':
            niters, nlines, npts = datas.shape
            nperseg_time = get_nperseg(niters)
            
            r = np.abs(coords[delta_x] - coords[0])

            frequency, psd = signal.welch(datas[:, :, :], fs=1/dt, window='hann', nperseg=nperseg_time//8, scaling='density', axis=0)
            psd = np.mean(np.mean(psd, axis=2), axis=1)

            avg_corr = np.zeros((niters))
            cnt = 0
            for line in range(nlines):
                for x in range(0, npts-delta_x-1, 2):
                    corr = signal.correlate(datas[:, line, x], datas[:, line, x+delta_x], mode='full')
                    avg_corr += corr[corr.size//2:]
                    cnt += 1
            avg_corr /= cnt
            n = frequency.shape[0]
            phi = np.real(fft.fft(avg_corr, n=n, workers=3))
            
            # phi, R = Correlation_2d(datas, geom = geom, axis = axis)
            # phi = np.mean(phi, axis=0)
            
            kc = 2*np.pi*frequency/Uc

            funct = np.abs(phi[:n]*np.exp(-1j * kc * r) / psd)
            funct = np.log(funct)
            

    return (frequency, phi, psd, funct, kc, r)


def Space_correation(data1, data2, geom = "plan", mode_corr = 'half', axis = "streamwise"):
    
    if geom == 'line':
            
        if axis == 'streamwise' or axis == 'spanwise' or axis == 'wallnormal':
            niters, npts = data1.shape
            full_corr = np.zeros((2*npts-1))
            R = np.zeros((npts))
            for iters in range(niters):
                full_corr[:] += signal.correlate(data1[iters,:], data2[iters,:], mode='full', method='auto')
            full_corr /= (niters)
            full_corr /= max(full_corr)
            R[:] = full_corr[full_corr.size //2 :]
            
        else:
            print('Error: The argument axis is incorrect. Please choose between \'time\', \'streamwise\', \'spanwise\' and \'wallnormal\'.')
            exit(-1)
        
        if mode_corr == 'full':
            return(full_corr)
        elif mode_corr == 'half':
            return(R)
        else:
            print('Error: The argument mode_corr is incorrect. Please choose between \'full\' and \'half\'.')
            exit(-1)
        
    elif geom == 'plan':
            
        if axis == 'streamwise' or axis == 'wallnormal':
            niters, npts, nlines = data1.shape
            full_corr = np.zeros((2*npts-1))
            R = np.zeros((npts))
            for lines in range(nlines):
                for iters in range(niters):
                    full_corr[:] += signal.correlate(data1[iters,:,lines], data2[iters,:,lines], mode='full', method='auto')
            full_corr /= (niters*nlines)
            full_corr /= max(full_corr)
            R[:] = full_corr[full_corr.size //2 : ]
            
        elif axis == 'spanwise':
            niters, nlines, npts = data1.shape
            full_corr = np.zeros((2*npts-1))
            R = np.zeros((npts))
            for lines in range(nlines):
                for iters in range(niters):
                    full_corr[:] += signal.correlate(data1[iters,lines,:], data2[iters,lines,:], mode='full', method='auto')
            full_corr /= (niters*nlines)
            full_corr /= max(full_corr)
            R[:] = full_corr[full_corr.size //2 :]
            
        else:
            print('Error: The argument axis is incorrect. Please choose between \'time\', \'streamwise\', \'spanwise\' and \'wallnormal\'.')
            exit(-1)
        
        if mode_corr == 'full':
            return(full_corr)
        elif mode_corr == 'half':
            return(R)
        else:
            print('Error: The argument geom is incorrect. Please choose between \'line\' and \'plan\'.')
            exit(-1)
            
    else:
        print('Error: The argument mode_corr is incorrect. Please choose between \'full\' and \'half\'.')
        exit(-1)
        
        
def von_karman_spectra(data1, data2, geom = 'plan', axis = "streamwise"):
    
    if geom == "plan":
        
        if axis == 'streamwise' or axis == 'wallnormal':
            niters, npts, nlines = data1.shape
            split_t = int(2**10)
            num_split_t = niters // split_t
            phi = np.zeros((split_t, npts, nlines))
            
            for n in tqdm(range(1,num_split_t), desc=f'psd', colour= 'GREEN'):
                f = fft.fft(data1[(n-1)*split_t:n*split_t,:,:],axis=1)
                ft = fft.fft(f,axis=0)
                fourier1 = fft.fft(ft,axis=2)
                f = fft.fft(data2[(n-1)*split_t:n*split_t,:,:],axis=1)
                ft = fft.fft(f,axis=0)
                fourier2 = fft.fft(ft,axis=2)
                phi += np.real(fourier1*np.conj(fourier2))
            phi /= (num_split_t-1)
            
        elif axis == "spanwise":
            niters, nlines, npts = data1.shape
            split_t = int(2**10)
            num_split_t = niters // split_t
            phi = np.zeros((split_t, nlines, npts))
            
            for n in tqdm(range(1,num_split_t), desc=f'psd', colour= 'GREEN'):
                f = fft.fft(data1[(n-1)*split_t:n*split_t,:,:],axis=1)
                fourier1 = fft.fft(f,axis=2)
                f = fft.fft(data2[(n-1)*split_t:n*split_t,:,:],axis=1)
                fourier2 = fft.fft(f,axis=2)
                phi += np.real(fourier1*np.conj(fourier2))
            phi /= (num_split_t-1)
            
        else:
            print('Error: The argument axis is incorrect. Please choose between \'streamwise\', \'spanwise\' and \'wallnormal\'.')
            exit(-1)
            
    elif geom == "line":
        
        niters, npts = data1.shape
        split_t = int(2**10)
        num_split_t = niters // split_t
        phi = np.zeros((split_t, npts, nlines))
        
        for n in tqdm(range(1,num_split_t), desc=f'psd', colour= 'GREEN'):
            fourier1 = fft.fft(data1[(n-1)*split_t:n*split_t,:],axis=1)
            fourier2 = fft.fft(data2[(n-1)*split_t:n*split_t,:],axis=1)
            phi += np.real(fourier1*np.conj(fourier2))
        phi /= (num_split_t-1)
        
    else:
        print('Error: The argument geom is incorrect. Please choose between \'plan\' and \'line\'.')
        exit(-1)
        
        
    return(phi)
        