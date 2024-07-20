import numpy as np
import time
from scipy import integrate
import glob
from tqdm import tqdm

from read_fpar import *
from parameters import *
from tools import *
from engine import *

# from statistics import *
# from van-karman import *

RED = '\033[91m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RESET = '\033[0m'

class ChannelFlow(object):
    
    def __init__(self, xlen, ylen, zlen, rho, uinf, re, ret):
        self.xlen = xlen #length of the channel
        self.ylen = ylen #width of the channel
        self.zlen = zlen #height of the channel
        self.rho  = rho #volumic mass
        self.uinf = uinf #inlet velocity
        self.re   = re #Reynolds number
        self.ret  = ret #friction Reynolds nummber
        
        self.h    = self.zlen/2 #half height of the channel
        self.nu   = self.uinf * self.zlen / self.re #kinetic viscoity 
        self.mu   = self.nu * self.rho #dynamic viscossity
        self.ut   = self.ret * self.nu / self.h #friction velocity
        
        
    def get_geometry(self):
        return(self.xlen, self.ylen, self.zlen, self.h)
    
    def get_fluid_proprety(self):
        return(self.rho, self.re, self.ret, self.nu, self.mu, self.ut, self.uinf)
    
    
def main():
    
    
    cflow = ChannelFlow(xlen=xlen, ylen=ylen, zlen=zlen, rho=rho, uinf=uinf, re=re, ret=ret)
    
    #Reading input files.
    zp = [5, 20, 40, 60, 80, 98, 151, 199, 251, 302, 392]
    z = []
    for he in zp:
        z.append(he * cflow.nu / cflow.ut)
    
    
    #### Streamwise files ####
    fpars_files_streamwise_p = sorted(glob.glob("/media/julien/Verbatim/julien/channel_wrles_retau395/streamwise/*fpar.p"))
    
    fpars_files_streamwise_u1 = sorted(glob.glob("/media/julien/Verbatim/julien/channel_wrles_retau395/streamwise/*fpar.u1"))
    
    fpars_files_streamwise_u2 = sorted(glob.glob("/media/julien/Verbatim/julien/channel_wrles_retau395/streamwise/*fpar.u2"))
    
    fpars_files_streamwise_u3 = sorted(glob.glob("/media/julien/Verbatim/julien/channel_wrles_retau395/streamwise/*fpar.u3"))
    
    #### Spanwise files ####
    fpars_files_spanwise_p = sorted(glob.glob("/media/julien/Verbatim/julien/channel_wrles_retau395/spanwise/*fpar.p"))
    
    fpars_files_spanwise_u1 = sorted(glob.glob("/media/julien/Verbatim/julien/channel_wrles_retau395/spanwise/*fpar.u1"))
    
    fpars_files_spanwise_u2 = sorted(glob.glob("/media/julien/Verbatim/julien/channel_wrles_retau395/spanwise/*fpar.u2"))
    
    fpars_files_spanwise_u3 = sorted(glob.glob("/media/julien/Verbatim/julien/channel_wrles_retau395/spanwise/*fpar.u3"))
    
    #### Wall normal files ####
    fpars_files_normal = sorted(glob.glob("/media/julien/Verbatim/julien/channel_wrles_retau395/normal/*fpar*"))
    
    if len(fpars_files_streamwise_u1) == 0:
        print("Error: Input fiels didn't find.")
        exit(0)
    
############################################################
    
    # Ask for splitting time
    print("Do you want to split the time serie. This action allow to use less RAM but the calculus will be on a smaller time series")
    split_time = input('(Y/n) : ')
    
    if split_time == 'Y':
        print("Time serie splited")
        
    elif split_time == 'n':
        print("Complete time serie kept")
    
    else:
        print("Please chose between 'Y' (Yes) or 'n' (no)")
        split_time = input('(Y/n) : ')
        
        if split_time == 'Y':
            print("Time serie splited")
        
        elif split_time == 'n':
            print("Complete time serie")
            
        else: 
            exit(0)

############################################################
    
    # Streamwise
    print("========================================")
    print("Reading input files ...")
    
    _, x1, x2, _, nt, n1, n2, _, tEnd, version, iprecision, shift, quaternion = read_fpar_extract_plane(fpars_files_streamwise_u1[0])
    
    nt = nt - 1
    t = np.linspace(0,(tEnd-tStart)*dt,nt)
    X = np.linspace(-np.pi, np.pi, n1)
    
    dx = cflow.xlen / n1
    dy = cflow.ylen / n2
    
    print('{YELLOW}Streamwise paraeters{RESET}')
    print('len x1:', len(x1))
    print('len x2:', len(x2))
    print('nt:', nt)
    print('n1:', n1)
    print('n2:', n2)
    print('iprecision:', iprecision)
    print('shift:', shift)
    print('quaternion', quaternion)
    

    ##### FROZEN TURBULENCE #####
    col = 1
    fig1, fig2 = init_figures(z, n2) #Figure initialization
    start_time = time.time()
    for zplan in np.arange(0, n2, n2//3, dtype=int):
        
        print("========================================")
        print(f'Frozen turbulence validation for {YELLOW}z={z[zplan]:.2f}{RESET}')
        
        print("Reading input files ...\n")
        _,_,_,var,_,_,_,_,_,_,_,_, _ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
        datas_u1 = var[1:,:,:]
        Uc = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
        print('Uc:', Uc)
        
        #### Spectra ####
        omega, k, time_spectra, space_spectra = frozen_turbulence(datas_u1, zplan, z, nt, split_time, dt, n1, dx=dx, ch="spectra")
        
        frozen_turbulence_plot(fig1, col, omega = omega, Uc = Uc, time_spectra = time_spectra, k = k, space_spectra = space_spectra, ch = "spectra")
        
        del time_spectra
        del space_spectra
        del omega
        del k
        
        #### Autocorrelation ####
        Dt, Dx, R_time, R_space = frozen_turbulence(datas_u1, zplan, z, nt, split_time, dt, n1, tEnd=tEnd, tStart=tStart, ch="corr")

        frozen_turbulence_plot(fig2, col, Uc=Uc, R_time = R_time, R_space = R_space, Dt = Dt, Dx = Dx, ch = "corr")
        
        del Dt
        del Dx
        del R_time
        del R_space
        
        col +=1
        del datas_u1
        
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
        
    # Update layout properties
    fig1.update_layout(height=600, width=900, title_text="Power spectra streawise", font=font, legend=dict(y=1.2, x=0.9))
    fig2.update_layout(height=600, width=900,title_text='Autocorrelation comparison Streamwise', legend=dict(y=1.2, x=0.9), font=font)
    
    if split_time == 'Y':
        save_figures(fig1, "output/split_time/frozen_turbulence/PowerSpectraComparison.png")
        save_figures(fig2, "output/split_time/frozen_turbulence/CorrelationComparison.png")
    
    print(f'\n Frozen turbulence valiation done in : {int(minutes)}m {seconds:.2f}s \n')
            
        # _,_,_, datas_u2,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
        # del datas_u2
        # _,_,_, datas_u3,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
        # del datas_u3
        # _,_,_, datas_p,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_p[zplan])
        # del datas_p
        
        
    
        
        
        
        
############################################################

    # # Spanwise
    # print("========================================")
    # print("Reading input files ...")
    # _, x1, x2, _, nt, n1, n2, _, tEnd, version, iprecision, shift, quaternion = read_fpar_extract_plane(fpars_files_spanwise_u1[0])

    # Y = np.linspace(-np.pi//2, np.pi//2, n2)
    
    # dx = cflow.xlen / n1
    # dy = cflow.ylen / n2
    
    # print('{YELLOW}Spanwise paraeters{RESET}')
    # print('len x1:', len(x1))
    # print('len x2:', len(x2))
    # print("z:", z)
    # print('nt:', nt)
    # print('n1:', n1)
    # print('n2:', n2)
    # print('iprecision:', iprecision)
    # print('shift:', shift)
    # print('quaternion', quaternion)
    # print("========================================")
    
    # for zplan in range(10):
        
    #     print(f'Spanwise computation for {YELLOW}z={z[zplan]}{RESET}')
        
    #     _,_,_,datas_u1,_,_,_,_,_,_,_,_, _ = read_fpar_extract_plane(fpars_files_spanwise_u1[zplan])
        
    #     if zplan == 0:
    #         print('datas.shape:', datas_u1.shape)
            
    #     _,_,_, datas_u2,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_spanwise_u2[zplan])
    #     _,_,_, datas_u3,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_spanwise_u3[zplan])
    #     _,_,_, datas_p,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_spanwise_p[zplan])
        
        
    # del datas_u1
    # del datas_u2
    # del datas_u3
    # del datas_p
        
    
    
    
if __name__ == "__main__":
    main()