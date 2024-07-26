import numpy as np
import time
from scipy import integrate
import glob
from tqdm import tqdm

from read_fpar import *
from parameters import *
from tools import *
from engine import *
from von_karman import*
from read_rans import*

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
    zp = [5, 20, 40, 60, 80, 98, 151, 199, 251, 302]
    z = []
    for he in zp:
        z.append(he * cflow.nu / cflow.ut)
    
    # path = 'input/' #change it as you want
    path = '/media/julien/Verbatim/julien/channel_wrles_retau395/' #on the external disc
    #### Streamwise files ####
    fpars_files_streamwise_p = sorted(glob.glob(path + "streamwise/*fpar.p"))
    
    fpars_files_streamwise_u1 = sorted(glob.glob(path + "streamwise/*fpar.u1"))
    
    fpars_files_streamwise_u2 = sorted(glob.glob(path + "streamwise/*fpar.u2"))
    
    fpars_files_streamwise_u3 = sorted(glob.glob(path + "streamwise/*fpar.u3"))
    
    #### Spanwise files ####
    fpars_files_spanwise_p = sorted(glob.glob(path + "spanwise/*fpar.p"))
    
    fpars_files_spanwise_u1 = sorted(glob.glob(path + "spanwise/*fpar.u1"))
    
    fpars_files_spanwise_u2 = sorted(glob.glob(path + "spanwise/*fpar.u2"))
    
    fpars_files_spanwise_u3 = sorted(glob.glob(path + "spanwise/*fpar.u3"))
    
    #### Wall normal files ####
    fpars_files_normal = sorted(glob.glob(path + "normal/*fpar*"))
    
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
    # print("========================================")
    print("\nReading input files ...")
    
    _, x1, x2, _, nt, n1, n2, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane(fpars_files_streamwise_u1[0])
    nt = nt - 1
    # t = np.linspace(0,(tEnd-tStart)*dt,nt)
    # X = np.linspace(-np.pi, np.pi, n1)
    
    dx = cflow.xlen / n1
    
    print(f'{YELLOW}Streamwise paraeters{RESET}')
    print('len x1:', len(x1))
    print('len x2:', len(x2))
    print('nt:', nt)
    print('n1:', n1)
    print('n2:', n2)
    print('iprecision:\n', iprecision)
    

    ##### FROZEN TURBULENCE #####
    # col = 1
    # fig1, fig2, fig3, fig4 = init_figures(z, n2) #Figure initialization
    # start_time = time.time()
    
    # for zplan in np.arange(0, n2, n2//3, dtype=int):
        
    #     print("========================================")
    #     print(f'Frozen turbulence validation for {YELLOW}z={z[zplan]:.2f}{RESET}')
    #     print('Plan number:', zplan)
    #     print("Reading input files ...\n")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
    #     nt = nt - 1
    #     dx = cflow.xlen / n1
    #     datas_u1 = var[1:,:,:]
        
    #     Uc = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
    #     print('Uc:', Uc)
        
    #     #### Spectra ####
    #     # omega, k, time_spectra, space_spectra = frozen_turbulence(datas_u1, zplan, z, nt, split_time, dt, n1, dx=dx, ch="spectra")
        
    #     # frozen_turbulence_plot(fig1, col, omega = omega, Uc = Uc, time_spectra = time_spectra, k = k, space_spectra = space_spectra, ch = "spectra")
        
    #     # del time_spectra
    #     # del space_spectra
    #     # del omega
    #     # del k
        
    #     #### Autocorrelation ####
    #     # Dt, Dx, R_time, R_space = frozen_turbulence(datas_u1, zplan, z, nt, split_time, dt, n1, tEnd=tEnd, tStart=tStart, ch="corr")

    #     # frozen_turbulence_plot(fig2, col, Uc=Uc, R_time = R_time, R_space = R_space, Dt = Dt, Dx = Dx, ch = "corr")
        
    #     # del Dt
    #     # del Dx
    #     # del R_time
    #     # del R_space
        
    #     #### Autocorrelation 2D ####
        
    #     # Dt, Dx, R2d, coef = frozen_turbulence(datas_u1, zplan, z, nt, split_time, dt, n1, tEnd = tEnd, tStart = tStart, ch = "corr2d")
        
    #     # frozen_turbulence_plot(fig3, col, Dt = Dt, Dx = Dx, R2d=R2d, coef=coef, ch = "corr2d")
        
    #     # del Dt
    #     # del Dx
    #     # del R2d
        
    #     #### Gamma ####
        
    #     for delta_x in range(1, 20, 5):
    #         kc, omega, funct = frozen_turbulence(datas_u1, zplan, z, nt, split_time, dt, n1, dx = delta_x, x1 = x1, Uc = Uc, ch = "gamma")
            
    #         r = np.abs(x1[delta_x] - x1[0])
            
            
    #         frozen_turbulence_plot(fig4, col, omega = omega, k = kc, delta_x=delta_x, funct=funct, r=r, ch = "w_gamma")
            
    #     del omega
    #     del funct
    #     del kc
        
    #     col +=1
    #     del datas_u1
        
    # elapsed_time = time.time() - start_time
    # minutes, seconds = divmod(elapsed_time, 60)
        
    # # Update layout properties
    # fig1.update_layout(height=600, width=900, title_text="Power spectra streawise", font=font, legend=dict(y=1.2, x=0.9))
    # fig2.update_layout(height=600, width=900, title_text='Autocorrelation comparison Streamwise', legend=dict(y=1.2, x=0.9), font=font)
    # fig3.update_layout(height=600, width=1100, title_text='Correlation 2D Streamwise', legend=dict(y=1.2, x=0.9), font=font)
    # fig4.update_layout(height=600, width=1100, title_text='Gamma determmination', legend=dict(y=1.2, x=1.0), font=font)
    
    # if split_time == 'Y':
    #     # save_figures(fig1, "output/split_time/frozen_turbulence/PowerSpectraComparison.png")
    #     # save_figures(fig2, "output/split_time/frozen_turbulence/CorrelationComparison.png")
    #     # save_figures(fig3, "output/split_time/frozen_turbulence/Correlation2D.png")
    #     save_figures(fig4, "output/split_time/frozen_turbulence/Crosscorrelation_omega.png")
    
    # print(f'\n Frozen turbulence valiation done in : {int(minutes)}m {seconds:.2f}s \n')
    
    #### MEAN VELOCITY PROFILE ####
    # Uc_list = []
    # for zplan, zvalue in enumerate(z):
    #     print("========================================")
    #     print(f'Veocity profile for {YELLOW}z={zvalue:.2f}{RESET}')
    #     print('Plan number:', zplan)
    #     print("Reading input files ...\n")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
        
    #     datas_u1 = var[1:,:,:]
    #     Uc = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
    #     Uc_list.append(Uc)
        
    # print('Uc:',Uc_list)
    # print('z:',z)
    # Uc_list = [0, *Uc_list]
    # z = [0, *z]
    # print('Uc:',Uc_list)
    # print('z:',z)
    # fig = go.Figure()
    # fig.add_trace(go.Scatter(x=Uc_list, y=z, mode= 'lines+markers', line=dict(color='firebrick', width=3), name='$Uc$'))
    # fig.add_trace(go.Scatter(x=np.linspace(0, 1.5, 10), y=np.zeros((10)), mode= 'lines', line=dict(color='black', width=3), name='$z=0$'))
    # fig.update_xaxes(title='$x(m)$')#, range=[-np.pi,np.pi])
    # fig.update_yaxes(title='$z(m)$')#, range=[0,1])
    # fig.update_layout(height=600, width=800, title="Mean velocity profile", font=font, showlegend=True)
    # save_figures(fig, "output/velocity_profile.png")
    
    
    #### VON KARMAN MODEL ####
    print("========================================")
    print(f'Van Karman theory')
    start_time = time.time()
    print(f"Reading input files (u1 velocity) for {YELLOW}z={z[0]:.2f}{RESET} ...\n")
    _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[0])
    nt = nt - 1
    datas_u1 = var[1:,:,:]
    print('datas_u1.shape:',datas_u1.shape)
    
    Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
    print('Ux:', Ux)
    
    if split_time == 'Y':
        split_t = 1000
        num_split_t = nt // split_t
        datas_fluct_u1 = np.zeros((split_t, n1, n2))
        for n in range(1,num_split_t):
            datas_fluct_u1 += np.abs(datas_u1[(n-1)*split_t:n*split_t,:,:] - Ux)
            
        datas_fluct_u1 /= (num_split_t-1)    
        sigma_u_squared = sigma_power2(datas_fluct_u1, axis = "streamwise")
        # sigma_u_squared = sigma(datas_fluct_u1, ko2_omega.shape[0], axis = "streamwise")
        # sigma_u_squared = sigma_u_squared**2
        
    if split_time == 'n':
        datas_fluct_u1 = np.abs(datas_u1 - Ux)
        sigma_u_squared = sigma_power2(datas_fluct_u1, axis = "streamwise")
        # sigma_u_squared = sigma(datas_fluct_u1, ko2_omega.shape[0], axis = "streamwise")
        # sigma_u_squared = sigma_u_squared**2
        
    del datas_u1
    del datas_fluct_u1
    print('sigma_u_squared:', sigma_u_squared)
    
    print("========================================")
    print(f"Reading input files (u2 velocity) for {YELLOW}z={z[0]:.2f}{RESET} ...\n")
    _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[0])
    nt = nt - 1
    datas_u2 = var[1:,:,:]
    print('datas_u1.shape:',datas_u2.shape)
    
    Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
    print('Uy:', Uy)
    
    if split_time == 'Y':
        split_t = 1000
        num_split_t = nt // split_t
        datas_fluct_u2 = np.zeros((split_t, n1, n2))
        for n in range(1,num_split_t):
            datas_fluct_u2 += np.abs(datas_u2[(n-1)*split_t:n*split_t,:,:] - Uy)
            
        datas_fluct_u2 /= (num_split_t-1)    
        sigma_v_squared = sigma_power2(datas_fluct_u2, axis = "streamwise")
        # sigma_v_squared = sigma(datas_fluct_u2, ko2_omega.shape[0], axis = "streamwise")
        # sigma_v_squared = sigma_v_squared**2
            
        
    if split_time == 'n':
        datas_fluct_u2 = np.abs(datas_u2 - Uy)
        sigma_v_squared = sigma_power2(datas_fluct_u2, axis = "streamwise")
        # sigma_v_squared = sigma(datas_fluct_u2, ko2_omega.shape[0], axis = "streamwise")
        # sigma_v_squared = sigma_v_squared**2
        
    del datas_u2
    del datas_fluct_u2
    print('sigma_v_squared:', sigma_v_squared)
    
    print("========================================")
    print(f"Reading input files (u2 velocity) for {YELLOW}z={z[0]:.2f}{RESET} ...\n")
    _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[0])
    nt = nt - 1
    datas_u3 = var[1:,:,:]
    print('datas_u1.shape:',datas_u3.shape)
    
    Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
    print('Uz:', Uz)
    
    if split_time == 'Y':
        split_t = 1000
        num_split_t = nt // split_t
        datas_fluct_u3 = np.zeros((split_t, n1, n2))
        for n in range(1,num_split_t):
            datas_fluct_u3 += np.abs(datas_u3[(n-1)*split_t:n*split_t,:,:] - Uz)
            
        datas_fluct_u3 /= (num_split_t-1)    
        sigma_w_squared = sigma_power2(datas_fluct_u3, axis = "streamwise")
        # sigma_w_squared = sigma(datas_fluct_u3, ko2_omega.shape[0], axis = "streamwise")
        # sigma_w_squared = sigma_w_squared**2
            
        
    if split_time == 'n':
        datas_fluct_u3 = np.abs(datas_u3 - Uz)
        sigma_w_squared = sigma_power2(datas_fluct_u3, axis = "streamwise")
        # sigma_w_squared = sigma(datas_fluct_u3, ko2_omega.shape[0], axis = "streamwise")
        # sigma_w_squared = sigma_w_squared**2
        
    del datas_u3
    del datas_fluct_u3
    print('sigma_w_squared:', sigma_w_squared)
    
    
    print("========================================")
    print(f"Computing sigma1c for {YELLOW}z={z[0]:.2f}{RESET} ...\n")
    
    sigma1c = sigma_1c(sigma_u_squared, 1/3., sigma_v_squared, 1/3., sigma_w_squared, 1/3.)
    
    print("========================================")
    print(f"Computing L for {YELLOW}z={z[0]:.2f}{RESET} ...\n")
    
    Le = L(0.09*0.519, ko2_tke, tens = 'omega', omega=ko2_omega) #list length omega
    
    R11 = R_11(ko2_omega, 1.0, Ux, sigma1c, Le)
    R22 = R_22(ko2_omega, 1.0, Ux, sigma1c, Le) #=R33
    
    print('R11.shape:',R11.shape)
    
    fig_vanK = go.Figure()
    fig_vanK.add_trace(go.Scatter(x=ko2_omega, y=R11, mode= 'lines+markers', line=dict(color='firebrick', width=3), name='$R^{0}_{11}$'))
    fig_vanK.add_trace(go.Scatter(x=ko2_omega, y=R22, mode= 'lines+markers', line=dict(color='midnightblue', width=3), name='$R^{0}_{22}$'))
    fig_vanK.update_xaxes(title='$\omega(s^{-1})$', type="log", exponentformat='power')
    fig_vanK.update_yaxes(title='$R_{ij}(\omega)$', type="log", exponentformat='power')
    fig_vanK.update_layout(height=600, width=800, title=f"Von Karman Correlation for z={z[0]:.2f}", font=font, showlegend=True)
    save_figures(fig_vanK, "output/von_karman/von_karman_correlation.png")
    
        
    
        
        
        
        
############################################################
    
            
        # _,_,_, datas_u2,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
        # del datas_u2
        # _,_,_, datas_u3,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
        # del datas_u3
        # _,_,_, datas_p,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_p[zplan])
        # del datas_p
        
        
    
        
        
        
        
############################################################
        
    
    
    
if __name__ == "__main__":
    main()