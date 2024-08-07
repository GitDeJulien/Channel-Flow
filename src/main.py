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
    fpars_files_normal_p = sorted(glob.glob(path + "normal/*fpar.p"))
    
    fpars_files_normal_u1 = sorted(glob.glob(path + "normal/*fpar.u1"))
    
    fpars_files_normal_u2 = sorted(glob.glob(path + "normal/*fpar.u2"))
    
    fpars_files_normal_u3 = sorted(glob.glob(path + "normal/*fpar.u3"))
    
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
    #=======================================================
            ##### FROZEN TURBULENCE #####
    #=======================================================
    # print("========================================")
    # print("\nReading input files ...")
    
    # _, x1, x2, _, nt, n1, n2, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane(fpars_files_streamwise_u1[0])
    # nt = nt - 1
    # # t = np.linspace(0,(tEnd-tStart)*dt,nt)
    # # X = np.linspace(-np.pi, np.pi, n1)
    
    # print(f'{YELLOW}Streamwise parameters{RESET}')
    # print('len x1:', len(x1))
    # print('len x2:', len(x2))
    # print('nt:', nt)
    # print('n1:', n1)
    # print('n2:', n2)
    # print('iprecision:\n', iprecision)
    

    # col = 1
    # row = 1
    # ch_plot ="all"
    # fig1u1, fig2u1, fig3u1 = init_figures_ft(zp, n2, ch=ch_plot) #Figure initialization
    # fig1u2, fig2u2, fig3u2 = init_figures_ft(zp, n2, ch=ch_plot)
    # fig1u3, fig2u3, fig3u3 = init_figures_ft(zp, n2, ch=ch_plot)
    
    # figU = go.Figure()
    # U_ratio = []
    # X_ratio = []
    
    # #figU1c = make_subplots(rows=1, cols=2, shared_yaxes= True, y_title='$z^+$')
    # figU1c = go.Figure()
    # U1_list = []
    # Uc_list = []
    
    # start_time = time.time()
    
    # # for zplan in np.arange(0, n2, n2//3, dtype=int):
    # for zplan, zpvalue in enumerate(zp):
        
    #     print("========================================")
    #     print(f'Frozen turbulence validation for {YELLOW}zp={zp[zplan]:.2f}{RESET}')
    #     print('Plan number:', zplan)
    #     print("\nReading input files u1 streamwise...")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
    #     nt = nt - 1
    #     dx = cflow.xlen / n1
    #     datas_u1 = var[1:,:,:]
        
    #     U1 = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
    #     print('U1:', U1)
        
    #     #datas_fluct_u1 = np.zeros_like(datas_u1)
    #     datas_u1 = datas_u1 - U1
    #     #del datas_u1
        
    #     #### Spectra ####
    #     omega, k, time_spectra, space_spectra = frozen_turbulence(datas_u1, zplan, z, nt, split_time, dt, n1, dx=dx, ch="spectra")
        
    #     frozen_turbulence_plot(fig1u1, col, row, omega = omega, Uc = U1, time_spectra = time_spectra, k = k, space_spectra = space_spectra, ch = "spectra")
        
    #     del time_spectra
    #     del space_spectra
    #     del omega
    #     del k
        
    #     #### Autocorrelation ####
    #     # ind1, ind2, Dt, Dx, R_time, R_space = frozen_turbulence(datas_u1, zplan, z, nt, split_time, dt, n1, tEnd=tEnd, tStart=tStart, ch="corr", dx=dx, Uc=U1)

    #     # frozen_turbulence_plot(fig2u1, col, row, Uc=U1, R_time = R_time[:ind1], R_space = R_space[:ind2], Dt = Dt[:ind1], Dx = Dx[:ind2], ch = "corr")
        
    #     # del Dt
    #     # del Dx
    #     # del R_time
    #     # del R_space
        
    #     #### Autocorrelation 2D ####
        
    #     # Dt, Dx, R2d, coef = frozen_turbulence(datas_u1, zplan, z, nt, split_time, dt, n1, tEnd = tEnd, tStart = tStart, ch = "corr2d")
        
    #     # frozen_turbulence_plot(fig3u1, col, row, Dt = Dt, Dx = Dx, R2d=R2d, coef=coef, ch = "corr2d")
        
    #     # U_ratio.append((1./coef[0])/U1)
    #     # X_ratio.append(zp[zplan])
        
    #     # U1_list.append(U1)
    #     # Uc_list.append(1./coef[0])
        
    #     # del Dt
    #     # del Dx
    #     # del R2d
        
        
    #     del datas_u1
        
    #     #### u2 ####
    #     print("\nReading input files u2 streamwise ...")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
    #     nt = nt - 1
    #     dx = cflow.xlen / n1
    #     datas_u2 = var[1:,:,:]
    #     Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
    #     #datas_fluct_u2 = np.zeros_like(datas_u2)
    #     datas_u2 = datas_u2 - Uy
    #     # del datas_u2
        
    #     #### Spectra ####
    #     omega, k, time_spectra, space_spectra = frozen_turbulence(datas_u2, zplan, z, nt, split_time, dt, n1, dx=dx, ch="spectra")
        
    #     frozen_turbulence_plot(fig1u2, col, row, omega = omega, Uc = U1, time_spectra = time_spectra, k = k, space_spectra = space_spectra, ch = "spectra")
        
    #     del time_spectra
    #     del space_spectra
    #     del omega
    #     del k
        
    #     #### Autocorrelation ####
    #     # ind1, ind2, Dt, Dx, R_time, R_space = frozen_turbulence(datas_u2, zplan, z, nt, split_time, dt, n1, tEnd=tEnd, tStart=tStart, ch="corr", dx=dx, Uc=U1)

    #     # frozen_turbulence_plot(fig2u2, col, row, Uc=U1, R_time = R_time[:ind1], R_space = R_space[:ind2], Dt = Dt[:ind1], Dx = Dx[:ind2], ch = "corr")
        
    #     # del Dt
    #     # del Dx
    #     # del R_time
    #     # del R_space
        
    #     #### Autocorrelation 2D ####
        
    #     # Dt, Dx, R2d, coef = frozen_turbulence(datas_u2, zplan, z, nt, split_time, dt, n1, tEnd = tEnd, tStart = tStart, ch = "corr2d")
        
    #     # frozen_turbulence_plot(fig3u2, col, row, Dt = Dt, Dx = Dx, R2d=R2d, coef=coef, ch = "corr2d")
        
    #     # del Dt
    #     # del Dx
    #     # del R2d

        
    #     del datas_u2
        
    #     #### u3 ####
    #     print("\nReading input files u3 streamwise ...")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
    #     nt = nt - 1
    #     dx = cflow.xlen / n1
    #     datas_u3 = var[1:,:,:]
    #     Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
    #     #datas_fluct_u3 = np.zeros_like(datas_u3)
    #     datas_u3 = datas_u3 - Uz
    #     #del datas_u3
        
    #     #### Spectra ####
    #     omega, k, time_spectra, space_spectra = frozen_turbulence(datas_u3, zplan, z, nt, split_time, dt, n1, dx=dx, ch="spectra")
        
    #     frozen_turbulence_plot(fig1u3, col, row, omega = omega, Uc = U1, time_spectra = time_spectra, k = k, space_spectra = space_spectra, ch = "spectra")
        
    #     del time_spectra
    #     del space_spectra
    #     del omega
    #     del k
        
    #     #### Autocorrelation ####
    #     # ind1, ind2, Dt, Dx, R_time, R_space = frozen_turbulence(datas_u3, zplan, z, nt, split_time, dt, n1, tEnd=tEnd, tStart=tStart, ch="corr", dx=dx, Uc=U1)

    #     # frozen_turbulence_plot(fig2u3, col, row, Uc=U1, R_time = R_time[:ind1], R_space = R_space[:ind2], Dt = Dt[:ind1], Dx = Dx[:ind2], ch = "corr")
        
    #     # del Dt
    #     # del Dx
    #     # del R_time
    #     # del R_space
        
    #     #### Autocorrelation 2D ####
        
    #     # Dt, Dx, R2d, coef = frozen_turbulence(datas_u3, zplan, z, nt, split_time, dt, n1, tEnd = tEnd, tStart = tStart, ch = "corr2d")
        
    #     # frozen_turbulence_plot(fig3u3, col, row, Dt = Dt, Dx = Dx, R2d=R2d, coef=coef, ch = "corr2d")
        
    #     # del Dt
    #     # del Dx
    #     # del R2d
        
    #     del datas_u3
        
    #     col +=1
    #     if zplan == 4:
    #         row +=1
    #         col = 1
            
        
    # # Update layout properties for 4 plots
    # if ch_plot == "normal":
    #     fig1u1.update_layout(height=600, width=900, title_text="Power spectra streawise u1", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.04, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
    #     fig2u1.update_layout(height=600, width=1100, title_text='Autocorrelation comparison Streamwise u1', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
    #     fig3u1.update_layout(height=600, width=1100, title_text='Correlation 2D Streamwise u1', legend=dict(y=1.2, x=0.9), font=font)
        
    #     fig1u2.update_layout(height=600, width=900, title_text="Power spectra Streawise u2", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.04, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
    #     fig2u2.update_layout(height=600, width=1100, title_text='Autocorrelation comparison Streamwise u2', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
    #     fig3u2.update_layout(height=600, width=1100, title_text='Correlation 2D Streamwise u2', legend=dict(y=1.2, x=0.9), font=font)
        
    #     fig1u3.update_layout(height=600, width=900, title_text="Power spectra Streawise u3", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
    #     fig2u3.update_layout(height=600, width=1100, title_text='Autocorrelation comparison Streamwise u3', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
    #     fig3u3.update_layout(height=600, width=1100, title_text='Correlation 2D Streamwise u3', legend=dict(y=1.2, x=0.9), font=font)
        
    #     figU.add_trace(go.Scatter(x=U_ratio, y=X_ratio, line=dict(color='midnightblue')))
    #     figU.update_xaxes(title='$U_c/U_1$')
    #     figU.update_yaxes(title='$z^+$')
    #     figU.update_layout(title='Velocity ratio', font=font, showlegend=False)
        
    #     figU1c.add_trace(go.Scatter(x=U1_list, y=X_ratio, name='$U_1$', line=dict(color='midnightblue')))
    #     figU1c.add_trace(go.Scatter(x=Uc_list, y=X_ratio, name='$U_c$', line=dict(color='firebrick')))
    #     figU1c.update_xaxes(title='$\\text{mean velocity}~(m.s^{-1})$')
    #     figU1c.update_yaxes(title='$z^+$')
    #     figU1c.update_layout(title='Velocity comparison', font=font, legend=dict(yanchor="bottom", xanchor="right"))
        
    #     if split_time == 'Y':
    #         save_figures(fig1u1, "output/split_time/frozen_turbulence/power_spectra/u1.png")
    #         #save_figures(fig2u1, "output/split_time/frozen_turbulence/correlation_st/u1.png")
    #         #save_figures(fig3u1, "output/split_time/frozen_turbulence/correlation2D/u1.png")
            
    #         save_figures(fig1u2, "output/split_time/frozen_turbulence/power_spectra/u2.png")
    #         #save_figures(fig2u2, "output/split_time/frozen_turbulence/correlation_st/u2.png")
    #         #save_figures(fig3u2, "output/split_time/frozen_turbulence/correlation2D/u2.png")
        
    #         save_figures(fig1u3, "output/split_time/frozen_turbulence/power_spectra/u3.png")
    #         #save_figures(fig2u3, "output/split_time/frozen_turbulence/correlation_st/u3.png")
    #         #save_figures(fig3u3, "output/split_time/frozen_turbulence/correlation2D/u3.png")
            
    #         #save_figures(figU, "output/split_time/frozen_turbulence/correlation2D/u_ratio.png")
            
    #         #save_figures(figU1c, "output/split_time/frozen_turbulence/correlation2D/u_1c.png")
        
    # # # Update layout properties for full plots.
    # if ch_plot == "all":
    #     fig1u1.update_layout(height=900, width=900, title_text="Power spectra streawise u1", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
    #     fig2u1.update_layout(height=600, width=1100, title_text='Autocorrelation comparison Streamwise u1', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
    #     fig3u1.update_layout(height=900, width=1100, title_text='Correlation 2D Streamwise u1', legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40), font=font)
        
    #     fig1u2.update_layout(height=900, width=900, title_text="Power spectra Streawise u2", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
    #     fig2u2.update_layout(height=600, width=1100, title_text='Autocorrelation comparison Streamwise u2', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
    #     fig3u2.update_layout(height=900, width=1100, title_text='Correlation 2D Streamwise u2', legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40), font=font)

    #     fig1u3.update_layout(height=900, width=900, title_text="Power spectra Streawise u3", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
    #     fig2u3.update_layout(height=600, width=1100, title_text='Autocorrelation comparison Streamwise u3', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
    #     fig3u3.update_layout(height=900, width=1100, title_text='Correlation 2D Streamwise u3', legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40), font=font)
        
    #     figU.add_trace(go.Scatter(x=U_ratio, y=X_ratio, line=dict(color='midnightblue')))
    #     figU.update_xaxes(title='$U_c/U_1$')
    #     figU.update_yaxes(title='$z^+$')
    #     figU.update_layout(title='Velocity ratio', font=font, showlegend=False)
        
    #     figU1c.add_trace(go.Scatter(x=U1_list, y=X_ratio, name='$U_1$', line=dict(color='midnightblue')))
    #     figU1c.add_trace(go.Scatter(x=Uc_list, y=X_ratio, name='$U_c$', line=dict(color='firebrick')))
    #     figU1c.update_xaxes(title='$\\text{mean velocity}~(m.s^{-1})$')
    #     figU1c.update_yaxes(title='$z^+$')
    #     figU1c.update_layout(title='Velocity comparison', font=font, legend=dict(yanchor="bottom", xanchor="right"))
        
    #     if split_time == 'Y':
    #         save_figures(fig1u1, "output/split_time/frozen_turbulence/power_spectra/u1_all.png")
    #         #save_figures(fig2u1, "output/split_time/frozen_turbulence/correlation_st/u1_all.png")
    #         #save_figures(fig3u1, "output/split_time/frozen_turbulence/correlation2D/u1_all.png")
            
    #         save_figures(fig1u2, "output/split_time/frozen_turbulence/power_spectra/u2_all.png")
    #         #save_figures(fig2u2, "output/split_time/frozen_turbulence/correlation_st/u2_all.png")
    #         #save_figures(fig3u2, "output/split_time/frozen_turbulence/correlation2D/u2_all.png")
            
    #         save_figures(fig1u3, "output/split_time/frozen_turbulence/power_spectra/u3_all.png")
    #         #save_figures(fig2u3, "output/split_time/frozen_turbulence/correlation_st/u3_all.png")
    #         #save_figures(fig3u3, "output/split_time/frozen_turbulence/correlation2D/u3_all.png")
            
    #         #save_figures(figU, "output/split_time/frozen_turbulence/correlation2D/u_ratio.png")
            
    #         #save_figures(figU1c, "output/split_time/frozen_turbulence/correlation2D/u_1c.png")
            
    # elapsed_time = time.time() - start_time
    # minutes, seconds = divmod(elapsed_time, 60)
    
    # print(f'\n Frozen turbulence valiation done in : {int(minutes)}m {seconds:.2f}s \n')
    
    
    #=======================================================
                    #### GAMMA ####
    #=======================================================
    # print("\n========================================")
    # print(f"{YELLOW}Gamma determination{RESET}")
    # start_time = time.time()
    # print("========================================")
    # print(f"{YELLOW}Streamwise{RESET}")
    # col = 1
    # row = 1
    
    # print("\nReading input files ...")
    # _, x1, x2, _, nt, n1, n2, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane(fpars_files_streamwise_u1[0])
    # nt = nt - 1
    # dx = cflow.xlen / n1
    
    # ch_plot ="normal" 
    # figu1g, figu2g, figu3g = init_figures_gamma(zp, n2, ch=ch_plot) #Figure initialization
    
    # figgamma = make_subplots(rows=1, cols=2, shared_yaxes=True, y_title='$\gamma$', subplot_titles=("First slop", "Second slop"))
    # gamma_u1_1 = []
    # gamma_u2_1 = []
    # gamma_u3_1 = []
    # gamma_u1_2 = []
    # gamma_u2_2 = []
    # gamma_u3_2 = []
    
    # z_axis = []
    
    # for zplan in np.arange(0, n2, n2//3, dtype=int):
    # #for zplan, zpvalue in enumerate(zp):
    
        
    #     print("========================================")
    #     print(f'Gamma determination for {YELLOW}zp={zp[zplan]:.2f}{RESET}')
    #     print('Plan number:', zplan)
        
    #     z_axis.append(zp[zplan])
        
    #     #### u1 ####
    #     print("\nReading input files u1 streamwise...")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
    #     datas_u1 = var[1:,:,:]
        
    #     U1 = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
    
    #     cpt = 0
    #     for delta_x in range(1, 20, 5):
    #         ind1, ind2, kc, omega, funct = frozen_turbulence(datas_u1, zplan, z, nt, split_time, dt, n1, dx = dx, delta_x=delta_x, x1 = x1, Uc = U1, ch = "gamma")
            
    #         r = np.abs(x1[delta_x] - x1[0])
            

    #         slop1, slop2 = frozen_turbulence_plot(figu1g, col, row, omega = omega, k = kc, delta_x=delta_x, funct=funct, r=r, z=zp[zplan], ind1=ind1, ind2=ind2, cpt=cpt, ch = "w_gamma")
    #         cpt += 1
            
    #     gamma_u1_1.append(slop1[0]*U1/r)
    #     gamma_u1_2.append(slop2[0]*U1/r)
            
    #     del omega
    #     del funct
    #     del kc
    #     del slop1
    #     del slop2
        
    #     del datas_u1
        
    #     #### u2 ####
    #     print("\nReading input files u2 streamwise ...")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
    #     datas_u2 = var[1:,:,:]
        
    #     cpt = 0
    #     for delta_x in range(1, 20, 5):
    #         ind1, ind2, kc, omega, funct = frozen_turbulence(datas_u2, zplan, z, nt, split_time, dt, n1, dx = dx, delta_x=delta_x, x1 = x1, Uc = U1, ch = "gamma")
            
    #         r = np.abs(x1[delta_x] - x1[0])
            

    #         slop1, slop2 = frozen_turbulence_plot(figu2g, col, row, omega = omega, k = kc, delta_x=delta_x, funct=funct, r=r, z=zp[zplan], ind1=ind1, ind2=ind2, cpt=cpt, ch = "w_gamma")
    #         cpt += 1
            
    #     gamma_u2_1.append(slop1[0]*U1/r)
    #     gamma_u2_2.append(slop2[0]*U1/r)
            
    #     del omega
    #     del funct
    #     del kc
    #     del slop1
    #     del slop2
        
    #     del datas_u2
        
    #     #### u3 ####
    #     print("\nReading input files u3 streamwise ...")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
    #     datas_u3 = var[1:,:,:]
        
    #     cpt = 0
    #     for delta_x in range(1, 20, 5):
    #         ind1, ind2, kc, omega, funct = frozen_turbulence(datas_u3, zplan, z, nt, split_time, dt, n1, dx = dx, delta_x=delta_x, x1 = x1, Uc = U1, ch = "gamma")
            
    #         r = np.abs(x1[delta_x] - x1[0])
        

    #         slop1, slop2 = frozen_turbulence_plot(figu3g, col, row, omega = omega, k = kc, delta_x=delta_x, funct=funct, r=r, z=zp[zplan], ind1=ind1, ind2=ind2, cpt=cpt, ch = "w_gamma")
    #         cpt += 1
            
    #     gamma_u3_1.append(slop1[0]*U1/r)
    #     gamma_u3_2.append(slop2[0]*U1/r)
            
    #     del omega
    #     del funct
    #     del kc
    #     del slop1
    #     del slop2
        
        
    #     del datas_u3
        
    #     col +=1
    #     if zplan == 4:
    #         row +=1
    #         col = 1
            
            
    # # Update layout properties for 4 plots
    # if ch_plot == "normal":
    #     figu1g.update_layout(height=600, width=1100, title_text='Gamma determination u1', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu2g.update_layout(height=600, width=1100, title_text='Gamma determination u2', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu3g.update_layout(height=600, width=1100, title_text='Gamma determination u3', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_1, name='$\gamma_{u1}$', line=dict(color='midnightblue')), row=1, col=1)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_1, name='$\gamma_{u2}$', line=dict(color='firebrick')), row=1, col=1)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_1, name='$\gamma_{u3}$', line=dict(color='darkgreen')), row=1, col=1)
        
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_2, showlegend=False, line=dict(color='midnightblue')), row=1, col=2)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_2, showlegend=False, line=dict(color='firebrick')), row=1, col=2)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_2, showlegend=False, line=dict(color='darkgreen')), row=1, col=2)
        
    #     figgamma.update_xaxes(title_text='$z^+$', row=1, col=1)
    #     figgamma.update_xaxes(title_text='$z^+$', row=1, col=2)
    #     figgamma.update_layout(height=600, width=800, title_text='Gamma evolution', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     if split_time == 'Y':

    #         save_figures(figu1g, "output/split_time/gamma/gamma_u1.png")
            
    #         save_figures(figu2g, "output/split_time/gamma/gamma_u2.png")
        
    #         save_figures(figu3g, "output/split_time/gamma/gamma_u3.png")
            
    #         save_figures(figgamma, "output/split_time/gamma/gamma_view.png")
        
    # # # Update layout properties for full plots.
    # if ch_plot == "all":

    #     figu1g.update_layout(height=900, width=1100, title_text='Gamma determination u1', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu2g.update_layout(height=900, width=1100, title_text='Gamma determination u2', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu3g.update_layout(height=900, width=1100, title_text='Gamma determination u3', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_1, name='$\gamma_{u1}$', line=dict(color='midnightblue')), row=1, col=1)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_1, name='$\gamma_{u2}$', line=dict(color='firebrick')), row=1, col=1)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_1, name='$\gamma_{u3}$', line=dict(color='darkgreen')), row=1, col=1)
        
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_2, showlegend=False, line=dict(color='midnightblue')), row=1, col=2)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_2, showlegend=False, line=dict(color='firebrick')), row=1, col=2)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_2, showlegend=False, line=dict(color='darkgreen')), row=1, col=2)
        
    #     figgamma.update_xaxes(title_text='$z^+$', row=1, col=1)
    #     figgamma.update_xaxes(title_text='$z^+$', row=1, col=2)
    #     figgamma.update_layout(height=600, width=800, title_text='Gamma evolution', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     if split_time == 'Y':

    #         save_figures(figu1g, "output/split_time/gamma/gamma_u1.png")

    #         save_figures(figu2g, "output/split_time/gamma/gamma_u2_all.png")

    #         save_figures(figu3g, "output/split_time/gamma/gamma_u3_all.png")
            
    #         save_figures(figgamma, "output/split_time/gamma/gamma_view_all.png")
            
    # elapsed_time = time.time() - start_time
    # minutes, seconds = divmod(elapsed_time, 60)
    
    # print(f'\n Gamma determination done in : {int(minutes)}m {seconds:.2f}s \n')
    
    
    #=======================================================
                #### SPACE CORRELATION ####
    #=======================================================
    print("\n========================================")
    print(f"{YELLOW}Space correlation computation{RESET}")
    start_time = time.time()
    print("========================================")
    print(f"{YELLOW}Streamwise{RESET}")
    col = 1
    row = 1
    print("\nReading input files ...")
    
    _, x1, x2, _, nt, n1, n2, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane(fpars_files_streamwise_u1[0])
    nt = nt - 1
    dx = cflow.xlen / n1
    
    ch_plot ="normal"
    fig1sc= init_figures_sc(zp, n2, ch=ch_plot) #Figure initialization
    
    for zplan in np.arange(0, n2, n2//3, dtype=int):
    #for zplan, zpvalue in enumerate(zp):
        
        print("========================================")
        print(f'Space correlation for {YELLOW}zp={zp[zplan]:.2f}{RESET}')
        print('Plan number:', zplan)
        print("Reading input files streamwise...")
        
        split_t = int(2**10)
        num_split_t = nt // split_t
        Dx = np.linspace(0,xlen//2,n1)
        
        ## uu ##
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
        datas_u1 = var[1:,:,:]
        
        Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
        print('Ux:', Ux)
        
        datas_u1 = datas_u1 - Ux
        
        if split_time == 'Y':
            
            # datas_fluct_u1 = np.zeros((split_t, n1, n2))
            # for n in range(1,num_split_t):
            #     datas_fluct_u1 += datas_u1[(n-1)*split_t:n*split_t,:,:] - Ux
            # del datas_u1
            # datas_fluct_u1 /= (num_split_t-1)
            
            R_uu = np.zeros((n1))
            for n in tqdm(range(1,num_split_t), desc=f'PSD (z={z[zplan]:.2f})', colour= 'GREEN'):
                R_uu += Space_correation(datas_u1[(n-1)*split_t:n*split_t,:,:], datas_u1[(n-1)*split_t:n*split_t,:,:], geom = "plan", mode_corr = 'half', axis = "streamwise")
            R_uu /= (num_split_t-1)
            
        if split_time == 'n':
            R_uu = Space_correation(datas_u1, datas_u1, geom = "plan", mode_corr = 'half', axis = "streamwise")
            
        space_correlation_plot(fig1sc, col, row, Dx, R_uu, name = '$R_{UU}$', color='midnightblue', axis='streamwise')
        del R_uu
        del datas_u1
        
        ## vv ##
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
        datas_u2 = var[1:,:,:]
        
        Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
        
        datas_u2 = datas_u2 - Uy
        
        if split_time == 'Y':
            
            # datas_fluct_u2 = np.zeros((split_t, n1, n2))
            # for n in range(1,num_split_t):
            #     datas_fluct_u2 += datas_u2[(n-1)*split_t:n*split_t,:,:] - Uy
            # del datas_u2
            # datas_fluct_u2 /= (num_split_t-1)
            
            R_vv = np.zeros((n1))
            for n in tqdm(range(1,num_split_t), desc=f'PSD (z={z[zplan]:.2f})', colour= 'GREEN'):
                R_vv += Space_correation(datas_u2[(n-1)*split_t:n*split_t,:,:], datas_u2[(n-1)*split_t:n*split_t,:,:], geom = "plan", mode_corr = 'half', axis = "streamwise")
            R_vv /= (num_split_t-1)
        
        if split_time == 'n':
            R_vv = Space_correation(datas_u2, datas_u2, geom = "plan", mode_corr = 'half', axis = "streamwise")
            
        space_correlation_plot(fig1sc, col, row, Dx, R_vv, name = '$R_{VV}$', color='firebrick', axis='streamwise')
        del R_vv
        del datas_u2
        
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
        datas_u3 = var[1:,:,:]
        
        Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
        
        datas_u3 = datas_u3 - Uz
        
        if split_time == 'Y':
            
            # datas_fluct_u3 = np.zeros((split_t, n1, n2))
            # for n in range(1,num_split_t):
            #     datas_fluct_u3 += datas_u3[(n-1)*split_t:n*split_t,:,:] - Uz
            # del datas_u3
            # datas_fluct_u3 /= (num_split_t-1)
            
            R_ww = np.zeros((n1))
            for n in tqdm(range(1,num_split_t), desc=f'PSD (z={z[zplan]:.2f})', colour= 'GREEN'):
                R_ww += Space_correation(datas_u3[(n-1)*split_t:n*split_t,:,:], datas_u3[(n-1)*split_t:n*split_t,:,:], geom = "plan", mode_corr = 'half', axis = "streamwise")
            R_ww /= (num_split_t-1)
        
        if split_time == 'n':
            R_ww = Space_correation(datas_u3, datas_u3, geom = "plan", mode_corr = 'half', axis = "streamwise")
            
        space_correlation_plot(fig1sc, col, row, Dx, R_ww, name = '$R_{WW}$', color='darkgreen', axis='streamwise')
        del R_ww
        del datas_u3
        
        col +=1
        if zplan == 4:
            row +=1
            col = 1
            
        
    if ch_plot == 'normal':
        fig1sc.update_layout(height=600, width=900, title_text='Space Correlation Streamwise', font=font,  legend=dict(yanchor="bottom", y=1.03, xanchor="right", x=1, orientation='h'))
    if ch_plot == 'all':
        fig1sc.update_layout(height=900, width=900, title_text='Space Correlation Streamwise', font=font,  legend=dict(yanchor="bottom", y=1.03, xanchor="right", x=1, orientation='h'))
        
    if split_time == 'Y':
        
        if ch_plot == 'normal':
            save_figures(fig1sc, "output/split_time/space_correlation/streamwise.png")
        if ch_plot == 'all':
            save_figures(fig1sc, "output/split_time/space_correlation/streamwise_all.png")
        
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    print(f'\n Space correlation done in : {int(minutes)}m {seconds:.2f}s \n')
        
    
    
    #=======================================================
                #### RANS AND LES ANALYSIS ####
    #=======================================================
    # print("\n========================================")
    # print(f"{YELLOW}Von Karman theorical computation{RESET}")
    # start_time = time.time()
    # print("========================================")
    # zp_RANS = cflow.ut/cflow.nu * normal
    # fig_RANS = make_subplots(rows=1, cols=2, shared_yaxes= True, y_title='$z^+$')
    # fig_vel_profil = make_subplots(rows=1, cols=3, shared_yaxes= True, y_title='$z^+$')
    # fig_var = make_subplots(rows=1, cols=3, shared_yaxes= True, y_title='$z^+$')
    # fig_vel_ratio_profil = go.Figure()
    # fig_RANS.add_trace(go.Scatter(x=u_velocity, y=zp_RANS, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'), name='$U\\text{(RANS)}$'), row=1, col=1)
    # fig_RANS.add_trace(go.Scatter(x=ko2_tke, y=zp_RANS, mode= 'lines+markers', line=dict(color='darkgreen', width=2), marker=dict(symbol='circle'), name='$k_T\\text{(RANS)}$'), row=1, col=2)
    
    # U = []
    # V = []
    # W = []
    # var1 = []
    # var2 = []
    # var3 = []
    # kt = []
    # for zplan, zpvalue in enumerate(zp):
    #     print('zp:', zpvalue)
        
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
    #     datas_u1 = var[1:,:,:]
    #     del var
    #     Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
    #     varUx = np.mean(np.mean(np.mean((datas_u1[:,:,:]-Ux)**2, axis=-1), axis=-1))
    #     U.append(Ux)
    #     var1.append(varUx)
    #     del datas_u1
        
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
    #     datas_u2 = var[1:,:,:]
    #     del var
    #     Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
    #     varUy = np.mean(np.mean(np.mean((datas_u2[:,:,:]-Uy)**2, axis=-1), axis=-1))
    #     V.append(Uy)
    #     var2.append(varUy)
    #     del datas_u2
        
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
    #     datas_u3 = var[1:,:,:]
    #     del var
    #     Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
    #     varUz = np.mean(np.mean(np.mean((datas_u3[:,:,:]-Uy)**2, axis=-1), axis=-1))
    #     W.append(Uz)
    #     var3.append(varUz)
    #     del datas_u3
        
    #     kt.append(varUx + varUy + varUz)
        
    # kt = np.array(kt)
    # U = np.array(U)
    # V = np.array(V)
    # W = np.array(W)
    # var1 = np.array(var1)
    # var2 = np.array(var2)
    # var3 = np.array(var3)
    # kt /= 2.
    # ratio2 = var2 / var1
    # ratio3 = var3 / var1
    
    # fig_RANS.add_trace(go.Scatter(x=U, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='x'), name='$U\\text{(LES)}$'), row=1, col=1)
    # fig_RANS.add_trace(go.Scatter(x=kt, y=zp, mode= 'lines+markers', line=dict(color='darkgreen', width=2), marker=dict(symbol='x'), name='$k_T\\text{(LES)}$'), row=1, col=2)
    
    # fig_RANS.update_xaxes(title='$\\text{velocity}~(m.s^{-1})$', row=1, col=1)
    # fig_RANS.update_xaxes(title='$\\text{kinetic energy}~(m^2.s^{-2})$', row=1, col=2)
    # fig_RANS.update_layout(height=600, width=800, title="RANS and LES data profiles", font=font, showlegend=True, legend=dict(yanchor='bottom', y=1.03, xanchor='left', x=0.9))
    
    # save_figures(fig_RANS, "output/RANS/RANS_profiles.png")
    
    # fig_vel_profil.add_trace(go.Scatter(x=U, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'), name='$U\\text{(LES)}$'), row=1, col=1)
    # fig_vel_profil.add_trace(go.Scatter(x=V, y=zp, mode= 'lines+markers', line=dict(color='midnightblue', width=2), marker=dict(symbol='diamond'), name='$V\\text{(LES)}$'), row=1, col=2)
    # fig_vel_profil.add_trace(go.Scatter(x=W, y=zp, mode= 'lines+markers', line=dict(color='darkgreen', width=2), marker=dict(symbol='x'), name='$W\\text{(LES)}$'),  row=1, col=3)
    
    # fig_vel_profil.update_xaxes(title='$U~(m.s^{-1})$', row=1, col=1)
    # fig_vel_profil.update_xaxes(title='$V~(m.s^{-1})$', row=1, col=2, exponentformat = 'e')
    # fig_vel_profil.update_xaxes(title='$W~(m.s^{-1})$', row=1, col=3, exponentformat = 'e')
    # fig_vel_profil.update_layout(height=600, width=800, title="Velocity profile (LES data)", font=font, showlegend=True, legend=dict(yanchor='bottom', y=1.03, xanchor='left', x=0.9))
    
    # save_figures(fig_vel_profil, "output/RANS/mean_velocity_profiles.png")
    
    # fig_var.add_trace(go.Scatter(x=var1, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'), name='$u_1$'), row=1, col=1)
    # fig_var.add_trace(go.Scatter(x=var2, y=zp, mode= 'lines+markers', line=dict(color='midnightblue', width=2), marker=dict(symbol='diamond'), name='$u_2$'), row=1, col=2)
    # fig_var.add_trace(go.Scatter(x=var3, y=zp, mode= 'lines+markers', line=dict(color='darkgreen', width=2), marker=dict(symbol='x'), name='$u_3$'),  row=1, col=3)
    
    # fig_var.update_xaxes(title='$u_1~(m.s^{-1})$', row=1, col=1)
    # fig_var.update_xaxes(title='$u_2~(m.s^{-1})$', row=1, col=2, exponentformat = 'e')
    # fig_var.update_xaxes(title='$u_3~(m.s^{-1})$', row=1, col=3, exponentformat = 'e')
    # fig_var.update_layout(height=600, width=800, title="Variance velocity fluctuation profile", font=font, showlegend=True, legend=dict(yanchor='bottom', y=1.03, xanchor='left', x=0.9))
    
    # save_figures(fig_var, "output/RANS/var_velocity_profiles.png")
    
    # fig_vel_ratio_profil.add_trace(go.Scatter(x=ratio2, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'), name='$u_2^2 / u_1^2$'))
    # fig_vel_ratio_profil.add_trace(go.Scatter(x=ratio3, y=zp, mode= 'lines+markers', line=dict(color='midnightblue', width=2), marker=dict(symbol='diamond'), name='$u_3^2 / u_1^2$'))
    
    # fig_vel_ratio_profil.update_yaxes(title='$z^+$')
    # fig_vel_ratio_profil.update_xaxes(title='velocity ratio')
    # fig_vel_ratio_profil.update_layout(height=600, width=800, title="Variance velocity ratio profile (LES data)", font=font, showlegend=True, legend=dict(yanchor='bottom', y=1.03, xanchor='left', x=0.9))
    
    # save_figures(fig_vel_ratio_profil, "output/RANS/velocity_ratio_profiles.png")
    
    # elapsed_time = time.time() - start_time
    # minutes, seconds = divmod(elapsed_time, 60)
    
    # print(f'\n RANS study done in : {int(minutes)}m {seconds:.2f}s \n')
    
    
    #=======================================================
                #### VON KARMAN MODEL ####
    #=======================================================
    # print("\n========================================")
    # print(f"{YELLOW}Von Karman theorical computation{RESET}")
    # start_time = time.time()
    # print("========================================")
    
    # col = 1
    # row = 1
    # print("\nReading input files ...")
    
    # _, x1, x2, _, nt, n1, n2, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane(fpars_files_streamwise_u1[0])
    # nt = nt - 1
    # dx = cflow.xlen / n1
    # dy = cflow.ylen / n2
    
    # ch_plot ="normal"
    # fig_vanK = init_figures_vk(zp, n2, ch=ch_plot)
    
    # for zplan in np.arange(0, n2, n2//3, dtype=int):
    # #for zplan, zpvalue in enumerate(zp):
        
    #     # print("========================================")
    #     # print(f'Von Karman theory for {YELLOW}zp={zp[zplan]:.2f}{RESET}')
    #     # print('Plan number:', zplan)
    
    #     # print(f"Reading input files (u1 velocity) for {YELLOW}z={zp[zplan]:.2f}{RESET} ...\n")
    #     # _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
    #     # nt = nt - 1
    #     # datas_u1 = var[1:,:,:]
    #     # print('datas_u1.shape:',datas_u1.shape)
        
    #     # Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
    #     # print('Ux:', Ux)
        
    #     split_t = int(2**10)
    #     num_split_t = nt // split_t
        
    #     # if split_time == 'Y':
    #     #     datas_fluct_u1 = np.zeros((split_t, n1, n2))
    #     #     for n in range(1,num_split_t):
    #     #         datas_fluct_u1 += np.abs(datas_u1[(n-1)*split_t:n*split_t,:,:] - Ux)
                
    #     #     datas_fluct_u1 /= (num_split_t-1)    
    #     #     sigma_u_squared = sigma_power2(datas_fluct_u1, axis = "streamwise")
    #     #     # sigma_u_squared = sigma(datas_fluct_u1, ko2_omega.shape[0], axis = "streamwise")
    #     #     # sigma_u_squared = sigma_u_squared**2
            
    #     # if split_time == 'n':
    #     #     datas_fluct_u1 = np.abs(datas_u1 - Ux)
    #     #     sigma_u_squared = sigma_power2(datas_fluct_u1, axis = "streamwise")
    #     #     # sigma_u_squared = sigma(datas_fluct_u1, ko2_omega.shape[0], axis = "streamwise")
    #     #     # sigma_u_squared = sigma_u_squared**2
            
    #     # del datas_u1
    #     # del datas_fluct_u1
    #     # print('sigma_u_squared:', sigma_u_squared)
        
    #     # print("========================================")
    #     # print(f"Reading input files (u2 velocity) for {YELLOW}z={zp[zplan]:.2f}{RESET} ...\n")
    #     # _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
    #     # nt = nt - 1
    #     # datas_u2 = var[1:,:,:]
    #     # print('datas_u1.shape:',datas_u2.shape)
        
    #     # Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
    #     # print('Uy:', Uy)
        
    #     # if split_time == 'Y':

    #     #     datas_fluct_u2 = np.zeros((split_t, n1, n2))
    #     #     for n in range(1,num_split_t):
    #     #         datas_fluct_u2 += np.abs(datas_u2[(n-1)*split_t:n*split_t,:,:] - Uy)
                
    #     #     datas_fluct_u2 /= (num_split_t-1)    
    #     #     sigma_v_squared = sigma_power2(datas_fluct_u2, axis = "streamwise")
    #     #     # sigma_v_squared = sigma(datas_fluct_u2, ko2_omega.shape[0], axis = "streamwise")
    #     #     # sigma_v_squared = sigma_v_squared**2
                
            
    #     # if split_time == 'n':
    #     #     datas_fluct_u2 = np.abs(datas_u2 - Uy)
    #     #     sigma_v_squared = sigma_power2(datas_fluct_u2, axis = "streamwise")
    #     #     # sigma_v_squared = sigma(datas_fluct_u2, ko2_omega.shape[0], axis = "streamwise")
    #     #     # sigma_v_squared = sigma_v_squared**2
            
    #     # del datas_u2
    #     # del datas_fluct_u2
    #     # print('sigma_v_squared:', sigma_v_squared)
        
    #     # print("========================================")
    #     # print(f"Reading input files (u2 velocity) for {YELLOW}z={zp[zplan]:.2f}{RESET} ...\n")
    #     # _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
    #     # nt = nt - 1
    #     # datas_u3 = var[1:,:,:]
    #     # print('datas_u1.shape:',datas_u3.shape)
        
    #     # Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
    #     # print('Uz:', Uz)
        
    #     # if split_time == 'Y':
    #     #     datas_fluct_u3 = np.zeros((split_t, n1, n2))
    #     #     for n in range(1,num_split_t):
    #     #         datas_fluct_u3 += np.abs(datas_u3[(n-1)*split_t:n*split_t,:,:] - Uz)
                
    #     #     datas_fluct_u3 /= (num_split_t-1)    
    #     #     sigma_w_squared = sigma_power2(datas_fluct_u3, axis = "streamwise")
    #     #     # sigma_w_squared = sigma(datas_fluct_u3, ko2_omega.shape[0], axis = "streamwise")
    #     #     # sigma_w_squared = sigma_w_squared**2
                
            
    #     # if split_time == 'n':
    #     #     datas_fluct_u3 = np.abs(datas_u3 - Uz)
    #     #     sigma_w_squared = sigma_power2(datas_fluct_u3, axis = "streamwise")
    #     #     # sigma_w_squared = sigma(datas_fluct_u3, ko2_omega.shape[0], axis = "streamwise")
    #     #     # sigma_w_squared = sigma_w_squared**2
            
    #     # del datas_u3
    #     # del datas_fluct_u3
    #     # print('sigma_w_squared:', sigma_w_squared)
        
        
    #     # print("========================================")
    #     # print(f"Computing sigma1c for {YELLOW}z={zp[zplan]:.2f}{RESET} ...\n")
        
    #     # sigma1c = sigma_1c(sigma_u_squared, 1/3., sigma_v_squared, 1/3., sigma_w_squared, 1/3.)
        
    #     # print("========================================")
    #     # print(f"Computing L for {YELLOW}z={zp[zplan]:.2f}{RESET} ...\n")
        
    #     # Le = L(0.09*0.519, ko2_tke, tens = 'omega', omega=ko2_omega) #list length omega
        
    #     # kc, phi11 = phi_11(ko2_omega, 1.0, Ux, sigma1c, Le)
    #     # kc, phi22 = phi_22(ko2_omega, 1.0, Ux, sigma1c, Le) #=R33
        
    #     # von_karman_plot(fig_vanK, col, row, kc[:-2], kc[:-2]*phi11[:-2], name = '$\phi^{0}_{11}$', color = 'firebrick')
    #     # von_karman_plot(fig_vanK, col, row, kc[:-2], kc[:-2]*phi22[:-2], name = '$\phi^{0}_{22}$', color = 'midnightblue')
        
    #     # del kc
    #     # del phi11
    #     # del phi22
        
        
    #     ### experience ###
    #     ## u1 ##
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
    #     datas_u1 = var[1:,:,:]
    #     Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
        
    #     phi11_exp__ = von_karman_spectra(datas_u1, datas_u1, geom='plan', axis='streamwise')
    #     del datas_u1
    #     phi11_exp = np.mean(np.mean(phi11_exp__, axis=-1), axis=-1)
    #     # freq1 = 2*np.pi*fft.fftfreq(n1, d=dx)
    #     # freq2 = 2*np.pi*fft.fftfreq(n2, d=dy)
        
    #     # phi11_exp_ = np.zeros((phi11_exp__.shape[0], phi11_exp__.shape[2]))
    #     # for w in range (phi11_exp__.shape[0]):
    #     #     for line in range (phi11_exp__.shape[2]):
    #     #         phi11_exp_[w,line] = integrate.simpson(phi11_exp__[w,:,line], freq1)
                
    #     del phi11_exp__
                
    #     # phi11_exp = np.zeros((phi11_exp_.shape[0]))
    #     # for w in range (phi11_exp_.shape[0]):
    #     #     phi11_exp[w] = integrate.simpson(phi11_exp_[w,:], freq2)
            
    #     # del phi11_exp_
        
    #     freq = fft.fftfreq(2*split_t, d=dt)
    #     freq = freq[:freq.shape[0]//2]
    #     om = 2*np.pi*freq
    #     k1 = Ux*om
        
    #     von_karman_plot(fig_vanK, col, row, k1, k1*phi11_exp, name = '$\phi_{11}exp$', color = 'firebrick', symbols='x')
        
    #     del phi11_exp
    #     del om
    #     del freq
        
    #     ## u2 ##
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
    #     nt = nt - 1
    #     datas_u2 = var[1:,:,:]
        
    #     phi22_exp__ = von_karman_spectra(datas_u2, datas_u2, geom='plan', axis='streamwise')
    #     del datas_u2
    #     phi22_exp = np.mean(np.mean(phi22_exp__, axis=-1), axis=-1)
        
    #     # phi22_exp_ = np.zeros((phi22_exp__.shape[0], phi22_exp__.shape[2]))
    #     # for w in range (phi22_exp__.shape[0]):
    #     #     for line in range (phi22_exp__.shape[2]):
    #     #         phi22_exp_[w,line] = integrate.simpson(phi22_exp__[w,:,line], freq1)
                
    #     del phi22_exp__
                
    #     # phi22_exp = np.zeros((phi22_exp_.shape[0]))
    #     # for w in range (phi22_exp_.shape[0]):
    #     #     phi22_exp[w] = integrate.simpson(phi22_exp_[w,:], freq2)

    #     # del phi22_exp_
        
    #     von_karman_plot(fig_vanK, col, row, k1, k1*phi22_exp, name = '$\phi_{22}exp$', color = 'midnightblue', symbols='x')
            
    #     del phi22_exp
    #     del k1
    #     # del freq2
    #     # del freq1
        
        
    #     col +=1
    #     if zplan == 4:
    #         row +=1
    #         col = 1
            
            
    # if ch_plot == 'normal':
    #     fig_vanK.update_layout(height=600, width=900, title=f"Von Karman Spectra", font=font, showlegend=True, legend=dict(yanchor="bottom", y=1.03, xanchor="right", x=1.03))
    #     save_figures(fig_vanK, "output/von_karman/von_karman_spectra_exp.png")
    # if ch_plot == 'all':
    #     fig_vanK.update_layout(height=900, width=900, title=f"Von Karman Spectra", font=font, showlegend=True, legend=dict(yanchor="bottom", y=1.03, xanchor="right", x=1.03))
    #     save_figures(fig_vanK, "output/von_karman/von_karman_spectra_all.png")
        
    # elapsed_time = time.time() - start_time
    # minutes, seconds = divmod(elapsed_time, 60)
    
    # print(f'\n Von Karman study done in : {int(minutes)}m {seconds:.2f}s \n')
        
        
        
        
############################################################
        
    
    
    
if __name__ == "__main__":
    main()