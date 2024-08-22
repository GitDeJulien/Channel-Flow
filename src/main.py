import numpy as np
import time
from scipy import integrate
import glob
from tqdm import tqdm

from read_fpars import *
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
    
    start_time_all = time.time()
    
    cflow = ChannelFlow(xlen=xlen, ylen=ylen, zlen=zlen, rho=rho, uinf=uinf, re=re, ret=ret)
    
    chplot = 'normal'
    
    ## Choose number of plans
    if chplot == "all":
        zp_ind = np.arange(0, len(zplus), 1, dtype=int)
        zp = zplus
        print('zp:',zp)
        print("mode 'all' chosen\n")
    elif chplot == "normal":
        zp_ind = np.arange(0, len(zplus), len(zplus)//3, dtype=int)
        zp = zplus[::3]
        print('zp:',zp)
        print('zp_ind:', zp_ind)
        print("mode 'normal' chosen\n")
    else:
        print("Error: 'chplot' in parameters.py have to be 'normal' or 'all' ")
        
    z = []
    for he in zp:
        z.append(he * cflow.nu / cflow.ut)
        
    #### Streamwise files ####
    #fpars_files_streamwise_p = sorted(glob.glob(path + "streamwise/*fpar.p"))
    
    fpars_files_streamwise_u1 = sorted(glob.glob(in_path + "streamwise/*fpar.u1"))
    
    fpars_files_streamwise_u2 = sorted(glob.glob(in_path + "streamwise/*fpar.u2"))
    
    fpars_files_streamwise_u3 = sorted(glob.glob(in_path + "streamwise/*fpar.u3"))
    
    #### Spanwise files ####
    #fpars_files_spanwise_p = sorted(glob.glob(path + "spanwise/*fpar.p"))
    
    fpars_files_spanwise_u1 = sorted(glob.glob(in_path + "spanwise/*fpar.u1"))
    
    fpars_files_spanwise_u2 = sorted(glob.glob(in_path + "spanwise/*fpar.u2"))
    
    fpars_files_spanwise_u3 = sorted(glob.glob(in_path + "spanwise/*fpar.u3"))
    
    #### Wall normal files ####
    
    fpars_files_normal_u1 = sorted(glob.glob(in_path + "normal_line/*fpar.u1"))
    
    fpars_files_normal_u2 = sorted(glob.glob(in_path + "normal_line/*fpar.u2"))
    
    fpars_files_normal_u3 = sorted(glob.glob(in_path + "normal_line/*fpar.u3"))
    
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
    print("========================================")
    print("\nReading input files ...")
    
    _, x1, x2, _, nt, n1, n2, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane(fpars_files_streamwise_u1[0])
    nt = nt - 1
    # t = np.linspace(0,(tEnd-tStart)*dt,nt)
    # X = np.linspace(-np.pi, np.pi, n1)
    
    print(f'{YELLOW}Streamwise parameters{RESET}')
    print('len x1:', len(x1))
    print('len x2:', len(x2))
    print('nt:', nt)
    print('n1:', n1)
    print('n2:', n2)
    print('iprecision:\n', iprecision)
    

    col = 1
    row = 1
    fig1u1, fig2u1, fig3u1 = init_figures_ft(zp, ch=chplot) #Figure initialization
    fig1u2, fig2u2, fig3u2 = init_figures_ft(zp, ch=chplot)
    fig1u3, fig2u3, fig3u3 = init_figures_ft(zp, ch=chplot)
    
    figU = go.Figure()
    U_ratio = []
    X_ratio = []
    
    #figU1c = make_subplots(rows=1, cols=2, shared_yaxes= True, y_title='$z^+$')
    figU1c = go.Figure()
    U1_list = []
    Uc_list = []
    cpt = 0
    start_time = time.time()
    

    for ind, zplan in enumerate(zp_ind):
        
        print("========================================")
        print(f'Frozen turbulence validation for {YELLOW}zp={zp[ind]:.2f}{RESET}')
        print('Plan number:', zplan)
        print("\nReading input files u1 streamwise...")
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
        nt = nt - 1
        dx = cflow.xlen / n1
        datas_u1 = var[1:,:,:]
        
        U1 = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
        print('U1:', U1)
        
        datas_u1 = datas_u1 - U1
        
        #### Autocorrelation ####
        ind1, ind2, Dt, Dx, R_time, R_space = frozen_turbulence(datas_u1, ind, zp, nt, split_time, dt, n1, tEnd=tEnd, tStart=tStart, ch="corr", dx=dx, Uc=U1)

        frozen_turbulence_plot(fig2u1, col, row, Uc=U1, R_time = R_time[:ind1], R_space = R_space[:ind2], Dt = Dt[:ind1], Dx = Dx[:ind2], ch = "corr")
        
        del Dt
        del Dx    #for zplan in np.arange(0, n2, n2//3, dtype=int):
        del R_time
        del R_space
        
        #### Autocorrelation 2D ####
        
        Dt, Dx, R2d, coef = frozen_turbulence(datas_u1, ind, zp, nt, split_time, dt, n1, tEnd = tEnd, tStart = tStart, ch = "corr2d")
        
        frozen_turbulence_plot(fig3u1, col, row, Dt = Dt, Dx = Dx, R2d=R2d, coef=coef, ch = "corr2d")
        
        U_ratio.append((1./coef[0])/U1)
        X_ratio.append(zp[ind])
        
        U1_list.append(U1)
        Uc_list.append(1./coef[0])
        
        del Dt
        del Dx
        del R2d
    
        #### Spectra ####
        omega, k, time_spectra, space_spectra = frozen_turbulence(datas_u1, ind, zp, nt, split_time, dt, n1, dx=dx, ch="spectra")
        
        frozen_turbulence_plot(fig1u1, col, row, omega = omega, Uc = U1, time_spectra = time_spectra, k = k, space_spectra = space_spectra, ch = "spectra")
        
        del time_spectra
        del space_spectra
        del omega
        del k
        
        
        del datas_u1
        
        #### u2 ####
        print("\nReading input files u2 streamwise ...")
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
        nt = nt - 1
        dx = cflow.xlen / n1
        datas_u2 = var[1:,:,:]
        Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))

        datas_u2 = datas_u2 - Uy
    
        
        #### Autocorrelation ####
        ind1, ind2, Dt, Dx, R_time, R_space = frozen_turbulence(datas_u2, ind, zp, nt, split_time, dt, n1, tEnd=tEnd, tStart=tStart, ch="corr", dx=dx, Uc=U1)

        frozen_turbulence_plot(fig2u2, col, row, Uc=U1, R_time = R_time[:ind1], R_space = R_space[:ind2], Dt = Dt[:ind1], Dx = Dx[:ind2], ch = "corr")
        
        del Dt
        del Dx
        del R_time
        del R_space
        
        #### Autocorrelation 2D ####
        
        Dt, Dx, R2d, coef = frozen_turbulence(datas_u2, ind, zp, nt, split_time, dt, n1, tEnd = tEnd, tStart = tStart, ch = "corr2d")
        
        frozen_turbulence_plot(fig3u2, col, row, Dt = Dt, Dx = Dx, R2d=R2d, coef=coef, ch = "corr2d")
        
        del Dt
        del Dx
        del R2d
    
        #### Spectra ####
        omega, k, time_spectra, space_spectra = frozen_turbulence(datas_u2, ind, zp, nt, split_time, dt, n1, dx=dx, ch="spectra")
        
        frozen_turbulence_plot(fig1u2, col, row, omega = omega, Uc = U1, time_spectra = time_spectra, k = k, space_spectra = space_spectra, ch = "spectra")
        
        del time_spectra
        del space_spectra
        del omega
        del k
        
        del datas_u2
        
        #### u3 ####
        print("\nReading input files u3 streamwise ...")
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
        nt = nt - 1
        dx = cflow.xlen / n1
        datas_u3 = var[1:,:,:]
        Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))

        datas_u3 = datas_u3 - Uz
        
        #### Autocorrelation ####
        ind1, ind2, Dt, Dx, R_time, R_space = frozen_turbulence(datas_u3, ind, zp, nt, split_time, dt, n1, tEnd=tEnd, tStart=tStart, ch="corr", dx=dx, Uc=U1)

        frozen_turbulence_plot(fig2u3, col, row, Uc=U1, R_time = R_time[:ind1], R_space = R_space[:ind2], Dt = Dt[:ind1], Dx = Dx[:ind2], ch = "corr")
        
        del Dt
        del Dx
        del R_time
        del R_space
        
        #### Autocorrelation 2D ####
        
        Dt, Dx, R2d, coef = frozen_turbulence(datas_u3, ind, zp, nt, split_time, dt, n1, tEnd = tEnd, tStart = tStart, ch = "corr2d")
        
        frozen_turbulence_plot(fig3u3, col, row, Dt = Dt, Dx = Dx, R2d=R2d, coef=coef, ch = "corr2d")
        
        del Dt
        del Dx
        del R2d
    
        #### Spectra ####
        omega, k, time_spectra, space_spectra = frozen_turbulence(datas_u3, ind, zp, nt, split_time, dt, n1, dx=dx, ch="spectra")
        
        frozen_turbulence_plot(fig1u3, col, row, omega = omega, Uc = U1, time_spectra = time_spectra, k = k, space_spectra = space_spectra, ch = "spectra")
        
        del time_spectra
        del space_spectra
        del omega
        del k
        
        del datas_u3
        
        col +=1
        cpt +=1
        if zplan == 4:
            row +=1
            col = 1
            
        
    # Update layout properties for 4 plots
    if chplot == "normal":
        fig1u1.update_layout(height=600, width=900, title_text="Power spectra streawise u1", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.04, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
        fig2u1.update_layout(height=600, width=1100, title_text='Autocorrelation comparison Streamwise u1', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
        fig3u1.update_layout(height=600, width=1100, title_text='Correlation 2D Streamwise u1', legend=dict(y=1.2, x=0.9), font=font)
        
        fig1u2.update_layout(height=600, width=900, title_text="Power spectra Streawise u2", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.04, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
        fig2u2.update_layout(height=600, width=1100, title_text='Autocorrelation comparison Streamwise u2', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
        fig3u2.update_layout(height=600, width=1100, title_text='Correlation 2D Streamwise u2', legend=dict(y=1.2, x=0.9), font=font)
        
        fig1u3.update_layout(height=600, width=900, title_text="Power spectra Streawise u3", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
        fig2u3.update_layout(height=600, width=1100, title_text='Autocorrelation comparison Streamwise u3', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
        fig3u3.update_layout(height=600, width=1100, title_text='Correlation 2D Streamwise u3', legend=dict(y=1.2, x=0.9), font=font)
        
        figU.add_trace(go.Scatter(x=U_ratio, y=X_ratio, line=dict(color='midnightblue')))
        figU.update_xaxes(title='$U_c/U_1$')
        figU.update_yaxes(title='$z^+$')
        figU.update_layout(title='Velocity ratio', font=font, showlegend=False)
        
        figU1c.add_trace(go.Scatter(x=U1_list, y=X_ratio, name='$U_1$', line=dict(color='midnightblue')))
        figU1c.add_trace(go.Scatter(x=Uc_list, y=X_ratio, name='$U_c$', line=dict(color='firebrick')))
        figU1c.update_xaxes(title='$\\text{mean velocity}~(m.s^{-1})$')
        figU1c.update_yaxes(title='$z^+$')
        figU1c.update_layout(title='Velocity comparison', font=font, legend=dict(yanchor="bottom", xanchor="right"))
        
        if split_time == 'Y':
            save_figures(fig1u1, "split_time/frozen_turbulence/power_spectra/u1.png")
            #save_figures(fig2u1, "split_time/frozen_turbulence/correlation_st/u1.png")
            save_figures(fig3u1, "split_time/frozen_turbulence/correlation2D/u1.png")
            
            save_figures(fig1u2, "split_time/frozen_turbulence/power_spectra/u2.png")
            #save_figures(fig2u2, "split_time/frozen_turbulence/correlation_st/u2.png")
            save_figures(fig3u2, "split_time/frozen_turbulence/correlation2D/u2.png")
        
            save_figures(fig1u3, "split_time/frozen_turbulence/power_spectra/u3.png")
            #save_figures(fig2u3, "split_time/frozen_turbulence/correlation_st/u3.png")
            save_figures(fig3u3, "split_time/frozen_turbulence/correlation2D/u3.png")
            
            save_figures(figU, "split_time/frozen_turbulence/correlation2D/u_ratio.png")
            
            save_figures(figU1c, "split_time/frozen_turbulence/correlation2D/u_1c.png")
        
    # # Update layout properties for full plots.
    if chplot == "all":
        fig1u1.update_layout(height=900, width=900, title_text="Power spectra streawise u1", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
        fig2u1.update_layout(height=600, width=1100, title_text='Autocorrelation comparison Streamwise u1', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
        fig3u1.update_layout(height=900, width=1100, title_text='Correlation 2D Streamwise u1', legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40), font=font)
        
        fig1u2.update_layout(height=900, width=900, title_text="Power spectra Streawise u2", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
        fig2u2.update_layout(height=600, width=1100, title_text='Autocorrelation comparison Streamwise u2', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
        fig3u2.update_layout(height=900, width=1100, title_text='Correlation 2D Streamwise u2', legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40), font=font)

        fig1u3.update_layout(height=900, width=900, title_text="Power spectra Streawise u3", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
        fig2u3.update_layout(height=600, width=1100, title_text='Autocorrelation comparison Streamwise u3', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
        fig3u3.update_layout(height=900, width=1100, title_text='Correlation 2D Streamwise u3', legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40), font=font)
        
        figU.add_trace(go.Scatter(x=U_ratio, y=X_ratio, line=dict(color='midnightblue')))
        figU.update_xaxes(title='$U_c/U_1$')
        figU.update_yaxes(title='$z^+$')
        figU.update_layout(title='Velocity ratio', font=font, showlegend=False)
        
        figU1c.add_trace(go.Scatter(x=U1_list, y=X_ratio, name='$U_1$', line=dict(color='midnightblue')))
        figU1c.add_trace(go.Scatter(x=Uc_list, y=X_ratio, name='$U_c$', line=dict(color='firebrick')))
        figU1c.update_xaxes(title='$\\text{mean velocity}~(m.s^{-1})$')
        figU1c.update_yaxes(title='$z^+$')
        figU1c.update_layout(title='Velocity comparison', font=font, legend=dict(yanchor="bottom", xanchor="right"))
        
        if split_time == 'Y':
            save_figures(fig1u1, "split_time/frozen_turbulence/power_spectra/u1_all.png")
            #save_figures(fig2u1, "split_time/frozen_turbulence/correlation_st/u1_all.png")
            save_figures(fig3u1, "split_time/frozen_turbulence/correlation2D/u1_all.png")
            
            save_figures(fig1u2, "split_time/frozen_turbulence/power_spectra/u2_all.png")
            #save_figures(fig2u2, "split_time/frozen_turbulence/correlation_st/u2_all.png")
            save_figures(fig3u2, "split_time/frozen_turbulence/correlation2D/u2_all.png")
            
            save_figures(fig1u3, "split_time/frozen_turbulence/power_spectra/u3_all.png")
            #save_figures(fig2u3, "split_time/frozen_turbulence/correlation_st/u3_all.png")
            save_figures(fig3u3, "split_time/frozen_turbulence/correlation2D/u3_all.png")
            
            save_figures(figU, "split_time/frozen_turbulence/correlation2D/u_ratio_all.png")
            
            save_figures(figU1c, "split_time/frozen_turbulence/correlation2D/u_1c_all.png")
            
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    print(f'\n Frozen turbulence valiation done in : {int(minutes)}m {seconds:.2f}s \n')
    
    
    #=======================================================
                    #### GAMMA ####
    #=======================================================
    print("\n========================================")
    print(f"{YELLOW}Gamma determination{RESET}")
    start_time = time.time()
    print("========================================")
    print(f"{YELLOW}Streamwise{RESET}")
    col = 1
    row = 1
    
    print("\nReading input files ...")
    _, x1, x2, _, nt, n1, n2, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane(fpars_files_streamwise_u1[0])
    nt = nt - 1
    dx = cflow.xlen / n1
    num_split_t = nt // split_t
    
    figu1g, figu2g, figu3g = init_figures_gamma(zp, ch=chplot) #Figure initialization
    figu1g_r, figu2g_r, figu3g_r = init_figures_gamma(zp, ch=chplot)
    
    figgamma = go.Figure()
    gamma_u1_1 = []
    gamma_u2_1 = []
    gamma_u3_1 = []
    
    figgamma_r = make_subplots(rows=1, cols=3, shared_yaxes=True, y_title='$\gamma$', subplot_titles=("$u_1$)", "$u_2$)", "$u_3$"))
    gamma_u1_1_r = []
    gamma_u2_1_r = []
    gamma_u3_1_r = []
    gamma_u1_2_r = []
    gamma_u2_2_r = []
    gamma_u3_2_r = []
    gamma_u1_3_r = []
    gamma_u2_3_r = []
    gamma_u3_3_r = []
    gamma_u1_4_r = []
    gamma_u2_4_r = []
    gamma_u3_4_r = []
    
    z_axis = zp
    
    # for zplan in np.arange(0, n2, n2//3, dtype=int):
    for ind, zplan in enumerate(zp_ind):
    
        
        print("========================================")
        print(f'Gamma determination for {YELLOW}zp={zp[ind]:.2f}{RESET}')
        print('Plan number:', zplan)
        
        #### u1 ####
        print("\nReading input files u1 streamwise...")
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
        datas_u1 = var[1:,:,:]
        
        U1 = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=0), axis=-1))
        
        datas_u1 = datas_u1 - U1
    
        #cpt = 0
        # for delta_x in range(1, 20, 5):
        #     ind1, ind2, kc, omega, funct = frozen_turbulence(datas_u1, zplan, z, nt, split_time, dt, n1, dx = dx, delta_x=delta_x, x1 = x1, Uc = U1, ch = "gamma")
            
        #     r = np.abs(x1[delta_x] - x1[0])
            

        #     slop1, slop2 = frozen_turbulence_plot(figu1g, col, row, omega = omega, k = kc, delta_x=delta_x, funct=funct, r=r, z=zp[zplan], ind1=ind1, ind2=ind2, cpt=cpt, ch = "w_gamma")
        #     cpt += 1
        
        _, omega, _ = Gamma_function(datas_u1[0:split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")
        
        funct = np.zeros((omega.shape[0], n1))
        omega = np.zeros((omega.shape[0]))
        Dx = np.zeros((n1))
        for n in tqdm(range(1,num_split_t), desc=f'Gamma', colour= 'GREEN'):
            
            funct += Gamma_function(datas_u1[(n-1)*split_t:n*split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")[0]
            omega += Gamma_function(datas_u1[(n-1)*split_t:n*split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")[1]
            Dx += Gamma_function(datas_u1[(n-1)*split_t:n*split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")[2]
            
        funct /= (num_split_t-1)
        omega /= (num_split_t-1)
        Dx /= (num_split_t-1)
        
        del datas_u1
        
        omega_lim = 3 ##have to be computed (mesh size?)
        for i in range(omega.shape[0]):
            if omega[i]>omega_lim:
                ind0 = i
                break
            
        omega_lim = 40 ##have to be computed (mesh size?)
        for i in range(omega.shape[0]):
            if omega[i]>omega_lim:
                ind1 = i
                break
                
        omega_lim = 70
        for i in range(omega.shape[0]):
            if omega[i]>omega_lim:
                ind2 = i
                break
            
        slop1 = gamma_plot(figu1g, col, row, funct, omega, Dx, ind0, ind1, ind2, ch = 'w')
            
        ind0 = Dx.shape[0]//2
            
        r_lim = 4
        for i in range(Dx.shape[0]):
            if Dx[i]>r_lim:
                ind1 = i
                break
            
        print("ind0",ind0) ; print("ind1",ind1)
        # print("Dx shape:", Dx[n1//2:].shape)
        
        moy1, moy2, moy3, moy4 = gamma_plot(figu1g_r, col, row, funct, omega, Dx, ind0, ind1, ind2, ch = 'x')
        
        r_lim = 0.4
        for i in range(Dx.shape[0]):
            if Dx[i]>r_lim:
                ind2 = i
                break
        
        gamma_u1_1.append(slop1[0]*U1/Dx[ind2])
        #gamma_u1_2.append(slop2[0]*U1/Dx[15])
        
        gamma_u1_1_r.append(moy1*U1/(omega[1]))
        gamma_u1_2_r.append(moy2*U1/(omega[2]))
        gamma_u1_3_r.append(moy3*U1/(omega[3]))
        gamma_u1_4_r.append(moy4*U1/(omega[4]))
            
        del omega
        del funct
        del Dx
        del slop1
        # del slop2
        
        #### u2 ####
        print("\nReading input files u2 streamwise ...")
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
        datas_u2 = var[1:,:,:]
        
        U2 = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=0), axis=-1))
        
        datas_u2 = datas_u2 - U2
        
        # cpt = 0
        # for delta_x in range(1, 20, 5):
        #     ind1, ind2, kc, omega, funct = frozen_turbulence(datas_u2, zplan, z, nt, split_time, dt, n1, dx = dx, delta_x=delta_x, x1 = x1, Uc = U1, ch = "gamma")
            
        #     r = np.abs(x1[delta_x] - x1[0])
            

        #     slop1, slop2 = frozen_turbulence_plot(figu2g, col, row, omega = omega, k = kc, delta_x=delta_x, funct=funct, r=r, z=zp[zplan], ind1=ind1, ind2=ind2, cpt=cpt, ch = "w_gamma")
        #     cpt += 1
        
        _, omega, _ = Gamma_function(datas_u2[0:split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")
        
        funct = np.zeros((omega.shape[0], n1))
        omega = np.zeros((omega.shape[0]))
        Dx = np.zeros((n1))
        for n in tqdm(range(1,num_split_t), desc=f'Gamma', colour= 'GREEN'):
            
            funct += Gamma_function(datas_u2[(n-1)*split_t:n*split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")[0]
            omega += Gamma_function(datas_u2[(n-1)*split_t:n*split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")[1]
            Dx += Gamma_function(datas_u2[(n-1)*split_t:n*split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")[2]
            
        funct /= (num_split_t-1)
        omega /= (num_split_t-1)
        Dx /= (num_split_t-1)
        
        del datas_u2
        
        omega_lim = 3 ##have to be computed (mesh size?)
        for i in range(omega.shape[0]):
            if omega[i]>omega_lim:
                ind0 = i
                break
            
        omega_lim = 40 ##have to be computed (mesh size?)
        for i in range(omega.shape[0]):
            if omega[i]>omega_lim:
                ind1 = i
                break
                
        omega_lim = 70
        for i in range(omega.shape[0]):
            if omega[i]>omega_lim:
                ind2 = i
                break
            
        slop1 = gamma_plot(figu2g, col, row, funct, omega, Dx, ind0, ind1, ind2, ch = 'w')
            
        ind0 = Dx.shape[0]//2
            
        r_lim = 4
        for i in range(Dx.shape[0]):
            if Dx[i]>r_lim:
                ind1 = i
                break
        
        moy1, moy2, moy3, moy4 = gamma_plot(figu2g_r, col, row, funct, omega, Dx, ind0, ind1, ind2, ch = 'x')
        
        r_lim = 0.4
        for i in range(Dx.shape[0]):
            if Dx[i]>r_lim:
                ind2 = i
                break
        
        gamma_u2_1.append(slop1[0]*U1/Dx[ind2])
        #gamma_u2_2.append(slop2[0]*U1/Dx[ind1])
        
        gamma_u2_1_r.append(moy1*U1/(omega[1]))
        gamma_u2_2_r.append(moy2*U1/(omega[2]))
        gamma_u2_3_r.append(moy3*U1/(omega[3]))
        gamma_u2_4_r.append(moy4*U1/(omega[4]))
            
        del omega
        del funct
        del Dx
        del slop1
        # del slop2
        
        #### u3 ####
        print("\nReading input files u3 streamwise ...")
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
        datas_u3 = var[1:,:,:]
        
        U3 = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=0), axis=-1))
        
        datas_u3 = datas_u3 - U3
        
        # cpt = 0
        # for delta_x in range(1, 20, 5):
        #     ind1, ind2, kc, omega, funct = frozen_turbulence(datas_u3, zplan, z, nt, split_time, dt, n1, dx = dx, delta_x=delta_x, x1 = x1, Uc = U1, ch = "gamma")
            
        #     r = np.abs(x1[delta_x] - x1[0])
        

        #     slop1, slop2 = frozen_turbulence_plot(figu3g, col, row, omega = omega, k = kc, delta_x=delta_x, funct=funct, r=r, z=zp[zplan], ind1=ind1, ind2=ind2, cpt=cpt, ch = "w_gamma")
        #     cpt += 1
        
        _, omega, _ = Gamma_function(datas_u3[0:split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")
        
        funct = np.zeros((omega.shape[0], n1))
        omega = np.zeros((omega.shape[0]))
        Dx = np.zeros((n1))
        for n in tqdm(range(1,num_split_t), desc=f'Gamma', colour= 'GREEN'):
            
            funct += Gamma_function(datas_u3[(n-1)*split_t:n*split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")[0] #I hould take Uc here
            omega += Gamma_function(datas_u3[(n-1)*split_t:n*split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")[1]
            Dx += Gamma_function(datas_u3[(n-1)*split_t:n*split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")[2]
            
        funct /= (num_split_t-1)
        omega /= (num_split_t-1)
        Dx /= (num_split_t-1)
        
        del datas_u3
        
        omega_lim = 3 ##have to be computed (mesh size?)
        for i in range(omega.shape[0]):
            if omega[i]>omega_lim:
                ind0 = i
                break
            
        omega_lim = 40 ##have to be computed (mesh size?)
        for i in range(omega.shape[0]):
            if omega[i]>omega_lim:
                ind1 = i
                break
                
        omega_lim = 70
        for i in range(omega.shape[0]):
            if omega[i]>omega_lim:
                ind2 = i
                break
            
        slop1 = gamma_plot(figu3g, col, row, funct, omega, Dx, ind0, ind1, ind2, ch = 'w')
            
        ind0 = Dx.shape[0]//2
            
        r_lim = 4
        for i in range(Dx.shape[0]):
            if Dx[i]>r_lim:
                ind1 = i
                break
        
        moy1, moy2, moy3, moy4 = gamma_plot(figu3g_r, col, row, funct, omega, Dx, ind0, ind1, ind2, ch = 'x')
        
        r_lim = 0.4
        for i in range(Dx.shape[0]):
            if Dx[i]>r_lim:
                ind2 = i
                break
        
        gamma_u3_1.append(slop1[0]*U1/Dx[ind2])
        #gamma_u3_2.append(slop2[0]*U1/Dx[ind1])
        
        gamma_u3_1_r.append(moy1*U1/(omega[1]))
        gamma_u3_2_r.append(moy2*U1/(omega[2]))
        gamma_u3_3_r.append(moy3*U1/(omega[3]))
        gamma_u3_4_r.append(moy4*U1/(omega[4]))
            
        del funct
        del Dx
        del slop1
        # del slop2
        
        
        col +=1
        if zplan == 4:
            row +=1
            col = 1
            
            
    # Update layout properties for 4 plots
    if chplot == "normal":
        figu1g.update_layout(height=600, width=1100, title_text='Gamma determination u1', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        figu2g.update_layout(height=600, width=1100, title_text='Gamma determination u2', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        figu3g.update_layout(height=600, width=1100, title_text='Gamma determination u3', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        figu1g_r.update_layout(height=600, width=1100, title_text='Gamma determination u1', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        figu2g_r.update_layout(height=600, width=1100, title_text='Gamma determination u2', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        figu3g_r.update_layout(height=600, width=1100, title_text='Gamma determination u3', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        ## Gamma for omega ##
        figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_1, name='$\gamma_{u1}$', line=dict(color='midnightblue')))
        figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_1, name='$\gamma_{u2}$', line=dict(color='firebrick')))
        figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_1, name='$\gamma_{u3}$', line=dict(color='darkgreen')))
        
        figgamma.update_xaxes(title_text='$z^+$')
        figgamma.update_yaxes(title_text='$\gamma$')
        figgamma.update_layout(height=600, width=800, title_text='Gamma evolution', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        ## Gamme for r ##
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_1_r, showlegend=False, line=dict(color='midnightblue')), row=1, col=1)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_2_r, showlegend=False, line=dict(color='firebrick')), row=1, col=1)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_3_r, showlegend=False, line=dict(color='green')), row=1, col=1)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_4_r, showlegend=False, line=dict(color='purple')), row=1, col=1)
        
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_1_r, showlegend=False, line=dict(color='midnightblue')), row=1, col=2)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_2_r, showlegend=False, line=dict(color='firebrick')), row=1, col=2)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_3_r, showlegend=False, line=dict(color='green')), row=1, col=2)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_4_r, showlegend=False, line=dict(color='purple')), row=1, col=2)
        
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_1_r, name=f'$\omega={omega[1]:.2f}$', line=dict(color='midnightblue')), row=1, col=3)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_2_r, name=f'$\omega={omega[2]:.2f}$', line=dict(color='firebrick')), row=1, col=3)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_3_r, name=f'$\omega={omega[3]:.2f}$', line=dict(color='green')), row=1, col=3)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_4_r, name=f'$\omega={omega[4]:.2f}$', line=dict(color='purple')), row=1, col=3)
        
        del omega
        
        figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=1)
        figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=2)
        figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=3)
        figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=4)
        figgamma_r.update_layout(height=600, width=1100, title_text='Gamma evolution', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        
        if split_time == 'Y':

            save_figures(figu1g, "split_time/gamma/gamma_u1_w.png")
            
            save_figures(figu2g, "split_time/gamma/gamma_u2_w.png")
        
            save_figures(figu3g, "split_time/gamma/gamma_u3_w.png")
            
            save_figures(figu1g_r, "split_time/gamma/gamma_u1_r.png")
            
            save_figures(figu2g_r, "split_time/gamma/gamma_u2_r.png")
        
            save_figures(figu3g_r, "split_time/gamma/gamma_u3_r.png")
            
            save_figures(figgamma, "split_time/gamma/gamma_view_w.png")
            
            save_figures(figgamma_r, "split_time/gamma/gamma_view_r.png")
        
    # # Update layout properties for full plots.
    if chplot == "all":

        figu1g.update_layout(height=900, width=1100, title_text='Gamma determination u1', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        figu2g.update_layout(height=900, width=1100, title_text='Gamma determination u2', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        figu3g.update_layout(height=900, width=1100, title_text='Gamma determination u3', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        figu1g_r.update_layout(height=900, width=1100, title_text='Gamma determination u1', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        figu2g_r.update_layout(height=900, width=1100, title_text='Gamma determination u2', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        figu3g_r.update_layout(height=900, width=1100, title_text='Gamma determination u3', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        ## Gamma for omega ##
        figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_1, name='$\gamma_{u1}$', line=dict(color='midnightblue')))
        figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_1, name='$\gamma_{u2}$', line=dict(color='firebrick')))
        figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_1, name='$\gamma_{u3}$', line=dict(color='darkgreen')))
        
        figgamma.update_xaxes(title_text='$z^+$')
        figgamma.update_yaxes(title_text='$\gamma$')
        figgamma.update_layout(height=600, width=800, title_text='Gamma evolution', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        ## Gamme for r ##
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_1_r, showlegend=False, line=dict(color='midnightblue')), row=1, col=1)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_2_r, showlegend=False, line=dict(color='firebrick')), row=1, col=1)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_3_r, showlegend=False, line=dict(color='green')), row=1, col=1)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_4_r, showlegend=False, line=dict(color='purple')), row=1, col=1)
        
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_1_r, showlegend=False, line=dict(color='midnightblue')), row=1, col=2)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_2_r, showlegend=False, line=dict(color='firebrick')), row=1, col=2)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_3_r, showlegend=False, line=dict(color='green')), row=1, col=2)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_4_r, showlegend=False, line=dict(color='purple')), row=1, col=2)
        
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_1_r, name=f'$\omega={omega[1]:.2f}$', line=dict(color='midnightblue')), row=1, col=3)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_2_r, name=f'$\omega={omega[2]:.2f}$', line=dict(color='firebrick')), row=1, col=3)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_3_r, name=f'$\omega={omega[3]:.2f}$', line=dict(color='green')), row=1, col=3)
        figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_4_r, name=f'$\omega={omega[4]:.2f}$', line=dict(color='purple')), row=1, col=3)
        
        del omega
        
        figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=1)
        figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=2)
        figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=3)
        figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=4)
        figgamma_r.update_layout(height=600, width=1100, title_text='Gamma evolution', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        if split_time == 'Y':

            save_figures(figu1g, "split_time/gamma/gamma_u1_w_all.png")

            save_figures(figu2g, "split_time/gamma/gamma_u2_w_all.png")

            save_figures(figu3g, "split_time/gamma/gamma_u3_w_all.png")
            
            save_figures(figu1g_r, "split_time/gamma/gamma_u1_r_all.png")
            
            save_figures(figu2g_r, "split_time/gamma/gamma_u2_r_all.png")
        
            save_figures(figu3g_r, "split_time/gamma/gamma_u3_r_all.png")
            
            save_figures(figgamma, "split_time/gamma/gamma_view_w_all.png")
            
            save_figures(figgamma_r, "split_time/gamma/gamma_view_r_all.png")
            
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    print(f'\n Gamma determination done in : {int(minutes)}m {seconds:.2f}s \n')
    
    
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
    
    fig1sc= init_figures_sc(zp, ch=chplot) #Figure initialization
    
    #for zplan in np.arange(0, n2, n2//3, dtype=int):
    for ind, zplan in enumerate(zp_ind):
        
        print("========================================")
        print(f'Space correlation for {YELLOW}zp={zp[ind]:.2f}{RESET}')
        print('Plan number:', zplan)
        print("Reading input files streamwise...")
        
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
            for n in tqdm(range(1,num_split_t), desc=f'PSD (zp={zp[ind]:.2f})', colour= 'GREEN'):
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
            for n in tqdm(range(1,num_split_t), desc=f'PSD (zp={zp[ind]:.2f})', colour= 'GREEN'):
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
            for n in tqdm(range(1,num_split_t), desc=f'PSD (zp={zp[ind]:.2f})', colour= 'GREEN'):
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
            
        
    if chplot == 'normal':
        fig1sc.update_layout(height=600, width=900, title_text='Space Correlation Streamwise', font=font,  legend=dict(yanchor="bottom", y=1.03, xanchor="right", x=1, orientation='h'))
    if chplot == 'all':
        fig1sc.update_layout(height=900, width=900, title_text='Space Correlation Streamwise', font=font,  legend=dict(yanchor="bottom", y=1.03, xanchor="right", x=1, orientation='h'))
        
    if split_time == 'Y':
        
        if chplot == 'normal':
            save_figures(fig1sc, "split_time/space_correlation/streamwise.png")
        if chplot == 'all':
            save_figures(fig1sc, "split_time/space_correlation/streamwise_all.png")
            
    print("\n========================================")
    print(f"{YELLOW}Spanwise{RESET}")
    col = 1
    row = 1
    print("Reading input files ...")
    
    _, x1, x2, _, nt, n2, n1, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane(fpars_files_spanwise_u1[0])
    nt = nt - 1
    print('n1:', n1)
    print('n2:', n2)
    dx = cflow.xlen / n1
    
    fig2sc= init_figures_sc(zp, ch=chplot) #Figure initialization
    
    #for zplan in np.arange(0, n2, n2//3, dtype=int):
    for ind, zplan in enumerate(zp_ind):
        
        print("========================================")
        print(f'Space correlation for {YELLOW}zp={zp[ind]:.2f}{RESET}')
        print('Plan number:', zplan)
        print("Reading input files spanwise...")
        
        num_split_t = nt // split_t
        Dx = np.linspace(0,ylen//2,n1)
        
        ## uu ##
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_spanwise_u1[zplan])
        datas_u1 = var[1:,:,:]
        
        Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
        
        datas_u1 = datas_u1 - Ux
        
        if split_time == 'Y':
            
            R_uu = np.zeros((n1))
            for n in tqdm(range(1,num_split_t), desc=f'PSD (z={zp[ind]:.2f})', colour= 'GREEN'):
                R_uu += Space_correation(datas_u1[(n-1)*split_t:n*split_t,:,:], datas_u1[(n-1)*split_t:n*split_t,:,:], geom = "plan", mode_corr = 'half', axis = "spanwise")
            R_uu /= (num_split_t-1)
            
        if split_time == 'n':
            R_uu = Space_correation(datas_u1, datas_u1, geom = "plan", mode_corr = 'half', axis = "spanwise")
            
        space_correlation_plot(fig2sc, col, row, Dx, R_uu, name = '$R_{UU}$', color='midnightblue', axis='spanwise')
        del R_uu
        del datas_u1
        
        ## vv ##
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_spanwise_u2[zplan])
        datas_u2 = var[1:,:,:]
        
        Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
        
        datas_u2 = datas_u2 - Uy
        
        if split_time == 'Y':
            
            R_vv = np.zeros((n1))
            for n in tqdm(range(1,num_split_t), desc=f'PSD (z={zp[ind]:.2f})', colour= 'GREEN'):
                R_vv += Space_correation(datas_u2[(n-1)*split_t:n*split_t,:,:], datas_u2[(n-1)*split_t:n*split_t,:,:], geom = "plan", mode_corr = 'half', axis = "spanwise")
            R_vv /= (num_split_t-1)
        
        if split_time == 'n':
            R_vv = Space_correation(datas_u2, datas_u2, geom = "plan", mode_corr = 'half', axis = "spanwise")
            
        space_correlation_plot(fig2sc, col, row, Dx, R_vv, name = '$R_{VV}$', color='firebrick', axis='spanwise')
        del R_vv
        del datas_u2
        
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_spanwise_u3[zplan])
        datas_u3 = var[1:,:,:]
        
        Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
        
        datas_u3 = datas_u3 - Uz
        
        if split_time == 'Y':
            
            R_ww = np.zeros((n1))
            for n in tqdm(range(1,num_split_t), desc=f'PSD (z={zp[ind]:.2f})', colour= 'GREEN'):
                R_ww += Space_correation(datas_u3[(n-1)*split_t:n*split_t,:,:], datas_u3[(n-1)*split_t:n*split_t,:,:], geom = "plan", mode_corr = 'half', axis = "spanwise")
            R_ww /= (num_split_t-1)
        
        if split_time == 'n':
            R_ww = Space_correation(datas_u3, datas_u3, geom = "plan", mode_corr = 'half', axis = "spanwise")
            
        space_correlation_plot(fig2sc, col, row, Dx, R_ww, name = '$R_{WW}$', color='darkgreen', axis='spanwise')
        del R_ww
        del datas_u3
        
        col +=1
        if zplan == 4:
            row +=1
            col = 1
            
        
    if chplot == 'normal':
        fig2sc.update_layout(height=600, width=900, title_text='Space Correlation Spanwise', font=font,  legend=dict(yanchor="bottom", y=1.03, xanchor="right", x=1, orientation='h'))
    if chplot == 'all':
        fig2sc.update_layout(height=900, width=900, title_text='Space Correlation Spanwise', font=font,  legend=dict(yanchor="bottom", y=1.03, xanchor="right", x=1, orientation='h'))
        
    if split_time == 'Y':
        
        if chplot == 'normal':
            save_figures(fig2sc, "split_time/space_correlation/spanwise.png")
        if chplot == 'all':
            save_figures(fig2sc, "split_time/space_correlation/spanwise_all.png")
            
        
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    print(f'\n Space correlation done in : {int(minutes)}m {seconds:.2f}s \n')
        
    
    
    #=======================================================
                #### RANS AND LES ANALYSIS ####
    #=======================================================
    # print("\n========================================")
    # print(f"{YELLOW}RANS AND LES COMPARISON{RESET}")
    # start_time = time.time()
    # print("========================================")
    # zp_RANS = cflow.ut/cflow.nu * normal
    # fig_RANS = make_subplots(rows=1, cols=2, shared_yaxes= True, y_title='$z^+$')
    # fig_vel_profil = make_subplots(rows=1, cols=3, shared_yaxes= True, y_title='$z^+$')
    # fig_var = make_subplots(rows=1, cols=3, shared_yaxes= True, y_title='$z^+$')
    # fig_vel_ratio_profil = go.Figure()
    # fig_up = go.Figure()
    # fig_RANS.add_trace(go.Scatter(x=u_velocity, y=zp_RANS, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'), name='$U\\text{(RANS)}$'), row=1, col=1)
    # fig_RANS.add_trace(go.Scatter(x=ko2_tke, y=zp_RANS, mode= 'lines+markers', line=dict(color='darkgreen', width=2), marker=dict(symbol='circle'), name='$k_T\\text{(RANS)}$'), row=1, col=2)
    
    # U = []
    # V = []
    # W = []
    # var1 = []
    # var2 = []
    # var3 = []
    # kt = []
    # for zplan in zp_ind:
    #     print('zp:', zp[zplan])
        
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
    #     datas_u1 = var[1:,:,:]
    #     del var
    #     #Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
    #     Ux = np.mean(datas_u1[:,:,:], axis=0)
    #     datas_u1 = datas_u1 - Ux[np.newaxis,:,:]
    #     varUx = np.mean(np.mean(np.mean((datas_u1)**2, axis=0), axis=-1))
    #     U.append(np.mean(np.mean(Ux, axis=-1)))
    #     var1.append(varUx)
    #     del datas_u1
        
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
    #     datas_u2 = var[1:,:,:]
    #     del var
    #     #Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
    #     Uy = np.mean(datas_u2[:,:,:], axis=0)
    #     datas_u2 = datas_u2 - Uy[np.newaxis,:,:]
    #     varUy = np.mean(np.mean(np.mean((datas_u2)**2, axis=0), axis=-1))
    #     V.append(np.mean(np.mean(Uy, axis=-1)))
    #     var2.append(varUy)
    #     del datas_u2
        
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
    #     datas_u3 = var[1:,:,:]
    #     del var
    #     #Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
    #     Uz = np.mean(datas_u3[:,:,:], axis=0)
    #     datas_u3 = datas_u3 - Uz[np.newaxis,:,:]
    #     varUz = np.mean(np.mean(np.mean((datas_u3)**2, axis=0), axis=-1))
    #     W.append(np.mean(np.mean(Uz, axis=-1)))
    #     var3.append(varUz)
    #     del datas_u3
        
    #     kt.append(varUx + varUy + varUz)
        
    # kt = np.array(kt)
    # U = np.array(U)
    # V = np.array(V)
    # W = np.array(W)
    # up = U / cflow.ut
    # var1 = np.array(var1)
    # var2 = np.array(var2)
    # var3 = np.array(var3)
    # kt /= 2.
    # ratio2 = var2 / var1
    # ratio3 = var3 / var1
    # ratio1 = var3 / var2
    
    # fig_RANS.add_trace(go.Scatter(x=U, y=z, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='x'), name='$U\\text{(LES)}$'), row=1, col=1)
    # fig_RANS.add_trace(go.Scatter(x=kt, y=zp, mode= 'lines+markers', line=dict(color='darkgreen', width=2), marker=dict(symbol='x'), name='$k_T\\text{(LES)}$'), row=1, col=2)
    
    # fig_RANS.update_xaxes(title='$\\text{velocity}~(m.s^{-1})$', row=1, col=1)
    # fig_RANS.update_xaxes(title='$\\text{kinetic energy}~(m^2.s^{-2})$', row=1, col=2)
    # fig_RANS.update_layout(height=600, width=800, title="RANS and LES data profiles", font=font, showlegend=True, legend=dict(yanchor='bottom', y=1.03, xanchor='left', x=0.9))
    
    # save_figures(fig_RANS, "RANS/RANS_profiles.png")
    
    # fig_vel_profil.add_trace(go.Scatter(x=U, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'), name='$U\\text{(LES)}$'), row=1, col=1)
    # fig_vel_profil.add_trace(go.Scatter(x=V, y=zp, mode= 'lines+markers', line=dict(color='midnightblue', width=2), marker=dict(symbol='diamond'), name='$V\\text{(LES)}$'), row=1, col=2)
    # fig_vel_profil.add_trace(go.Scatter(x=W, y=zp, mode= 'lines+markers', line=dict(color='darkgreen', width=2), marker=dict(symbol='x'), name='$W\\text{(LES)}$'),  row=1, col=3)
    
    # fig_vel_profil.update_xaxes(title='$U~(m.s^{-1})$', row=1, col=1)
    # fig_vel_profil.update_xaxes(title='$V~(m.s^{-1})$', row=1, col=2, exponentformat = 'e')
    # fig_vel_profil.update_xaxes(title='$W~(m.s^{-1})$', row=1, col=3, exponentformat = 'e')
    # fig_vel_profil.update_layout(height=600, width=800, title="Velocity profile (LES data)", font=font, showlegend=True, legend=dict(yanchor='bottom', y=1.03, xanchor='left', x=0.9))
    
    # save_figures(fig_vel_profil, "RANS/mean_velocity_profiles.png")
    
    # fig_var.add_trace(go.Scatter(x=var1, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'), name='$\overline{u_1u_1}$'), row=1, col=1)
    # fig_var.add_trace(go.Scatter(x=var2, y=zp, mode= 'lines+markers', line=dict(color='midnightblue', width=2), marker=dict(symbol='diamond'), name='$\overline{u_2u_2}$'), row=1, col=2)
    # fig_var.add_trace(go.Scatter(x=var3, y=zp, mode= 'lines+markers', line=dict(color='darkgreen', width=2), marker=dict(symbol='x'), name='$\overline{u_3u_3}$'),  row=1, col=3)
    
    # fig_var.update_xaxes(title='$\overline{u_1u_1}$', row=1, col=1)
    # fig_var.update_xaxes(title='$\overline{u_2u_2}$', row=1, col=2, exponentformat = 'e')
    # fig_var.update_xaxes(title='$\overline{u_3u_3}$', row=1, col=3, exponentformat = 'e')
    # fig_var.update_layout(height=600, width=800, title="Variance velocity fluctuation profile", font=font, showlegend=True, legend=dict(yanchor='bottom', y=1.03, xanchor='left', x=0.9))
    
    # save_figures(fig_var, "RANS/var_velocity_profiles.png")
    
    # fig_vel_ratio_profil.add_trace(go.Scatter(x=ratio2, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'), name='$\overline{u_2u_2} / \overline{u_1u_1}$'))
    # fig_vel_ratio_profil.add_trace(go.Scatter(x=ratio3, y=zp, mode= 'lines+markers', line=dict(color='midnightblue', width=2), marker=dict(symbol='diamond'), name='$\overline{u_3u_3} / \overline{u_1u_1}$'))
    # fig_vel_ratio_profil.add_trace(go.Scatter(x=ratio1, y=zp, mode= 'lines+markers', line=dict(color='green', width=2), marker=dict(symbol='x'), name='$\overline{u_3u_3} / \overline{u_2u_2}$'))
    
    # fig_vel_ratio_profil.update_yaxes(title='$z^+$')
    # fig_vel_ratio_profil.update_xaxes(title='velocity ratio')
    # fig_vel_ratio_profil.update_layout(height=600, width=800, title="Variance velocity ratio profile (LES data)", font=font, showlegend=True, legend=dict(yanchor='bottom', y=1.03, xanchor='left', x=0.9))
    
    # save_figures(fig_vel_ratio_profil, "RANS/velocity_ratio_profiles.png")
    
    # fig_up.add_trace(go.Scatter(x=zp, y=up, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle-open'), showlegend=False))
    
    # fig_up.update_yaxes(title='$u^+$')
    # fig_up.update_xaxes(title='$z^+$')
    # fig_up.update_layout(height=600, width=800, title="Adimentionelized velocity profile", font=font, showlegend=True, legend=dict(yanchor='bottom', y=1.03, xanchor='left', x=0.9))
    
    # save_figures(fig_up, "RANS/velocity_ad.png")
    
    # elapsed_time = time.time() - start_time
    # minutes, seconds = divmod(elapsed_time, 60)
    
    # print(f'\n RANS study done in : {int(minutes)}m {seconds:.2f}s \n')
    
    
    #=======================================================
                #### NORMAL PLAN COMPUTATION ####
    #=======================================================
    # print("\n========================================")
    # print(f"{YELLOW}NORMAL PLAN COMPUTATION{RESET}")
    # start_time = time.time()
    # print("========================================")
    # nlines = len(fpars_files_normal_u1)
    # _, x1, x2, _, nt, n2, n1, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane_line(fpars_files_normal_u1[0])
    # nt = nt - 1
    # print('n1:', n1)
    # print('n2:', n2)
    # print('x2:', x2)
    # x2 = np.array(x2)
    # x2 *= cflow.ut / cflow.nu
    # zp_RANS = cflow.ut/cflow.nu * normal
    # print('x2[0]:',x2[0])
    # print('x2[len(x2)]:',x2[len(x2)-1])
    
    # #Figure initiaization##
    # fig_u_z = make_subplots(rows=1, cols=3, shared_yaxes= True, y_title='$z^+$')
    # #fig_corr_z = go.Figure()
    
    # print('U1 profile')
    # data_u1 = np.zeros((nt, n1), dtype=float)
    # for line in range(nlines):
    #     print('line number:', line)
        
    #     var = read_fpar_extract_plane_line(fpars_files_normal_u3[line])[3]
    #     data_u1 += var[1:,:]
    #     del var
        
    # data_u1 /= nlines
    # data_u1 = np.mean(data_u1, axis=0)
    # #data_u1 /= cflow.ut
    # fig_u_z.add_trace(go.Scatter(x=data_u1, y=x2, mode= 'lines', line=dict(color='firebrick', width=2), name='$U_1(LES)$'), row=1, col=1)
    # fig_u_z.add_trace(go.Scatter(x=u_velocity, y=zp_RANS, mode= 'markers', marker=dict(color='firebrick', symbol='circle-open'), name='$U_1(RANS)$'), row=1, col=1)
    
    # del data_u1
    
    # print('U2 profile')
    # data_u2 = np.zeros((nt, n1), dtype=float)
    # for line in range(nlines):
    #     print('line number:', line)
        
    #     var = read_fpar_extract_plane_line(fpars_files_normal_u2[line])[3]
    #     data_u2 += var[1:,:]
    #     del var
        
    # data_u2 /= nlines
    # data_u2 = np.mean(data_u2, axis=0)
    # #data_u2 /= cflow.ut
    # fig_u_z.add_trace(go.Scatter(x=data_u2, y=x2, mode= 'lines', line=dict(color='midnightblue', width=2), name='$U_2(LES)$'), row=1, col=2)
    # #fig_u_z.add_trace(go.Scatter(x=v_velocity, y=zp_RANS, mode= 'markers', marker=dict(color='midnightblue', symbol='circle-open'), name='$U_2(RANS)$'), row=1, col=2)
    
    # del data_u2
    
    # print('U3 profile')
    # data_u3 = np.zeros((nt, n1), dtype=float)
    # for line in range(nlines):
    #     print('line number:', line)
        
    #     var = read_fpar_extract_plane_line(fpars_files_normal_u1[line])[3]
    #     data_u3 += var[1:,:]
    #     del var
        
    # data_u3 /= nlines
    # data_u3 = np.mean(data_u3, axis=0)
    # #data_u3 /= cflow.ut
    # fig_u_z.add_trace(go.Scatter(x=data_u3, y=x2, mode= 'lines', line=dict(color='darkgreen', width=2), name='$U_3(LES)$'), row=1, col=3)
    # #fig_u_z.add_trace(go.Scatter(x=w_velocity, y=zp_RANS, mode= 'markers', marker=dict(color='darkgreen', symbol='circle-open'), name='$U_3(RANS)$'), row=1, col=3)
    
    # del data_u3
    
    # #fig_u1_z.update_yaxes(title='$z^+$')
    # fig_u_z.update_xaxes(title='velocity', row=1, col=1)
    # fig_u_z.update_xaxes(title='velocity', row=1, col=2)
    # fig_u_z.update_xaxes(title='velocity', row=1, col=3)
    # fig_u_z.update_layout(height=600, width=800, title="Wall-normal velocity profile", font=font, showlegend=True, legend=dict(yanchor='bottom', xanchor='right'))
    
    # save_figures(fig_u_z, "Normal_plan/velocity_profiles.png")
    
    # ###############
    # # Correation ##
    # ###############
    
    # if len(zp) == 4:
    #     fig_corr_z = make_subplots(rows=1, cols=4, shared_yaxes= True, y_title='$\delta z^+$', subplot_titles=(f"$z_0={zp[0]:.2f}$", f"$z_0={zp[1]:.2f}$", f"$z_0={zp[2]:.2f}$", f"$z_0={zp[3]:.2f}$"))
    # else:
    #     fig_corr_z = make_subplots(rows=1, cols=4, shared_yaxes= True, y_title='$\delta z^+$', subplot_titles=(f"$z_0={z[0]:.2f}$", f"$z_0={z[1]:.2f}$", f"$z_0={z[2]:.2f}$", f"$z_0={z[3]:.2f}$", f"$z_0={z[4]:.2f}$", f"$z_0={z[5]:.2f}$", f"$z_0={z[6]:.2f}$", f"$z_0={z[7]:.2f}$", f"$z_0={z[8]:.2f}$", f"$z_0={z[9]:.2f}$"))
    
    # print('\nUU correlation')
    # data_u1 = np.zeros((nt, n1), dtype=float)
    # for line in range(nlines):
    #     print('line number:', line)
        
    #     var = read_fpar_extract_plane_line(fpars_files_normal_u3[line])[3]
    #     Ux = np.mean(var[1:,:], axis=0)
    #     data_u1 += var[1:,:] - Ux[np.newaxis, :]
    #     del var
        
    # data_u1 /= nlines
    # corr = np.zeros((n1), dtype=float)
    
    # col = 1
    # row = 1
    # for v, lim in enumerate(zp):
    #     for i in range(x2.shape[0]):
    #         if x2[i] > lim:
    #             ind = i
    #             break
            
    #     for t in range (nt):
    #         full_corr = signal.correlate(data_u1[t,ind:], data_u1[t,ind:], mode='full', method='auto')
    #         corrp = full_corr[full_corr.shape[0]//2:]
    #         corrp /= max(corrp)
    #         full_corr = signal.correlate(data_u1[t,0:ind], data_u1[t,0:ind], mode='full', method='auto')
    #         corrm = full_corr[full_corr.shape[0]//2:]
    #         corrm /= max(corrm)
    #         corr[:] += np.concatenate((corrm[::-1], corrp))
            
    #     corr /= nt
        
    #     print('corr shape:', corr.shape)
    #     print('x2 shape:', x2.shape)
        
    #     if row == 1 and col == 4:
    #         fig_corr_z.add_trace(go.Scatter(x=corr, y=x2, mode= 'lines', line=dict(color='firebrick', width=2), name='$R_{UU}^{(3)}$'), row=row, col=col)
    #     else:
    #         fig_corr_z.add_trace(go.Scatter(x=corr, y=x2, mode= 'lines', line=dict(color='firebrick', width=2), showlegend=False), row=row, col=col)
        
    #     col +=1
    #     if v == 4:
    #         row +=1
    #         col = 1
    
    # del data_u1
    # del corr
    
    # print('\nVV correlation')
    # data_u2 = np.zeros((nt, n1), dtype=float)
    # for line in range(nlines):
    #     print('line number:', line)
        
    #     var = read_fpar_extract_plane_line(fpars_files_normal_u2[line])[3]
    #     Uy = np.mean(var[1:,:], axis=0)
    #     data_u2 += var[1:,:] - Uy[np.newaxis, :]
    #     del var
        
    # data_u2 /= nlines
    # corr = np.zeros((n1), dtype=float)
    
    # col = 1
    # row = 1
    # for v, lim in enumerate(zp):
    #     for i in range(x2.shape[0]):
    #         if x2[i] > lim:
    #             ind = i
    #             break
            
    #     for t in range (nt):
    #         full_corr = signal.correlate(data_u2[t,ind:], data_u2[t,ind:], mode='full', method='auto')
    #         corrp = full_corr[full_corr.shape[0]//2:]
    #         corrp /= max(corrp)
    #         full_corr = signal.correlate(data_u2[t,0:ind], data_u2[t,0:ind], mode='full', method='auto')
    #         corrm = full_corr[full_corr.shape[0]//2:]
    #         corrm /= max(corrm)
    #         corr[:] += np.concatenate((corrm[::-1], corrp))
            
    #     corr /= nt
        
    #     print('corr shape:', corr.shape)
    #     print('x2 shape:', x2.shape)
        
    #     if row == 1 and col == 4:
    #         fig_corr_z.add_trace(go.Scatter(x=corr, y=x2, mode= 'lines', line=dict(color='midnightblue', width=2), name='$R_{VV}^{(3)}$'), row=row, col=col)
    #     else:
    #         fig_corr_z.add_trace(go.Scatter(x=corr, y=x2, mode= 'lines', line=dict(color='midnightblue', width=2), showlegend=False), row=row, col=col)
        
    #     col +=1
    #     if v == 4:
    #         row +=1
    #         col = 1
    
    # del data_u2
    # del corr
    
    # print('\nWW correlation')
    # data_u3 = np.zeros((nt, n1), dtype=float)
    # for line in range(nlines):
    #     print('line number:', line)
        
    #     var = read_fpar_extract_plane_line(fpars_files_normal_u1[line])[3]
    #     Uz = np.mean(var[1:,:], axis=0)
    #     data_u3 += var[1:,:] - Uz[np.newaxis, :]
    #     del var
        
    # data_u3 /= nlines
    # corr = np.zeros((n1), dtype=float)
    
    # col = 1
    # row = 1
    # for v, lim in enumerate(zp):
    #     for i in range(x2.shape[0]):
    #         if x2[i] > lim:
    #             ind = i
    #             break
            
    #     for t in range (nt):
    #         full_corr = signal.correlate(data_u3[t,ind:], data_u3[t,ind:], mode='full', method='auto')
    #         corrp = full_corr[full_corr.shape[0]//2:]
    #         corrp /= max(corrp)
    #         full_corr = signal.correlate(data_u3[t,0:ind], data_u3[t,0:ind], mode='full', method='auto')
    #         corrm = full_corr[full_corr.shape[0]//2:]
    #         corrm /= max(corrm)
    #         corr[:] += np.concatenate((corrm[::-1], corrp))
            
    #     corr /= nt
        
    #     print('corr shape:', corr.shape)
    #     print('x2 shape:', x2.shape)
        
    #     if row == 1 and col == 4:
    #         fig_corr_z.add_trace(go.Scatter(x=corr, y=x2, mode= 'lines', line=dict(color='darkgreen', width=2), name='$R_{WW}^{(3)}$'), row=row, col=col)
    #     else:
    #         fig_corr_z.add_trace(go.Scatter(x=corr, y=x2, mode= 'lines', line=dict(color='darkgreen', width=2), showlegend=False), row=row, col=col)
        
    #     col +=1
    #     if v == 4:
    #         row +=1
    #         col = 1
    
    # del data_u3
    # del corr
    
    # row = 1
    # col = 1
    # for i in range(len(zp)):
    #     fig_corr_z.update_xaxes(title='$R_{ii}^{(3)}$', row=row, col=col)
    #     col +=1
    #     if i == 4:
    #         row +=1
    #         col = 1
        
    # fig_corr_z.update_layout(height=600, width=900, title="Wall-normal auto-correlation", font=font, showlegend=True, legend=dict(yanchor='bottom', xanchor='right'))
    
    # save_figures(fig_corr_z, "Normal_plan/autocorrelation_z.png")
    
    # #########################
    # # Decorrelation Length ##
    # #########################
    
    # file_int = open(out_path+'Normal_plan/Lii.dat', 'w')
    
    # if len(zp) == 4:
    #     fig_liip_z = make_subplots(rows=1, cols=4, shared_yaxes= True, y_title='$R_{ii}^+(\omega)$', subplot_titles=(f"$z_0={zp[0]:.2f}$", f"$z_0={zp[1]:.2f}$", f"$z_0={zp[2]:.2f}$", f"$z_0={zp[3]:.2f}$"))
    #     fig_liim_z = make_subplots(rows=1, cols=4, shared_yaxes= True, y_title='$R_{ii}^-(\omega)$', subplot_titles=(f"$z_0={zp[0]:.2f}$", f"$z_0={zp[1]:.2f}$", f"$z_0={zp[2]:.2f}$", f"$z_0={zp[3]:.2f}$"))
    # else:
    #     fig_liip_z = make_subplots(rows=1, cols=4, shared_yaxes= True, y_title=('$R_{ii}^+(\omega)$','$R_{ii}^+(\omega)$'), subplot_titles=(f"$z_0={z[0]:.2f}$", f"$z_0={z[1]:.2f}$", f"$z_0={z[2]:.2f}$", f"$z_0={z[3]:.2f}$", f"$z_0={z[4]:.2f}$", f"$z_0={z[5]:.2f}$", f"$z_0={z[6]:.2f}$", f"$z_0={z[7]:.2f}$", f"$z_0={z[8]:.2f}$", f"$z_0={z[9]:.2f}$"))
    #     fig_liim_z = make_subplots(rows=1, cols=4, shared_yaxes= True, y_title=('$R_{ii}^-(\omega)$','$R_{ii}^-(\omega)$'), subplot_titles=(f"$z_0={z[0]:.2f}$", f"$z_0={z[1]:.2f}$", f"$z_0={z[2]:.2f}$", f"$z_0={z[3]:.2f}$", f"$z_0={z[4]:.2f}$", f"$z_0={z[5]:.2f}$", f"$z_0={z[6]:.2f}$", f"$z_0={z[7]:.2f}$", f"$z_0={z[8]:.2f}$", f"$z_0={z[9]:.2f}$"))
    
    # print('UU decorrelation length')
    # data_u1 = np.zeros((nt, n1), dtype=float)
    # for line in range(nlines):
    #     print('line number:', line)
        
    #     var = read_fpar_extract_plane_line(fpars_files_normal_u3[line])[3]
    #     Ux = np.mean(var[1:,:], axis=0)
    #     data_u1 += var[1:,:] - Ux[np.newaxis, :]
    #     del var
        
    # data_u1 /= nlines

    # file_int.write("u1\n")
    # file_int.write("z+ \t L+ \t L-\n")
    # col = 1
    # row = 1
    # for v, lim in enumerate(zp):
    #     for i in range(x2.shape[0]):
    #         if x2[i] > lim:
    #             ind = i
    #             break
            
    #     corrp = np.zeros((split_t, n1 - ind), dtype=float)
    #     corrm = np.zeros((split_t, ind), dtype=float)
    
    #     num_split_t = nt // split_t
    #     for n in tqdm(range(1,num_split_t), desc=f'PSD normal plan', colour= 'GREEN'):
    #         fourier = fft.fft(data_u1[(n-1)*split_t:n*split_t,:], axis=0, workers=3)
    #         for w in range(split_t):
    #             full_corr = np.real(signal.correlate(fourier[w,ind:], fourier[w,ind:], mode='full', method='auto'))
    #             corrp[w,:] += full_corr[full_corr.shape[0]//2:]

    #             full_corr = np.real(signal.correlate(fourier[w,0:ind], fourier[w,0:ind], mode='full', method='auto'))
    #             corrm[w,:] += full_corr[full_corr.shape[0]//2:]

        
    #     for w in range(split_t):
    #         corrp[w,:] /= (num_split_t-1)
    #         corrm[w,:] /= (num_split_t-1)

    #     omega = 2*np.pi*np.linspace(0, 1./(dt), split_t)
    #     Liip = np.mean(integrate.simps(corrp, omega, axis=0))
    #     Liim = np.mean(integrate.simps(corrm, omega, axis=0))
        
    #     file_int.write(f"{lim} \t {Liip} \t {Liim}\n")
        
    #     corrp = np.mean(corrp, axis=1)
    #     corrm = np.mean(corrm, axis=1)
        
    #     corrp /= max(corrp)
    #     corrm /= max(corrm)
        
    #     if row == 1 and col == 4:
    #         fig_liip_z.add_trace(go.Scatter(x=omega[:split_t//2+1], y=corrp[:split_t//2+1], mode= 'lines', line=dict(color='firebrick', width=2), name='$R_{UU}^+(\omega)$'), row=row, col=col)
    #         fig_liim_z.add_trace(go.Scatter(x=omega[:split_t//2+1], y=corrm[:split_t//2+1], mode= 'lines', line=dict(color='firebrick', width=2), name='$R_{UU}^-(\omega)$'), row=row, col=col)
            
    #     else:
    #         fig_liip_z.add_trace(go.Scatter(x=omega[:split_t//2+1], y=corrp[:split_t//2+1], mode= 'lines', line=dict(color='firebrick', width=2), showlegend=False), row=row, col=col)
    #         fig_liim_z.add_trace(go.Scatter(x=omega[:split_t//2+1], y=corrm[:split_t//2+1], mode= 'lines', line=dict(color='firebrick', width=2), showlegend=False), row=row, col=col)
            
        
    #     col +=1
    #     if v == 4:
    #         row +=1
    #         col = 1
    # #del corr
    # del corrp
    # del corrm
    # del data_u1

    
    
    # print('VV decorrelation length')
    # data_u2 = np.zeros((nt, n1), dtype=float)
    # for line in range(nlines):
    #     print('line number:', line)
        
    #     var = read_fpar_extract_plane_line(fpars_files_normal_u2[line])[3]
    #     Uy = np.mean(var[1:,:], axis=0)
    #     data_u2 += var[1:,:] - Uy[np.newaxis, :]
    #     del var
        
    # data_u2 /= nlines
    
    # file_int.write("\nu2\n")
    # file_int.write("z+ \t L+ \t L-\n")
    # col = 1
    # row = 1
    # for v, lim in enumerate(zp):
    #     for i in range(x2.shape[0]):
    #         if x2[i] > lim:
    #             ind = i
    #             break
            
    #     corrp = np.zeros((split_t, n1 - ind), dtype=float)
    #     corrm = np.zeros((split_t, ind), dtype=float)

    
    #     num_split_t = nt // split_t
    #     for n in tqdm(range(1,num_split_t), desc=f'PSD normal plan', colour= 'GREEN'):
    #         fourier = fft.fft(data_u2[(n-1)*split_t:n*split_t,:], axis=0, workers=3)
    #         for w in range(split_t):
    #             full_corr = np.real(signal.correlate(fourier[w,ind:], fourier[w,ind:], mode='full', method='auto'))
    #             corrp[w,:] += full_corr[full_corr.shape[0]//2:]

    #             full_corr = np.real(signal.correlate(fourier[w,0:ind], fourier[w,0:ind], mode='full', method='auto'))
    #             corrm[w,:] += full_corr[full_corr.shape[0]//2:]

        
    #     for w in range(split_t):
    #         corrp[w,:] /= (num_split_t-1)
    #         corrm[w,:] /= (num_split_t-1)

            
    #     omega = 2*np.pi*np.linspace(0, 1./(dt), split_t)
    #     Liip = np.mean(integrate.simps(corrp, omega, axis=0))
    #     Liim = np.mean(integrate.simps(corrm, omega, axis=0))
        
    #     file_int.write(f"{lim} \t {Liip} \t {Liim}\n")
        
    #     corrp = np.mean(corrp, axis=1)
    #     corrm = np.mean(corrm, axis=1)
        
    #     corrp /= max(corrp)
    #     corrm /= max(corrm)
        
    #     if row == 1 and col == 4:
    #         fig_liip_z.add_trace(go.Scatter(x=omega[:split_t//2+1], y=corrp[:split_t//2+1], mode= 'lines', line=dict(color='midnightblue', width=2), name='$R_{VV}^+(\omega)$'), row=row, col=col)
    #         fig_liim_z.add_trace(go.Scatter(x=omega[:split_t//2+1], y=corrp[:split_t//2+1], mode= 'lines', line=dict(color='midnightblue', width=2), name='$R_{VV}^-(\omega)$'), row=row, col=col)
            
    #     else:
    #         fig_liip_z.add_trace(go.Scatter(x=omega[:split_t//2+1], y=corrp[:split_t//2+1], mode= 'lines', line=dict(color='midnightblue', width=2), showlegend=False), row=row, col=col)
    #         fig_liim_z.add_trace(go.Scatter(x=omega[:split_t//2+1], y=corrm[:split_t//2+1], mode= 'lines', line=dict(color='midnightblue', width=2), showlegend=False), row=row, col=col)
            
        
    #     col +=1
    #     if v == 4:
    #         row +=1
    #         col = 1
            
    # del corrp
    # del corrm
    # del data_u2
    
    # print('WW decorrelation length')
    # data_u3 = np.zeros((nt, n1), dtype=float)
    # for line in range(nlines):
    #     print('line number:', line)
        
    #     var = read_fpar_extract_plane_line(fpars_files_normal_u1[line])[3]
    #     Uz = np.mean(var[1:,:], axis=0)
    #     data_u3 += var[1:,:] - Uz[np.newaxis, :]
    #     del var
        
    # data_u3 /= nlines
    
    # file_int.write("\nu3\n")
    # file_int.write("z+ \t L+ \t L-\n")
    # col = 1
    # row = 1
    # for v, lim in enumerate(zp):
    #     for i in range(x2.shape[0]):
    #         if x2[i] > lim:
    #             ind = i
    #             break
            
    #     corrp = np.zeros((split_t, n1 - ind), dtype=float)
    #     corrm = np.zeros((split_t, ind), dtype=float)
    
    #     num_split_t = nt // split_t
    #     for n in tqdm(range(1,num_split_t), desc=f'PSD normal plan', colour= 'GREEN'):
    #         fourier = fft.fft(data_u3[(n-1)*split_t:n*split_t,:], axis=0, workers=3)
    #         for w in range(split_t):
    #             full_corr = np.real(signal.correlate(fourier[w,ind:], fourier[w,ind:], mode='full', method='auto'))
    #             corrp[w,:] += full_corr[full_corr.shape[0]//2:]
                
    #             full_corr = np.real(signal.correlate(fourier[w,0:ind], fourier[w,0:ind], mode='full', method='auto'))
    #             corrm[w,:] += full_corr[full_corr.shape[0]//2:]
                
                
        
    #     for w in range(split_t):
    #         corrp[w,:] /= (num_split_t-1)
    #         corrm[w,:] /= (num_split_t-1)

            
    #     omega = 2*np.pi*np.linspace(0, 1./dt, split_t)
    #     Liip = np.mean(integrate.simps(corrp, omega, axis=0))
    #     Liim = np.mean(integrate.simps(corrm, omega, axis=0))
        
    #     file_int.write(f"{lim} \t {Liip} \t {Liim} \n")
        
    #     corrp = np.mean(corrp, axis=1)
    #     corrm = np.mean(corrm, axis=1)
        
    #     corrp /= max(corrp)
    #     corrm /= max(corrm)
        
    #     if row == 1 and col == 4:
    #         fig_liip_z.add_trace(go.Scatter(x=omega[:split_t//2+1], y=corrp[:split_t//2+1], mode= 'lines', line=dict(color='darkgreen', width=2), name='$R_{WW}^+(\omega)$'), row=row, col=col)
    #         fig_liim_z.add_trace(go.Scatter(x=omega[:split_t//2+1], y=corrm[:split_t//2+1], mode= 'lines', line=dict(color='darkgreen', width=2), name='$R_{WW}^-(\omega)$'), row=row, col=col)
            
    #     else:
    #         fig_liip_z.add_trace(go.Scatter(x=omega[:split_t//2+1], y=corrp[:split_t//2+1], mode= 'lines', line=dict(color='darkgreen', width=2), showlegend=False), row=row, col=col)
    #         fig_liim_z.add_trace(go.Scatter(x=omega[:split_t//2+1], y=corrm[:split_t//2+1], mode= 'lines', line=dict(color='darkgreen', width=2), showlegend=False), row=row, col=col)
            
        
    #     col +=1
    #     if v == 4:
    #         row +=1
    #         col = 1
            
    # del corrp
    # del corrm
    # del data_u3
    
    # row = 1
    # col = 1
    # for i in range(len(zp)):
    #     fig_liip_z.update_xaxes(title='$\omega$', row=row, col=col, range=[0,100])
    #     fig_liim_z.update_xaxes(title='$\omega$', row=row, col=col, range=[0,100])
    #     col +=1
    #     if i == 4:
    #         row +=1
    #         col = 1
    
    # fig_liip_z.update_layout(height=600, width=800, title="Wall-normal correlation length (+)", font=font, showlegend=True, legend=dict(yanchor='bottom', xanchor='right'))
    # fig_liim_z.update_layout(height=600, width=800, title="Wall-normal correlation length (-)", font=font, showlegend=True, legend=dict(yanchor='bottom', xanchor='right'))
    
    # save_figures(fig_liip_z, "Normal_plan/correlation_length_plus.png")
    # save_figures(fig_liim_z, "Normal_plan/correlation_length_minus.png")
    
    
    # file_int.close()
    
    # elapsed_time = time.time() - start_time
    # minutes, seconds = divmod(elapsed_time, 60)
    
    # print(f'\n Normal plan study done in : {int(minutes)}m {seconds:.2f}s \n')

    
    
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
    
    # fig_vanK = init_figures_vk(zp, ch=chplot)
    
    # #for zplan in np.arange(0, n2, n2//3, dtype=int):
    # for ind, zplan in enumerate(zp_ind):
        
    #     print("========================================")
    #     print(f'Von Karman theory for {YELLOW}zp={zp[ind]:.2f}{RESET}')
    #     print('Plan number:', zplan)
    
    #     ### experience ###
    #     ## u1 ##
    #     print(f"Computing spectra from LES datas for {YELLOW}z={zp[ind]:.2f}{RESET} ...\n")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
    #     datas_u1 = var[1:,:,:]
    #     Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
    #     datas_u1 = datas_u1 - Ux
        
    #     omega, k, _, phi11_exp = frozen_turbulence(datas_u1, ind, zp, nt, split_time, dt, n1, dx=dx, ch="spectra")
    #     k /= Ux
        
    #     lim = 400
    #     for i in range(k.shape[0]):
    #         if k[i] > lim:
    #             tick_w = i
    #             break
            
    #     del datas_u1
        
    #     von_karman_plot(fig_vanK, col, row, k[1:tick_w], k[1:tick_w]*phi11_exp[1:tick_w], name = '$\phi_{11}exp$', color = 'firebrick', symbols='x')
        
    #     int_exp1 = integrate.simpson(phi11_exp[1:], k[1:])
        
    #     del phi11_exp
    #     del k

        
    #     ## u2 ##
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
    #     nt = nt - 1
    #     datas_u2 = var[1:,:,:]
    #     Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
    #     datas_u2 = datas_u2 - Uy

    #     omega, k, _, phi22_exp = frozen_turbulence(datas_u2, ind, zp, nt, split_time, dt, n1, dx=dx, ch="spectra")
    #     k /= Ux
        
    #     del datas_u2
                
    #     von_karman_plot(fig_vanK, col, row, k[1:tick_w], k[1:tick_w]*phi22_exp[1:tick_w], name = '$\phi_{22}exp$', color = 'midnightblue', symbols='x')
        
    #     int_exp2 = integrate.simpson(phi22_exp[1:], k[1:])
            
    #     del phi22_exp
    #     del k

        
    #     ## u3 ##
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
    #     nt = nt - 1
    #     datas_u3 = var[1:,:,:]
    #     Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
    #     datas_u3 = datas_u3 - Uz

    #     omega, k, _, phi33_exp = frozen_turbulence(datas_u3, ind, zp, nt, split_time, dt, n1, dx=dx, ch="spectra")
    #     k /= Ux
    #     print('k:', k)
        
    #     del datas_u3
        
    #     von_karman_plot(fig_vanK, col, row, k[1:tick_w], k[1:tick_w]*phi33_exp[1:tick_w], name = '$\phi_{33}exp$', color = 'darkgreen', symbols='x')
            
    #     del phi33_exp
    #     #del k
    
    #     ## Theoric VK ##
    #     print(f"Reading input files (u1 velocity) for {YELLOW}z={zp[ind]:.2f}{RESET} ...\n")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
    #     nt = nt - 1
    #     datas_u1 = var[1:,:,:]
    #     del var
    #     print('datas_u1.shape:',datas_u1.shape)
        
    #     Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
    #     print('Ux:', Ux)

    #     sigma_u_squared = np.mean(np.mean(np.mean((datas_u1[:,:,:]-Ux)**2, axis=0), axis=-1))
        
            
    #     del datas_u1
    #     print('sigma_u_squared:', sigma_u_squared)
        
    #     print("========================================")
    #     print(f"Reading input files (u2 velocity) for {YELLOW}z={zp[ind]:.2f}{RESET} ...\n")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
    #     nt = nt - 1
    #     datas_u2 = var[1:,:,:]
    #     del var
    #     print('datas_u1.shape:',datas_u2.shape)
        
    #     Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
    #     print('Uy:', Uy)

    #     sigma_v_squared = np.mean(np.mean(np.mean((datas_u2[:,:,:]-Uy)**2, axis=0), axis=-1))
        
            
    #     del datas_u2
    #     print('sigma_v_squared:', sigma_v_squared)
        
    #     print("========================================")
    #     print(f"Reading input files (u3 velocity) for {YELLOW}z={zp[ind]:.2f}{RESET} ...\n")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
    #     nt = nt - 1
    #     datas_u3 = var[1:,:,:]
    #     del var
    #     print('datas_u1.shape:',datas_u3.shape)
        
    #     Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
    #     print('Uz:', Uz)

    #     sigma_w_squared = np.mean(np.mean(np.mean((datas_u3[:,:,:]-Uz)**2, axis=0), axis=-1))
        
    #     kt = 1./2 * (sigma_u_squared + sigma_v_squared + sigma_w_squared)
            
    #     del datas_u3
    #     print('sigma_w_squared:', sigma_w_squared)
        
        
    #     print("========================================")
    #     print(f"Computing sigma1c for {YELLOW}z={zp[ind]:.2f}{RESET} ...\n")
        
    #     sigma11 = sigma_1c(sigma_u_squared, 1., sigma_v_squared, 0., sigma_w_squared, 0.)
    #     sigma12 = sigma_1c(sigma_u_squared, 0., sigma_v_squared, 1., sigma_w_squared, 0.)
        
    #     print("========================================")
    #     print(f"Computing L for {YELLOW}z={zp[ind]:.2f}{RESET} ...\n")
        
    #     dis = ko2_omega.shape[0]
    #     Le = L(0.519, ko2_tke[dis//2:], ko2_omega[dis//2:]) #list length omega
    #     #Le = L(0.519, kt, k*Ux)
        
    #     kc, phi11 = phi_11(ko2_omega[dis//2:], 1.0, Ux, sigma11, Le)
    #     kc, phi22 = phi_22(ko2_omega[dis//2:], 1.0, Ux, sigma12, Le) #=R33
    #     # kc, phi11 = phi_11(k, 1.0, Ux, sigma11, Le)
    #     # kc, phi22 = phi_22(k, 1.0, Ux, sigma12, Le)

        
    #     lim = 400
    #     for i in range(kc.shape[0]):
    #         if kc[i] > lim:
    #             tick_w = i
    #             break
            
    #     print('tick_w',tick_w)
        
    #     # ratio11 = sigma11**2*Le
    #     # ratio12 = sigma12**2*Le
    #     # ratio21 = (3 + 8*(1*kc*Le)**2) / (1 + (1*kc * Le)**2)**(11/6.)
    #     # ratio22 = 1. / (1 + (1*kc*Le)**2)**(5/6.)
        
    #     # print('sum21:', sum(ratio21))
    #     # print('sum22:', sum(ratio22))
        
    #     # ratio21 = integrate.simps(ratio21, ko2_omega)
    #     # ratio22 = integrate.simps(ratio22, ko2_omega)
        
    #     int_theo1 = integrate.simpson(phi11[:tick_w], kc[:tick_w])
    #     int_theo2 = integrate.simpson(phi22[:tick_w], kc[:tick_w])
        
    #     # int_theo1 = sum(phi11)
    #     # int_theo2 = sum(phi22)
        
    #     # print('ratio11:', ratio11)
    #     # print('ratio12:', ratio12)
    #     # print('ratio21:', ratio21)
    #     # print('ratio22:', ratio22)
    #     print('int_theo1:', int_theo1)
    #     print('int_theo2:', int_theo2)
        
    #     new_phi11 = kc[:tick_w]*phi11[:tick_w] * int_exp1**2 / int_theo1**2
    #     new_phi22 = kc[:tick_w]*phi22[:tick_w] * int_exp2**2 / int_theo2**2
        
    #     von_karman_plot(fig_vanK, col, row, kc[:tick_w], new_phi11, name = '$\phi^{0}_{11}$', color = 'firebrick', symbols='circle-open')
    #     von_karman_plot(fig_vanK, col, row, kc[:tick_w], new_phi22, name = '$\phi^{0}_{22}$', color = 'midnightblue', symbols='circle-open')
        
    #     del kc
    #     del phi11
    #     del phi22
        
        
    #     col +=1
    #     if zplan == 4:
    #         row +=1
    #         col = 1
            
            
    # if chplot == 'normal':
    #     fig_vanK.update_layout(height=600, width=900, title=f"Von Karman Spectra", font=font, showlegend=True, legend=dict(yanchor="top", xanchor="right"))
    #     save_figures(fig_vanK, "von_karman/von_karman_spectra_.png")
    # if chplot == 'all':
    #     fig_vanK.update_layout(height=900, width=900, title=f"Von Karman Spectra", font=font, showlegend=True, legend=dict(yanchor="top", xanchor="right"))
    #     save_figures(fig_vanK, "von_karman/von_karman_spectra_all.png")
        
    # elapsed_time = time.time() - start_time
    # minutes, seconds = divmod(elapsed_time, 60)
    
    # print(f'\n Von Karman study done in : {int(minutes)}m {seconds:.2f}s \n')
        
        
    del fpars_files_streamwise_u1
    del fpars_files_streamwise_u2
    del fpars_files_streamwise_u3
    del fpars_files_spanwise_u1
    del fpars_files_spanwise_u2
    del fpars_files_spanwise_u3
    del fpars_files_normal_u1
    del fpars_files_normal_u2
    del fpars_files_normal_u3
    
    elapsed_time_all = time.time() - start_time_all
    minutes, seconds = divmod(elapsed_time_all, 60)
    
    print(f'\n The whole computation time is : {int(minutes)}m {seconds:.2f}s \n')
        
############################################################
######################## END ###############################
############################################################
    
if __name__ == "__main__":
    main()