import numpy as np
import time
import glob
from tqdm import tqdm


from read_fpars import *
from tools import *
from engine import *
from von_karman import*
from parameters import *
from plot_figures import *
from read_moser import read_moser

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
    
    print(" Starting the code ...\n")
    
    start_time_all = time.time()
    
    cflow = ChannelFlow(xlen=xlen, ylen=ylen, zlen=zlen, rho=rho, uinf=uinf, re=re, ret=ret)
    
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
    
    print('nu:', cflow.nu)
    print('mu:', cflow.mu)
    print('ut:', cflow.ut)

############################################################

    ##=======================================================
            ##### FROZEN TURBULENCE #####
    ##=======================================================
    print("========================================")
    print("\nReading input files ...")
    
    _, x1, x2, _, nt, n1, n2, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane(fpars_files_streamwise_u1[0])
    nt = nt - 1
    #t = np.linspace(0, nt*dt, nt)
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
        
        #### Autocorrelation 2D ####
        
        Dt, Dx, R2d, coef = frozen_turbulence(datas_u1, ind, zp, nt, split_time, dt, n1, ch = "corr2d")
        
        frozen_turbulence_plot(fig3u1, col, row, Dt = Dt, Dx = Dx, R2d=R2d, coef=coef, ch = "corr2d")
        
        if split_time == 'Y':
            if chplot == 'all':
                save_datas([Dt, Dx, R2d[0], R2d[1]], ['dt', 'dx', 'Rdt', 'Rdx'], f'split_time/frozen_turbulence/correlation2D/u1_z{zp[ind]}.dat', '2D correlation ellipsies')
        if split_time == 'n':
            if chplot == 'all':
                save_datas([Dt, Dx, R2d[0], R2d[1]], ['dt', 'dx', 'Rdt', 'Rdx'], f'whole_time/frozen_turbulence/correlation2D/u1_z{zp[ind]}.dat', '2D correlation ellipsies')
        
        U_ratio.append((1./coef[0])/U1)
        X_ratio.append(zp[ind])
        
        U1_list.append(U1)
        Uc_list.append(1./coef[0])
        
        
        print('Uc:', Uc_list[cpt])
        
        del Dt
        del Dx
        del R2d
        
        #### Autocorrelation ####
        # ind1, ind2, Dt, Dx, R_time, R_space = frozen_turbulence(datas_u1, ind, zp, nt, split_time, dt, n1, ch="corr", dx=dx, Uc=U1)

        # frozen_turbulence_plot(fig2u1, col, row, Uc=Uc_list[cpt], R_time = R_time[:ind1], R_space = R_space[:ind2], Dt = Dt[:ind1], Dx = Dx[:ind2], ch = "corr")
        # # frozen_turbulence_plot(fig2u1, col, row, Uc=Uc_list[cpt], R_time = R_time[:], R_space = R_space[:], Dt = Dt[:], Dx = Dx[:], ch = "corr")
        
        # if split_time == 'Y':
        #     if chplot == 'all':
        #         save_datas([Dt[:ind1], Dx[:ind2], R_time[:ind1], R_space[:ind2]], ['dt', 'dx', 'Rtime', 'Rspace'], f'split_time/frozen_turbulence/correlation_st/u1_z{zp[ind]}.dat', 'Comparison space/time correlation')
        # if split_time == 'n':
        #     if chplot == 'all':
        #         save_datas([Dt[:ind1], Dx[:ind2], R_time[:ind1], R_space[:ind2]], ['dt', 'dx', 'Rtime', 'Rspace'], f'whole_time/frozen_turbulence/correlation_st/u1_z{zp[ind]}.dat', 'Comparison space/time correlation')
        
        # del Dt
        # del Dx
        # del R_time
        # del R_space
    
        #### Spectra ####
        omega, k, time_spectra, space_spectra = frozen_turbulence(datas_u1, ind, zp, nt, split_time, dt, n1, dx=dx, ch="spectra", scaling='density')
        
        frozen_turbulence_plot(fig1u1, col, row, omega = omega, Uc = Uc_list[cpt], time_spectra = time_spectra, k = k, space_spectra = space_spectra, ch = "spectra")
        
        if split_time == 'Y':
            if chplot == 'all':
                save_datas([omega[1:]/Uc_list[cpt], k[1:], time_spectra[1:], space_spectra[1:]/Uc_list[cpt]], ['omega', 'kx', 'time spectra', 'space spectra'], f'split_time/frozen_turbulence/power_spectra/u1_z{zp[ind]}.dat', 'Power spectra space/time comparison')
        if split_time == 'n':
            if chplot == 'all':
                save_datas([omega[1:]/Uc_list[cpt], k[1:], time_spectra[1:], space_spectra[1:]/Uc_list[cpt]], ['omega', 'kx', 'time spectra', 'space spectra'], f'whole_time/frozen_turbulence/power_spectra/u1_z{zp[ind]}.dat', 'Power spectra space/time comparison')
        
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
    
        
        #### Autocorrelation 2D ####
        
        Dt, Dx, R2d, coef = frozen_turbulence(datas_u2, ind, zp, nt, split_time, dt, n1, ch = "corr2d")
        
        frozen_turbulence_plot(fig3u2, col, row, Dt = Dt, Dx = Dx, R2d=R2d, coef=coef, ch = "corr2d")
        
        if split_time == 'Y':
            if chplot == 'all':
                save_datas([Dt, Dx, R2d[0], R2d[1]], ['dt', 'dx', 'Rdt', 'Rdx'], f'split_time/frozen_turbulence/correlation2D/u2_z{zp[ind]}.dat', '2D correlation ellipsies')
        if split_time == 'n':
            if chplot == 'all':
                save_datas([Dt, Dx, R2d[0], R2d[1]], ['dt', 'dx', 'Rdt', 'Rdx'], f'whole_time/frozen_turbulence/correlation2D/u2_z{zp[ind]}.dat', '2D correlation ellipsies')
        
        del Dt
        del Dx
        del R2d
        
        #### Autocorrelation ####
        # ind1, ind2, Dt, Dx, R_time, R_space = frozen_turbulence(datas_u2, ind, zp, nt, split_time, dt, n1, ch="corr", dx=dx, Uc=U1)

        # frozen_turbulence_plot(fig2u2, col, row, Uc=Uc_list[cpt], R_time = R_time[:ind1], R_space = R_space[:ind2], Dt = Dt[:ind1], Dx = Dx[:ind2], ch = "corr")
        # # frozen_turbulence_plot(fig2u2, col, row, Uc=Uc_list[cpt], R_time = R_time[:], R_space = R_space[:], Dt = Dt[:], Dx = Dx[:], ch = "corr")
        
        # if split_time == 'Y':
        #     if chplot == 'all':
        #         save_datas([Dt[:ind1], Dx[:ind2], R_time[:ind1], R_space[:ind2]], ['dt', 'dx', 'Rtime', 'Rspace'], f'split_time/frozen_turbulence/correlation_st/u2_z{zp[ind]}.dat', 'Comparison space/time correlation')
        # if split_time == 'n':
        #     if chplot == 'all':
        #         save_datas([Dt[:ind1], Dx[:ind2], R_time[:ind1], R_space[:ind2]], ['dt', 'dx', 'Rtime', 'Rspace'], f'whole_time/frozen_turbulence/correlation_st/u2_z{zp[ind]}.dat', 'Comparison space/time correlation')
        
        
        # del Dt
        # del Dx
        # del R_time
        # del R_space
    
        #### Spectra ####
        omega, k, time_spectra, space_spectra = frozen_turbulence(datas_u2, ind, zp, nt, split_time, dt, n1, dx=dx, ch="spectra", scaling='density')
        
        frozen_turbulence_plot(fig1u2, col, row, omega = omega, Uc = Uc_list[cpt], time_spectra = time_spectra, k = k, space_spectra = space_spectra, ch = "spectra")
        
        if split_time == 'Y':
            if chplot == 'all':
                save_datas([omega[1:]/Uc_list[cpt], k[1:], time_spectra[1:], space_spectra[1:]/Uc_list[cpt]], ['omega', 'kx', 'time spectra', 'space spectra'], f'split_time/frozen_turbulence/power_spectra/u2_z{zp[ind]}.dat', 'Power spectra space/time comparison')
        if split_time == 'n':
            if chplot == 'all':
                save_datas([omega[1:]/Uc_list[cpt], k[1:], time_spectra[1:], space_spectra[1:]/Uc_list[cpt]], ['omega', 'kx', 'time spectra', 'space spectra'], f'whole_time/frozen_turbulence/power_spectra/u2_z{zp[ind]}.dat', 'Power spectra space/time comparison')
        
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
        
        #### Autocorrelation 2D ####
        
        Dt, Dx, R2d, coef = frozen_turbulence(datas_u3, ind, zp, nt, split_time, dt, n1, ch = "corr2d")
        
        frozen_turbulence_plot(fig3u3, col, row, Dt = Dt, Dx = Dx, R2d=R2d, coef=coef, ch = "corr2d")
        
        if split_time == 'Y':
            if chplot == 'all':
                save_datas([Dt, Dx, R2d[0], R2d[1]], ['dt', 'dx', 'Rdt', 'Rdx'], f'split_time/frozen_turbulence/correlation2D/u3_z{zp[ind]}.dat', '2D correlation ellipsies')
        if split_time == 'n':
            if chplot == 'all':
                save_datas([Dt, Dx, R2d[0], R2d[1]], ['dt', 'dx', 'Rdt', 'Rdx'], f'whole_time/frozen_turbulence/correlation2D/u3_z{zp[ind]}.dat', '2D correlation ellipsies')
        
        del Dt
        del Dx
        del R2d
        
        
        #### Autocorrelation ####
        # ind1, ind2, Dt, Dx, R_time, R_space = frozen_turbulence(datas_u3, ind, zp, nt, split_time, dt, n1, ch="corr", dx=dx, Uc=U1)

        # frozen_turbulence_plot(fig2u3, col, row, Uc=Uc_list[cpt], R_time = R_time[:ind1], R_space = R_space[:ind2], Dt = Dt[:ind1], Dx = Dx[:ind2], ch = "corr")
        # # frozen_turbulence_plot(fig2u3, col, row, Uc=Uc_list[cpt], R_time = R_time[:], R_space = R_space[:], Dt = Dt[:], Dx = Dx[:], ch = "corr")
        
        # if split_time == 'Y':
        #     if chplot == 'all':
        #         save_datas([Dt[:ind1], Dx[:ind2], R_time[:ind1], R_space[:ind2]], ['dt', 'dx', 'Rtime', 'Rspace'], f'split_time/frozen_turbulence/correlation_st/u3_z{zp[ind]}.dat', 'Comparison space/time correlation')
        # if split_time == 'n':
        #     if chplot == 'all':
        #         save_datas([Dt[:ind1], Dx[:ind2], R_time[:ind1], R_space[:ind2]], ['dt', 'dx', 'Rtime', 'Rspace'], f'whole_time/frozen_turbulence/correlation_st/u3_z{zp[ind]}.dat', 'Comparison space/time correlation')
        
        
        # del Dt
        # del Dx
        # del R_time
        # del R_space
    
    
        #### Spectra ####
        omega, k, time_spectra, space_spectra = frozen_turbulence(datas_u3, ind, zp, nt, split_time, dt, n1, dx=dx, ch="spectra", scaling='density')
        
        frozen_turbulence_plot(fig1u3, col, row, omega = omega, Uc = Uc_list[cpt], time_spectra = time_spectra, k = k, space_spectra = space_spectra, ch = "spectra")
        
        if split_time == 'Y':
            if chplot == 'all':
                save_datas([omega[1:]/Uc_list[cpt], k[1:], time_spectra[1:], space_spectra[1:]/Uc_list[cpt]], ['omega', 'kx', 'time spectra', 'space spectra'], f'split_time/frozen_turbulence/power_spectra/u3_z{zp[ind]}.dat', 'Power spectra space/time comparison')
        if split_time == 'n':
            if chplot == 'all':
                save_datas([omega[1:]/Uc_list[cpt], k[1:], time_spectra[1:], space_spectra[1:]/Uc_list[cpt]], ['omega', 'kx', 'time spectra', 'space spectra'], f'whole_time/frozen_turbulence/power_spectra/u3_z{zp[ind]}.dat', 'Power spectra space/time comparison')
        
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
            
    X_ratio = np.array(X_ratio)
    U1_list = np.array(U1_list)
    Uc_list = np.array(Uc_list)
    U_ratio = np.array(U_ratio)
            
    ## save data ##
    if split_time == 'Y':
        if chplot == 'all':
            save_datas([X_ratio, U1_list, Uc_list, U_ratio], ['X', 'U1', 'Uc', 'U_ratio'], 'split_time/frozen_turbulence/correlation2D/U_compare.dat', 'Vellocity comparison')
    if split_time == 'n':
        if chplot == 'all':
            save_datas([X_ratio, U1_list, Uc_list, U_ratio], ['X', 'U1', 'Uc', 'U_ratio'], 'whole_time/frozen_turbulence/correlation2D/U_compare.dat', 'Vellocity comparison')
            
        
    # Update layout properties for 4 plots
    if chplot == "normal":
        fig1u1.update_layout(height=600, width=1100, title_text="Power spectra streawise u1", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.04, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
        fig2u1.update_layout(height=600, width=1100, title_text='Autocorrelation comparison streamwise u1', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
        fig3u1.update_layout(height=600, width=1100, title_text='Correlation 2D streamwise u1', legend=dict(y=1.2, x=0.9), font=font)
        
        fig1u2.update_layout(height=600, width=1100, title_text="Power spectra streawise u2", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.04, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
        fig2u2.update_layout(height=600, width=1100, title_text='Autocorrelation comparison streamwise u2', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
        fig3u2.update_layout(height=600, width=1100, title_text='Correlation 2D streamwise u2', legend=dict(y=1.2, x=0.9), font=font)
        
        fig1u3.update_layout(height=600, width=1100, title_text="Power spectra streawise u3", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
        fig2u3.update_layout(height=600, width=1100, title_text='Autocorrelation comparison streamwise u3', font=font,  legend=dict(yanchor="bottom", y=1.01, xanchor="left", x=0.93))
        fig3u3.update_layout(height=600, width=1100, title_text='Correlation 2D streamwise u3', legend=dict(y=1.2, x=0.9), font=font)
        
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
            # save_figures(fig2u1, "split_time/frozen_turbulence/correlation_st/u1.png")
            save_figures(fig3u1, "split_time/frozen_turbulence/correlation2D/u1.png")
            
            save_figures(fig1u2, "split_time/frozen_turbulence/power_spectra/u2.png")
            # save_figures(fig2u2, "split_time/frozen_turbulence/correlation_st/u2.png")
            save_figures(fig3u2, "split_time/frozen_turbulence/correlation2D/u2.png")
        
            save_figures(fig1u3, "split_time/frozen_turbulence/power_spectra/u3.png")
            # save_figures(fig2u3, "split_time/frozen_turbulence/correlation_st/u3.png")
            save_figures(fig3u3, "split_time/frozen_turbulence/correlation2D/u3.png")
            
            save_figures(figU, "split_time/frozen_turbulence/correlation2D/u_ratio.png")
            
            save_figures(figU1c, "split_time/frozen_turbulence/correlation2D/u_1c.png")
            
        if split_time == 'n':
            save_figures(fig1u1, "whole_time/frozen_turbulence/power_spectra/u1.png")
            save_figures(fig2u1, "whole_time/frozen_turbulence/correlation_st/u1.png")
            save_figures(fig3u1, "whole_time/frozen_turbulence/correlation2D/u1.png")
            
            save_figures(fig1u2, "whole_time/frozen_turbulence/power_spectra/u2.png")
            save_figures(fig2u2, "whole_time/frozen_turbulence/correlation_st/u2.png")
            save_figures(fig3u2, "whole_time/frozen_turbulence/correlation2D/u2.png")
        
            save_figures(fig1u3, "whole_time/frozen_turbulence/power_spectra/u3.png")
            save_figures(fig2u3, "whole_time/frozen_turbulence/correlation_st/u3.png")
            save_figures(fig3u3, "whole_time/frozen_turbulence/correlation2D/u3.png")
            
            save_figures(figU, "whole_time/frozen_turbulence/correlation2D/u_ratio.png")
            
            save_figures(figU1c, "whole_time/frozen_turbulence/correlation2D/u_1c.png")
        
    # # Update layout properties for full plots.
    if chplot == "all":
        fig1u1.update_layout(height=900, width=900, title_text="Power spectra streawise u1", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
        fig2u1.update_layout(height=900, width=1100, title_text='Autocorrelation comparison Streamwise u1', font=font,  legend=dict(yanchor="top", xanchor="right"))
        fig3u1.update_layout(height=900, width=1100, title_text='Correlation 2D Streamwise u1', legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40), font=font)
        
        fig1u2.update_layout(height=900, width=900, title_text="Power spectra Streawise u2", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
        fig2u2.update_layout(height=900, width=1100, title_text='Autocorrelation comparison Streamwise u2', font=font,  legend=dict(yanchor="top", xanchor="right"))
        fig3u2.update_layout(height=900, width=1100, title_text='Correlation 2D Streamwise u2', legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40), font=font)

        fig1u3.update_layout(height=900, width=900, title_text="Power spectra Streawise u3", font=font, legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40))
        fig2u3.update_layout(height=900, width=1100, title_text='Autocorrelation comparison Streamwise u3', font=font,  legend=dict(yanchor="top", xanchor="right"))
        fig3u3.update_layout(height=900, width=1100, title_text='Correlation 2D Streamwise u3', legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="right", x=1), margin=dict(l=40, r=40, t=60, b=40), font=font)
        
        figU.add_trace(go.Scatter(x=U_ratio, y=X_ratio, line=dict(color='midnightblue')))
        figU.update_xaxes(title='$U_c/U_1$', dtick=0.1)
        figU.update_yaxes(title='$z^+$')
        figU.update_layout(title='Velocity ratio', font=font, showlegend=False)
        
        figU1c.add_trace(go.Scatter(x=U1_list, y=X_ratio, name='$U_1$', line=dict(color='midnightblue')))
        figU1c.add_trace(go.Scatter(x=Uc_list, y=X_ratio, name='$U_c$', line=dict(color='firebrick')))
        figU1c.update_xaxes(title='$\\text{mean velocity}~(m.s^{-1})$')
        figU1c.update_yaxes(title='$z^+$')
        figU1c.update_layout(title='Velocity comparison', font=font, legend=dict(yanchor="bottom", xanchor="right"))
        
        if split_time == 'Y':
            save_figures(fig1u1, "split_time/frozen_turbulence/power_spectra/u1_all.png")
            save_figures(fig2u1, "split_time/frozen_turbulence/correlation_st/u1_all.png")
            save_figures(fig3u1, "split_time/frozen_turbulence/correlation2D/u1_all.png")
            
            save_figures(fig1u2, "split_time/frozen_turbulence/power_spectra/u2_all.png")
            save_figures(fig2u2, "split_time/frozen_turbulence/correlation_st/u2_all.png")
            save_figures(fig3u2, "split_time/frozen_turbulence/correlation2D/u2_all.png")
            
            save_figures(fig1u3, "split_time/frozen_turbulence/power_spectra/u3_all.png")
            save_figures(fig2u3, "split_time/frozen_turbulence/correlation_st/u3_all.png")
            save_figures(fig3u3, "split_time/frozen_turbulence/correlation2D/u3_all.png")
            
            save_figures(figU, "split_time/frozen_turbulence/correlation2D/u_ratio_all.png")
            
            save_figures(figU1c, "split_time/frozen_turbulence/correlation2D/u_1c_all.png")
            
        if split_time == 'n':
            save_figures(fig1u1, "whole_time/frozen_turbulence/power_spectra/u1_all.png")
            save_figures(fig2u1, "whole_time/frozen_turbulence/correlation_st/u1_all.png")
            save_figures(fig3u1, "whole_time/frozen_turbulence/correlation2D/u1_all.png")
            
            save_figures(fig1u2, "whole_time/frozen_turbulence/power_spectra/u2_all.png")
            save_figures(fig2u2, "whole_time/frozen_turbulence/correlation_st/u2_all.png")
            save_figures(fig3u2, "whole_time/frozen_turbulence/correlation2D/u2_all.png")
        
            save_figures(fig1u3, "whole_time/frozen_turbulence/power_spectra/u3_all.png")
            save_figures(fig2u3, "whole_time/frozen_turbulence/correlation_st/u3_all.png")
            save_figures(fig3u3, "whole_time/frozen_turbulence/correlation2D/u3_all.png")
            
            save_figures(figU, "whole_time/frozen_turbulence/correlation2D/u_ratio_all.png")
            
            save_figures(figU1c, "whole_time/frozen_turbulence/correlation2D/u_1c_all.png")
            
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    print(f'\n Frozen turbulence valiation done in : {int(minutes)}m {seconds:.2f}s \n')
    
    
    
    ##=======================================================
                    #### GAMMA ####
    ##=======================================================
    
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
    # num_split_t = nt // split_t
    
    # figu1g, figu2g, figu3g = init_figures_gamma(zp, ch=chplot) #Figure initialization
    # figu1g_r, figu2g_r, figu3g_r = init_figures_gamma(zp, ch=chplot)
    
    # figgamma = make_subplots(rows=1, cols=3, shared_yaxes=True, y_title='$\gamma$', subplot_titles=("$u_1$", "$u_2$", "$u_3$"))
    # gamma_u1_1 = []
    # gamma_u2_1 = []
    # gamma_u3_1 = []
    # gamma_u1_2 = []
    # gamma_u2_2 = []
    # gamma_u3_2 = []
    # gamma_u1_3 = []
    # gamma_u2_3 = []
    # gamma_u3_3 = []
    # gamma_u1_4 = []
    # gamma_u2_4 = []
    # gamma_u3_4 = []
    
    # figgamma_r = make_subplots(rows=1, cols=3, shared_yaxes=True, y_title='$\gamma$', subplot_titles=("$u_1$", "$u_2$", "$u_3$"))
    # gamma_u1_1_r = []
    # gamma_u2_1_r = []
    # gamma_u3_1_r = []
    # gamma_u1_2_r = []
    # gamma_u2_2_r = []
    # gamma_u3_2_r = []
    # gamma_u1_3_r = []
    # gamma_u2_3_r = []
    # gamma_u3_3_r = []
    # gamma_u1_4_r = []
    # gamma_u2_4_r = []
    # gamma_u3_4_r = []
    
    # z_axis = zp
    # cpt = 0
    
    # for ind, zplan in enumerate(zp_ind):
    
        
    #     print("========================================")
    #     print(f'Gamma determination for {YELLOW}zp={zp[ind]:.2f}{RESET}')
    #     print('Plan number:', zplan)
        
    #     #### u1 ####
    #     print("\nReading input files u1 streamwise...")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
    #     datas_u1 = var[1:,:,:]
        
    #     U1 = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=0), axis=-1))
        
    #     datas_u1 = datas_u1 - U1
        
    #     _, omega, _ = Gamma_function(datas_u1[0:split_t,:,:], dt, Uc_list[cpt], geom = "plan", axis = "streamwise")
        
    #     funct = np.zeros((omega.shape[0], n1))
    #     omega = np.zeros((omega.shape[0]))
    #     Dx = np.zeros((n1))
    #     for n in tqdm(range(1,num_split_t), desc=f'Gamma', colour= 'GREEN'):
            
    #         funct += Gamma_function(datas_u1[(n-1)*split_t:n*split_t,:,:], dt, Uc_list[cpt], geom = "plan", axis = "streamwise")[0]
    #         omega += Gamma_function(datas_u1[(n-1)*split_t:n*split_t,:,:], dt, Uc_list[cpt], geom = "plan", axis = "streamwise")[1]
    #         Dx += Gamma_function(datas_u1[(n-1)*split_t:n*split_t,:,:], dt, Uc_list[cpt], geom = "plan", axis = "streamwise")[2]
            
    #     funct /= (num_split_t-1)
    #     omega /= (num_split_t-1)
    #     Dx /= (num_split_t-1)
        
    #     del datas_u1
        

    #     for i in range(omega.shape[0]):
    #         if omega[i]>g_limit0:
    #             ind0 = i
    #             break
            
    #     for i in range(omega.shape[0]):
    #         if omega[i]>g_limit1:
    #             ind1 = i
    #             break
            
    #     slop1 = gamma_plot(figu1g, col, row, funct, omega, Dx, ind0, ind1, ch = 'w')
        
    #     if split_time == 'Y':
    #         if chplot == 'all':
    #             save_datas([omega[ind0:ind1], funct[ind0:ind1,1], funct[ind0:ind1,5], funct[ind0:ind1,10], funct[ind0:ind1,15]], ['omega', 'fct(r=1)', 'fct(r=5)', 'fct(r=10)', 'fct(r=15)'], f'split_time/gamma/u1_w_{zp[ind]}.dat', 'Gamma in function of w (u1)')
    #     if split_time == 'n':
    #         if chplot == 'all':
    #             save_datas([omega[ind0:ind1], funct[ind0:ind1,1], funct[ind0:ind1,5], funct[ind0:ind1,10], funct[ind0:ind1,15]], ['omega', 'fct(r=1)', 'fct(r=5)', 'fct(r=10)', 'fct(r=15)'], f'whole_time/gamma/u1_w_{zp[ind]}.dat', 'Gamma in function of w (u1)')
            
    #     ind0 = Dx.shape[0]//2
            
    #     r_lim = 4
    #     for i in range(Dx.shape[0]):
    #         if Dx[i]>r_lim:
    #             ind1 = i
    #             break
            
    #     print("ind0",ind0) ; print("ind1",ind1)
        
    #     moy1, moy2, moy3, moy4 = gamma_plot(figu1g_r, col, row, funct, omega, Dx, ind0, ind1, ch = 'x')
        
    #     if split_time == 'Y':
    #         if chplot == 'all':
    #             save_datas([Dx[ind0:ind1]-xlen//2, funct[1,ind0:ind1], funct[2,ind0:ind1], funct[3,ind0:ind1], funct[4,ind0:ind1]], ['dx', 'fct(w=1)', 'fct(w=2)', 'fct(w=3)', 'fct(w=4)'], f'split_time/gamma/u1_r_{zp[ind]}.dat', 'Gamma in function of r (u1)')
    #     if split_time == 'n':
    #         if chplot == 'all':
    #             save_datas([Dx[ind0:ind1]-xlen//2, funct[1,ind0:ind1], funct[2,ind0:ind1], funct[3,ind0:ind1], funct[4,ind0:ind1]], ['dx', 'fct(w=1)', 'fct(w=2)', 'fct(w=3)', 'fct(w=4)'], f'whole_time/gamma/u1_r_{zp[ind]}.dat', 'Gamma in function of r (u1)')
        
        
    #     gamma_u1_1.append(slop1[0]*Uc_list[cpt]/Dx[1])
    #     gamma_u1_2.append(slop1[0]*Uc_list[cpt]/Dx[5])
    #     gamma_u1_3.append(slop1[0]*Uc_list[cpt]/Dx[10])
    #     gamma_u1_4.append(slop1[0]*Uc_list[cpt]/Dx[15])
        
    #     gamma_u1_1_r.append(moy1*Uc_list[cpt]/(omega[1]))
    #     gamma_u1_2_r.append(moy2*Uc_list[cpt]/(omega[2]))
    #     gamma_u1_3_r.append(moy3*Uc_list[cpt]/(omega[3]))
    #     gamma_u1_4_r.append(moy4*Uc_list[cpt]/(omega[4]))
        
            
    #     del omega
    #     del funct
    #     del Dx
    #     del slop1
    #     # del slop2
        
    #     #### u2 ####
    #     print("\nReading input files u2 streamwise ...")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
    #     datas_u2 = var[1:,:,:]
        
    #     U2 = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=0), axis=-1))
        
    #     datas_u2 = datas_u2 - U2
        
        
    #     _, omega, _ = Gamma_function(datas_u2[0:split_t,:,:], dt, Uc_list[cpt], geom = "plan", axis = "streamwise")
        
    #     funct = np.zeros((omega.shape[0], n1))
    #     omega = np.zeros((omega.shape[0]))
    #     Dx = np.zeros((n1))
    #     for n in tqdm(range(1,num_split_t), desc=f'Gamma', colour= 'GREEN'):
            
    #         funct += Gamma_function(datas_u2[(n-1)*split_t:n*split_t,:,:], dt, Uc_list[cpt], geom = "plan", axis = "streamwise")[0]
    #         omega += Gamma_function(datas_u2[(n-1)*split_t:n*split_t,:,:], dt, Uc_list[cpt], geom = "plan", axis = "streamwise")[1]
    #         Dx += Gamma_function(datas_u2[(n-1)*split_t:n*split_t,:,:], dt, Uc_list[cpt], geom = "plan", axis = "streamwise")[2]
            
    #     funct /= (num_split_t-1)
    #     omega /= (num_split_t-1)
    #     Dx /= (num_split_t-1)
        
    #     del datas_u2
        
    #     for i in range(omega.shape[0]):
    #         if omega[i]>g_limit0:
    #             ind0 = i
    #             break
            
    #     for i in range(omega.shape[0]):
    #         if omega[i]>g_limit1:
    #             ind1 = i
    #             break
                
    #     slop1 = gamma_plot(figu2g, col, row, funct, omega, Dx, ind0, ind1, ch = 'w')
        
    #     if split_time == 'Y':
    #         if chplot == 'all':
    #             save_datas([omega[ind0:ind1], funct[ind0:ind1,1], funct[ind0:ind1,5], funct[ind0:ind1,10], funct[ind0:ind1,15]], ['omega', 'fct(r=1)', 'fct(r=5)', 'fct(r=10)', 'fct(r=15)'], f'split_time/gamma/u2_w_{zp[ind]}.dat', 'Gamma in function of w (u2)')
    #     if split_time == 'n':
    #         if chplot == 'all':
    #             save_datas([omega[ind0:ind1], funct[ind0:ind1,1], funct[ind0:ind1,5], funct[ind0:ind1,10], funct[ind0:ind1,15]], ['omega', 'fct(r=1)', 'fct(r=5)', 'fct(r=10)', 'fct(r=15)'], f'whole_time/gamma/u2_w_{zp[ind]}.dat', 'Gamma in function of w (u2)')
            
    #     ind0 = Dx.shape[0]//2
            
    #     r_lim = 4
    #     for i in range(Dx.shape[0]):
    #         if Dx[i]>r_lim:
    #             ind1 = i
    #             break
        
    #     moy1, moy2, moy3, moy4 = gamma_plot(figu2g_r, col, row, funct, omega, Dx, ind0, ind1, ch = 'x')
        
    #     if split_time == 'Y':
    #         if chplot == 'all':
    #             save_datas([Dx[ind0:ind1]-xlen//2, funct[1,ind0:ind1], funct[2,ind0:ind1], funct[3,ind0:ind1], funct[4,ind0:ind1]], ['dx', 'fct(w=1)', 'fct(w=2)', 'fct(w=3)', 'fct(w=4)'], f'split_time/gamma/u2_r_{zp[ind]}.dat', 'Gamma in function of r (u2)')
    #     if split_time == 'n':
    #         if chplot == 'all':
    #             save_datas([Dx[ind0:ind1]-xlen//2, funct[1,ind0:ind1], funct[2,ind0:ind1], funct[3,ind0:ind1], funct[4,ind0:ind1]], ['dx', 'fct(w=1)', 'fct(w=2)', 'fct(w=3)', 'fct(w=4)'], f'whole_time/gamma/u2_r_{zp[ind]}.dat', 'Gamma in function of r (u2)')
        
    #     gamma_u2_1.append(slop1[0]*Uc_list[cpt]/Dx[1])
    #     gamma_u2_2.append(slop1[0]*Uc_list[cpt]/Dx[5])
    #     gamma_u2_3.append(slop1[0]*Uc_list[cpt]/Dx[10])
    #     gamma_u2_4.append(slop1[0]*Uc_list[cpt]/Dx[15])
        
    #     gamma_u2_1_r.append(moy1*Uc_list[cpt]/(omega[1]))
    #     gamma_u2_2_r.append(moy2*Uc_list[cpt]/(omega[2]))
    #     gamma_u2_3_r.append(moy3*Uc_list[cpt]/(omega[3]))
    #     gamma_u2_4_r.append(moy4*Uc_list[cpt]/(omega[4]))
            
    #     del omega
    #     del funct
    #     del Dx
    #     del slop1
    #     # del slop2
        
    #     #### u3 ####
    #     print("\nReading input files u3 streamwise ...")
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
    #     datas_u3 = var[1:,:,:]
        
    #     U3 = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=0), axis=-1))
        
    #     datas_u3 = datas_u3 - U3
        
    #     _, omega, _ = Gamma_function(datas_u3[0:split_t,:,:], dt, U1, geom = "plan", axis = "streamwise")
        
    #     funct = np.zeros((omega.shape[0], n1))
    #     omega = np.zeros((omega.shape[0]))
    #     Dx = np.zeros((n1))
    #     for n in tqdm(range(1,num_split_t), desc=f'Gamma', colour= 'GREEN'):
            
    #         funct += Gamma_function(datas_u3[(n-1)*split_t:n*split_t,:,:], dt, Uc_list[cpt], geom = "plan", axis = "streamwise")[0] #I hould take Uc here
    #         omega += Gamma_function(datas_u3[(n-1)*split_t:n*split_t,:,:], dt, Uc_list[cpt], geom = "plan", axis = "streamwise")[1]
    #         Dx += Gamma_function(datas_u3[(n-1)*split_t:n*split_t,:,:], dt, Uc_list[cpt], geom = "plan", axis = "streamwise")[2]
            
    #     funct /= (num_split_t-1)
    #     omega /= (num_split_t-1)
    #     Dx /= (num_split_t-1)
        
    #     del datas_u3
        
    #     for i in range(omega.shape[0]):
    #         if omega[i]>g_limit0:
    #             ind0 = i
    #             break
            
    #     for i in range(omega.shape[0]):
    #         if omega[i]>g_limit1:
    #             ind1 = i
    #             break
            
    #     slop1 = gamma_plot(figu3g, col, row, funct, omega, Dx, ind0, ind1, ch = 'w')
        
    #     if split_time == 'Y':
    #         if chplot == 'all':
    #             save_datas([omega[ind0:ind1], funct[ind0:ind1,1], funct[ind0:ind1,5], funct[ind0:ind1,10], funct[ind0:ind1,15]], ['omega', 'fct(r=1)', 'fct(r=5)', 'fct(r=10)', 'fct(r=15)'], f'split_time/gamma/u3_w_{zp[ind]}.dat', 'Gamma in function of w (u3)')
    #     if split_time == 'n':
    #         if chplot == 'all':
    #             save_datas([omega[ind0:ind1], funct[ind0:ind1,1], funct[ind0:ind1,5], funct[ind0:ind1,10], funct[ind0:ind1,15]], ['omega', 'fct(r=1)', 'fct(r=5)', 'fct(r=10)', 'fct(r=15)'], f'whole_time/gamma/u3_w_{zp[ind]}.dat', 'Gamma in function of w (u3)')
            
    #     ind0 = Dx.shape[0]//2
            
    #     r_lim = 4
    #     for i in range(Dx.shape[0]):
    #         if Dx[i]>r_lim:
    #             ind1 = i
    #             break
        
    #     moy1, moy2, moy3, moy4 = gamma_plot(figu3g_r, col, row, funct, omega, Dx, ind0, ind1, ch = 'x')
        
    #     if split_time == 'Y':
    #         if chplot == 'all':
    #             save_datas([Dx[ind0:ind1]-xlen//2, funct[1,ind0:ind1], funct[2,ind0:ind1], funct[3,ind0:ind1], funct[4,ind0:ind1]], ['dx', 'fct(w=1)', 'fct(w=2)', 'fct(w=3)', 'fct(w=4)'], f'split_time/gamma/u3_r_{zp[ind]}.dat', 'Gamma in function of r (u3)')
    #     if split_time == 'n':
    #         if chplot == 'all':
    #             save_datas([Dx[ind0:ind1]-xlen//2, funct[1,ind0:ind1], funct[2,ind0:ind1], funct[3,ind0:ind1], funct[4,ind0:ind1]], ['dx', 'fct(w=1)', 'fct(w=2)', 'fct(w=3)', 'fct(w=4)'], f'whole_time/gamma/u3_r_{zp[ind]}.dat', 'Gamma in function of r (u3)')
        
    #     gamma_u3_1.append(slop1[0]*Uc_list[cpt]/Dx[1])
    #     gamma_u3_2.append(slop1[0]*Uc_list[cpt]/Dx[5])
    #     gamma_u3_3.append(slop1[0]*Uc_list[cpt]/Dx[10])
    #     gamma_u3_4.append(slop1[0]*Uc_list[cpt]/Dx[15])
        
    #     gamma_u3_1_r.append(moy1*Uc_list[cpt]/(omega[1]))
    #     gamma_u3_2_r.append(moy2*Uc_list[cpt]/(omega[2]))
    #     gamma_u3_3_r.append(moy3*Uc_list[cpt]/(omega[3]))
    #     gamma_u3_4_r.append(moy4*Uc_list[cpt]/(omega[4]))
            
    #     del funct
    #     del slop1
    #     # del slop2
        
        
    #     col +=1
    #     cpt +=1
    #     if zplan == 4:
    #         row +=1
    #         col = 1
            
            
    # # Update layout properties for 4 plots
    # if chplot == "normal":
    #     figu1g.update_layout(height=600, width=1100, title_text='Gamma determination u1', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu2g.update_layout(height=600, width=1100, title_text='Gamma determination u2', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu3g.update_layout(height=600, width=1100, title_text='Gamma determination u3', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu1g_r.update_layout(height=600, width=1100, title_text='Gamma determination u1', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu2g_r.update_layout(height=600, width=1100, title_text='Gamma determination u2', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu3g_r.update_layout(height=600, width=1100, title_text='Gamma determination u3', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     ## Gamma for omega ##
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_1, name=f'$r={Dx[1]:.2f}$', line=dict(color='midnightblue')), row=1, col=1)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_2, name=f'$r={Dx[5]:.2f}$', line=dict(color='firebrick')), row=1, col=1)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_3, name=f'$r={Dx[10]:.2f}$', line=dict(color='darkgreen')), row=1, col=1)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_4, name=f'$r={Dx[15]:.2f}$', line=dict(color='purple')), row=1, col=1)
        
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_1, showlegend=False, line=dict(color='midnightblue')), row=1, col=2)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_2, showlegend=False, line=dict(color='firebrick')), row=1, col=2)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_3, showlegend=False, line=dict(color='darkgreen')), row=1, col=2)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_4, showlegend=False, line=dict(color='purple')), row=1, col=2)
        
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_1, showlegend=False, line=dict(color='midnightblue')), row=1, col=3)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_2, showlegend=False, line=dict(color='firebrick')), row=1, col=3)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_3, showlegend=False, line=dict(color='darkgreen')), row=1, col=3)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_4, showlegend=False, line=dict(color='purple')), row=1, col=3)
        
    #     del Dx
        
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=1)
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=2)
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=3)
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=4)
    #     figgamma.update_layout(height=600, width=900, title_text='Gamma evolution', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        
    #     ## Gamme for r ##
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_1_r, showlegend=False, line=dict(color='midnightblue')), row=1, col=1)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_2_r, showlegend=False, line=dict(color='firebrick')), row=1, col=1)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_3_r, showlegend=False, line=dict(color='green')), row=1, col=1)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_4_r, showlegend=False, line=dict(color='purple')), row=1, col=1)
        
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_1_r, showlegend=False, line=dict(color='midnightblue')), row=1, col=2)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_2_r, showlegend=False, line=dict(color='firebrick')), row=1, col=2)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_3_r, showlegend=False, line=dict(color='green')), row=1, col=2)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_4_r, showlegend=False, line=dict(color='purple')), row=1, col=2)
        
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_1_r, name=f'$\omega={omega[1]:.2f}$', line=dict(color='midnightblue')), row=1, col=3)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_2_r, name=f'$\omega={omega[2]:.2f}$', line=dict(color='firebrick')), row=1, col=3)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_3_r, name=f'$\omega={omega[3]:.2f}$', line=dict(color='green')), row=1, col=3)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_4_r, name=f'$\omega={omega[4]:.2f}$', line=dict(color='purple')), row=1, col=3)
        
    #     del omega
        
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=1)
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=2)
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=3)
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=4)
    #     figgamma_r.update_layout(height=600, width=900, title_text='Gamma evolution', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
        
    #     if split_time == 'Y':

    #         save_figures(figu1g, "split_time/gamma/gamma_u1_w.png")
            
    #         save_figures(figu2g, "split_time/gamma/gamma_u2_w.png")
        
    #         save_figures(figu3g, "split_time/gamma/gamma_u3_w.png")
            
    #         save_figures(figu1g_r, "split_time/gamma/gamma_u1_r.png")
            
    #         save_figures(figu2g_r, "split_time/gamma/gamma_u2_r.png")
        
    #         save_figures(figu3g_r, "split_time/gamma/gamma_u3_r.png")
            
    #         save_figures(figgamma, "split_time/gamma/gamma_view_w.png")
            
    #         save_figures(figgamma_r, "split_time/gamma/gamma_view_r.png")
            
    #     if split_time == 'n':

    #         save_figures(figu1g, "whole_time/gamma/gamma_u1_w.png")
            
    #         save_figures(figu2g, "whole_time/gamma/gamma_u2_w.png")
        
    #         save_figures(figu3g, "whole_time/gamma/gamma_u3_w.png")
            
    #         save_figures(figu1g_r, "whole_time/gamma/gamma_u1_r.png")
            
    #         save_figures(figu2g_r, "whole_time/gamma/gamma_u2_r.png")
        
    #         save_figures(figu3g_r, "whole_time/gamma/gamma_u3_r.png")
            
    #         save_figures(figgamma, "whole_time/gamma/gamma_view_w.png")
            
    #         save_figures(figgamma_r, "whole_time/gamma/gamma_view_r.png")
        
    # # # Update layout properties for full plots.
    # if chplot == "all":

    #     figu1g.update_layout(height=900, width=1100, title_text='Gamma determination u1', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu2g.update_layout(height=900, width=1100, title_text='Gamma determination u2', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu3g.update_layout(height=900, width=1100, title_text='Gamma determination u3', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu1g_r.update_layout(height=900, width=1100, title_text='Gamma determination u1', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu2g_r.update_layout(height=900, width=1100, title_text='Gamma determination u2', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     figu3g_r.update_layout(height=900, width=1100, title_text='Gamma determination u3', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     ## Gamma for omega ##
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_1, name=f'$r={Dx[1]:.2f}$', line=dict(color='midnightblue')), row=1, col=1)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_2, name=f'$r={Dx[5]:.2f}$', line=dict(color='firebrick')), row=1, col=1)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_3, name=f'$r={Dx[10]:.2f}$', line=dict(color='darkgreen')), row=1, col=1)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u1_4, name=f'$r={Dx[15]:.2f}$', line=dict(color='purple')), row=1, col=1)
        
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_1, showlegend=False, line=dict(color='midnightblue')), row=1, col=2)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_2, showlegend=False, line=dict(color='firebrick')), row=1, col=2)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_3, showlegend=False, line=dict(color='darkgreen')), row=1, col=2)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u2_4, showlegend=False, line=dict(color='purple')), row=1, col=2)
        
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_1, showlegend=False, line=dict(color='midnightblue')), row=1, col=3)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_2, showlegend=False, line=dict(color='firebrick')), row=1, col=3)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_3, showlegend=False, line=dict(color='darkgreen')), row=1, col=3)
    #     figgamma.add_trace(go.Scatter(x=z_axis, y=gamma_u3_4, showlegend=False, line=dict(color='purple')), row=1, col=3)
        
    #     del Dx
        
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=1)
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=2)
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=3)
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=4)
    #     figgamma.update_layout(height=600, width=900, title_text='Gamma evolution', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     ## Gamme for r ##
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_1_r, showlegend=False, line=dict(color='midnightblue')), row=1, col=1)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_2_r, showlegend=False, line=dict(color='firebrick')), row=1, col=1)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_3_r, showlegend=False, line=dict(color='green')), row=1, col=1)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u1_4_r, showlegend=False, line=dict(color='purple')), row=1, col=1)
        
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_1_r, showlegend=False, line=dict(color='midnightblue')), row=1, col=2)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_2_r, showlegend=False, line=dict(color='firebrick')), row=1, col=2)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_3_r, showlegend=False, line=dict(color='green')), row=1, col=2)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u2_4_r, showlegend=False, line=dict(color='purple')), row=1, col=2)
        
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_1_r, name=f'$\omega={omega[1]:.2f}$', line=dict(color='midnightblue')), row=1, col=3)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_2_r, name=f'$\omega={omega[2]:.2f}$', line=dict(color='firebrick')), row=1, col=3)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_3_r, name=f'$\omega={omega[3]:.2f}$', line=dict(color='green')), row=1, col=3)
    #     figgamma_r.add_trace(go.Scatter(x=z_axis, y=gamma_u3_4_r, name=f'$\omega={omega[4]:.2f}$', line=dict(color='purple')), row=1, col=3)
        
    #     del omega
        
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=1)
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=2)
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=3)
    #     figgamma_r.update_xaxes(title_text='$z^+$', row=1, col=4)
    #     figgamma_r.update_layout(height=600, width=900, title_text='Gamma evolution', legend=dict(yanchor="bottom", xanchor="right"), font=font)
        
    #     if split_time == 'Y':

    #         save_figures(figu1g, "split_time/gamma/gamma_u1_w_all.png")

    #         save_figures(figu2g, "split_time/gamma/gamma_u2_w_all.png")

    #         save_figures(figu3g, "split_time/gamma/gamma_u3_w_all.png")
            
    #         save_figures(figu1g_r, "split_time/gamma/gamma_u1_r_all.png")
            
    #         save_figures(figu2g_r, "split_time/gamma/gamma_u2_r_all.png")
        
    #         save_figures(figu3g_r, "split_time/gamma/gamma_u3_r_all.png")
            
    #         save_figures(figgamma, "split_time/gamma/gamma_view_w_all.png")
            
    #         save_figures(figgamma_r, "split_time/gamma/gamma_view_r_all.png")
            
    #     if split_time == 'n':

    #         save_figures(figu1g, "whole_time/gamma/gamma_u1_w_all.png")

    #         save_figures(figu2g, "whole_time/gamma/gamma_u2_w_all.png")

    #         save_figures(figu3g, "whole_time/gamma/gamma_u3_w_all.png")
            
    #         save_figures(figu1g_r, "whole_time/gamma/gamma_u1_r_all.png")
            
    #         save_figures(figu2g_r, "whole_time/gamma/gamma_u2_r_all.png")
        
    #         save_figures(figu3g_r, "whole_time/gamma/gamma_u3_r_all.png")
            
    #         save_figures(figgamma, "whole_time/gamma/gamma_view_w_all.png")
            
    #         save_figures(figgamma_r, "whole_time/gamma/gamma_view_r_all.png")
            
    # z_axis = np.array(z_axis)
    # gamma_u1_1 = np.array(gamma_u1_1)
    # gamma_u2_1 = np.array(gamma_u2_1)
    # gamma_u3_1 = np.array(gamma_u3_1)
    # gamma_u1_2 = np.array(gamma_u1_2)
    # gamma_u2_2 = np.array(gamma_u2_2)
    # gamma_u3_2 = np.array(gamma_u3_2)
    # gamma_u1_3 = np.array(gamma_u1_3)
    # gamma_u2_3 = np.array(gamma_u2_3)
    # gamma_u3_3 = np.array(gamma_u3_3)
    # gamma_u1_4 = np.array(gamma_u1_4)
    # gamma_u2_4 = np.array(gamma_u2_4)
    # gamma_u3_4 = np.array(gamma_u3_4)
    # gamma_u1_1_r = np.array(gamma_u1_1_r)
    # gamma_u2_1_r = np.array(gamma_u2_1_r)
    # gamma_u3_1_r = np.array(gamma_u3_1_r)
    # gamma_u1_2_r = np.array(gamma_u1_2_r)
    # gamma_u2_2_r = np.array(gamma_u2_2_r)
    # gamma_u3_2_r = np.array(gamma_u3_2_r)
    # gamma_u1_3_r = np.array(gamma_u1_3_r)
    # gamma_u2_3_r = np.array(gamma_u2_3_r)
    # gamma_u3_3_r = np.array(gamma_u3_3_r)
    # gamma_u1_4_r = np.array(gamma_u1_4_r)
    # gamma_u2_4_r = np.array(gamma_u2_4_r)
    # gamma_u3_4_r = np.array(gamma_u3_4_r)        
    
    # if split_time == 'Y':
    #     if chplot == 'all':
    #         save_datas([z_axis, gamma_u1_1, gamma_u1_2, gamma_u1_3, gamma_u1_4], ['zp', 'gamma(r=1)', 'gamma(r=5)', 'gamma(r=10)', 'gamma(r=15)'], f'split_time/gamma/u1_evolve_w.dat', 'Gamma evolve w (u1)')
    #         save_datas([z_axis, gamma_u2_1, gamma_u2_2, gamma_u2_3, gamma_u2_4], ['zp', 'gamma(r=1)', 'gamma(r=5)', 'gamma(r=10)', 'gamma(r=15)'], f'split_time/gamma/u2_evolve_w.dat', 'Gamma evolve w (u2)')
    #         save_datas([z_axis, gamma_u3_1, gamma_u3_2, gamma_u3_3, gamma_u3_4], ['zp', 'gamma(r=1)', 'gamma(r=5)', 'gamma(r=10)', 'gamma(r=15)'], f'split_time/gamma/u3_evolve_w.dat', 'Gamma evolve w (u3)')
            
    #         save_datas([z_axis, gamma_u1_1_r, gamma_u1_2_r, gamma_u1_3_r, gamma_u1_4_r], ['zp', 'gamma(w=1)', 'gamma(w=2)', 'gamma(w=3)', 'gamma(w=4)'], f'split_time/gamma/u1_evolve_r.dat', 'Gamma evolve r (u1)')
    #         save_datas([z_axis, gamma_u2_1_r, gamma_u2_2_r, gamma_u2_3_r, gamma_u2_4_r], ['zp', 'gamma(w=1)', 'gamma(w=2)', 'gamma(w=3)', 'gamma(w=4)'], f'split_time/gamma/u2_evolve_r.dat', 'Gamma evolve r (u2)')
    #         save_datas([z_axis, gamma_u3_1_r, gamma_u3_2_r, gamma_u3_3_r, gamma_u3_4_r], ['zp', 'gamma(w=1)', 'gamma(w=2)', 'gamma(w=3)', 'gamma(w=4)'], f'split_time/gamma/u3_evolve_r.dat', 'Gamma evolve r (u3)')
            
    # if split_time == 'n':
    #     if chplot == 'all':
    #         save_datas([z_axis, gamma_u1_1, gamma_u1_2, gamma_u1_3, gamma_u1_4], ['zp', 'gamma(r=1)', 'gamma(r=5)', 'gamma(r=10)', 'gamma(r=15)'], f'whole_time/gamma/u1_evolve_w.dat', 'Gamma evolve w (u1)')
    #         save_datas([z_axis, gamma_u2_1, gamma_u2_2, gamma_u2_3, gamma_u2_4], ['zp', 'gamma(r=1)', 'gamma(r=5)', 'gamma(r=10)', 'gamma(r=15)'], f'whole_time/gamma/u2_evolve_w.dat', 'Gamma evolve w (u2)')
    #         save_datas([z_axis, gamma_u3_1, gamma_u3_2, gamma_u3_3, gamma_u3_4], ['zp', 'gamma(r=1)', 'gamma(r=5)', 'gamma(r=10)', 'gamma(r=15)'], f'whole_time/gamma/u3_evolve_w.dat', 'Gamma evolve w (u3)')
            
    #         save_datas([z_axis, gamma_u1_1_r, gamma_u1_2_r, gamma_u1_3_r, gamma_u1_4_r], ['zp', 'gamma(w=1)', 'gamma(w=2)', 'gamma(w=3)', 'gamma(w=4)'], f'whole_time/gamma/u1_evolve_r.dat', 'Gamma evolve r (u1)')
    #         save_datas([z_axis, gamma_u2_1_r, gamma_u2_2_r, gamma_u2_3_r, gamma_u2_4_r], ['zp', 'gamma(w=1)', 'gamma(w=2)', 'gamma(w=3)', 'gamma(w=4)'], f'whole_time/gamma/u2_evolve_r.dat', 'Gamma evolve r (u2)')
    #         save_datas([z_axis, gamma_u3_1_r, gamma_u3_2_r, gamma_u3_3_r, gamma_u3_4_r], ['zp', 'gamma(w=1)', 'gamma(w=2)', 'gamma(w=3)', 'gamma(w=4)'], f'whole_time/gamma/u3_evolve_r.dat', 'Gamma evolve r (u3)')
            
            
    # elapsed_time = time.time() - start_time
    # minutes, seconds = divmod(elapsed_time, 60)
    
    # print(f'\n Gamma determination done in : {int(minutes)}m {seconds:.2f}s \n')
    
    
    
    ##=======================================================
                #### SPACE CORRELATION ####
    ##=======================================================
    
    # print("\n========================================")
    # print(f"{YELLOW}Space correlation computation{RESET}")
    # start_time = time.time()
    # print("========================================")
    # print(f"{YELLOW}Streamwise{RESET}")
    # col = 1
    # row = 1
    # print("\nReading input files ...")
    
    # _, x1, x2, _, nt, n1, n2, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane(fpars_files_streamwise_u1[0])
    # nt = nt - 1
    # dx = cflow.xlen / n1
    
    # fig1sc= init_figures_sc(zp, ch=chplot) #Figure initialization
    
    # #for zplan in np.arange(0, n2, n2//3, dtype=int):
    # for ind, zplan in enumerate(zp_ind):
        
    #     print("========================================")
    #     print(f'Space correlation for {YELLOW}zp={zp[ind]:.2f}{RESET}')
    #     print('Plan number:', zplan)
    #     print("Reading input files streamwise...")
        
    #     num_split_t = nt // split_t
    #     Dx = np.linspace(0,xlen//2,n1)
        
    #     ## uu ##
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
    #     datas_u1 = var[1:,:,:]
        
    #     Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
    #     print('Ux:', Ux)
        
    #     datas_u1 = datas_u1 - Ux
        
    #     if split_time == 'Y':
            
    #         R_uu = np.zeros((n1))
    #         for n in tqdm(range(1,num_split_t), desc=f'PSD (zp={zp[ind]:.2f})', colour= 'GREEN'):
    #             R_uu += Space_correation(datas_u1[(n-1)*split_t:n*split_t,:,:], datas_u1[(n-1)*split_t:n*split_t,:,:], geom = "plan", mode_corr = 'half', axis = "streamwise")
    #         R_uu /= (num_split_t-1)
            
    #     if split_time == 'n':
    #         R_uu = Space_correation(datas_u1, datas_u1, geom = "plan", mode_corr = 'half', axis = "streamwise")
            
    #     space_correlation_plot(fig1sc, col, row, Dx, R_uu, name = '$R_{UU}$', color='midnightblue', axis='streamwise')
        
    #     if chplot == 'all':
    #         if split_time == 'Y':
    #             save_datas([Dx, R_uu], ['Dx', 'R_uu'], f'split_time/space_correlation/uu_streamwise_{zp[ind]}.dat', 'Space correlation streamwise (uu)')
    #         if split_time == 'n':
    #             save_datas([Dx, R_uu], ['Dx', 'R_uu'], f'whole_time/space_correlation/uu_streamwise_{zp[ind]}.dat', 'Space correlation streamwise (uu)')
        
    #     del R_uu
    #     del datas_u1
        
    #     ## vv ##
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
    #     datas_u2 = var[1:,:,:]
        
    #     Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
        
    #     datas_u2 = datas_u2 - Uy
        
    #     if split_time == 'Y':
            
    #         R_vv = np.zeros((n1))
    #         for n in tqdm(range(1,num_split_t), desc=f'PSD (zp={zp[ind]:.2f})', colour= 'GREEN'):
    #             R_vv += Space_correation(datas_u2[(n-1)*split_t:n*split_t,:,:], datas_u2[(n-1)*split_t:n*split_t,:,:], geom = "plan", mode_corr = 'half', axis = "streamwise")
    #         R_vv /= (num_split_t-1)
        
    #     if split_time == 'n':
    #         R_vv = Space_correation(datas_u2, datas_u2, geom = "plan", mode_corr = 'half', axis = "streamwise")
            
    #     space_correlation_plot(fig1sc, col, row, Dx, R_vv, name = '$R_{VV}$', color='firebrick', axis='streamwise')
        
    #     if chplot == 'all':
    #         if split_time == 'Y':
    #             save_datas([Dx, R_vv], ['Dx', 'R_vv'], f'split_time/space_correlation/vv_streamwise_{zp[ind]}.dat', 'Space correlation streamwise (vv)')
    #         if split_time == 'n':
    #             save_datas([Dx, R_vv], ['Dx', 'R_vv'], f'whole_time/space_correlation/vv_streamwise_{zp[ind]}.dat', 'Space correlation streamwise (vv)')
        
    #     del R_vv
    #     del datas_u2
        
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
    #     datas_u3 = var[1:,:,:]
        
    #     Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
        
    #     datas_u3 = datas_u3 - Uz
        
    #     if split_time == 'Y':
            
    #         R_ww = np.zeros((n1))
    #         for n in tqdm(range(1,num_split_t), desc=f'PSD (zp={zp[ind]:.2f})', colour= 'GREEN'):
    #             R_ww += Space_correation(datas_u3[(n-1)*split_t:n*split_t,:,:], datas_u3[(n-1)*split_t:n*split_t,:,:], geom = "plan", mode_corr = 'half', axis = "streamwise")
    #         R_ww /= (num_split_t-1)
        
    #     if split_time == 'n':
    #         R_ww = Space_correation(datas_u3, datas_u3, geom = "plan", mode_corr = 'half', axis = "streamwise")
            
    #     space_correlation_plot(fig1sc, col, row, Dx, R_ww, name = '$R_{WW}$', color='darkgreen', axis='streamwise')
        
    #     if chplot == 'all':
    #         if split_time == 'Y':
    #             save_datas([Dx, R_ww], ['Dx', 'R_ww'], f'split_time/space_correlation/ww_streamwise_{zp[ind]}.dat', 'Space correlation streamwise (ww)')
    #         if split_time == 'n':
    #             save_datas([Dx, R_ww], ['Dx', 'R_ww'], f'whole_time/space_correlation/ww_streamwise_{zp[ind]}.dat', 'Space correlation streamwise (ww)')
        
    #     del R_ww
    #     del datas_u3
        
    #     col +=1
    #     if zplan == 4:
    #         row +=1
    #         col = 1
            
        
    # if chplot == 'normal':
    #     fig1sc.update_layout(height=600, width=900, title_text='Space Correlation Streamwise', font=font,  legend=dict(yanchor="bottom", y=1.03, xanchor="right", x=1, orientation='h'))
    # if chplot == 'all':
    #     fig1sc.update_layout(height=900, width=900, title_text='Space Correlation Streamwise', font=font,  legend=dict(yanchor="bottom", y=1.03, xanchor="right", x=1, orientation='h'))
        
    # if split_time == 'Y':
        
    #     if chplot == 'normal':
    #         save_figures(fig1sc, "split_time/space_correlation/streamwise.png")
    #     if chplot == 'all':
    #         save_figures(fig1sc, "split_time/space_correlation/streamwise_all.png")
            
    # if split_time == 'n':
        
    #     if chplot == 'normal':
    #         save_figures(fig1sc, "whole_time/space_correlation/streamwise.png")
    #     if chplot == 'all':
    #         save_figures(fig1sc, "whole_time/space_correlation/streamwise_all.png")
            
    # print("\n========================================")
    # print(f"{YELLOW}Spanwise{RESET}")
    # col = 1
    # row = 1
    # print("Reading input files ...")
    
    # _, x1, x2, _, nt, n2, n1, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane(fpars_files_spanwise_u1[0])
    # nt = nt - 1
    # print('n1:', n1)
    # print('n2:', n2)
    # dx = cflow.xlen / n1
    
    # fig2sc= init_figures_sc(zp, ch=chplot) #Figure initialization
    
    # #for zplan in np.arange(0, n2, n2//3, dtype=int):
    # for ind, zplan in enumerate(zp_ind):
        
    #     print("========================================")
    #     print(f'Space correlation for {YELLOW}zp={zp[ind]:.2f}{RESET}')
    #     print('Plan number:', zplan)
    #     print("Reading input files spanwise...")
        
    #     num_split_t = nt // split_t
    #     Dx = np.linspace(0,ylen//2,n1)
        
    #     ## uu ##
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_spanwise_u1[zplan])
    #     datas_u1 = var[1:,:,:]
        
    #     Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
        
    #     datas_u1 = datas_u1 - Ux
        
    #     if split_time == 'Y':
            
    #         R_uu = np.zeros((n1))
    #         for n in tqdm(range(1,num_split_t), desc=f'PSD (z={zp[ind]:.2f})', colour= 'GREEN'):
    #             R_uu += Space_correation(datas_u1[(n-1)*split_t:n*split_t,:,:], datas_u1[(n-1)*split_t:n*split_t,:,:], geom = "plan", mode_corr = 'half', axis = "spanwise")
    #         R_uu /= (num_split_t-1)
            
    #     if split_time == 'n':
    #         R_uu = Space_correation(datas_u1, datas_u1, geom = "plan", mode_corr = 'half', axis = "spanwise")
            
    #     space_correlation_plot(fig2sc, col, row, Dx, R_uu, name = '$R_{UU}$', color='midnightblue', axis='spanwise')
        
    #     if chplot == 'all':
    #         if split_time == 'Y':
    #             save_datas([Dx, R_uu], ['Dx', 'R_uu'], f'split_time/space_correlation/uu_streamwise_{zp[ind]}.dat', 'Space correlation spanwise (uu)')
    #         if split_time == 'n':
    #             save_datas([Dx, R_uu], ['Dx', 'R_uu'], f'whole_time/space_correlation/uu_streamwise_{zp[ind]}.dat', 'Space correlation spanwise (uu)')
        
    #     del R_uu
    #     del datas_u1
        
    #     ## vv ##
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_spanwise_u2[zplan])
    #     datas_u2 = var[1:,:,:]
        
    #     Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
        
    #     datas_u2 = datas_u2 - Uy
        
    #     if split_time == 'Y':
            
    #         R_vv = np.zeros((n1))
    #         for n in tqdm(range(1,num_split_t), desc=f'PSD (z={zp[ind]:.2f})', colour= 'GREEN'):
    #             R_vv += Space_correation(datas_u2[(n-1)*split_t:n*split_t,:,:], datas_u2[(n-1)*split_t:n*split_t,:,:], geom = "plan", mode_corr = 'half', axis = "spanwise")
    #         R_vv /= (num_split_t-1)
        
    #     if split_time == 'n':
    #         R_vv = Space_correation(datas_u2, datas_u2, geom = "plan", mode_corr = 'half', axis = "spanwise")
            
    #     space_correlation_plot(fig2sc, col, row, Dx, R_vv, name = '$R_{VV}$', color='firebrick', axis='spanwise')
        
    #     if chplot == 'all':
    #         if split_time == 'Y':
    #             save_datas([Dx, R_vv], ['Dx', 'R_vv'], f'split_time/space_correlation/vv_streamwise_{zp[ind]}.dat', 'Space correlation spanwise (vv)')
    #         if split_time == 'n':
    #             save_datas([Dx, R_vv], ['Dx', 'R_vv'], f'whole_time/space_correlation/vv_streamwise_{zp[ind]}.dat', 'Space correlation spanwise (vv)')
        
    #     del R_vv
    #     del datas_u2
        
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_spanwise_u3[zplan])
    #     datas_u3 = var[1:,:,:]
        
    #     Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
        
    #     datas_u3 = datas_u3 - Uz
        
    #     if split_time == 'Y':
            
    #         R_ww = np.zeros((n1))
    #         for n in tqdm(range(1,num_split_t), desc=f'PSD (z={zp[ind]:.2f})', colour= 'GREEN'):
    #             R_ww += Space_correation(datas_u3[(n-1)*split_t:n*split_t,:,:], datas_u3[(n-1)*split_t:n*split_t,:,:], geom = "plan", mode_corr = 'half', axis = "spanwise")
    #         R_ww /= (num_split_t-1)
        
    #     if split_time == 'n':
    #         R_ww = Space_correation(datas_u3, datas_u3, geom = "plan", mode_corr = 'half', axis = "spanwise")
            
    #     space_correlation_plot(fig2sc, col, row, Dx, R_ww, name = '$R_{WW}$', color='darkgreen', axis='spanwise')
        
    #     if chplot == 'all':
    #         if split_time == 'Y':
    #             save_datas([Dx, R_ww], ['Dx', 'R_ww'], f'split_time/space_correlation/ww_streamwise_{zp[ind]}.dat', 'Space correlation spanwise (ww)')
    #         if split_time == 'n':
    #             save_datas([Dx, R_ww], ['Dx', 'R_ww'], f'whole_time/space_correlation/ww_streamwise_{zp[ind]}.dat', 'Space correlation spanwise (ww)')
        
    #     del R_ww
    #     del datas_u3
        
    #     col +=1
    #     if zplan == 4:
    #         row +=1
    #         col = 1
            
        
    # if chplot == 'normal':
    #     fig2sc.update_layout(height=600, width=900, title_text='Space Correlation Spanwise', font=font,  legend=dict(yanchor="bottom", y=1.03, xanchor="right", x=1, orientation='h'))
    # if chplot == 'all':
    #     fig2sc.update_layout(height=900, width=900, title_text='Space Correlation Spanwise', font=font,  legend=dict(yanchor="bottom", y=1.03, xanchor="right", x=1, orientation='h'))
        
    # if split_time == 'Y':
        
    #     if chplot == 'normal':
    #         save_figures(fig2sc, "split_time/space_correlation/spanwise.png")
    #     if chplot == 'all':
    #         save_figures(fig2sc, "split_time/space_correlation/spanwise_all.png")
            
    # if split_time == 'n':
        
    #     if chplot == 'normal':
    #         save_figures(fig2sc, "whole_time/space_correlation/spanwise.png")
    #     if chplot == 'all':
    #         save_figures(fig2sc, "whole_time/space_correlation/spanwise_all.png")
            
        
    # elapsed_time = time.time() - start_time
    # minutes, seconds = divmod(elapsed_time, 60)
    
    # print(f'\n Space correlation done in : {int(minutes)}m {seconds:.2f}s \n')
        
    
    
    ##=======================================================
                #### RANS AND LES ANALYSIS ####
    ##=======================================================
    
    # print("\n========================================")
    # print(f"{YELLOW}RANS AND LES COMPARISON{RESET}")
    # start_time = time.time()
    # print("========================================")
    # fig_RANS = make_subplots(rows=1, cols=4, shared_yaxes= True, y_title='$z^+$')
    # fig_var = make_subplots(rows=2, cols=3, shared_yaxes='rows', row_titles=('$z^+$','$z^+$'), vertical_spacing=0.15)
    # fig_vel_ratio_profil = go.Figure()
    
    
    # base_path = os.path.join(os.path.dirname(__file__), "../")
    # _, u_velocity, v_velocity, w_velocity, _, _, _, normal, ko2_tke, _ = read_rans(os.path.join(base_path,rans_path))
    # u_velocity = u_velocity / cflow.ut
    # v_velocity /= cflow.ut
    # w_velocity /= cflow.ut
    # ko2_tke /= cflow.ut**2
    # zp_RANS = cflow.ut/cflow.nu * normal
    
    # U = []
    # V = []
    # W = []
    # uu = []
    # vv = []
    # ww = []
    # uv = []
    # uw = []
    # vw = []
    # kt = []
    # for ind, zplan in enumerate(zp_ind):
    #     print('zp:', zp[ind])
        
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
    #     datas_u1 = var[1:,:,:]
    #     del var
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
    #     datas_u2 = var[1:,:,:]
    #     del var
    #     _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
    #     datas_u3 = var[1:,:,:]
    #     del var
        
    #     Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
    #     Uy = - np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
    #     Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
        
    #     datas_u1 = datas_u1 - Ux
    #     datas_u2 = -(datas_u2 - Uy)
    #     datas_u3 = datas_u3 - Uz
        
    #     varUx = np.mean(np.mean(np.mean((datas_u1)**2, axis=0), axis=-1))
    #     varUy = np.mean(np.mean(np.mean((datas_u2)**2, axis=0), axis=-1))
    #     varUz = np.mean(np.mean(np.mean((datas_u3)**2, axis=0), axis=-1))
    #     varUV = np.mean(np.mean(np.mean(datas_u1*datas_u2, axis=0), axis=-1))
    #     varUW = np.mean(np.mean(np.mean(datas_u1*datas_u3, axis=0), axis=-1))
    #     varVW = np.mean(np.mean(np.mean(datas_u3*datas_u2, axis=0), axis=-1))
        
    #     del datas_u1
    #     del datas_u2
    #     del datas_u3
        
    #     Ux /= cflow.ut
    #     Uy /= cflow.ut
    #     Uz /= cflow.ut
    #     varUx /= cflow.ut**2
    #     varUy /= cflow.ut**2
    #     varUz /= cflow.ut**2
    #     varUV /= cflow.ut**2
    #     varUW /= cflow.ut**2
    #     varVW /= cflow.ut**2
        
    #     U.append(Ux)
    #     V.append(Uy)
    #     W.append(Uz)
        
    #     uu.append(varUx)
    #     vv.append(varUy)
    #     ww.append(varUz)
    #     uv.append(varUV)
    #     uw.append(varUW)
    #     vw.append(varVW)
        
    #     kt.append(varUx + varUy + varUz)
        
    # kt = np.array(kt)
    # U = np.array(U)
    # V = np.array(V)
    # W = np.array(W)
    # uu = np.array(uu)
    # vv = np.array(vv)
    # ww = np.array(ww)
    # uv = np.array(uv)
    # uw = np.array(uw)
    # vw = np.array(vw)
    # kt /= 2.
    # ratio2 = vv / uu
    # ratio3 = ww / uu
    # ratio1 = ww / vv
    
    # if model_name == "WRLES_Retau395":
    #     zp_moser = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.means.txt')[1]
    #     U_moser = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.means.txt')[2]
    #     V_moser = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.means.txt')[4]
        
    
    #     fig_RANS.add_trace(go.Scatter(x=u_velocity, y=zp_RANS, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle-open'), name='RANS'), row=1, col=1)
    #     fig_RANS.add_trace(go.Scatter(x=v_velocity, y=zp_RANS, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle-open'), showlegend=False), row=1, col=2)
    #     fig_RANS.add_trace(go.Scatter(x=w_velocity, y=zp_RANS, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle-open'), showlegend=False), row=1, col=3)
    #     fig_RANS.add_trace(go.Scatter(x=ko2_tke, y=zp_RANS, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle-open'), showlegend=False), row=1, col=4)
    #     fig_RANS.add_trace(go.Scatter(x=U, y=zp, mode= 'lines+markers', line=dict(color='forestgreen', width=2), marker=dict(symbol='x'), name='LES'), row=1, col=1)
    #     fig_RANS.add_trace(go.Scatter(x=V, y=zp, mode= 'lines+markers', line=dict(color='forestgreen', width=2), marker=dict(symbol='x'), showlegend=False), row=1, col=2)
    #     fig_RANS.add_trace(go.Scatter(x=W, y=zp, mode= 'lines+markers', line=dict(color='forestgreen', width=2), marker=dict(symbol='x'), showlegend=False), row=1, col=3)
    #     fig_RANS.add_trace(go.Scatter(x=kt, y=zp, mode= 'lines+markers', line=dict(color='forestgreen', width=2), marker=dict(symbol='x'), showlegend=False), row=1, col=4)
    #     fig_RANS.add_trace(go.Scatter(x=U_moser, y=zp_moser, mode= 'lines+markers', line=dict(color='midnightblue', width=2), marker=dict(symbol='diamond-open'), name='MOSER'), row=1, col=1)
    #     fig_RANS.add_trace(go.Scatter(x=V_moser, y=zp_moser, mode= 'lines+markers', line=dict(color='midnightblue', width=2), marker=dict(symbol='diamond-open'), showlegend=False), row=1, col=2)
        
        
    #     fig_RANS.update_xaxes(title='$U$', row=1, col=1)
    #     fig_RANS.update_xaxes(title='$V$', row=1, col=2)
    #     fig_RANS.update_xaxes(title='$W$', row=1, col=3)
    #     fig_RANS.update_xaxes(title='$k_T$', row=1, col=4)
    #     fig_RANS.update_layout(height=600, width=900, title=f"Mean velocity profile LES, DNS and RANS comparison (Retau={ret})", font=font, showlegend=True, legend=dict(yanchor='bottom', xanchor='left'))
        
    # else:
    #     zp_moser = read_moser('input/Moser/Re1000/RANS_comparison/Moser_chan1000_mean_prof.txt')[1]
    #     U_moser = read_moser('input/Moser/Re1000/RANS_comparison/Moser_chan1000_mean_prof.txt')[2]
    #     V_moser = read_moser('input/Moser/Re1000/RANS_comparison/Moser_chan1000_mean_prof.txt')[4]
    #     kt_moser = read_moser('input/Moser/Re1000/RANS_comparison/Moser_chan1000_vel_fluc_prof.txt')[8]
        
    #     fig_RANS.add_trace(go.Scatter(x=U, y=zp, mode= 'lines+markers', line=dict(color='forestgreen', width=2), marker=dict(symbol='x'), name='LES'), row=1, col=1)
    #     fig_RANS.add_trace(go.Scatter(x=V, y=zp, mode= 'lines+markers', line=dict(color='forestgreen', width=2), marker=dict(symbol='x'), showlegend=False), row=1, col=2)
    #     fig_RANS.add_trace(go.Scatter(x=W, y=zp, mode= 'lines+markers', line=dict(color='forestgreen', width=2), marker=dict(symbol='x'), showlegend=False), row=1, col=3)
    #     fig_RANS.add_trace(go.Scatter(x=kt, y=zp, mode= 'lines+markers', line=dict(color='forestgreen', width=2), marker=dict(symbol='x'), showlegend=False), row=1, col=4)
    #     fig_RANS.add_trace(go.Scatter(x=U_moser, y=zp_moser, mode= 'lines+markers', line=dict(color='black', width=2), marker=dict(symbol='diamond-open'), name='MOSER$'), row=1, col=1)
    #     fig_RANS.add_trace(go.Scatter(x=V_moser, y=zp_moser, mode= 'lines+markers', line=dict(color='black', width=2), marker=dict(symbol='diamond-open'), showlegend=False), row=1, col=2)
    #     fig_RANS.add_trace(go.Scatter(x=kt_moser, y=zp_moser, mode= 'lines+markers', line=dict(color='black', width=2), marker=dict(symbol='diamond-open'), showlegend=False), row=1, col=4)
        
    #     fig_RANS.update_yaxes(range=[0,400])
    #     fig_RANS.update_xaxes(title='$U$', row=1, col=1)
    #     fig_RANS.update_xaxes(title='$V$', row=1, col=2)
    #     fig_RANS.update_xaxes(title='$W$', row=1, col=3)
    #     fig_RANS.update_xaxes(title='$k_T$', row=1, col=4)
    #     fig_RANS.update_layout(height=600, width=900, title=f"Mean velocity profile LES and DNS comparison (Retau={ret})", font=font, showlegend=True, legend=dict(yanchor='bottom', xanchor='left'))
        

    # if split_time == 'Y':
    #     if chplot == 'all':
    #         save_figures(fig_RANS, "split_time/RANS/RANS_LES_MOSER_profiles_all.png")
    # if split_time == 'n':
    #     if chplot == 'all':
    #         save_figures(fig_RANS, "whole_time/RANS/RANS_LES_MOSER_profiles_all.png")
            
    
    # if model_name == "WRLES_Retau395":
    #     zp_moser = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.reystress.txt')[1]
    #     uu_moser = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.reystress.txt')[2]
    #     ww_moser = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.reystress.txt')[3]
    #     vv_moser = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.reystress.txt')[4]
    #     uw_moser = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.reystress.txt')[5]
    #     uv_moser = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.reystress.txt')[6]
    #     vw_moser = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.reystress.txt')[7]
    # else:
    #     zp_moser = read_moser('input/Moser/Re1000/RANS_comparison/Moser_chan1000_vel_fluc_prof.txt')[1]
    #     uu_moser = read_moser('input/Moser/Re1000/RANS_comparison/Moser_chan1000_vel_fluc_prof.txt')[2]
    #     ww_moser = read_moser('input/Moser/Re1000/RANS_comparison/Moser_chan1000_vel_fluc_prof.txt')[3]
    #     vv_moser = read_moser('input/Moser/Re1000/RANS_comparison/Moser_chan1000_vel_fluc_prof.txt')[4]
    #     uw_moser = read_moser('input/Moser/Re1000/RANS_comparison/Moser_chan1000_vel_fluc_prof.txt')[5]
    #     uv_moser = read_moser('input/Moser/Re1000/RANS_comparison/Moser_chan1000_vel_fluc_prof.txt')[6]
    #     vw_moser = read_moser('input/Moser/Re1000/RANS_comparison/Moser_chan1000_vel_fluc_prof.txt')[7]
        
        
    # fig_var.add_trace(go.Scatter(x=uu, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'), name='LES'), row=1, col=1)
    # fig_var.add_trace(go.Scatter(x=vv, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'),showlegend=False), row=1, col=2)
    # fig_var.add_trace(go.Scatter(x=ww, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'),showlegend=False),  row=1, col=3)
    # fig_var.add_trace(go.Scatter(x=uv, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'),showlegend=False),  row=2, col=1)
    # fig_var.add_trace(go.Scatter(x=uw, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'),showlegend=False),  row=2, col=2)
    # fig_var.add_trace(go.Scatter(x=vw, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'),showlegend=False),  row=2, col=3)
    
    # fig_var.add_trace(go.Scatter(x=uu_moser, y=zp_moser, mode= 'lines+markers', line=dict(color='black', width=2), marker=dict(symbol='circle-open'), name='MOSER'), row=1, col=1)
    # fig_var.add_trace(go.Scatter(x=vv_moser, y=zp_moser, mode= 'lines+markers', line=dict(color='black', width=2), marker=dict(symbol='circle-open'),showlegend=False), row=1, col=2)
    # fig_var.add_trace(go.Scatter(x=ww_moser, y=zp_moser, mode= 'lines+markers', line=dict(color='black', width=2), marker=dict(symbol='circle-open'),showlegend=False),  row=1, col=3)
    # fig_var.add_trace(go.Scatter(x=uv_moser, y=zp_moser, mode= 'lines+markers', line=dict(color='black', width=2), marker=dict(symbol='circle-open'),showlegend=False),  row=2, col=1)
    # fig_var.add_trace(go.Scatter(x=uw_moser, y=zp_moser, mode= 'lines+markers', line=dict(color='black', width=2), marker=dict(symbol='circle-open'),showlegend=False),  row=2, col=2)
    # fig_var.add_trace(go.Scatter(x=vw_moser, y=zp_moser, mode= 'lines+markers', line=dict(color='black', width=2), marker=dict(symbol='circle-open'),showlegend=False),  row=2, col=3)
    
    # fig_RANS.update_yaxes(range=[0,400])
    # fig_var.update_xaxes(title='$\overline{u_1u_1}$', row=1, col=1)
    # fig_var.update_xaxes(title='$\overline{u_2u_2}$', row=1, col=2)
    # fig_var.update_xaxes(title='$\overline{u_3u_3}$', row=1, col=3)
    # fig_var.update_xaxes(title='$\overline{u_1u_2}$', row=2, col=1)
    # fig_var.update_xaxes(title='$\overline{u_1u_3}$', row=2, col=2)
    # fig_var.update_xaxes(title='$\overline{u_2u_3}$', row=2, col=3)
    # fig_var.update_layout(height=800, width=800, title=f"Raynolds stress component profiles (Retau={ret})", font=font, showlegend=True, legend=dict(yanchor='bottom', xanchor='left'))
    
    # if split_time == 'Y':
    #     if chplot == 'all':
    #         save_figures(fig_var, "split_time/RANS/var_velocity_profiles_all.png")
    # if split_time == 'n':
    #     if chplot == 'all':
    #         save_figures(fig_var, "whole_time/RANS/var_velocity_profiles_all.png")
    
    # fig_vel_ratio_profil.add_trace(go.Scatter(x=ratio2, y=zp, mode= 'lines+markers', line=dict(color='firebrick', width=2), marker=dict(symbol='circle'), name='$\overline{u_2u_2} / \overline{u_1u_1}$'))
    # fig_vel_ratio_profil.add_trace(go.Scatter(x=ratio3, y=zp, mode= 'lines+markers', line=dict(color='midnightblue', width=2), marker=dict(symbol='diamond'), name='$\overline{u_3u_3} / \overline{u_1u_1}$'))
    # fig_vel_ratio_profil.add_trace(go.Scatter(x=ratio1, y=zp, mode= 'lines+markers', line=dict(color='green', width=2), marker=dict(symbol='x'), name='$\overline{u_3u_3} / \overline{u_2u_2}$'))
    
    # fig_vel_ratio_profil.update_yaxes(title='$z^+$')
    # fig_vel_ratio_profil.update_xaxes(title='velocity ratio')
    # fig_vel_ratio_profil.update_layout(height=600, width=800, title="Variance velocity ratio profile (LES data)", font=font, showlegend=True, legend=dict(yanchor='bottom', xanchor='left'))
    
    # zp_array = np.array(zp)
    
    # if split_time == 'Y':
    #     if chplot == 'all':
    #         save_figures(fig_vel_ratio_profil, "split_time/RANS/velocity_ratio_profiles_all.png")
    #         save_datas([zp_array, U, V, W, uu, vv, ww, uv, uw, vw, ratio1, ratio2, ratio3], ['zp', 'U', 'V', 'W', 'uu', 'vv', 'ww', 'uv', 'uw', 'vw', 'ww/vv', 'vv/uu', 'ww/uu'], f'split_time/RANS/all.dat', 'First analysis (vellocity, variance, ratio proiles)')
    # if split_time == 'n':
    #     if chplot == 'all':
    #         save_figures(fig_vel_ratio_profil, "whole_time/RANS/velocity_ratio_profiles_all.png")
    #         save_datas([zp_array, U, V, W, uu, vv, ww, uv, uw, vw, ratio1, ratio2, ratio3], ['zp', 'U', 'V', 'W', 'uu', 'vv', 'ww', 'uv', 'uw', 'vw', 'ww/vv', 'vv/uu', 'ww/uu'], f'whole_time/RANS/all.dat', 'First analysis (vellocity, variance, ratio proiles)')
    
    
    # elapsed_time = time.time() - start_time
    # minutes, seconds = divmod(elapsed_time, 60)
    
    # print(f'\n RANS study done in : {int(minutes)}m {seconds:.2f}s \n')
    
    
    ##=======================================================
                #### NORMAL PLAN COMPUTATION ####
    ##=======================================================
    # if model_name == 'WRLES_Retau395':
    #     print("\n========================================")
    #     print(f"{YELLOW}NORMAL PLAN COMPUTATION{RESET}")
    #     start_time = time.time()
    #     print("========================================")
    #     nlines = len(fpars_files_normal_u1)
    #     _, x1, x2, _, nt, n2, n1, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane_line(fpars_files_normal_u1[0])
    #     nt = nt - 1
        
    #     base_path = os.path.join(os.path.dirname(__file__), "../")
    #     _, u_velocity, v_velocity, w_velocity, _, _, _, normal, ko2_tke, _ = read_rans(os.path.join(base_path,rans_path))
    #     u_velocity = u_velocity / cflow.ut
    #     v_velocity /= cflow.ut
    #     w_velocity /= cflow.ut
    #     ko2_tke /= cflow.ut**2
    #     zp_RANS = cflow.ut/cflow.nu * normal
        
    #     print('n1:', n1)
    #     print('n2:', n2)
    #     print('x2:', x2)
    #     x2 = np.array(x2)
    #     x2 *= cflow.ut / cflow.nu
    #     zp_RANS = cflow.ut/cflow.nu * normal
    #     print('x2[0]:',x2[0])
    #     print('x2[len(x2)]:',x2[len(x2)-1])
        
    #     #Figure initiaization##
    #     fig_u_z = make_subplots(rows=1, cols=3, shared_yaxes= True, y_title='$z^+$')
    #     fig_up_prof = go.Figure()
    #     #fig_corr_z = go.Figure()
        
    #     print('U1 profile')
    #     data_u1 = np.zeros((nt, n1), dtype=float)
    #     for line in range(nlines):
    #         print('line number:', line)
            
    #         var = read_fpar_extract_plane_line(fpars_files_normal_u3[line])[3]
    #         data_u1 += var[1:,:]
    #         del var
            
    #     data_u1 /= nlines
    #     data_u1 = np.mean(data_u1, axis=0)
        
    #     zp_moser = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.means.txt')[1]
    #     U_moser = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.means.txt')[2]
        
    #     fig_u_z.add_trace(go.Scatter(x=data_u1/cflow.ut, y=x2, mode= 'lines', line=dict(color='firebrick', width=2), name='$U_1(LES)$'), row=1, col=1)
    #     fig_u_z.add_trace(go.Scatter(x=u_velocity, y=zp_RANS, mode= 'markers', marker=dict(color='firebrick', symbol='circle-open'), name='$U_1(RANS)$'), row=1, col=1)
        
    #     fig_up_prof.add_trace(go.Scatter(x=x2, y=data_u1/cflow.ut, mode= 'lines + markers', line=dict(color='firebrick', width=2), marker=dict(color='firebrick', symbol='circle-open'), name='$WRLES$'))
    #     fig_up_prof.add_trace(go.Scatter(x=zp_RANS, y=u_velocity, mode= 'markers', marker=dict(color='firebrick', symbol='circle-open'), name='$RANS$'))
    #     fig_up_prof.add_trace(go.Scatter(x=zp_moser, y=U_moser, mode= 'lines', marker=dict(color='black'), name='$DNS$'))
    #     fig_up_prof.update_xaxes(type='log', exponentformat = 'power', title_text = 'z^+')
    #     fig_up_prof.update_yaxes(title_text = 'U_1^+')
        
    #     fig_up_prof.update_layout(height=450, width=700, title="Stramwise velocity", font=font, showlegend=True, legend=dict(yanchor='bottom', xanchor='right'))
        
    #     if split_time == 'Y':
    #         save_figures(fig_u_z, "split_time/Normal_plan/u1_ad_comp.png")
    #     if split_time == 'n':
    #         save_figures(fig_u_z, "whole_time/Normal_plan/u1_ad_comp.png")
        
    #     del data_u1
        
    #     print('U2 profile')
    #     data_u2 = np.zeros((nt, n1), dtype=float)
    #     for line in range(nlines):
    #         print('line number:', line)
            
    #         var = read_fpar_extract_plane_line(fpars_files_normal_u2[line])[3]
    #         data_u2 += var[1:,:]
    #         del var
            
    #     data_u2 /= nlines
    #     data_u2 = np.mean(data_u2, axis=0)
    #     #data_u2 /= cflow.ut
    #     fig_u_z.add_trace(go.Scatter(x=data_u2, y=x2, mode= 'lines', line=dict(color='midnightblue', width=2), name='$U_2(LES)$'), row=1, col=2)
    #     #fig_u_z.add_trace(go.Scatter(x=v_velocity, y=zp_RANS, mode= 'markers', marker=dict(color='midnightblue', symbol='circle-open'), name='$U_2(RANS)$'), row=1, col=2)
        
    #     del data_u2
        
    #     print('U3 profile')
    #     data_u3 = np.zeros((nt, n1), dtype=float)
    #     for line in range(nlines):
    #         print('line number:', line)
            
    #         var = read_fpar_extract_plane_line(fpars_files_normal_u1[line])[3]
    #         data_u3 += var[1:,:]
    #         del var
            
    #     data_u3 /= nlines
    #     data_u3 = np.mean(data_u3, axis=0)
    #     #data_u3 /= cflow.ut
    #     fig_u_z.add_trace(go.Scatter(x=data_u3, y=x2, mode= 'lines', line=dict(color='darkgreen', width=2), name='$U_3(LES)$'), row=1, col=3)
    #     #fig_u_z.add_trace(go.Scatter(x=w_velocity, y=zp_RANS, mode= 'markers', marker=dict(color='darkgreen', symbol='circle-open'), name='$U_3(RANS)$'), row=1, col=3)
        
    #     del data_u3
        
    #     #fig_u1_z.update_yaxes(title='$z^+$')
    #     fig_u_z.update_xaxes(title='$U$', row=1, col=1)
    #     fig_u_z.update_xaxes(title='$V$', row=1, col=2)
    #     fig_u_z.update_xaxes(title='$W$', row=1, col=3)
    #     fig_u_z.update_layout(height=600, width=800, title="Wall-normal velocity profile", font=font, showlegend=True, legend=dict(yanchor='bottom', xanchor='right'))
        
    #     if split_time == 'Y':
    #         save_figures(fig_u_z, "split_time/Normal_plan/velocity_profiles.png")
    #     if split_time == 'n':
    #         save_figures(fig_u_z, "whole_time/Normal_plan/velocity_profiles.png")
        
        ###############
        # Correation ##
        ###############
        
        # if len(zp) == 4:
        #     fig_corr_z = make_subplots(rows=1, cols=4, shared_yaxes= True, y_title='$\delta z^+$', subplot_titles=(f"$z_0={zp[0]:.2f}$", f"$z_0={zp[1]:.2f}$", f"$z_0={zp[2]:.2f}$", f"$z_0={zp[3]:.2f}$"))
        # else:
        #     fig_corr_z = make_subplots(rows=2, cols=5, shared_yaxes= True, y_title='$\delta z^+$', subplot_titles=(f"$z_0={z[0]:.2f}$", f"$z_0={z[1]:.2f}$", f"$z_0={z[2]:.2f}$", f"$z_0={z[3]:.2f}$", f"$z_0={z[4]:.2f}$", f"$z_0={z[5]:.2f}$", f"$z_0={z[6]:.2f}$", f"$z_0={z[7]:.2f}$", f"$z_0={z[8]:.2f}$", f"$z_0={z[9]:.2f}$"))
        
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
        
        # if split_time == 'Y':
        #     save_figures(fig_corr_z, "split_time/Normal_plan/autocorrelation_z.png")
        # if split_time == 'n':
        #     save_figures(fig_corr_z, "whole_time/Normal_plan/autocorrelation_z.png")
        
        # #########################
        # # Decorrelation Length ##
        # #########################
        
        # file_int = open(out_figure_path+'split_time/Normal_plan/Lii.dat', 'w')
        
        # if len(zp) == 4:
        #     fig_liip_z = make_subplots(rows=1, cols=4, shared_yaxes= True, y_title='$R_{ii}^+(\omega)$', subplot_titles=(f"$z_0={zp[0]:.2f}$", f"$z_0={zp[1]:.2f}$", f"$z_0={zp[2]:.2f}$", f"$z_0={zp[3]:.2f}$"))
        #     fig_liim_z = make_subplots(rows=1, cols=4, shared_yaxes= True, y_title='$R_{ii}^-(\omega)$', subplot_titles=(f"$z_0={zp[0]:.2f}$", f"$z_0={zp[1]:.2f}$", f"$z_0={zp[2]:.2f}$", f"$z_0={zp[3]:.2f}$"))
        # else:
        #     fig_liip_z = make_subplots(rows=2, cols=5, shared_yaxes= True, y_title=('$R_{ii}^+(\omega)$','$R_{ii}^+(\omega)$'), subplot_titles=(f"$z_0={z[0]:.2f}$", f"$z_0={z[1]:.2f}$", f"$z_0={z[2]:.2f}$", f"$z_0={z[3]:.2f}$", f"$z_0={z[4]:.2f}$", f"$z_0={z[5]:.2f}$", f"$z_0={z[6]:.2f}$", f"$z_0={z[7]:.2f}$", f"$z_0={z[8]:.2f}$", f"$z_0={z[9]:.2f}$"))
        #     fig_liim_z = make_subplots(rows=2, cols=5, shared_yaxes= True, y_title=('$R_{ii}^-(\omega)$','$R_{ii}^-(\omega)$'), subplot_titles=(f"$z_0={z[0]:.2f}$", f"$z_0={z[1]:.2f}$", f"$z_0={z[2]:.2f}$", f"$z_0={z[3]:.2f}$", f"$z_0={z[4]:.2f}$", f"$z_0={z[5]:.2f}$", f"$z_0={z[6]:.2f}$", f"$z_0={z[7]:.2f}$", f"$z_0={z[8]:.2f}$", f"$z_0={z[9]:.2f}$"))
        
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
        
        # if split_time == 'Y':
        #     save_figures(fig_liip_z, "split_time/Normal_plan/correlation_length_plus.png")
        #     save_figures(fig_liim_z, "split_time/Normal_plan/correlation_length_minus.png")
        # if split_time == 'n':
        #     save_figures(fig_liip_z, "whole_time/Normal_plan/correlation_length_plus.png")
        #     save_figures(fig_liim_z, "whole_time/Normal_plan/correlation_length_minus.png")
        
        
        # file_int.close()
        
        # elapsed_time = time.time() - start_time
        # minutes, seconds = divmod(elapsed_time, 60)
        
        # print(f'\n Normal plan study done in : {int(minutes)}m {seconds:.2f}s \n')

    
    
    ##=======================================================
                #### VON KARMAN MODEL ####
    ##=======================================================
    
    print("\n========================================")
    print(f"{YELLOW}Von Karman theorical computation{RESET}")
    start_time = time.time()
    print("========================================")
    
    col = 1
    row = 1
    print("\nReading input files ...")
    
    _, x1, x2, _, nt, n1, n2, _, tEnd, _, iprecision, _, _ = read_fpar_extract_plane(fpars_files_streamwise_u1[0])
    nt = nt - 1
    dx = cflow.xlen / n1
    dy = cflow.ylen / n2
    
    fig_vanK = init_figures_vk(zp, ch=chplot)
    fig_int_scale = go.Figure()
    
    sigma_u = np.zeros((len(zp)))
    sigma_v = np.zeros((len(zp)))
    sigma_w = np.zeros((len(zp)))
    
    int_exp1 = np.zeros((len(zp)))
    int_exp2 = np.zeros((len(zp)))
    int_exp3 = np.zeros((len(zp)))
    
    #for zplan in np.arange(0, n2, n2//3, dtype=int):
    for ind, zplan in enumerate(zp_ind):
        
        print("========================================")
        print(f'Von Karman theory for {YELLOW}zp={zp[ind]:.2f}{RESET}')
        print('Plan number:', zplan)
    
        ### experience ###
        ## u1 ##
        print(f"Computing spectra from LES datas for {YELLOW}z={zp[ind]:.2f}{RESET} ...\n")
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u1[zplan])
        datas_u1 = var[1:,:,:]
        Ux = np.mean(np.mean(np.mean(datas_u1[:,:,:], axis=-1), axis=-1))
        datas_u1 = datas_u1 - Ux
        sigma_u[ind] = np.mean(np.mean(np.mean((datas_u1[:,:,:])**2, axis=0), axis=-1))
        
        omega, k, _, phi11_exp = frozen_turbulence(datas_u1, ind, zp, nt, split_time, dt, n1, dx=dx, ch="spectra", scaling='spectrum')
            
        del datas_u1
        
        von_karman_plot(fig_vanK, col, row, k, phi11_exp, name = '$\phi_{11}LES$', color = 'firebrick', symbols='x')
        
        if chplot == 'all':
            if split_time == 'Y':
                save_datas([k, phi11_exp], ['k', 'phi11_exp'], f'split_time/von_karman/spactra_exp_11_{zp[ind]}.dat', 'LES datas spectra 11 streamwise')
            if split_time == 'n':
                save_datas([k, phi11_exp], ['k', 'phi11_exp'], f'whole_time/von_karman/spactra_exp_11_{zp[ind]}.dat', 'LES datas spectra 11 streamwise')
        
        int_exp1[ind] = np.trapz(phi11_exp, k)
        
        del phi11_exp
        del k

        
        ## u2 ##
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u2[zplan])
        nt = nt - 1
        datas_u2 = var[1:,:,:]
        Uy = np.mean(np.mean(np.mean(datas_u2[:,:,:], axis=-1), axis=-1))
        datas_u2 = datas_u2 - Uy
        sigma_v[ind] = np.mean(np.mean(np.mean((datas_u2[:,:,:])**2, axis=0), axis=-1))

        omega, k, _, phi22_exp = frozen_turbulence(datas_u2, ind, zp, nt, split_time, dt, n1, dx=dx, ch="spectra", scaling='spectrum')
        
        del datas_u2
                
        von_karman_plot(fig_vanK, col, row, k, phi22_exp, name = '$\phi_{22}LES$', color = 'midnightblue', symbols='x')
        
        if chplot == 'all':
            if split_time == 'Y':
                save_datas([k, phi22_exp], ['k', 'phi22_exp'], f'split_time/von_karman/spactra_exp_22_{zp[ind]}.dat', 'LES datas spectra 22 streamwise')
            if split_time == 'n':
                save_datas([k, phi22_exp], ['k', 'phi22_exp'], f'whole_time/von_karman/spactra_exp_22_{zp[ind]}.dat', 'LES datas spectra 22 streamwise')
        
        int_exp2[ind] = np.trapz(phi22_exp, k)
            
        del phi22_exp
        del k

        
        ## u3 ##
        _,_,_,var,_,_,_,_,_,_,_,_,_ = read_fpar_extract_plane(fpars_files_streamwise_u3[zplan])
        nt = nt - 1
        datas_u3 = var[1:,:,:]
        Uz = np.mean(np.mean(np.mean(datas_u3[:,:,:], axis=-1), axis=-1))
        datas_u3 = datas_u3 - Uz
        sigma_w[ind] = np.mean(np.mean(np.mean((datas_u3[:,:,:])**2, axis=0), axis=-1))

        omega, k, _, phi33_exp = frozen_turbulence(datas_u3, ind, zp, nt, split_time, dt, n1, dx=dx, ch="spectra", scaling='spectrum')
        
        del datas_u3
        
        von_karman_plot(fig_vanK, col, row, k, phi33_exp, name = '$\phi_{33}LES$', color = 'darkgreen', symbols='x')
        
        if chplot == 'all':
            if split_time == 'Y':
                save_datas([k, phi33_exp], ['k', 'phi33_exp'], f'split_time/von_karman/spactra_exp_33_{zp[ind]}.dat', 'LES datas spectra 33 streamwise')
            if split_time == 'n':
                save_datas([k, phi33_exp], ['k', 'phi33_exp'], f'whole_time/von_karman/spactra_exp_33_{zp[ind]}.dat', 'LES datas spectra 33 streamwise')
        
        int_exp3[ind] = np.trapz(phi33_exp, k)
            
        del phi33_exp
        del k
        
        print('int_exp 1:', int_exp1[ind])
        print('int_exp 2:', int_exp2[ind])
        print('int_exp 3:', int_exp3[ind])
        
        col +=1
        if zplan == 4:
            row +=1
            col = 1
    
    
    
    ## Theoric VK ##
        
    base_path = os.path.join(os.path.dirname(__file__), "../")
    z_rans, _, _, _, _, _, _, _, ko2_tke, ko2_omega = read_rans(os.path.join(base_path,rans_path))
    zp_rans = z_rans * cflow.ut / cflow.nu
    
    print("========================================")
    print(f"Computing length scale ...\n")
    

    Le = L(0.519, ko2_tke, ko2_omega) #length scale function of z
    
    kmina = np.zeros((ko2_tke.shape[0]))
    kmaxa = np.zeros((ko2_tke.shape[0]))
    for i in range(ko2_tke.shape[0]):
        eps = ko2_omega[i] * 0.09 * ko2_tke[i]
        lmin = ((cflow.nu**3)/eps)**(1/4) ; lmax = lmin/(cflow.re**(-0.75))
        kmin = 2*np.pi/lmax ; kmax = 2*np.pi/lmin
        kmina[i] = kmin ; kmaxa[i] = kmax


    kstart = np.floor(np.log10(np.min(kmina[1:])))
    kend = np.ceil(np.log10(np.max(kmaxa)))

    kx = np.logspace(int(kstart),int(kend)+1, 1000*int(abs(kstart)+abs(kend)+1))
    kx = np.append(0,kx)
    norm_kx = 2*np.pi*np.linalg.norm(kx)
    # kx = read_moser('input/Moser/Re395/Spectra/Streamwise/chan395.xspec.10.txt')[0]
    
    print('Le shape:', Le.shape)
    print('kx shape:', kx.shape)
    
    int_theo1 = np.zeros((len(zp)))
    int_theo2 = np.zeros((len(zp)))
    int_theo3 = np.zeros((len(zp)))
    

    col = 1
    row = 1
    for ind, zplan in enumerate(zp_ind):
        
        ind_z = np.argsort(np.abs(zp_rans - zp[ind]))[0]
        print('ind_z:', ind_z)
        
        phi11 = phi_11(kx, 1.0, sigma_u[ind], Le[ind_z])
        phi22 = phi_22(kx, 1.0, sigma_v[ind], Le[ind_z])
        phi33 = phi_22(kx, 1.0, sigma_w[ind], Le[ind_z])
        
        print('phi11 shape:', phi11.shape)
        
        int_theo1[ind] = np.trapz(phi11, kx)
        int_theo2[ind] = np.trapz(phi22, kx)
        int_theo3[ind] = np.trapz(phi33, kx)
        
        
        print('int_theo1:', int_theo1[ind])
        print('int_theo2:', int_theo2[ind])
        print('int_theo3:', int_theo3[ind])
        
        if model_name == "WRLES_Retau395":
            kx_moser = read_moser(f'input/Moser/Re395/Spectra/Streamwise/chan395.xspec.{zp[ind]}.txt')[0]
            E_uu_moser = read_moser(f'input/Moser/Re395/Spectra/Streamwise/chan395.xspec.{zp[ind]}.txt')[1]
            E_ww_moser = read_moser(f'input/Moser/Re395/Spectra/Streamwise/chan395.xspec.{zp[ind]}.txt')[2]
            E_vv_moser = read_moser(f'input/Moser/Re395/Spectra/Streamwise/chan395.xspec.{zp[ind]}.txt')[3]
            
        u_tau_moser = 0.0572
        E_uu_moser *= u_tau_moser**2
        E_vv_moser *= u_tau_moser**2
        E_ww_moser *= u_tau_moser**2
        
        
        von_karman_plot(fig_vanK, col, row, kx, phi11, name = '$\phi^{1}_{11}\\text{VK}$', color = 'firebrick', line='yes' )
        von_karman_plot(fig_vanK, col, row, kx, phi22, name = '$\phi^{1}_{22}\\text{VK}$', color = 'midnightblue', line='yes' )
        von_karman_plot(fig_vanK, col, row, kx, phi33, name = '$\phi^{1}_{33}\\text{VK}$', color = 'darkgreen', line='yes' )
        
        von_karman_plot(fig_vanK, col, row, kx_moser, E_uu_moser, name = '$\phi^{1}_{11}\\text{DNS}$', color = 'firebrick', symbols='diamond-open')
        von_karman_plot(fig_vanK, col, row, kx_moser, E_vv_moser, name = '$\phi^{1}_{22}\\text{DNS}$', color = 'midnightblue', symbols='diamond-open')
        von_karman_plot(fig_vanK, col, row, kx_moser, E_ww_moser, name = '$\phi^{1}_{33}\\text{DNS}$', color = 'darkgreen', symbols='diamond-open')
    
        # new_phi11 = kc*phi11 # * int_exp1**2 / int_theo1**2
        # new_phi22 = kc*phi22 # * int_exp2**2 / int_theo2**2
        
        # slop_test = np.log(np.abs(phi11[-6] - phi11[-4])) / np.log(np.abs(kx[-6] - kx[-4]))
        # slop_test_2 = np.log(np.abs(phi11[-6] - phi11[-4]) / np.abs(kx[-6] - kx[-4]))
        # print('slop_test:',slop_test)
        # print('slop_test_2:',slop_test_2)
        # print('-5/3:', -5./3)
    
        
        if chplot == 'all':
            if split_time == 'Y':
                save_datas([kx, phi11, phi22, phi33], ['kx', 'phi11', 'phi22', 'phi33'], f'split_time/von_karman/vk_spactra_{zp[ind]}.dat', 'von Karman spectra')
            if split_time == 'n':
                save_datas([kx, phi11, phi22, phi33], ['kx', 'phi11', 'phi22', 'phi33'], f'whole_time/von_karman/vk_spactra_{zp[ind]}.dat', 'von Karman spectra')
        
        col +=1
        if zplan == 4:
            row +=1
            col = 1
            
    fig_vk_moser = make_subplots(rows=1, cols=3, subplot_titles=('$E_{UU}(k_x)$', '$E_{VV}(k_x)$', '$E_{WW}(k_x)$'), y_title='$E_{ii}(k_x)$')
            
    ind_z_mid = np.argsort(np.abs(zp_rans - 395))[0]
    print('ind_z_mid:', ind_z_mid)
    
    phi_mid_1 = phi_11(kx, 1.0, 1.1*sigma_u[-1], Le[ind_z_mid])
    phi_mid_2 = phi_22(kx, 1.0, 1.1*sigma_v[-1], Le[ind_z_mid])
    phi_mid_3 = phi_22(kx, 1.0, 1.1*sigma_w[-1], Le[ind_z_mid])
    
    if model_name == "WRLES_Retau395":
        kx_moser = read_moser('input/Moser/Re395/Spectra/Streamwise/chan395.xspec.392.txt')[0]
        E_uu_moser = read_moser('input/Moser/Re395/Spectra/Streamwise/chan395.xspec.392.txt')[1]
        E_ww_moser = read_moser('input/Moser/Re395/Spectra/Streamwise/chan395.xspec.392.txt')[2]
        E_vv_moser = read_moser('input/Moser/Re395/Spectra/Streamwise/chan395.xspec.392.txt')[3]
        
        u_tau_moser = 0.0572
        E_uu_moser *= u_tau_moser**2
        E_vv_moser *= u_tau_moser**2
        E_ww_moser *= u_tau_moser**2
        
        fig_vk_moser.add_trace(go.Scatter(x=kx_moser, y=E_uu_moser, name='DNS', mode='lines', line=dict(color='black', width=1)), row=1, col=1)
        fig_vk_moser.add_trace(go.Scatter(x=kx_moser, y=E_vv_moser, showlegend=False, mode='lines', line=dict(color='black', width=1)), row=1, col=2)
        fig_vk_moser.add_trace(go.Scatter(x=kx_moser, y=E_ww_moser, showlegend=False, mode='lines', line=dict(color='black', width=1)), row=1, col=3)
        
        fig_vk_moser.add_trace(go.Scatter(x=kx, y=phi_mid_1, name='von K.', mode='lines+markers', line=dict(color='firebrick', width=1), marker=dict(symbol='circle-open')), row=1, col=1)
        fig_vk_moser.add_trace(go.Scatter(x=kx, y=phi_mid_2, showlegend=False, mode='lines+markers', line=dict(color='firebrick', width=1), marker=dict(symbol='circle-open')), row=1, col=2)
        fig_vk_moser.add_trace(go.Scatter(x=kx, y=phi_mid_3, showlegend=False, mode='lines+markers', line=dict(color='firebrick', width=1), marker=dict(symbol='circle-open')), row=1, col=3)
        
        fig_vk_moser.update_xaxes(title='$k_x$', row=1, col=1, type="log", exponentformat='power')
        fig_vk_moser.update_xaxes(title='$k_x$', row=1, col=2, type="log", exponentformat='power')
        fig_vk_moser.update_xaxes(title='$k_x$', row=1, col=3, type="log", exponentformat='power')
        
        fig_vk_moser.update_layout(height=600, width=900, title=f"Von Karman and LES spectra comparison", font=font, showlegend=True, legend=dict(yanchor="bottom", xanchor="left"))
        
        if split_time == 'Y': 
            save_figures(fig_vk_moser, "split_time/von_karman/vk_DNS_mid.png")
        if split_time == 'n': 
            save_figures(fig_vk_moser, "whole_time/von_karman/vk_DNS_mid.png")
    
    
        
    del kx
    del phi11
    del phi22
    del phi33
    del E_uu_moser
    del E_vv_moser
    del E_ww_moser
    del kx_moser
    
    
    if chplot == 'all':
        if split_time == 'Y':
            save_datas([int_theo1, int_theo2, int_theo3], ['int1', 'int2', 'int3'], f'split_time/von_karman/vk_integral_scale.dat', 'von Karman integral scale')
            save_datas([int_exp1, int_exp2, int_exp3], ['int1', 'int2', 'int3'], f'split_time/von_karman/LES_integral_scale.dat', 'LES integral scale')
        if split_time == 'n':
            save_datas([int_theo1, int_theo2, int_theo3], ['int1', 'int2', 'int3'], f'whole_time/von_karman/vk_integral_scale.dat', 'von Karman integral scale')
            save_datas([int_exp1, int_exp2, int_exp3], ['int1', 'int2', 'int3'], f'whole_time/von_karman/LES_integral_scale.dat', 'LES integral scale')
        
    if split_time == 'Y':    
        if chplot == 'normal':
            fig_vanK.update_layout(height=600, width=900, title=f"Von Karman and LES spectra comparison", font=font, showlegend=True, legend=dict(yanchor="top", xanchor="right"))
            save_figures(fig_vanK, "split_time/von_karman/von_karman_spectra.png")
        if chplot == 'all':
            fig_vanK.update_layout(height=900, width=900, title=f"Von Karman and LES spectra comparison", font=font, showlegend=True, legend=dict(yanchor="top", xanchor="right"))
            save_figures(fig_vanK, "split_time/von_karman/von_karman_spectra_all.png")
    if split_time == 'n':
        if chplot == 'normal':
            fig_vanK.update_layout(height=600, width=900, title=f"Von Karman and LES spectra comparison", font=font, showlegend=True, legend=dict(yanchor="top", xanchor="right"))
            save_figures(fig_vanK, "whole_time/von_karman/von_karman_spectra.png")
        if chplot == 'all':
            fig_vanK.update_layout(height=900, width=900, title=f"Von Karman and LES spectra comparison", font=font, showlegend=True, legend=dict(yanchor="top", xanchor="right"))
            save_figures(fig_vanK, "whole_time/von_karman/von_karman_spectra_all.png")
            
    fig_int_scale.add_trace(go.Scatter(x=int_exp1/int_theo1, y=zp, name="$l_{11}^{LES}/l_{11}^{VK}$", mode= 'lines+markers', line=dict(color='midnightblue', width=3), marker=dict(symbol='circle')))
    fig_int_scale.add_trace(go.Scatter(x=int_exp2/int_theo2, y=zp, name="$l_{22}^{LES}/l_{22}^{VK}$", mode= 'lines+markers', line=dict(color='midnightblue', width=3), marker=dict(symbol='circle')))
    fig_int_scale.add_trace(go.Scatter(x=int_exp3/int_theo3, y=zp, name="$l_{33}^{LES}/l_{33}^{VK}$", mode= 'lines+markers', line=dict(color='midnightblue', width=3), marker=dict(symbol='circle')))
    fig_int_scale.update_xaxes(title='$l_{ii}^{LES}/l_{ii}^{VK}$')
    fig_int_scale.update_yaxes(title='$z^+$')
    
    if split_time == 'Y':
        if chplot == 'normal':
            fig_int_scale.update_layout(height=400, width=600, title=f"Integral length scale ratio", font=font, showlegend=True, legend=dict(yanchor="bottom", xanchor="left"))
            save_figures(fig_int_scale, "split_time/von_karman/integral_lenght_scale.png")
        if chplot == 'all':
            fig_int_scale.update_layout(height=400, width=600, title=f"Integral length scale ratio", font=font, showlegend=True, legend=dict(yanchor="bottom", xanchor="left"))
            save_figures(fig_int_scale, "split_time/von_karman/integral_lenght_scale_all.png")
            
    if split_time == 'n':
        if chplot == 'normal':
            fig_int_scale.update_layout(height=400, width=600, title=f"Integral length scale ratio", font=font, showlegend=True, legend=dict(yanchor="bottom", xanchor="left"))
            save_figures(fig_int_scale, "whole_time/von_karman/integral_lenght_scale.png")
        if chplot == 'all':
            fig_int_scale.update_layout(height=400, width=600, title=f"Integral length scale ratio", font=font, showlegend=True, legend=dict(yanchor="bottom", xanchor="left"))
            save_figures(fig_int_scale, "whole_time/von_karman/integral_lenght_scale_all.png")
        
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    print(f'\n Von Karman study done in : {int(minutes)}m {seconds:.2f}s \n')
        
        
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