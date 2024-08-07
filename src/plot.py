import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from math import *

font = dict(family="serif",
        size=16,
        color="black")

def init_figures_ft(z, n2, ch='normal'):
    if ch == 'normal':
        fig1 = make_subplots(rows=1, cols=4, shared_yaxes=True, y_title='$\phi(k_x)$', subplot_titles=(f"$z^+={z[0]}$", f"$z^+={z[n2//3]}$", f"$z^+={z[2*n2//3]}$", f"$z^+={z[3*floor(n2//3)]}$"))
        
        fig2 = make_subplots(rows=1, cols=4, shared_yaxes=True, y_title='$R_{UU}(\delta t)$', horizontal_spacing=0.02, subplot_titles=(f"$z^+={z[0]}$", f"$z^+={z[n2//3]}$", f"$z^+={z[2*n2//3]}$", f"$z^+={z[3*floor(n2//3)]}$"))
        
        fig3 = make_subplots(rows=1, cols=4, shared_yaxes=True, y_title='$\delta t$', subplot_titles=(f"$z^+={z[0]}$", f"$z^+={z[n2//3]}$", f"$z^+={z[2*n2//3]}$", f"$z^+={z[3*floor(n2//3)]}$"))

        
    if ch == "all":
        fig1 = make_subplots(rows=2, cols=5, shared_yaxes='rows', row_titles=('$\phi(k_x)$', '$\phi(k_x)$'), vertical_spacing=0.15, subplot_titles=(f"$z^+={z[0]}$", f"$z^+={z[1]}$", f"$z^+={z[2]}$", f"$z^+={z[3]}$", f"$z^+={z[4]}$", f"$z^+={z[5]}$", f"$z^+={z[6]}$", f"$z^+={z[7]}$", f"$z^+={z[8]}$", f"$z^+={z[9]}$"))
        
        fig2 = make_subplots(rows=2, cols=5, shared_yaxes='rows', row_titles=('$R_{UU}(\delta t)$', '$R_{UU}(\delta t)$'), horizontal_spacing=0.02, subplot_titles=(f"$z^+={z[0]}$", f"$z^+={z[1]}$", f"$z^+={z[2]}$", f"$z^+={z[3]}$", f"$z^+={z[4]}$", f"$z^+={z[5]}$", f"$z^+={z[6]}$", f"$z^+={z[7]}$", f"$z^+={z[8]}$", f"$z^+={z[9]}$"))
        
        fig3 = make_subplots(rows=2, cols=5, shared_yaxes='rows', row_titles=('$\delta t$', '$\delta t$'), vertical_spacing=0.15, subplot_titles=(f"$z^+={z[0]}$", f"$z^+={z[1]}$", f"$z^+={z[2]}$", f"$z^+={z[3]}$", f"$z^+={z[4]}$", f"$z^+={z[5]}$", f"$z^+={z[6]}$", f"$z^+={z[7]}$", f"$z^+={z[8]}$", f"$z^+={z[9]}$"))
        
    
    return(fig1, fig2, fig3)


def frozen_turbulence_plot(fig, col, row, omega = None, Uc = None, time_spectra = None, k = None, space_spectra = None, R_time = None, R_space = None, Dt = None, Dx = None, R2d=None, coef=None, delta_x=None, funct=None, r=None, z=None, ind1=None, ind2=None, cpt=None, ch = "spectra"):
    
    if ch == "spectra":
        
        if col == 4 and row == 1:
            fig.add_trace(go.Scatter(x=omega[1:]/Uc, y=time_spectra[1:], mode='lines', name='$F(\omega/Uc)_{UU}$', line=dict(color='midnightblue', width=3)), row=row, col=col)
        
            fig.add_trace(go.Scatter(x=k[1:], y=space_spectra[1:]/Uc, mode='lines', name='$P(k_1)_{UU}/Uc$', line=dict(color='firebrick', width=3)), row=row, col=col)
            
            lin1 = np.logspace(0,3)
            fig.add_trace(go.Scatter(x=lin1, y=lin1**(-5./3), line=dict(color='darkgreen', dash='dash', width=2), name='slop: -5/3'), row=row, col=col)
            
            # Update axis properties
            fig.update_xaxes(range=[0,2], title_text="$k_x$", type="log", exponentformat='power', minexponent= 1, dtick=1, row=row, col=col)
            fig.update_yaxes(range=[-9,0], type="log", exponentformat='power')
            
            
        else:
            fig.add_trace(go.Scatter(x=omega[1:]/Uc, y=time_spectra[1:], mode='lines', line=dict(color='midnightblue', width=3), showlegend=False), row=row, col=col)
            
            fig.add_trace(go.Scatter(x=k[1:], y=space_spectra[1:]/Uc, line=dict(color='firebrick', width=3), showlegend=False), row=row, col=col)
            
            lin1 = np.logspace(0,3)
            fig.add_trace(go.Scatter(x=lin1, y=lin1**(-5./3), line=dict(color='darkgreen', dash='dash', width=2), showlegend=False), row=row, col=col)
            
            # Update axis properties
            fig.update_xaxes(range=[0,2],title_text="$k_x$", type="log", exponentformat='power', minexponent= 1, dtick=1, row=row, col=col)
            fig.update_yaxes(range=[-9,0], type="log", exponentformat='power')
            
        
    if ch == "corr":
        
        if col == 4 and row == 1:
            fig.add_trace(go.Scatter(x=Dt, y=R_time, name='time', line=dict(color='midnightblue', width=3)), row=row, col=col)
            # R_space_interpol = np.interp(Dt*Uc, Dx, R_space)
            # fig.add_trace(go.Scatter(x=Dt, y=R_space_interpol, name='space', line=dict(color='firebrick', width=3)), row=row, col=col)
            fig.add_trace(go.Scatter(x=Dx/Uc, y=R_space, name='space', line=dict(color='firebrick', width=3)), row=row, col=col)
            
            # fit_space = np.polyfit(Dt, np.log(np.abs(R_space_interpol)), 1)
            fit_space = np.polyfit(Dx/Uc, np.log(np.abs(R_space)), 1)
            fit_time = np.polyfit(Dt, np.log(np.abs(R_time)), 1)
            curve_fit_space = np.exp(fit_space[0]*Dt)
            curve_fit_time = np.exp(fit_time[0]*Dt)
            # fig.add_trace(go.Scatter(x=Dt, y=curve_fit_space, line=dict(color='firebrick', dash='dash', width=3), name='$\\text{fit space } (\Gamma)$'), row=row, col=col)
            fig.add_trace(go.Scatter(x=Dx/Uc, y=curve_fit_space, line=dict(color='firebrick', dash='dash', width=3), name='$\\text{fit space } (\Gamma)$'), row=row, col=col)
            fig.add_trace(go.Scatter(x=Dt, y=curve_fit_time, line=dict(color='midnightblue', dash='dash', width=3), name='$\\text{fit time } (\gamma)$'), row=row, col=col)
            fig.add_annotation(xanchor='left',x=1.5, yanchor='bottom', y=1, text=f'$\Gamma \simeq {np.abs(fit_space[0]):.2f}$', font=font, showarrow=False, row=row, col=col)
            fig.add_annotation(xanchor='left',x=1.5, yanchor='bottom', y=0.9, text=f'$\gamma \simeq {np.abs(fit_time[0]):.2f}$', font=font, showarrow=False, row=row, col=col)
            
            print('Gamma:', np.abs(fit_space[0]))
            print('gamma:', np.abs(fit_time[0]))
            
            # Update axis properties
            fig.update_xaxes(title_text="$\delta t$")
            #fig.update_yaxes(type="log", exponentformat='power')
            
        else:
            fig.add_trace(go.Scatter(x=Dt, y=R_time, line=dict(color='midnightblue', width=3), showlegend=False), row=row, col=col)
            # R_space_interpol = np.interp(Dt*Uc, Dx, R_space)
            # fig.add_trace(go.Scatter(x=Dt, y=R_space_interpol, line=dict(color='firebrick', width=3), showlegend=False), row=row, col=col)
            fig.add_trace(go.Scatter(x=Dx/Uc, y=R_space, line=dict(color='firebrick', width=3), showlegend=False), row=row, col=col)
            
            # fit_space = np.polyfit(Dt, np.log(np.abs(R_space_interpol)), 1)
            fit_space = np.polyfit(Dx/Uc, np.log(np.abs(R_space)), 1)
            fit_time = np.polyfit(Dt, np.log(np.abs(R_time)), 1)
            curve_fit_space = np.exp(fit_space[0]*Dt)
            curve_fit_time = np.exp(fit_time[0]*Dt)
            # fig.add_trace(go.Scatter(x=Dt, y=curve_fit_space, line=dict(color='firebrick', dash='dash', width=3), name='$y=\\frac{1}{U_c}e^{-\Gamma\delta t*U_c}$', showlegend=False), row=row, col=col)
            fig.add_trace(go.Scatter(x=Dx/Uc, y=curve_fit_space, line=dict(color='firebrick', dash='dash', width=3), name='$y=\\frac{1}{U_c}e^{-\Gamma\delta t*U_c}$', showlegend=False), row=row, col=col)
            fig.add_trace(go.Scatter(x=Dt, y=curve_fit_time, line=dict(color='midnightblue', dash='dash', width=3), name='$y=\\frac{1}{U_c}e^{-\gamma\delta t}$', showlegend=False), row=row, col=col)
            fig.add_annotation(xanchor='left',x=1.5, yanchor='bottom', y=1, text=f'$\Gamma \simeq {np.abs(fit_space[0]):.2f}$', font=font, showarrow=False, row=row, col=col)
            fig.add_annotation(xanchor='left',x=1.5, yanchor='bottom', y=0.9, text=f'$\gamma \simeq {np.abs(fit_time[0]):.2f}$', font=font, showarrow=False, row=row, col=col)
            
            print('Gamma:', np.abs(fit_space[0]))
            print('gamma:', np.abs(fit_time[0]))
            
            # Update axis properties
            fig.update_xaxes(title_text="$\delta t$", row=row, col=col)
            #fig.update_yaxes(type="log", exponentformat='power')
            

    if ch == 'corr2d':
        
        if col == 1 and row == 1:
            fig.add_trace(go.Contour(x=Dx, y=Dt, z=R2d, contours= dict(start = 0.2, end = 1, size=0.1), contours_coloring='lines', line_width=2, name='$R_{UU}(\delta t,\delta x)$'), row=row, col=col)
            fig.add_trace(go.Scatter(x=Dx, y=coef[0]*Dx, line=dict(color='royalblue', dash='dash'), showlegend=False), row=row, col=col)
            
            fig.update_yaxes(range=[-7,7], row=row, col=col)
            fig.update_xaxes(title_text="$\delta x$", row=row, col=col, range=[-np.pi,np.pi])
            
            fig.add_annotation(x=-1.5, y=4, text=f'$Uc\simeq{1./coef[0]:.2f}$', showarrow = False, row=row, col=col)
            
            
            
        else:
            fig.add_trace(go.Contour(x=Dx, y=Dt, z=R2d, contours= dict(start = 0.2, end = 1, size=0.1), contours_coloring='lines', line_width=2, showlegend=False, showscale=False), row=row, col=col)
            fig.add_trace(go.Scatter(x=Dx, y=coef[0]*Dx, line=dict(color='royalblue', dash='dash'), showlegend=False), row=row, col=col)
            
            fig.update_yaxes(range=[-7,7], row=row, col=col)
            fig.update_xaxes(row=row, col=col, title_text="$\delta x$", range=[-np.pi,np.pi])
            
            fig.add_annotation(x=-1.5, y=4, text=f'$Uc\simeq{1./coef[0]:.2f}$', showarrow = False, row=row, col=col)
            
            
    if ch == 'w_gamma':
        
        if col == 4 and row == 1:
            slop1 = np.polyfit(omega[2:ind2], funct[2:ind2], 1)
            slop2 = np.polyfit(omega[ind2:ind1], funct[ind2:ind1], 1)
            if delta_x == 1:
                
                fig.add_trace(go.Scatter(x=omega[2:ind2], y=slop1[0]*omega[2:ind2]+5, mode='lines', showlegend=False, line=dict(color='firebrick', dash='dash', width=2)), row=1, col=col)
                
                fig.add_annotation(x=60, y=10, text=f'$slop:~{slop1[0]:.3f}$', showarrow = False, row=row, col=col)
                
                fig.add_trace(go.Scatter(x=omega[ind2:ind1], y=slop2[0]*omega[ind2:ind1]+10, mode='lines', showlegend=False, line=dict(color='midnightblue', dash='dash', width=2)), row=1, col=col)
                
                fig.add_annotation(x=40, y=30, text=f'$slop:~{slop2[0]:.3f}$', showarrow = False, row=row, col=col)
                
            if cpt == 0:
                fig.add_trace(go.Scatter(x=omega[2:ind1], y=funct[2:ind1], mode='lines', name=f'$r={r:.2f}$', line=dict(color='olive')), row=row, col=col)
            if cpt == 1:
                fig.add_trace(go.Scatter(x=omega[2:ind1], y=funct[2:ind1], mode='lines', name=f'$r={r:.2f}$', line=dict(color='orangered')), row=row, col=col)
            if cpt == 2:
                fig.add_trace(go.Scatter(x=omega[2:ind1], y=funct[2:ind1], mode='lines', name=f'$r={r:.2f}$', line=dict(color='purple')), row=row, col=col)
            if cpt == 3:
                fig.add_trace(go.Scatter(x=omega[2:ind1], y=funct[2:ind1], mode='lines', name=f'$r={r:.2f}$', line=dict(color='mediumvioletred')), row=row, col=col)
            
        else:
            slop1 = np.polyfit(omega[2:ind2], funct[2:ind2], 1)
            slop2 = np.polyfit(omega[ind2:ind1], funct[ind2:ind1], 1)
            if delta_x == 1:
                
                fig.add_trace(go.Scatter(x=omega[2:ind2], y=slop1[0]*omega[2:ind2]+5, mode='lines', showlegend=False, line=dict(color='firebrick', dash='dash', width=2)), row=1, col=col)
                
                fig.add_annotation(x=60, y=10, text=f'$slop:~{slop1[0]:.3f}$', showarrow = False, row=row, col=col)
                
                fig.add_trace(go.Scatter(x=omega[ind2:ind1], y=slop2[0]*omega[ind2:ind1]+7, mode='lines', showlegend=False, line=dict(color='midnightblue', dash='dash', width=2)), row=1, col=col)
                
                fig.add_annotation(x=40, y=30, text=f'$slop:~{slop2[0]:.3f}$', showarrow = False, row=row, col=col)
            
            if cpt == 0:
                fig.add_trace(go.Scatter(x=omega[2:ind1], y=funct[2:ind1], mode='lines', showlegend = False, line=dict(color='olive')), row=row, col=col)
            if cpt == 1:
                fig.add_trace(go.Scatter(x=omega[2:ind1], y=funct[2:ind1], mode='lines', showlegend = False, line=dict(color='orangered')), row=row, col=col)
            if cpt == 2:
                fig.add_trace(go.Scatter(x=omega[2:ind1], y=funct[2:ind1], mode='lines', showlegend = False, line=dict(color='purple')), row=row, col=col)
            if cpt == 3:
                fig.add_trace(go.Scatter(x=omega[2:ind1], y=funct[2:ind1], mode='lines', showlegend = False, line=dict(color='mediumvioletred')), row=row, col=col)
            
        # Update axis properties
        fig.update_xaxes(title_text="$\omega(s^{-1})$", row=row, col=col)
        
        return(slop1, slop2)
            
            
            
    return(None)

def init_figures_gamma(z, n2, ch='normal'):
    
    if ch == 'normal':
        figu1 = make_subplots(rows=1, cols=4, shared_yaxes=True, y_title='$log(\phi(r_1,\omega)e^{-ik_cr_1} / \phi_{UU}(\omega))$', subplot_titles=(f"z={z[0]}", f"z={z[n2//3]}", f"z={z[2*n2//3]}", f"z={z[3*floor(n2//3)]}"))
        figu2 = make_subplots(rows=1, cols=4, shared_yaxes=True, y_title='$log(\phi(r_1,\omega)e^{-ik_cr_1} / \phi_{VV}(\omega))$', subplot_titles=(f"z={z[0]}", f"z={z[n2//3]}", f"z={z[2*n2//3]}", f"z={z[3*floor(n2//3)]}"))
        figu3 = make_subplots(rows=1, cols=4, shared_yaxes=True, y_title='$log(\phi(r_1,\omega)e^{-ik_cr_1} / \phi_{WW}(\omega))$', subplot_titles=(f"z={z[0]}", f"z={z[n2//3]}", f"z={z[2*n2//3]}", f"z={z[3*floor(n2//3)]}"))
        
    if ch == "all":
        figu1 = make_subplots(rows=2, cols=5, shared_yaxes='rows', row_titles=('$log(\phi(r_1,\omega)e^{-ik_cr_1} / \phi_{UU}(\omega))$', '$log(\phi(r_1,\omega)e^{-ik_cr_1} / \phi_{UU}(\omega))$'), subplot_titles=(f"$z^+={z[0]}$", f"$z^+={z[1]}$", f"$z^+={z[2]}$", f"$z^+={z[3]:.2f}$", f"$z^+={z[4]}$", f"$z^+={z[5]}$", f"$z^+={z[6]}$", f"$z^+={z[7]}$", f"$z^+={z[8]}$", f"$z^+={z[9]}$"))
        figu2 = make_subplots(rows=2, cols=5, shared_yaxes='rows', row_titles=('$log(\phi(r_1,\omega)e^{-ik_cr_1} / \phi_{UU}(\omega))$', '$log(\phi(r_1,\omega)e^{-ik_cr_1} / \phi_{UU}(\omega))$'), subplot_titles=(f"$z^+={z[0]}$", f"$z^+={z[1]}$", f"$z^+={z[2]}$", f"$z^+={z[3]}$", f"$z^+={z[4]}$", f"$z^+={z[5]}$", f"$z^+={z[6]}$", f"$z^+={z[7]}$", f"$z^+={z[8]}$", f"$z^+={z[9]}$"))
        figu3 = make_subplots(rows=2, cols=5, shared_yaxes='rows', row_titles=('$log(\phi(r_1,\omega)e^{-ik_cr_1} / \phi_{UU}(\omega))$', '$log(\phi(r_1,\omega)e^{-ik_cr_1} / \phi_{UU}(\omega))$'), subplot_titles=(f"$z^+={z[0]}$", f"$z^+={z[1]}$", f"$z^+={z[2]}$", f"$z^+={z[3]}$", f"$z^+={z[4]}$", f"$z^+={z[5]}$", f"$z^+={z[6]}$", f"$z^+={z[7]}$", f"$z^+={z[8]}$", f"$z^+={z[9]}$"))
        
    return(figu1, figu2, figu3)


def init_figures_sc(z, n2, ch='normal'):
    
    if ch == 'normal':
        fig = make_subplots(rows=1, cols=4, shared_yaxes=True, y_title='$R(\delta x)$', subplot_titles=(f"z={z[0]:.2f}", f"z={z[n2//3]:.2f}", f"z={z[2*n2//3]:.2f}", f"z={z[3*floor(n2//3)]:.2f}"))
        
    if ch == "all":
        fig = make_subplots(rows=2, cols=5, shared_yaxes='rows', row_titles=('$R(\delta x)$', '$\R(\delta x)$'), vertical_spacing=0.15, subplot_titles=(f"$z^+={z[0]:.2f}$", f"$z^+={z[1]:.2f}$", f"$z^+={z[2]:.2f}$", f"$z^+={z[3]:.2f}$", f"$z^+={z[4]:.2f}$", f"$z^+={z[5]:.2f}$", f"$z^+={z[6]:.2f}$", f"$z^+={z[7]:.2f}$", f"$z^+={z[8]:.2f}$", f"$z^+={z[9]:.2f}$"))
    
    return(fig)

def init_figures_vk(z, n2, ch='normal'):
    
    if ch == 'normal':
        fig = make_subplots(rows=1, cols=4, shared_yaxes=True, y_title='$k_c.\phi_{ij}(k_c)$', subplot_titles=(f"z={z[0]:.2f}", f"z={z[n2//3]:.2f}", f"z={z[2*n2//3]:.2f}", f"z={z[3*floor(n2//3)]:.2f}"))
        
    if ch == "all":
        fig = make_subplots(rows=2, cols=5, shared_yaxes='rows', row_titles=('$k_c.\phi_{ij}(k_c)$', '$k_c.\phi_{ij}(k_c)$'), vertical_spacing=0.15, subplot_titles=(f"$z^+={z[0]:.2f}$", f"$z^+={z[1]:.2f}$", f"$z^+={z[2]:.2f}$", f"$z^+={z[3]:.2f}$", f"$z^+={z[4]:.2f}$", f"$z^+={z[5]:.2f}$", f"$z^+={z[6]:.2f}$", f"$z^+={z[7]:.2f}$", f"$z^+={z[8]:.2f}$", f"$z^+={z[9]:.2f}$"))
    
    return(fig)

def von_karman_plot(fig, col, row, kc, phi, name = 'corr', color = 'firebrick', symbols='circle'):
    
    if col == 4 and row == 1:
        fig.add_trace(go.Scatter(x=kc, y=phi, name=name, mode= 'lines+markers', line=dict(color=color, width=3), marker=dict(symbol=symbols)), row=row, col=col)
        fig.update_xaxes(title='$k_c(m.s^{-2})$', type="log", exponentformat='power', row=row, col=col)
        fig.update_yaxes(type="log", exponentformat='power', row=row, col=col)
        
    else:
        fig.add_trace(go.Scatter(x=kc, y=phi, name=name, mode= 'lines+markers', line=dict(color=color, width=3), marker=dict(symbol=symbols), showlegend=False), row=row, col=col)
        fig.update_xaxes(title='$k_c(m.s^{-2})$', type="log", exponentformat='power', row=row, col=col)
        fig.update_yaxes(type="log", exponentformat='power', row=row, col=col)
        
    return(None)



def space_correlation_plot(fig, col, row, Dx, corr, name = 'corr', color = 'firebrick', axis = 'streamwise'):
    
    if col == 4 and row == 1:
        fig.add_trace(go.Scatter(x=Dx, y=corr, name=name, line=dict(color=color, width=3)), row=row, col=col)
        if axis == 'streamwise':
            fig.update_xaxes(row=row, col=col, title_text="$\delta x$")
        if axis == 'spanwise':
            fig.update_xaxes(row=row, col=col, title_text="$\delta y$")
        if axis == 'wallnormal':
            fig.update_xaxes(row=row, col=col, title_text="$\delta z$")
        
    else:
        fig.add_trace(go.Scatter(x=Dx, y=corr, name=name, line=dict(color=color, width=3), showlegend=False), row=row, col=col)
        if axis == 'streamwise':
            fig.update_xaxes(row=row, col=col, title_text="$\delta x$")
        if axis == 'spanwise':
            fig.update_xaxes(row=row, col=col, title_text="$\delta y$")
        if axis == 'wallnormal':
            fig.update_xaxes(row=row, col=col, title_text="$\delta z$")
        
    return(None)


def save_figures(fig, path):
    fig.write_image(path)
    return(None)