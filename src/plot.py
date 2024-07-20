import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

font = dict(family="serif",
        size=16,
        color="black")

def init_figures(z, n2):
    fig1 = make_subplots(rows=1, cols=4, shared_yaxes=True, y_title='$\phi(k_x)$', subplot_titles=(f"z={z[0]:.2f}", f"z={z[n2//3]:.2f}", f"z={z[2*n2//3]:.2f}", f"z={z[3*n2//3]:.2f}"))
    fig2 = make_subplots(rows=1, cols=4, shared_yaxes=True, y_title='$R_{UU}(\delta t)$', subplot_titles=(f"z={z[0]:.2f}", f"z={z[n2//3]:.2f}", f"z={z[2*n2//3]:.2f}", f"z={z[3*n2//3]:.2f}"))
    return(fig1, fig2)

def frozen_turbulence_plot(fig, col, omega = None, Uc = None, time_spectra = None, k = None, space_spectra = None, R_time = None, R_space = None, Dt = None, Dx = None, ch = "spectra"):
    
    if ch == "spectra":
        
        if col == 4:
            fig.add_trace(go.Scatter(x=omega[1:]/Uc, y=time_spectra[1:], mode='lines', name='$F(\omega/Uc)_{UU}$', line=dict(color='midnightblue', width=3)), row=1, col=col)
        
            fig.add_trace(go.Scatter(x=k[1:], y=space_spectra[1:]/Uc, mode='lines', name='$P(k_1)_{UU}/Uc$', line=dict(color='firebrick', width=3)), row=1, col=col)
            
            lin1 = np.logspace(0,3)
            fig.add_trace(go.Scatter(x=lin1, y=lin1**(-5./3), line=dict(color='darkgreen', dash='dash', width=2), name='slop: -5/3'), row=1, col=col)
            
            # Update axis properties
            fig.update_xaxes(title_text="$k_x$", type="log", exponentformat='power', minexponent= 1, dtick=1, row=1, col=col)
            fig.update_yaxes(type="log", exponentformat='power')
            
        else:
    
            fig.add_trace(go.Scatter(x=omega[1:]/Uc, y=time_spectra[1:], mode='lines', line=dict(color='midnightblue', width=3), showlegend=False), row=1, col=col)
            
            fig.add_trace(go.Scatter(x=k[1:], y=space_spectra[1:]/Uc, line=dict(color='firebrick', width=3), showlegend=False), row=1, col=col)
            
            lin1 = np.logspace(0,3)
            fig.add_trace(go.Scatter(x=lin1, y=lin1**(-5./3), line=dict(color='darkgreen', dash='dash', width=2), showlegend=False), row=1, col=col)
            
            # Update axis properties
            fig.update_xaxes(title_text="$k_x$", type="log", exponentformat='power', minexponent= 1, dtick=1, row=1, col=col)
            fig.update_yaxes(type="log", exponentformat='power')
            
        
    if ch == "corr":
        
        if col == 4:
            fig.add_trace(go.Scatter(x=Dt, y=R_time, name='$time$', line=dict(color='midnightblue', width=3)), row=1, col=col)
            fig.add_trace(go.Scatter(x=Dx/Uc, y=R_space, name='$space$', line=dict(color='firebrick', width=3)), row=1, col=col)
            
            fit_space = np.polyfit(Dx/Uc, np.log(R_space), 1)
            curve_fit_space = np.exp(fit_space[0]*Dx/Uc)
            fig.add_trace(go.Scatter(x=Dx/Uc, y=curve_fit_space, line=dict(color='darkgreen', dash='dash', width=3), name='$y=e^{-\delta x/(T*Uc)}$'), row=1, col=col)
            fig.add_annotation(xanchor='left',x=3.1, yanchor='bottom', y=0.1, text=f'$\gamma = 1/T \simeq {np.abs(fit_space[0]):.2f}$', font=font, showarrow=False)
            
            print('gamma:', np.abs(fit_space[0]))
            
            # Update axis properties
            fig.update_xaxes(title_text="$\delta t$")
            fig.update_yaxes(type="log", exponentformat='power')
            
        else:
            fig.add_trace(go.Scatter(x=Dt, y=R_time, line=dict(color='midnightblue', width=3), showlegend=False), row=1, col=col)
            fig.add_trace(go.Scatter(x=Dx/Uc, y=R_space, line=dict(color='firebrick', width=3), showlegend=False), row=1, col=col)
            
            fit_space = np.polyfit(Dx/Uc, np.log(R_space), 1)
            curve_fit_space = np.exp(fit_space[0]*Dx/Uc)
            fig.add_trace(go.Scatter(x=Dx/Uc, y=curve_fit_space, line=dict(color='darkgreen', dash='dash', width=3), showlegend=False), row=1, col=col)
            if col == 1:
                fig.add_annotation(xanchor='left',x=0.1, yanchor='bottom', y=0.1, text=f'$\gamma = 1/T \simeq {np.abs(fit_space[0]):.2f}$', font=font, showarrow=False)
            if col == 2:
                fig.add_annotation(xanchor='left',x=1.1, yanchor='bottom', y=0.1, text=f'$\gamma = 1/T \simeq {np.abs(fit_space[0]):.2f}$', font=font, showarrow=False)
            if col == 3:
                fig.add_annotation(xanchor='left',x=2.1, yanchor='bottom', y=0.1, text=f'$\gamma = 1/T \simeq {np.abs(fit_space[0]):.2f}$', font=font, showarrow=False)
            
            print('gamma:', np.abs(fit_space[0]))
            
            # Update axis properties
            fig.update_xaxes(title_text="$\delta t$", row=1, col=col)
            fig.update_yaxes(type="log", exponentformat='power')
    return(None)


def save_figures(fig, path):
    fig.write_image(path)
    return(None)