""" Parameters """

import os
from read_parameters import read_parameters
from read_rans import read_rans

base_path = os.path.join(os.path.dirname(__file__), "../")
    
## READ PARAMETERS ##
parameters_path = os.path.join(base_path, "parameters/parameters.dat")
params = read_parameters(parameters_path)
xlen = params['xlen']
ylen = params['ylen']
zlen = params['zlen']
rho = params['rho']
uinf = params['uinf']
re = params['re']
ret = params['ret']
dt = params['dt']
g_limit0 = params['g_limit0']
g_limit1 = params['g_limit1']
zplus = params['zplus']
in_path = params['in_path']
out_figure_path = params['out_figure_path']
rans_path = params['rans_path']
out_data_path = params['out_data_path']
split_time = params['split_time']
split_t = params['nb_split_t']
chplot = params['chplot']

## READ RANS ##
y, u_velocity, v_velocity, w_velocity, pressure, vis, vorticity, normal, ko2_tke, ko2_omega = read_rans(os.path.join(base_path,rans_path))

