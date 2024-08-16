""" Parameters """

import numpy as np

xlen = 2*np.pi #m
ylen = np.pi   #m
zlen = 2       #m
rho  = 1       #kg/m
uinf = 1       #m/s
re   = 13800   #_
ret  = 395     #_
tStart = 1.5500099999987353E+03 #starter time of the collected datas
dt = 0.01                       # time step of the simulation
split_t = int(2**10)            # size of time interval in time splited mode
ch_plot = "normal"              # calculation for all 10 plans "all" or only 4 plans "normal"
in_path = '/media/julien/Verbatim/julien/channel_wrles_retau395/' # path to the directory contening all the datas (on the external disc here)
out_path = 'output/channel_wrles_retau395/' # path to the directory that will contain all the results figures
zplus = [5, 20, 40, 60, 80, 98, 151, 199, 251, 302]  # adimentianalized height of plans (z+)