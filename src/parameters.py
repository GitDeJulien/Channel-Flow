""" Parameters """

import numpy as np

model = 'WMLES_retau1000'

###################################################################
#################### WRLES Ret = 395 ##############################
###################################################################

if model == 'WRLES_retau395':

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
    ch_plot = "all"              # calculation for all 10 plans "all" or only 4 plans "normal"
    zplus = [5, 20, 40, 60, 80, 98, 151, 199, 251, 302]  # adimentianalized height of plans (z+)

    ## in path ## path to the directory contening all the datas (on the external disc here)
    in_path = '/media/julien/Verbatim/julien/channel_wrles_retau395/' 

    ## out path ## path to the directory that will contain all the results figures
    out_path = 'output/channel_wrles_retau395/' 



###################################################################
#################### WMLES Ret = 1000 #############################
###################################################################

if model == 'WMLES_retau1000':
    
    xlen = 4*np.pi #m
    ylen = 1.5*np.pi   #m
    zlen = 2       #m
    rho  = 1       #kg/m
    uinf = 1       #m/s
    re   = 13800   #_
    ret  = 968.55     #_
    tStart = 1.5500099999987353E+03 #starter time of the collected datas
    dt = 0.01                       # time step of the simulation
    split_t = int(2**10)            # size of time interval in time splited mode
    ch_plot = "all"              # calculation for all 10 plans "all" or only 4 plans "normal"
    zplus = [5, 20, 40, 60, 80, 98, 151, 199, 251, 302]  # adimentianalized height of plans (z+)

    in_path = '/media/julien/Verbatim/julien/channel_wmles_retau1000/'

    out_path = 'output/channel_wmles_retau1000/'


###################################################################
#################### WRLES Ret = 1000 #############################
###################################################################

if model == 'WRLES_retau1000':
    
    xlen = 4*np.pi #m
    ylen = 1.5*np.pi   #m
    zlen = 2       #m
    rho  = 1       #kg/m
    uinf = 1       #m/s
    re   = 13800   #_
    ret  = 980.86     #_
    tStart = 1.5500099999987353E+03 #starter time of the collected datas
    dt = 0.01                       # time step of the simulation
    split_t = int(2**10)            # size of time interval in time splited mode
    ch_plot = "all"              # calculation for all 10 plans "all" or only 4 plans "normal"
    zplus = [5, 20, 40, 60, 80, 98, 151, 199, 251, 302]  # adimentianalized height of plans (z+)

    in_path = '/media/julien/Verbatim/julien/channel_wrles_retau1000/'

    out_path = 'output/channel_wrles_retau1000/'
