import numpy as np
import time 
import sys
import glob

from read_fpar import *
from parameters import *
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
    print("=========== Reading input files ===========")
    zp = [5, 20, 40, 60, 80, 98, 151, 199, 251, 302, 392]
    z = []
    for he in zp:
        z.append(he * cflow.nu / cflow.ut)
    
    print("z:", z)
    
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
    
    
    
    
    
    
if __name__ == "__main__":
    main()