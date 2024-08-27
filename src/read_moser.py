import numpy as np
import matplotlib.pyplot as plt

def read_moser(filename):
    
    """
    Read the Moser data :
    Ref : Myoungkyu Lee and Robert D. Moser,  
    Direct numerical simulation of turbulent channel flow up to Re_tau = 5200, 
    2015, Journal of Fluid Mechanics, vol. 774, pp. 395-415
    """
    
    data = np.loadtxt(filename, comments=('#','%'))
    
    # Split the data into a list of arrays, one for each column
    columns = [data[:, i] for i in range(data.shape[1])]
    
    return columns 

yp = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.reystress.txt')[1]
R_uu = read_moser('input/Moser/Re395/RANS_comparison/Moser_chan395.reystress.txt')[2]

