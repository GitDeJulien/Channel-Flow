import numpy as np

with open('input/chan_rans_mean.dat', 'r') as f:
    y = []
    ko2_omega = []
    ko2_tke = []
    normal = []
    pressure = []
    u_velocity = []
    v_velocity = []
    vis = []
    vorticity = []
    w_velocity = []
    
    for line in f:
        if "#" in line:
            continue
        data = line.split()
        y.append(float(data[0]))
        u_velocity.append(float(data[1]))
        v_velocity.append(float(data[2]))
        w_velocity.append(float(data[3]))
        pressure.append(float(data[4]))
        vis.append(float(data[5]))
        vorticity.append(float(data[6]))
        normal.append(float(data[7]))
        ko2_tke.append(float(data[8]))
        ko2_omega.append(float(data[9]))

y= np.array(y)
u_velocity = np.array(u_velocity)
v_velocity = np.array(v_velocity)
w_velocity = np.array(w_velocity)
pressure = np.array(pressure)
vis = np.array(vis)
vorticity = np.array(vorticity)
normal = np.array(normal)
ko2_tke = np.array(ko2_tke)
ko2_omega = np.array(ko2_omega)









