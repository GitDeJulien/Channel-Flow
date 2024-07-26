import numpy as np

with open('input/chan395_rans2d.dat', 'r') as f:
    ko2_omega = []
    ko2_tke = []
    normal = []
    Points_0 = []
    Points_1 = []
    Points_2 = []
    Points_Magnitude = []
    pressure = []
    u = []
    v = []
    vis = []
    vorticity = []
    w = []
    
    for line in f:
        if "#" in line:
            continue
        data = line.split()
        ko2_omega.append(float(data[1]))
        ko2_tke.append(float(data[2]))
        normal.append(float(data[3]))
        Points_0.append(float(data[4]))
        Points_1.append(float(data[5]))
        Points_2.append(float(data[6]))
        Points_Magnitude.append(float(data[7]))
        pressure.append(float(data[8]))
        u.append(float(data[9]))
        v.append(float(data[10]))
        vis.append(float(data[11]))
        vorticity.append(float(data[12]))
        w.append(float(data[13]))
        
ko2_omega = np.array(ko2_omega)
ko2_tke = np.array(ko2_tke)
normal = np.array(normal)
Points_0 = np.array(Points_0)
Points_1 = np.array(Points_1)
Points_2 = np.array(Points_2)
Points_Magnitude = np.array(Points_Magnitude)
pressure = np.array(pressure)
u = np.array(u)
v = np.array(v)
vis = np.array(vis)
vorticity = np.array(vorticity)
w = np.array(w) 



