""" 07/08/2024 Prakyath Pindi Nata """
""" Read binary data files from the LES code """

import numpy as np

def read_fpar_extract_plane(infile):
    # Read extract out plane saved in fparser binary format initially based on script from Mads B.
    with open(infile, 'rb') as fid:
        data_array = np.fromfile(fid, np.int32, 1)
        version = data_array[0]
        if version >= 0:
            # Old version
            print('Old version')
            nt = data_array[0]  # Number of timesteps
            data_array = np.fromfile(fid, np.int32, 3)
            n1 = data_array[0]
            n2 = data_array[1]
            n3 = 1
            iprecision = 2
            ndata = n1 * n2 * n3
            recordlength = np.maximum(32, 4 * ndata * iprecision)  # Choose 32 or ndata*8
            data_array2 = np.fromfile(fid, np.float64, 2)
            tStart = data_array2[0]
            tEnd = data_array2[1]
            shift = np.array([0.0, 0.0, 0.0])
            quaternion = np.array([1.0, 0.0, 0.0, 0.0])
        elif version == -1:
            # New version
            print('New version')
            data_array = np.fromfile(fid, np.int32, 5)
            nt = data_array[0]  # Number of timesteps
            n1 = data_array[1]
            n2 = data_array[2]
            n3 = 1
            iprecision = data_array[4]
            ndata = n1 * n2 * n3
            recordlength = np.maximum(96, 4 * ndata * iprecision)
            data_array2 = np.fromfile(fid, np.float64, 9)
            tStart = data_array2[0]
            tEnd = data_array2[1]
            shift = data_array2[2:5]
            quaternion = data_array2[5:9]
        if iprecision == 1:
            dtype = np.float32
        elif iprecision == 2:
            dtype = np.float64
        t = np.linspace(tStart, tEnd, nt)
        fid.seek(recordlength, 0)
        x1 = np.fromfile(fid, dtype, n1)
        fid.seek(recordlength * 2, 0)
        x2 = np.fromfile(fid, dtype, n2)
        # fid.seek(recordlength * 3,0)
        # x3 = np.fromfile(fid, np.float64, n3)
        var1d = np.zeros((nt * ndata))  # One very large 1d vector
        for i in range(nt):
            fid.seek(recordlength * (4 + i), 0)
            var1d[i * ndata:(i + 1) * ndata] = np.fromfile(fid, dtype, ndata)
    fid.close()
    var = np.reshape(var1d, (n2, n1, nt), order='F')
    # Change dimension to be compliant with CF/COARDS convention, i.e. have time first.
    var = np.moveaxis(var, -1, 0)
    var = np.moveaxis(var, 1, 2)
    return t, x1, x2, var, nt, n1, n2, tStart, tEnd, version, iprecision, shift, quaternion