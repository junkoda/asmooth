#
# Visualize density mesh data
#
# Usage:
#     python3 mesh.py <redshift>
#
# Inputs:
#     dir/<redshift>n_all.dat ntot_all.dat
#
# Output:
#     ./<redshift>.n
#

import sys, os
import numpy as np
import struct
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Configuration:
#     dir: filename = <dir>/<redshift>n_all.dat
#     nmean: mean density per grid = 8*mean number of particles per grid

nmean = 8*(4032/600)**3 # for 4032^3 particles in a 600^3 grid
dir = 'global/so/nc600'


if sys.argv.length == 1:
    print('python3 mesh.py 0.000')
    os.exit(1)

redshift = sys.argv[1]


#
# n_tot
#
filename = '%s/%sn_all.dat' % (dir, redshift)

with open(filename, 'rb') as f:
    nc = struct.unpack('iii', f.read(12))
    mesh = np.fromfile(f, np.float32, nc[0]*nc[1]*nc[2])
    mesh = np.reshape(mesh, nc, order='F')

mesh2d = np.log10(mesh[:,0,:]/nmean)

plt.figure(figsize = [12, 6])
plt.subplot(1, 2, 1)

plt.imshow(mesh2d, clim=(-0.5, 3.0))
plt.colorbar()

#
# n_tot
#
filename = '%s/%sntot_all.dat' % (dir, redshift)

with open(filename, 'rb') as f:
    nc = struct.unpack('iii', f.read(12))
    mesh = np.fromfile(f, np.float32, nc[0]*nc[1]*nc[2])
    mesh = np.reshape(mesh, nc, order='F')

mesh2d = np.log10(mesh[:,0,:]/nmean)

plt.subplot(1, 2, 2)
plt.imshow(mesh2d, clim=(-0.5, 3.0))
plt.colorbar()


plt.savefig('%sn.png' % redshift)
