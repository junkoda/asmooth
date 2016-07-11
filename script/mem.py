#
# Read np_local from headers and prints Mbytes necessary for particles
#
#

# Case 1: Estimate RAM necessary for -allocate <value>
#      -allocate option in clumping_tree is Mbytes allocated for particles
#      python3 mem.py --allocate=4096
#
# Case 2: Estimate minimum -allocate necessary for np_local_max
#      np_local_max is maximum number of particles in one MPI node
#
#      python3 mem.py --np-max=65638257 --buffer-factor=0.025
#
# Case 3: Read np_max from files
#      --redshift: redshift in filename
#      --nc-node:  number of MPI nodes per dim; n_files = nc_node**3
#
#      python3 mem.py --redshift=4.000 --nc-node=12 --bufer-factor=0.025
#      Files reading are: results/node<i>/<redshift>zip0_<i>.dat
#
from optparse import OptionParser
import struct

#
# Functions
#
def read_np_max(nfile):
    np_local = []
    for i in range(nfile):
        filename = 'results/node%d/%szip0_%d.dat' % (i, redshift, i)
        with open(filename, 'rb') as f:
            np_local.append(struct.unpack('i', f.read(4))[0])

            # For big endian data on little endian machine
            #np_local.append(struct.unpack('>L', f.read(4))[0])
    return max(np_local)

def mem_kdtree(np):
    nleaf = 1
    while 16*nleaf < np:
        nleaf *= 2

    nalloc = 2*nleaf - 1

    return 20*nalloc

#
# default values
#
nc_node = 12
redshift = '6.000'
buffer_factor = 0.05

#
# Commandline options
#
parser = OptionParser()
parser.add_option("--redshift", dest="redshift",
                  help="redshift of reading np_local")
parser.add_option("--nc-node", dest="nc_node",
                  help="number of MPI nodes per dim; nfile = nc-node**3")
parser.add_option("--np-max", dest="np_max",
                  help="np_local_max if you already know and need not to read from file")
parser.add_option("--buffer-factor", dest="buffer_factor",
                  help="buffer_factor in cluming_tree")
parser.add_option("--allocate", dest="allocate",
                  help="-allocate option in clumping_tree (Mbytes for particles)")


(opt, args) = parser.parse_args()


if opt.redshift:
    redshift = opt.redshift

if opt.allocate:
    allocate = int(opt.allocate)
else:
    # If --allocate is not given, compute 'allocate' from np_max
    # using --buffer-factor and np_max
    if opt.np_max:
        np_max = int(opt.np_max)
    else:
        # If --np_max is not give, read np_max from files
        # using --redshift and --nc-node
        if opt.nc_node:
            nc_node = int(opt.nc_node)

        nfile = nc_node**3
        print('nc_node = %d, %d files' % (nc_node, nfile))

        print('Reading np_locals at redshift %s' % redshift)
        np_max = read_np_max(nfile)

    print('computes -allcate from max np_local = %d' % np_max)

    if opt.buffer_factor:
        buffer_factor = float(opt.buffer_factor)
    print('buffer_factor = %f' % buffer_factor)

    mega = 1024**2
    size_particle = 4*9
    mem_particle = np_max*size_particle/mega*(1.0 + 2.0*buffer_factor)**3

    print('-alloate %d at least, for particles' % mem_particle)

    allocate = mem_particle

np = allocate*mega // size_particle
mem_tree = mem_kdtree(np)

print('RAM requirement for -allocate %d' % allocate)
print('RAM necessary for `temp` = %d' % (4*np))
print('RAM necessary for kdtree = %d' % mem_tree)
mem_total = np*(size_particle + 4) + mem_tree

print('RAM necessary per MPI node %d (%d mbytes)' % (mem_total, mem_total/1000**2))



