asmooth
=======

Adaptively Smoothed Fields for Cosmological Simulations

## Required inputs

### snapshot files of cubep3m simulation

#### For xv format

- `node_dir``<i>`/`<redshift>`xv`<i>`.dat

#### For zip format

- `node_dir``<i>`/`<redshift>`zip0_`<i>`.dat
- `node_dir``<i>`/`<redshift>`zip1_`<i>`.dat
- `node_dir``<i>`/`<redshift>`zip2_`<i>`.dat
- `node_dir``<i>`/`<redshift>`zip3_`<i>`.dat
- `node_dir``<i>`/`<redshift>`halo_`<i>`.dat

### redshift.txt

The list of redshifts should be in a text file "redshifts.txt", one redshift
per line, e.g.,

```
10.000
5.000
```

### command-line arguments

- `-node_dir results/node`: This is the directory name of the `node`
  in `snapshot files of cubep3m simulation` above.

- `-allocate <Mbytes>`: Amount memory allocated for particle in Mbytes
          about (`np_local*36*1.1^3`) byte is necessary (1.1^3 for 5% buffer).
          
#### For xv particle format

- `-xv <boxsize>`

`<boxsize>` is the boxsize of a local cube in cubep3m simulation in internal unit, in which mean particle spearation is 2.

For example, for 5472^3 particles in total with 12^3 MPI nodes, boxsize should be 5472/12*2.

#### For zip particle format
- `-nc_node_dim`: a parameter necessary to read snapshot, defined in a simulation code `cubepm.par`.

- `-nc <nc>`: This is the number of output cells per dimension *per node*.
          You can assign comma sperated numbers (no space in between).

- `-omegam 0.308`: This is Omega_m at z=0. This is used if Bryan-Normal virial mass is used in halo file.
          


### example

```bash
$ tree
.
├── asmooth
│   └── clumping_tree
├── redshifts.txt
├── results
│   ├── node0
│   │   ├── 0.000halo0.dat
│   │   ├── 0.000zip0_0.dat
...
│   ├── node7
...
```

```bash
$ mpirun -n 8 ./asmooth/clumping_tree -nc 32,64,128 -nc_node_dim 32 -node_dir results/node -allocate 128
```

## More options

### `-pm_redshift redshifts.txt`

Name of the file for the list of redshifts.
Each line should contain redshifts that is in the filename
e.g. 0.000 for 0.000zip0_0.dat. Default is redshifts.txt

### `-shift shift`

This is the name of the directory containing overall shift values (to correct for random shifts in cubep3m). Default is `shift`.

if you want to shift positions, create a file `<shift>/<redshift>shift.txt` and put one line of space spearated ascii numbers:

```
z shift[0] shift[1] shift[2]
```

The shift numbers are subtracted from particle and halo positions: x[i] = x[i] - shift[i];






