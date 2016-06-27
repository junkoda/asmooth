asmooth
=======

Adaptively Smoothed Fields for Cosmological Simulations

## Required inputs

### snapshot files of cubep3m simulation

- `node``<i>`/`<redshift>`xv`<i>`.dat
- `node``<i>`/`<redshift>`zip0_`<i>`.dat
- `node``<i>`/`<redshift>`zip1_`<i>`.dat
- `node``<i>`/`<redshift>`zip2_`<i>`.dat
- `node``<i>`/`<redshift>`zip3_`<i>`.dat
- `node``<i>`/`<redshift>`halo_`<i>`.dat

### redshift.txt

The list of redshifts should be in a text file "redshifts.txt", one redshift
per line.

```
100.000
10.000
```

### command-line arguments

- `-nc_node_dim`: a parameter necessary to read snapshot.

- `-node_dir results/node`: This is the directory name of the `node` in `snapshot files of cubep3m simulation` above.


- `-nc <nc>`: This is the number of output cells per dimension *per node*.
          You can assign comma sperated numbers (no space in between).
          
- `-allocate <Mbytes>`: Amount memory allocated for particle in Mbytes
          about (`np_local*36*1.1^3`) byte is necessary (1.1^3 for 5% buffer).
          
###