Optimizing CFD for memory locality with microdomains
====================================================

This repository contains demo code accompanying a talk for OpenFOAM workshop 2021, 09 June 2021

**Ilya Popov** and **Dmitry Astakhov**

**RnD-ISAN** and **ISTEQ**

Code installation and running
-----------------------------

### Get the code

- Either by downloading an archive
- Or cloning git repository

```
git clone https://github.com/ISTEQ-BV/OpenFOAM_Workshop_2021_Microdomains_Demo.git
```

### Compile the code:

```
cd OpenFOAM_Workshop_2021_Microdomains_Demo
mkdir build
cmake -B build/
cmake --build build/
```

### Run the demo case

Demo case is located in `test_cases/E02_cube_hierarchical_decomposition` subdirectory.

First, we have to create block mesh in the original case:

```
cd original_case
blockMesh
cd ..
```

then we run the script that creates microdomain decomposition

```
./run.sh
```

this can take a lot of time because it uses decomposePar to decompose into microdomains.

Then go into newly created case

```
cd new_case
```

and run the demo solver

```
mpirun -np 20 ../../../build/src/of_transport -parallel -method OpenFOAM
```

there is an additional argument `-method`
that tells which variant of the optimization to use. 
"OpenFOAM" is the baseline OpenFOAM transport euqtion. 
The final version is "microdomains_lsgrad_swapped".

Metric to be used is time of the third time step.

To change number of microdomains adjust `original_case/system/decomposeParDict_multilevel` and rerun run.sh

To change number of MPI domains adjust `original_case/system/decomposeParDict_multilevel`, `original_case/system/decomposeParDict_manual`, and `run.sh`
