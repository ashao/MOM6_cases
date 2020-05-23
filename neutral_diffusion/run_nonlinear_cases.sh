#!/bin/bash

cases='teos10_layer teos10_full teos10_int_drho teos10_linear'
#cases='teos10_layer teos10_linear'
cases='teos10_neutral'

module load intel_compiler hdf5/intel-1.8.20 netcdf/intel-4.6.1 netcdf-fortran/intel-4.4.4 openmpi/intel-4.0.0
ulimit -s unlimited

for c in $cases; do
    pushd $c
    for nlev in {10..200..10}; do
#    for nlev in 5; do
        echo "Case: $c Levels:$nlev"
        sed -i "s/^.*NK=.*$/#override NK=${nlev}/" MOM_override
        sed -i "s/^.*DAYMAX=.*$/#override DAYMAX=0.0625/" MOM_override
        sed -i 's/^.*DIAG_COORD_DEF_Z.*$/DIAG_COORD_DEF_Z="UNIFORM:100"/' MOM_override
        ./MOM6 &> run.$nlev.log
        mkdir nk$nlev
        mv prog*.nc Initial_state.nc nk$nlev/
    done
    popd
done
