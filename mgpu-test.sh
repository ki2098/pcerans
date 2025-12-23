#!/bin/sh

#PJM -L rscgrp=b-batch
#PJM -L gpu=4
#PJM -L elapse=00:30:00
#PJM -j

MPIEXEC=~/.julia/bin/mpiexecjl

module load julia

$MPIEXEC -n 4 -outfile-pattern mgpu-test.%r.log julia scripts/u-mc-mgpu.jl u-sampling-setup.json 100 --skip-history 2>&1