#!/usr/bin/bash

#PJM -L rscgrp=b-batch
#PJM -L gpu=1
#PJM -L elapse=00:30:00
#PJM -j
#PJM -o test-run.log

module load cuda nvidia julia

julia scripts/cuda-hpc-setup.jl 12.2.2

julia u-det-run.jl --shutup

julia u-nipce-sampling.jl --shutup

julia u-mc-sampling.jl 100 --shutup