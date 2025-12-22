#!/usr/bin/bash

julia scripts/u-det-run.jl u-sampling-setup.json --skip-history

julia scripts/u-nipce-sampling.jl u-sampling-setup.json 5 10 --skip-history

sleep 60

julia scripts/u-mc-sampling.jl u-sampling-setup.json 10 --skip-history

sleep 60

julia scripts/u-mc-sampling.jl u-sampling-setup.json 100 --skip-history

sleep 60

julia scripts/u-mc-sampling.jl u-sampling-setup.json 500 --skip-history

sleep 60

julia scripts/u-mc-sampling.jl u-sampling-setup.json 1000 --skip-history