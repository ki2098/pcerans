#!/usr/bin/bash

julia u-det-run.jl --shutup

julia u-nipce-sampling.jl --shutup

sleep 60

julia u-mc-sampling.jl 10 --shutup

sleep 60

julia u-mc-sampling.jl 100 --shutup

sleep 60

julia u-mc-sampling.jl 1000 --shutup