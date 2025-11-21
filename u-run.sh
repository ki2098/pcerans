#!/usr/bin/bash

julia u-det-run.jl
julia u-nipce-sampling.jl

sleep 60

julia u-mc-sampling.jl 10

sleep 60

julia u-mc-sampling.jl 100

sleep 60

julia u-mc-sampling.jl 1000