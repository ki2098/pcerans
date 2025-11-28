#!/usr/bin/bash

julia uti-det-run.jl --shutup

julia uti-nipce-sampling.jl --shutup

sleep 60

julia uti-mc-sampling.jl 10 --shutup

sleep 60

julia uti-mc-sampling.jl 50 --shutup

sleep 60

julia uti-mc-sampling.jl 100 --shutup