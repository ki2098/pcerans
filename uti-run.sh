#!/usr/bin/bash

julia scripts/uti-det-run.jl uti-sampling-setup.json --skip-history

julia scripts/uti-nipce-sampling.jl uti-sampling-setup.json 5 10 --skip-history

sleep 60

julia scripts/uti-mc-sampling.jl uti-sampling-setup.json 10 --skip-history

sleep 60

julia uti-mc-sampling.jl 20 --shutup

sleep 60

julia uti-mc-sampling.jl 100 --shutup

sleep 60

julia uti-mc-sampling.jl 200 --shutup