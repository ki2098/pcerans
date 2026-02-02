#!/bin/bash

julia --startup-file=no -e "using DaemonMode; serve()" > juliaserver.log 2>&1 &
pid=$!
echo "juliaserver's pid = $pid"

sleep 1

for i in {1..10}; do
    julia --startup-file=no -e "using DaemonMode; runargs()" dummy.jl
done

kill -TERM "$pid"
wait "$pid"