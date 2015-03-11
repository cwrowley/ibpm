#!/bin/bash

echo "Speed benchmarking tests"
echo "Stationary body, body fixed frame running...time taken is "
tfixedubf1="$(TIMEFORMAT='%lU'; (time ./build/ibpm -outdir benchmarking/benchout -tecplot 0 -restart 0 -nx 200 -ny 200 -ngrid 3 -length 4 -xoffset -2 -yoffset -2 -xshift 0.72 -geom benchmarking/cylinder2Pa.geom -ubf 1 -Re 100 -dt 0.01 -nsteps 20 -scheme euler  2>&1 1>/dev/null))"

rm benchmarking/benchout/*
echo "Stationary body, inertial frame running...time taken is "
tfixedubf0="$(TIMEFORMAT='%lU'; (time ./build/ibpm -outdir benchmarking/benchout -tecplot 0 -restart 0 -nx 200 -ny 200 -ngrid 3 -length 4 -xoffset -2 -yoffset -2 -xshift 0.72 -geom benchmarking/cylinder2Pa.geom -ubf 0 -Re 100 -dt 0.01 -nsteps 20 -scheme euler 2>&1 1>/dev/null))"

rm benchmarking/benchout/*
echo "Moving body, body fixed frame running...time taken is "
tplungeubf1="$(TIMEFORMAT='%lU'; (time ./build/ibpm -outdir benchmarking/benchout -tecplot 0 -restart 0 -nx 200 -ny 200 -ngrid 3 -length 4 -xoffset -2 -yoffset -2 -xshift 0.72 -geom benchmarking/cylinder2PaPlunge.geom -ubf 1 -Re 100 -dt 0.01 -nsteps 20 -scheme euler 2>&1 1>/dev/null))"

rm benchmarking/benchout/*
echo "Moving body, inertial frame running...time taken is"
tplungeubf0="$(TIMEFORMAT='%lU'; (time ./build/ibpm -outdir benchmarking/benchout -tecplot 0 -restart 0 -nx 200 -ny 200 -ngrid 3 -length 4 -xoffset -2 -yoffset -2 -xshift 0.72 -geom benchmarking/cylinder2PaPlunge.geom -ubf 0 -Re 100 -dt 0.01 -nsteps 20 -scheme euler 2>&1 1>/dev/null))"

rm benchmarking/benchout/*


