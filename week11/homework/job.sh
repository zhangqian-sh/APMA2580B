#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -n 8
#SBATCH --mem=32g
#SBATCH -p batch
#SBATCH -C 48core

#SBATCH -J 2D-Burgers

#SBATCH -o misc/output/2D-Burgers-%j.out
#SBATCH -e misc/output/2D-Burgers-%j.err

cd /users/qzhang64/Projects/Courses/APMA2580B/week11/homework/ || exit
module load julia/1.8.5

julia main.jl 0.3
julia main.jl 2.0
