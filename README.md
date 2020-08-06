# Kernel Density Estimation on Homomorphically Encrypted Data

This repo contains the codes for my MPhil thesis.

## Structure

- folder main contains the functions used to create the figures
- folder nonparametrics contains lower-level functions
- once main.jl is run, folder figures will contain all figures
- folders HEAAN-1.0 and ntl-11.4.3 might be required if HAAN.out needs to be recompiled (see below)

## How to Run?

Make sure that [julia](https://julialang.org) is installed with the following packages:
- Plots.jl
- Distributions.jl
- SpecialFunctions.jl
- QuadGK.jl
- JLD.jl
- LaTeXStrings.jl
- DelimitedFiles.jl

Once they are installed, clone this repo. If [git](https://git-scm.com/) is installed, this can be done from the git bash or terminal by the command `git clone https://github.com/kmmate/mphil_thesis_codes.git`.

Then from git bash or terminal `cd mphil_thesis_codes`, and then run `julia ./main/src/main.jl` to reproduce the thesis figures.

## Notes

- The executable  main/src/HEAAN.out assumes that the machine has 8 cores. If this is not the case, HEAAN.cpp must be recompiled by g++. See main/src/HEAAN.cpp for more details.

- HEAAN.out can be used in itself to compute HE-KDE on encrypted data. For this, replace the numbers in main/src/temp/x_data.txt with the non-encrypted *i.i.d.* sample and the numbers in main/src/temp/x_query.txt with the desired query points. Then run  `./main/src/HEAAN.out [args]`, replacing `[args]` with the arguments described in main/src/HEAAN.cpp. The results will be written to main/src/temp/result.txt.
