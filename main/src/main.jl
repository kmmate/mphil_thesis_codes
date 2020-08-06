#!/usr/local/bin/julia

#=

All figures in the thesis are reproducable by this script.

WARNING:
encrypted HEKDE examples need to be uncommented to run, because they can take
around 24 hours on
- OS: Ubuntu 18.04
- RAM: 16 GB
- Processor: Intel i7 vPro 8th gen, 8 cores
when run in parallel on 8 cores.

=#

include("../../nonparametrics/src/kernels.jl")
include("../../nonparametrics/src/kde.jl")

include("plot_kernels.jl")
include("plot_mise.jl")
include("plot_kde.jl")
include("plot_kde_bias.jl")

# Kernels
plot_kernel_hat()
plot_kernel_hat_sign()
plot_kernel_hat_properties()

# MISE
plot_h_integrals()
plot_optimal_parameters_bind()  # statistically optimal parameters
plot_optimal_parameters_bind_he()  # encryption-optimal parameters

# KDE
plot_kde_hat()  # examples on non-encrypted data
# plot_kde_hat_he() # examples on encrypted data
plot_kde_hat_bias_he()  # bias for encryption-optimal parameters
