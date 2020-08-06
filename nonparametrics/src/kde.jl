#!usr/local/bin/julia

#=

Kernel Density Estimation
- on non-encrypted data
- on encrypted data via C++ wrapper

=#

using Distributions, Plots, SpecialFunctions, DelimitedFiles


##########################################################
#   Non-Encypted Estimators
##########################################################

# Proper estimator

"""
	kde(x0::T where T<:Real, x::Array{<:Real, 1}, h::Float64)

Kernel density estimate at point `x0` based on i.i.d. data `x`.
Uses inverse quadratic kernel with bandwidth `h`.
"""
function kde(x0::T where T<:Real, x::Array{<:Real, 1}, h::Float64)
	n = length(x)
	fhat = sum(kernel((x0 - x_i) / h) for x_i in x)/(n * h)
end

# Approximation-based estimators

"""
	kde_hat(x0::T where T<:Real, x::Array{<:Real, 1}, h::Float64; d=4, τ=100)

Kernel density estimate at point `x0` based on i.i.d. data `x`.
Uses inverse approximation of the inverse quadratic kernel with degree `d`,
shift `τ`, and bandwidth `h`.
"""
function kde_hat(x0::T where T<:Real, x::Array{<:Real, 1}, h::Float64;
		d=6, τ=100)
	n = length(x)
	khat(u) = kernel_hat(u, d=d, τ=τ)
	fhat = sum(khat((x0 - x_i) / h) for x_i in x)/(n * h)
end

##########################################################
#   Approximation Estimator Properties
##########################################################

"""
Upper bound for Mean Integrated Squre Error (MISE) of ftilde_{`d`,`τ`},
evaluated at `h`.
"""
function _mise_bound(h::Float64, d::Integer, τ::Real, n::Integer, βf::Float64)
	 h_1 = h1(d, τ)
	 h_2 = h2(d, τ, βf)
	 h_1 / (n * h) + (h ^ 4) * h_2
 end

 """
 MISE upper bound evaluated at optimal bandwith with binding constraint.
 """
 function _mise_bound_bind(d::Integer, τ::Real, n::Integer, βf::Float64,
	  	a::Real, b::Real)
	 h_opt = _hoptimal_bind(τ, a, b)
	 _mise_bound(h_opt, d, τ, n, βf)
 end

 """
 Optimal bandwidth choice when the constraint is binding.
 """
 function _hoptimal_bind(τ::Real, a::Real, b::Real)
	h_opt = (b - a) / sqrt(τ)
 end

function round2(x::Real)
	2 * round(Int, x / 2, RoundNearestTiesAway)
end


##########################################################
#   Encypted Estimators
##########################################################

"""

Calls a wrapper to the HEAAN-1.0 C++ library

Returns
- `result`::Array{Float64,1} : `length(x0)`-long vector of KDE estimates
"""
function kde_hat_he_cpp(x0::Array{T, 1} where T<:Real, x::Array{T, 1} where T<:Real,
			h::Float64; d=6, τ=10)
	@assert (iseven(d) & (d > 0) & (τ > 0)) error()
	@assert d <= 14 error("degree `d` too large. HEKDE.cpp need to be adjusted.")
	n_query = length(x0)
	n_x = length(x)
	n = 2 * n_x
	α = try
		normal_constant_dict[(d, τ)]
	# if normal constant for (d,τ) is not precomputed, compute it now
	catch e
		if isa(e, KeyError)
			_update_normal_constant_dict(d, τ)
			normal_constant_dict[(d, τ)]
		else
			error()
		end
	end
	# write arrays to files to be passed to c++ code
	for (v, vname) in zip([x0, x], ["x_query_data", "x_data"])
		open("./main/src/temp/$(vname).txt", "w") do io
		        DelimitedFiles.writedlm(io, v)
		end
	end
	# construct and call c++ code
	cd("./main/src")
	run(`./HEKDE.out $(n_query) $(n) $(h) $(d) $(Int(τ)) $(α)`);
	cd("../..")
	# read result
	result = open("./main/src/temp/result.txt", "r") do io
	        readdlm(io)
	end
	return result
end
