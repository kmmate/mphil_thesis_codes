#!usr/local/bin/julia

#=

Kernels/kernel approximations for non-encrypted data

=#
using SpecialFunctions, QuadGK, JLD
using Plots, Plots.PlotMeasures, LaTeXStrings, Printf

##########################################################
#   Kernel (approximations)
##########################################################

const recompute_normal_constant = false
const recompute_h1_constant = false
const recompute_h2_constant = false
const normal_constant_dict_path = "./nonparametrics/src/_normalising_constants_dictionaries/norm_dict.jld"
const h1_constant_dict_path = "./nonparametrics/src/_mise_constants_dictionaries/h1_integral_dict.jld"
const h2_constant_dict_path = "./nonparametrics/src/_mise_constants_dictionaries/h2_integral_dict.jld"

################# Kernel approximation

_weight(k::Integer) = (-1) ^ (k + 1)/k + zeta(k + 1) * (-1) ^ k


function invhat(u::Float64; d=4, τ=6)
	# log(u) - digamma(u)
	γ = Base.MathConstants.eulergamma
	γ + sum(_weight(k) * (u-τ) ^ k for k in 1:d)
end

function kernel_hat_unnormalised(u::Float64; d=4, τ=6)
	invhat(1 + u ^ 2, d=d, τ=τ) / pi
end


#### Normalising factor for kernel approximation

"""
Compute "critical points" ±`u` where kernel approximation starts to grow.
I.e. where its "legs end".
"""
function _get_limits(d::Real, τ::Real)
	k(u::Float64) = kernel_hat_unnormalised(u, d=d, τ=τ)
	u_old = 0.0
	u_new = 0.0
	# increase u until value of kernel starts increasing
	while k(u_old) >= k(u_new)
		u_old = u_new
		u_new += 1e-2
	end
	(-u_new, u_new)
end

# Compute def integral between "critical points"
function _normal_constant(d::Integer, τ::Real)
	k(u::Float64) = kernel_hat_unnormalised(u, d=d, τ=τ)
	limits = _get_limits(d, τ)
	integral, err = quadgk(k, limits[1], limits[2])
	integral
end

function _normal_constant_dict(recompute_normal_constant::Bool)
	if recompute_normal_constant
		d_lst = 2:2:20
		τ_lst = 1:1:1e4
		dict = Dict((d, τ) => _normal_constant(d, τ) for τ in τ_lst, d in d_lst)
		println("_normal_constant_dict is computed. Saving file...")
		#JLD.save(normal_constant_dict_path, "dict", dict)
		jldopen(normal_constant_dict_path, "w") do file
    		write(file, "dict", dict)
		end
		println("_normal_constant_dict is saved.")
	else
		#dict = JLD.load(normal_constant_dict_path)["dict"]
		dict = jldopen(normal_constant_dict_path, "r") do file
    				read(file, "dict")
		end
	end
	dict
end
normal_constant_dict = _normal_constant_dict(recompute_normal_constant)

# singla (d,τ) updating
function _update_normal_constant_dict(d::Integer, τ::Real)
	# if normal constant for (d,τ) is not precomputed, compute it,
	# add to normal_constant_dict and write it to file
	error("Do batch update instead  involving d=$(d), τ=$(τ)")
	if !haskey(normal_constant_dict, (d, τ))
			α =  _normal_constant(d, τ)
			normal_constant_dict[(d, τ)] = α
			jldopen(normal_constant_dict_path, "w") do file
	    		write(file, "dict", normal_constant_dict)
			end
	end
end

# batch updating
function _update_normal_constant_dict(d::Union{Array{T, 1}, StepRange{T, T}} where T<:Integer,
	 	τ::Union{Array{T, 1}, UnitRange{T}, StepRangeLen{T, A, A}} where {T<:Real, A<:Base.TwicePrecision})
	updated = false
	for dval in d, τval in τ
		if !haskey(normal_constant_dict, (dval, τval))
				normalint =  _normal_constant(dval, τval)
				normal_constant_dict[(dval, τval)] = normalint
				updated = true
		end
	end
	if updated
		jldopen(normal_constant_dict_path, "w") do file
			write(file, "dict", normal_constant_dict)
		end
	end
end


function kernel_hat(u::Float64; d::Integer=4, τ::Real=6)
	@assert (iseven(d) & (d > 0) & (τ > 0)) error()
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
	kernel_hat_unnormalised(u, d=d, τ=τ) / α
end


"""
Find points ±`u` closest to zero such that kernel_hat(±`u`)≈0.
I.e. where the kernel approximation start flattening out, where the
"legs of the kernel approximation start".
"""
function _get_zeropoints(d::Real, τ::Real)
	k(u::Float64) = kernel_hat(u, d=d, τ=τ)
	u = 0.0
	itercount = 0
	itermax = 1e5
	stepsize = 1e-3 * sqrt(τ)
	# increase u until value of kernel starts increasing
	while !isapprox(k(u), 0, atol=1e-3)
		u += stepsize
		itercount +=1
		@assert itercount < itermax error("no convergence, kernel values do not near zero")
	end
	(-u, u)
end

#### terms in kernel approximation MISE

# h1

"""
h₁int = ∫_{-√τ}^{√τ} Ktilde_{d,τ}(z)^2 dz
"""
function _h1_integral(d::Integer, τ::Real)
	khat_sq(u::Float64) = kernel_hat(u, d=d, τ=τ) ^ 2
	limits = (-sqrt(τ), sqrt(τ))
	integral, err = quadgk(khat_sq, limits[1], limits[2])
	integral
end

function _h1_constant_dict(recompute_h1_constant::Bool)
	if recompute_h1_constant
		d_lst = 2:2:20
		τ_lst = 1:1:1e4
		dict = Dict((d, τ) => _h1_integral(d, τ) for τ in τ_lst, d in d_lst)
		println("_h1_constant_dict is computed. Saving file...")
		jldopen(h1_constant_dict_path, "w") do file
    		write(file, "dict", dict)
		end
		println("_h1_constant_dict is saved.")
	else
		dict = jldopen(h1_constant_dict_path, "r") do file
    				read(file, "dict")
		end
	end
	dict
end
h1_constant_dict = _h1_constant_dict(recompute_h1_constant)

function _update_h1_constant_dict(d::Integer, τ::Real)
	# if h1 constant for (d,τ) is not precomputed, compute it now,
	# add to h1_constant_dict and write it to file
	error("Do batch update instead involving d=$(d), τ=$(τ)")
	if !haskey(h1_constant_dict, (d, τ))
			h1int =  _h1_integral(d, τ)
			h1_constant_dict[(d, τ)] = h1int
			jldopen(h1_constant_dict_path, "w") do file
	    		write(file, "dict", h1_constant_dict)
			end
	end
end

function _update_h1_constant_dict(d::Union{Array{T, 1}, StepRange{T, T}} where T<:Integer,
	 	τ::Union{Array{T, 1}, UnitRange{T}, StepRangeLen{T, A, A}} where {T<:Real, A<:Base.TwicePrecision})
	_update_normal_constant_dict(d, τ)
	updated = false
	for dval in d, τval in τ
		if !haskey(h1_constant_dict, (dval, τval))
				h1int =  _h1_integral(dval, τval)
				h1_constant_dict[(dval, τval)] = h1int
				updated = true
		end
	end
	if updated
		jldopen(h1_constant_dict_path, "w") do file
			write(file, "dict", h1_constant_dict)
		end
	end
end

"""
h₁(d, τ) = h₁int
"""
function h1(d::Integer, τ::Real)
	h1int = try
		h1_constant_dict[(d, τ)]
	# if integral for (d,τ) is not precomputed, compute it now
	catch e
		if isa(e, KeyError)
			_update_h1_constant_dict(d, τ)
			h1_constant_dict[(d, τ)]
		else
			error()
		end
	end
	h1int
end


# h2

"""
h₂int = ∫_{-√τ}^{√τ} Ktilde_{d,τ}(z)*z^2 dz
"""
function _h2_integral(d::Integer, τ::Real)
	khatusq(u::Float64) = u^2 * kernel_hat(u, d=d, τ=τ)
	limits = (-sqrt(τ), sqrt(τ))
	integral, err = quadgk(khatusq, limits[1], limits[2])
	integral
end


function _h2_constant_dict(recompute_h2_constant::Bool)
	if recompute_h2_constant
		d_lst = 2:2:20
		τ_lst = 1:1:1e4
		dict = Dict((d, τ) => _h2_integral(d, τ) for τ in τ_lst, d in d_lst)
		println("_h2_constant_dict is computed. Saving file...")
		jldopen(h2_constant_dict_path, "w") do file
    		write(file, "dict", dict)
		end
		println("_h2_constant_dict is saved.")
	else
		dict = jldopen(h2_constant_dict_path, "r") do file
    				read(file, "dict")
		end
	end
	dict
end
h2_constant_dict = _h2_constant_dict(recompute_h2_constant)

function _update_h2_constant_dict(d::Integer, τ::Real)
	# if h2 constant for (d,τ) is not precomputed, compute it now,
	# add to h2_constant_dict and write it to file
	error("Do batch update instead  involving d=$(d), τ=$(τ)")
	if !haskey(h2_constant_dict, (d, τ))
			h2int =  _h2_integral(d, τ)
			h2_constant_dict[(d, τ)] = h2int
			jldopen(h2_constant_dict_path, "w") do file
	    		write(file, "dict", h2_constant_dict)
			end
	end
end


function _update_h2_constant_dict(d::Union{Array{T, 1}, StepRange{T, T}} where T<:Integer,
	 	τ::Union{Array{T, 1}, UnitRange{T}, StepRangeLen{T, A, A}} where {T<:Real, A<:Base.TwicePrecision})
	_update_normal_constant_dict(d, τ)
	updated = false
	for dval in d, τval in τ
		if !haskey(h2_constant_dict, (dval, τval))
				h2int =  _h2_integral(dval, τval)
				h2_constant_dict[(dval, τval)] = h2int
				updated = true
		end
	end
	if updated
		jldopen(h2_constant_dict_path, "w") do file
			write(file, "dict", h2_constant_dict)
		end
	end
end

"""
h₂(d, τ) = (h₂int ^ 2) * (βf / 3)
"""
function h2(d::Integer, τ::Real, βf::Real)
	h2int = try
		h2_constant_dict[(d, τ)]
	# if integral for (d,τ) is not precomputed, compute it now
	catch e
		if isa(e, KeyError)
			_update_h2_constant_dict(d, τ)
			h2_constant_dict[(d, τ)]
		else
			error()
		end
	end
	(h2int ^ 2) * (βf / 3)
end
