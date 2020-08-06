#!/usr/local/bin/julia

#=

Plots non-encrypted and encrypted examples of HE-KDE
=#

const fig_path = "./figures"

##########################################################
#   Non-Encypted Estimators
##########################################################


"""
KDE plot on non-encrypted data: single sample
"""
function _plot_kde_hat_dist_singlesample(dist::T where T<:UnivariateDistribution,
										 dist_name::String,
										 n_x::Integer,
										 use_legend::Bool)
	# parameters
	(low, high) = extrema(dist)
	optimaldeg(n::Integer) = round2(n^(4/5))
	optimalτ(n::Integer) = round(n^(2/5))
	deg_opt = optimaldeg(n_x)
	τ_opt = optimalτ(n_x)
	h_opt = _hoptimal_bind(optimalτ(n_x), low, high)
	δ = (high - low) / 100
	x0 = 1*low:δ:1*high  # query points
	# generate data and estimate
	xdata = rand(dist, n_x)
	fhat = zeros(length(x0))  # query values
	for (i, x) in enumerate(x0)
		fhat[i] = kde_hat(x, xdata, h_opt, d=deg_opt, τ=τ_opt)
	end
	f = pdf.(dist, x0)

	# plot
	ls(s::String) = LaTeXString(s)
	xlab_str = L"x"
	ylab_str = L"f(x),~\tilde f_{d_n, \tau_n}(x)"
	title_str = @sprintf("\$ %s,~ n = %d \$", dist_name, n_x) |> ls
	legendpos = use_legend ? (0.8, 0.75) : use_legend
	p = plot(x0,
			 [f fhat],
			 legend=legendpos,
			 label=[L"\textrm{True}" L"\textrm{Estimate}"],
			 color=[:blue :red],
			 linewidth=2,
			 top_margin=1mm,
			 bottom_margin=2mm,
			 left_margin=2mm,
			 right_margin=2mm,
			 xtickfontsize=6,
			 ytickfontsize=4,
			 guidefontsize=6,
			 legendfontsize=5,
			 title=title_str,
			 titlefontsize=8)
	scatter!(xdata,
			zeros(length(xdata)),
			label="",
			marker=:vline,
			color=:blue)
	xlabel!(xlab_str)
	ylabel!(ylab_str)
	return p
end

"""
KDE plot on non-encrypted data: monte carlo
"""
function _plot_kde_hat_dist_montecarlo(dist::T where T<:UnivariateDistribution,
										 dist_name::String,
										 n_x::Integer,
										 mc_reps::Integer,
										 use_legend::Bool)
    # parameters
 	(low, high) = extrema(dist)
 	optimaldeg(n::Integer) = round2(n^(4/5))
 	optimalτ(n::Integer) = round(n^(2/5))#round(n^(2/5))
 	deg_opt = optimaldeg(n_x)
 	τ_opt = optimalτ(n_x)
 	h_opt = _hoptimal_bind(optimalτ(n_x), low, high)
	δ = (high - low) / 100
	x0 = 1*low:δ:1*high  # query points
	f = pdf.(dist, x0)  # true density
	ls(s::String) = LaTeXString(s)
	fhat = zeros(Float64, length(x0), mc_reps)
	p = plot()
	# monte carlo reps
	for rep in 1:mc_reps
	 	# generate data and estimate
	 	xdata = rand(dist, n_x)
	 	#fhat = zeros(length(x0))  # query values
	 	for (i, x) in enumerate(x0)
	 		fhat[i, rep] = kde_hat(x, xdata, h_opt, d=deg_opt, τ=τ_opt)
	 	end
		# plot
		plot!(p, x0, fhat[:, rep], linewidth=0.5, label="")
	end
	fhat_mean = mean(fhat, dims=2)
	# plot true density
	xlab_str = L"x"
	# title_str = @sprintf("\$ \\textrm{MC:}~%s,~n = %d \$", dist_name, n_x) |> ls
	title_str = @sprintf("\$ \\textrm{Monte Carlo:} ~n = %d \$", n_x) |> ls
	legendpos = use_legend ? (0.8, 0.75) : use_legend
	plot!(p, x0,
			 [f fhat_mean],
			 legend=legendpos,
			 label=[L"\textrm{True}" L"\textrm{MC mean}"],
			 color=[:blue :red],
			 linewidth=2,
			 top_margin=1mm,
			 bottom_margin=2mm,
			 left_margin=2mm,
			 right_margin=4mm,
			 xtickfontsize=6,
			 ytickfontsize=4,
			 guidefontsize=6,
			 legendfontsize=5,
			 title=title_str,
			 titlefontsize=8)
 	xlabel!(xlab_str)
 	return p
end

"""
KDE plots on non-encrypted data: single sample and monte carlo
"""
function _plot_kde_hat_dist(dist_lst::Array{T, 1} where T<:UnivariateDistribution,
						distname_lst::Array{String, 1},
						figname::String)
	@assert length(dist_lst) == length(distname_lst) error("lists must have same length")
	print("_plot_kde_hat_dist is working on: ", figname, ". ")
	# plot parameters
	n_col = 3
	n_row = length(dist_lst)
	col_idx = 1
	p = Array{Plots.Plot{backend() |> typeof}, 2}(undef, n_col, n_row)
	# kde parameters
	mc_reps = 20
	n_lst = [2^7, 2^8]
	n_x_lst = Int.(n_lst ./ 2)
	for (row_idx, dist) in enumerate(dist_lst)
		println("- distribution: ", dist)
		# single sample plot
		legend = (row_idx == 1) ? true : false
		p[col_idx, row_idx] = _plot_kde_hat_dist_singlesample(dist,
		 										distname_lst[row_idx], n_x_lst[1],
												legend)
		col_idx += 1
		# monte carlo plots
		for n_x in n_x_lst
			legend = (row_idx == 1) ? true : false
			p[col_idx, row_idx] = _plot_kde_hat_dist_montecarlo(dist,
											distname_lst[row_idx], n_x, mc_reps,
											legend)
			col_idx += 1
		end
		col_idx = 1
	end
	p_main = plot(p..., layout=(n_row, n_col))
	display(p_main)
	#savefig(joinpath(fig_path, @sprintf("plot_kde_hat_%s.pdf", figname)))
end

"""
Multiple KDE plots on non-encrypted data for different distributions.
"""
function plot_kde_hat()
	# Cont: Beta, Cosine,GeneralizedPareto(0, 2, ξ=-2) support bounded only for ξ<0,
	# Semicircle, SymTriangular,
	# Triangular, VonMises
	# Discrete: BetaBinomial(high, 0.7, 2), Binomial(10, 0.5)
	# create lists
	dist_lst = [[Beta(1,2), Beta(2,2), Beta(3, 2), Beta(6, 2)],
				[Cosine(0, 0.5), Cosine(0, 2), Cosine(0, 4)],
				[Arcsine(-1, 1), Arcsine(-pi, pi)],
				[GeneralizedPareto(0, 2, -2), GeneralizedPareto(0, 1, -10), GeneralizedPareto(0, 4, -3), GeneralizedPareto(0, 1, -0.1)],
				[Semicircle(0.3), Semicircle(3)],
				[TriangularDist(-1, 1, -0.9), TriangularDist(-1, 1, 1/pi), TriangularDist(-10, 10, 0)],
				[VonMises(0, 0.1), VonMises(0, 0.5), VonMises(0, 2), VonMises(0, 3)],
				[Binomial(50, 0.1), Binomial(50, 0.5), Binomial(50, 0.8)],
				[BetaBinomial(50, 0.8, 1), BetaBinomial(10, 0.2, 0.2)]
				]
	distname_lst = [map(x->@sprintf("\\textrm{Beta}%s", x), ["(\\alpha=1,\\beta=2)", "(2,2)", "(3,2)", "(6,2)"]),
					map(x->@sprintf("\\textrm{Cosine}%s", x), ["(\\mu=0,\\sigma=0.5)", "(0,2)", "(0,4)"]),
					map(x->@sprintf("\\textrm{Arcsine}%s", x), ["(a=-1, b=1)", "(-\\pi,\\pi)"]),
					map(x->@sprintf("\\textrm{GPD}%s", x), ["(\\mu=0,\\sigma=2, \\xi=-2)", "(0,1,-10)", "(0,4,-3)", "(0, 1, -0.1)"]),
					map(x->@sprintf("\\textrm{Semicircle}%s", x), ["(r=0.3)", "(3)"]),
					map(x->@sprintf("\\textrm{Triangular}%s", x), ["(a=-1, b=1, c=-0.9)", "(-1, 1, 1/\\pi)", "(-10, 10, 0)"]),
					map(x->@sprintf("\\textrm{VM}%s", x), ["(\\mu=0, \\kappa=0.1)", "(0, 0.5)", "(0, 2)", "(0, 3)"]),
					map(x->@sprintf("\\textrm{Bin}%s", x), ["(n=50, p=0.1)", "(50, 0.5)", "(50, 0.8)"]),
					map(x->@sprintf("\\textrm{BB}%s", x), ["(n=50, \\alpha=0.8, \\beta=1)", "(10, 0.2, 0.2)"])
					]
	figname_lst = ["cont_beta",
				   "cont_cosine",
				   "cont_arcsine",
				   "cont_gpd",
				   "cont_semicircle",
				   "cont_triangular",
				   "cont_vonmises",
				   "disc_binomial",
				   "disc_betabinomial"
				   ]
	# plot
	for (dist, distname, figname) in zip(dist_lst, distname_lst, figname_lst)
		_plot_kde_hat_dist(dist, distname, figname)
	end
end


##########################################################
#   Encrypted Estimators
##########################################################


"""
KDE plot on encrypted data: single sample
"""
function _plot_kde_hat_dist_singlesample_he(dist::T where T<:UnivariateDistribution,
										 dist_name::String,
										 n_x::Integer,
										 use_legend::Bool)
	# kde parameters
	(low, high) = extrema(dist)
	optimaldeg(n::Integer) = 14
	optimaldeg_he(n::Integer) = 14
	optimalτ(n::Integer) = round(n^(1/5))
	deg_opt = optimaldeg(n_x)
	deg_opt_he = optimaldeg_he(n_x)
	τ_opt = optimalτ(n_x)
	h_opt = _hoptimal_bind(optimalτ(n_x), low, high)
	δ = (high - low) / 50
	x0 = (1*low):δ:(1*high) |> collect  # query points
	# generate data
	xdata = rand(dist, n_x)
	# estimate
	fhat = map(x->kde_hat(x, xdata, h_opt, d=deg_opt, τ=τ_opt), x0)
	fhat_he = kde_hat_he_cpp(x0, xdata, h_opt, d=deg_opt_he, τ=τ_opt)  # query values
	f = pdf.(dist, x0)

	# plot
	ls(s::String) = LaTeXString(s)
	xlab_str = L"x"
	ylab_str = L"f(x),~\tilde f_{d_n, \tau_n}(x)"
	title_str = @sprintf("\$ %s,~ n = %d \$", dist_name, n_x) |> ls
	legendpos = use_legend ? (0.8, 0.75) : use_legend
	topmargin = (dist_name == "cont_gpd") ? 3mm : 1mm
	p = plot(x0,
			 [f fhat fhat_he],
			 legend=legendpos,
			 label=[L"f(x)" L"\tilde f_{d_n, \tau_n}(x)" L"\tilde f_{d_n, \tau_n}(x) \textrm{ encrypted}"],
			 color=[:blue :red :cyan],
			 linewidth=2,
			 top_margin=topmargin,
			 bottom_margin=2mm,
			 left_margin=2mm,
			 right_margin=2mm,
			 xtickfontsize=6,
			 ytickfontsize=4,
			 guidefontsize=6,
			 legendfontsize=5,
			 title=title_str,
			 titlefontsize=8)
	scatter!(xdata,
			zeros(length(xdata)),
			label="",
			marker=:vline,
			color=:blue)
	xlabel!(xlab_str)
	ylabel!(ylab_str)
	return p
end

"""
KDE plot on encrypted data: monte carlo
"""
function _plot_kde_hat_dist_montecarlo_he(dist::T where T<:UnivariateDistribution,
										 dist_name::String,
										 n_x::Integer,
										 mc_reps::Integer,
										 use_legend::Bool)
	# parameters
 	(low, high) = extrema(dist)
	optimaldeg(n::Integer) = 14
 	optimaldeg_he(n::Integer) = 14
 	optimalτ(n::Integer) = round(n^(1/5))
	deg_opt = optimaldeg(n_x)
 	deg_opt_he = optimaldeg_he(n_x)
 	τ_opt = optimalτ(n_x)
 	h_opt = _hoptimal_bind(optimalτ(n_x), low, high)
	δ = (high - low) / 20
	x0 = (1*low):δ:(1*high) |> collect  # query points
	ls(s::String) = LaTeXString(s)
	fhat = zeros(Float64, length(x0), mc_reps)
	fhat_he = zeros(Float64, length(x0), mc_reps)
	p = plot()
	# monte carlo reps
	for rep in 1:mc_reps
		println("--- mc rep = ", rep)
	 	# generate data and estimate
	 	xdata = rand(dist, n_x)
		# non-encrypted
		fhat[:, rep] = map(x->kde_hat(x, xdata, h_opt, d=deg_opt, τ=τ_opt), x0)
		# encrypted
		fhat_he[:, rep] = kde_hat_he_cpp(x0, xdata, h_opt, d=deg_opt_he, τ=τ_opt)
		# plot
		plot!(p, x0, fhat_he[:, rep], linewidth=0.5, label="")
	end
	fhat_mean = mean(fhat, dims=2)  # non-encrypted MC mean
	fhat_mean_he = mean(fhat_he, dims=2)  # encrypted MC mean
	f = pdf.(dist, x0)  # true density
	# plot true density
	xlab_str = L"x"
	title_str = @sprintf("\$ \\textrm{Monte Carlo:} ~n = %d \$", n_x) |> ls
	legendpos = use_legend ? (0.8, 0.75) : use_legend
	plot!(p, x0,
			 [f fhat_mean fhat_mean_he],
			 legend=legendpos,
			 label=[L"\textrm{True}" L"\textrm{MC mean}" L"\textrm{MC mean enc.}"],
			 color=[:blue :red :cyan],
			 linewidth=2,
			 top_margin=1mm,
			 bottom_margin=2mm,
			 left_margin=2mm,
			 right_margin=4mm,
			 xtickfontsize=6,
			 ytickfontsize=4,
			 guidefontsize=6,
			 legendfontsize=5,
			 title=title_str,
			 titlefontsize=8
			 )
 	xlabel!(xlab_str)
 	return p
end

"""
KDE plots on encrypted data: single sample and monte carlo
"""
function _plot_kde_hat_dist_he(dist_lst::Array{T, 1} where T<:UnivariateDistribution,
						distname_lst::Array{String, 1},
						figname::String)
	@assert length(dist_lst) == length(distname_lst) error("lists must have same length")
	println("_plot_kde_hat_dist_he is working on: ", figname, ". ")
	# kde parameters
	mc_reps = 10
	n_lst = [2^7, 2^15]
	n_x_lst = Int.(n_lst ./ 2)
	# plot parameters
	n_col = length(n_lst) + 1
	n_row = length(dist_lst) + 1  # +1 for legend plot
	col_idx = 1
	p = Array{Plots.Plot{backend() |> typeof}, 2}(undef, n_col, n_row)
	for (row_idx, dist) in enumerate(dist_lst)
		println("- distribution: ", dist)
		# single sample plot
		legend = false#(row_idx == 1) ? true : false
		p[col_idx, row_idx] = _plot_kde_hat_dist_singlesample_he(dist,
		 										distname_lst[row_idx],
												n_x_lst[1],
												legend)
		col_idx += 1
		# monte carlo plots
		for n_x in n_x_lst
			legend = false#(row_idx == 1) ? true : false
			p[col_idx, row_idx] = _plot_kde_hat_dist_montecarlo_he(dist,
											distname_lst[row_idx], n_x, mc_reps,
											legend)
			col_idx += 1
		end
		col_idx = 1
	end

	# legends in separate subplots
	ls(s::String) = LaTeXString(s)
	heights = if n_row-1 == 3
		[9, 4.5, 0]
	elseif n_row-1 == 2
		[4, 2, 0]
	else
		[160, 80, 0]
	end
	fontsize = (n_row-1) > 3 ? 7 : 8
	# first column: single sample
	p[1, n_row] = plot(1:3,
					[heights[1] .* ones(3) heights[2] .* ones(3) heights[3] .* ones(3)],
					grid=false, showaxis=false, legend=false, xlims=(1, 11),
					color=[:blue :red :cyan], linewidth=2, top_margin=-8mm)
    annotate!(p[1, n_row], [(3.5, heights[1], Plots.text("\$f(x)\\textrm{ true}\$" |> ls, :left, pointsize=fontsize)),
                (3.5, heights[2], Plots.text("\$\\tilde{f}_{d_n,\\tau_n}(x)\\textrm{ non-encrypted}\$" |> ls, :left, pointsize=fontsize)),
				(3.5, heights[3], Plots.text("\$\\tilde{f}_{d_n,\\tau_n}(x)\\textrm{ encrypted}\$" |> ls, :left, pointsize=fontsize))])
	# other columns: monte carlo
	for col_idx in 2:n_col
		p[col_idx, n_row] = plot(1:3,
						[heights[1] .* ones(3) heights[2] .* ones(3) heights[3] .* ones(3)],
						grid=false, showaxis=false, legend=false, xlims=(1, 12.5),
						color=[:blue :red :cyan], linewidth=2, top_margin=-2mm)
	    annotate!(p[col_idx, n_row], [(3.5, heights[1], Plots.text("\$\\textrm{True}\$" |> ls, :left, pointsize=fontsize)),
	                (3.5, heights[2], Plots.text("\$\\textrm{MC mean non-encrypted}\$" |> ls, :left, pointsize=fontsize)),
					(3.5, heights[3], Plots.text("\$\\textrm{MC mean encrypted}\$" |> ls, :left, pointsize=fontsize))])
	end

	p_main = plot(p..., layout=(n_row, n_col))
	display(p_main)
	savefig(joinpath(fig_path, @sprintf("plot_kde_hat_he_%s.pdf", figname)))
end

"""
Multiple KDE plots on encrypted data for different distributions.
"""
function plot_kde_hat_he()
	# Cont: Beta, Cosine,GeneralizedPareto(0, 2, ξ=-2) support bounded only for ξ<0,
	# Semicircle, SymTriangular,
	# Triangular, VonMises
	# Discrete: BetaBinomial(high, 0.7, 2), Binomial(10, 0.5)
	# create lists
	dist_lst = [[Beta(1,2), Beta(2,2), Beta(3, 2), Beta(6, 2)],
				[Cosine(0, 0.5), Cosine(0, 2), Cosine(0, 4)],
				[Arcsine(-1, 1), Arcsine(-pi, pi)],
				[GeneralizedPareto(0, 2, -2), GeneralizedPareto(0, 1, -10), GeneralizedPareto(0, 4, -3), GeneralizedPareto(0, 1, -0.1)],
				[Semicircle(0.3), Semicircle(3)],
				[TriangularDist(-1, 1, -0.9), TriangularDist(-1, 1, 1/pi), TriangularDist(-10, 10, 0)],
				[VonMises(0, 0.1), VonMises(0, 0.5), VonMises(0, 2), VonMises(0, 3)],
				[Binomial(50, 0.1), Binomial(50, 0.5), Binomial(50, 0.8)],
				[BetaBinomial(50, 0.8, 1), BetaBinomial(10, 0.2, 0.2)]
				]
	distname_lst = [map(x->@sprintf("\\textrm{Beta}%s", x), ["(\\alpha=1,\\beta=2)", "(2,2)", "(3,2)", "(6,2)"]),
					map(x->@sprintf("\\textrm{Cosine}%s", x), ["(\\mu=0,\\sigma=0.5)", "(0,2)", "(0,4)"]),
					map(x->@sprintf("\\textrm{Arcsine}%s", x), ["(a=-1, b=1)", "(-\\pi,\\pi)"]),
					map(x->@sprintf("\\textrm{GPD}%s", x), ["(\\mu=0,\\sigma=2, \\xi=-2)", "(0,1,-10)", "(0,4,-3)", "(0, 1, -0.1)"]),
					map(x->@sprintf("\\textrm{Semicircle}%s", x), ["(r=0.3)", "(3)"]),
					map(x->@sprintf("\\textrm{Triangular}%s", x), ["(a=-1, b=1, c=-0.9)", "(-1, 1, 1/\\pi)", "(-10, 10, 0)"]),
					map(x->@sprintf("\\textrm{VM}%s", x), ["(\\mu=0, \\kappa=0.1)", "(0, 0.5)", "(0, 2)", "(0, 3)"]),
					map(x->@sprintf("\\textrm{Bin}%s", x), ["(n=50, p=0.1)", "(50, 0.5)", "(50, 0.8)"]),
					map(x->@sprintf("\\textrm{BB}%s", x), ["(n=50, \\alpha=0.8, \\beta=1)", "(10, 0.2, 0.2)"])
					]
	figname_lst = ["cont_beta",
				   "cont_cosine",
				   "cont_arcsine",
				   "cont_gpd",
				   "cont_semicircle",
				   "cont_triangular",
				   "cont_vonmises",
				   "disc_binomial",
				   "disc_betabinomial"
				   ]
	@assert (length(dist_lst) == length(distname_lst) == length(figname_lst)) error()
	# plot
	for (dist, distname, figname) in zip(dist_lst, distname_lst, figname_lst)
		_plot_kde_hat_dist_he(dist, distname, figname)
	end
end
