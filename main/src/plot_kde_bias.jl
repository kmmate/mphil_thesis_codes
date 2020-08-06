#!/usr/local/bin/julia

#=
Plots the empirical bias as a function of bandwidth h
=#

const fig_path = "./figures"

function plot_kde_hat_bias_he()
	ls(s::String) = LaTeXString(s)
    # lists
    dist_lst = [Beta(6, 2),
				Cosine(0, 0.5),
				Arcsine(-pi, pi),
				GeneralizedPareto(0, 4, -3),
				TriangularDist(-1, 1, 1/pi),
				VonMises(0, 2)
				]
	distname_lst = map(ls, ["\$\\textrm{Beta}(\\alpha=6, \\beta=2)\$",
					"\$\\textrm{Cosine}(\\mu=0,\\sigma=0.5)\$",
					"\$\\textrm{Arcsine}(a=\\pi, b=\\pi)\$",
					"\$\\textrm{GPD}(\\mu=0,\\sigma=4,\\xi=-3)\$",
					"\$\\textrm{Triangular}(a=-1, b=1, c=1/\\pi)\$",
					"\$\\textrm{VM}(\\mu=0,\\kappa=2)\$"
					]
					)
	x0_dict = Dict(1 => # Beta
					   ([mean(dist_lst[1]), median(dist_lst[1]), mode(dist_lst[1])],
					   ["\$x=\\textrm{Mean}\$" "\$x=\\textrm{Median}\$" "\$x=\\textrm{Mode}\$"],
						[:lime :crimson :coral]),
					2 => # Cosine
						([mean(dist_lst[2]), quantile(dist_lst[2], 0.1), quantile(dist_lst[2], 0.9)],
				   		 ["\$x=\\textrm{Mean}\$" "\$x=q_{0.1}\$" "\$x=q_{0.9}\$"],
						 [:lime :darkviolet :goldenrod]),
					3 => # Arcsine
						([mean(dist_lst[3]), quantile(dist_lst[3], 0.1), quantile(dist_lst[3], 0.9)],
						 ["\$x=\\textrm{Mean}\$" "\$x=q_{0.1}\$" "\$x=q_{0.9}\$"],
						 [:lime :darkviolet :goldenrod]),
					4 => # GPD
						([mean(dist_lst[4]), median(dist_lst[4])],
					   	 ["\$x=\\textrm{Mean}\$" "\$x=\\textrm{Median}\$"],
						 [:lime :crimson]),
					5 => # Triangular
						([mean(dist_lst[5]), median(dist_lst[5]), mode(dist_lst[5])],
						 ["\$x=\\textrm{Mean}\$" "\$x=\\textrm{Median}\$" "\$x=\\textrm{Mode}\$"],
						 [:lime :crimson :coral]),
					6 => # von Mises
						([mean(dist_lst[6]), mean(dist_lst[6])-2*std(dist_lst[6]), mean(dist_lst[6])+2*std(dist_lst[6])],
				   		 ["\$x=\\textrm{Mean}\$" "\$x=\\textrm{Mean}-2\\sigma\$" "\$x=\\textrm{Mean}+2\\sigma\$"],
						 [:lime :crimson :coral])
					)
	@assert (length(dist_lst) == length(distname_lst)) error()
	# !!!!!!!! NOTE: latex needs to change if h_lst changes !!!!!!!!!!!
	h_lst_fn(h_opt) = [c * h_opt for c in (0.7, 0.8, 0.9, 0.95, 1, 1.1, 1.3, 1.5, 1.7)]

	# parameters
	n = 2^15
	n_x = Int(n / 2)
	mc_reps = 20
	optimaldeg(n::Integer) = 14
	optimalτ(n::Integer) = round(n^(1/5))
	deg_opt = optimaldeg(n_x)
	τ_opt = optimalτ(n_x)

	# plot parameters
	n_col = 2
	n_row = 3
	p = Array{Plots.Plot{backend() |> typeof}, 1}(undef, n_col*n_row)
	# plot
	for (idx_dist, (dist, distname)) in enumerate(zip(dist_lst, distname_lst))
		# kde parameters
		(low, high) = extrema(dist)
	 	h_opt = _hoptimal_bind(optimalτ(n_x), low, high)
		h_lst = h_lst_fn(h_opt)
		(x0_lst, label, color) = x0_dict[idx_dist]
		fhat = zeros(Float64, mc_reps, length(h_lst), length(x0_lst))
		for rep in 1:mc_reps
			xdata = rand(dist, n_x)
			for (idx_h, h) in enumerate(h_lst)
				for (idx_x0, x0) in enumerate(x0_lst)
					fhat[rep, idx_h, idx_x0] = kde_hat(x0, xdata, h, d=deg_opt, τ=τ_opt)
				end
			end
		end
		f = zeros(length(h_lst), length(x0_lst))
		for idx_h in 1:length(h_lst)
			f[idx_h, :] .= pdf.(dist, x0_lst)
		end
		bias = mean(fhat, dims=1)[1,:,:] .- f
		xlab_str = L"h"
		ylab_str = "\$\\textrm{Bias}_\\textrm{MC}(\\tilde f_{d_n,\\tau_n}(x))\$" |> ls
		p[idx_dist] = plot(h_lst, bias,
							linewidth=2,
							marker=:cross,
							markersize=4,
							label=label,
							legend=(0.75, 0.8),
							color=color,
							top_margin=1mm,
				   			 bottom_margin=2mm,
				   			 left_margin=1mm,
				   			 right_margin=1mm,
				   			 xtickfontsize=6,
				   			 ytickfontsize=6,
				   			 guidefontsize=8,
				   			 legendfontsize=6,
			   			 	titlefontsize=9,
						 	title=distname)
		xlabel!(xlab_str)
		(idx_dist % n_col == 1) ? ylabel!(ylab_str) : nothing
		vline!([h_opt], linewidth=2, line=:dash, color=:black, label="")
	end
	p_main = plot(p..., layout=(n_row, n_col))
	display(p_main)
	savefig(joinpath(fig_path, "plot_kde_hat_bias_he.pdf"))
end
