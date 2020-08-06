#!/usr/local/bin/julia

#=

Plots
- approximate HE kernel for different d, τ
- approximate kernel properties as function of d, τ

=#

const fig_path = "./figures"

##########################################################
#   Plots
##########################################################


"""

Approximate kernel on non-encrypted data, with different degree and τ
"""
function plot_kernel_hat()
	ls(s::String) = LaTeXString(s)
	# kernel parameters
	degree_lst = (4, 8, 30, 60, 100)
	τ_lst = (6, 64, 100)
	n_row = length(degree_lst); n_col = length(τ_lst)
	p = Array{Plots.Plot{backend() |> typeof}, 2}(undef, n_col, n_row)
	# plots
	for (row_idx, degree) in enumerate(degree_lst)
		for (col_idx, τ) in enumerate(τ_lst)
			# query points
			factor = (degree >= 4 && τ >= 64) ? 1.4 : 1.28 #(degree >= 4 && τ >= 64) ? 15 : 4#
			low = -factor*sqrt(τ)
			high = -low
			x = collect(low:0.1:high)
			# computation
			y = kernel_hat.(x, d=degree, τ=τ)
			# ploting
			title_str =  @sprintf("\$d=%d, ~ \\tau=%d \$", degree, τ) |> ls
			p[col_idx, row_idx] = plot(x, y,
												legend=false,
												 linewidth=2,
												 top_margin=0mm,
												 bottom_margin=1mm,
												 left_margin=1mm,
												 right_margin=-1mm,
												 xtickfontsize=4,
												 ytickfontsize=3,
												 guidefontsize=7,
												 title=title_str,
												 titlefontsize=8)
			# plot and annotate critical/zeropoints
			fontsize = 8
			limit_clr = :red
			limits = _get_limits(degree, τ)
			vline!([limits...], color=limit_clr)
			annot_str = L"l_1"
			annot_x_location = maximum((0.5*(maximum(x)+limits[2]), limits[2] + 0.5))
			annot_y_location = 0.8 * maximum(y)
			annotate!([(annot_x_location, annot_y_location, Plots.text(annot_str, limit_clr, :center, pointsize=fontsize))])
			zeropoint_clr = :blue
			zeropoints = _get_zeropoints(degree, τ)
			vline!([zeropoints...], color=zeropoint_clr)
			annot_str = L"l_2"
			annot_x_location = (degree <= 8) ? 0.7*(zeropoints[2]) : 0.7*zeropoints[2]+0.3*limits[2]
			annot_y_location = 0.8 * maximum(y)
			annotate!([(annot_x_location, annot_y_location, Plots.text(annot_str, zeropoint_clr, :center, pointsize=fontsize))])

			xlabel!(L"u")
			ylabel!(L"\tilde{K}_{d,\tau}(u)")
		end
	end
	title_str =  "\$\\textrm{Kernel Approximations}\$" |> ls
	p_main = plot(p..., layout=(n_row, n_col), plot_title=title_str)
	display(p_main)
	savefig(joinpath(fig_path,"plot_kernel_hat.pdf"))
end


"""
Approximate kernel non-negativity
"""
function plot_kernel_hat_sign()
	ls(s::String) = LaTeXString(s)
	n_row = 1; n_col = 2
	p = Array{Plots.Plot{backend() |> typeof}, 2}(undef, n_col, n_row)
	# Col 1: weights
	k = 1:16
	title_str = L"k\mapsto w_k"
	p[1, 1] = plot(k, _weight.(k),
					legend=false,
					linewidth=2,
					marker=:x,
					top_margin=0mm,
					bottom_margin=1mm,
					left_margin=-1mm,
					right_margin=-1mm,
					xtickfontsize=6,
					ytickfontsize=6,
					guidefontsize=8,
					title=title_str,
					titlefontsize=8)
	xlabel!(L"k")
	# Col 2: non-negativity on subdomain
	t = 0:0.05:1
	degree_lst = 2:2:10000
	f(t) = minimum(invhat(t, d=d, τ=0) for d in degree_lst)
	title_str = L"t\mapsto \min_{d\in\mathcal{D}}\gamma+\sum_{k=1}^d w_kt^k"
	y = f.(t)
	p[2, 1] = plot(t, y,
					legend=false,
					linewidth=2,
					top_margin=0mm,
					bottom_margin=1mm,
					left_margin=-1mm,
					right_margin=-1mm,
					xtickfontsize=6,
					ytickfontsize=6,
					guidefontsize=8,
					title=title_str,
					titlefontsize=8)
	my = minimum(y)
	annot_str = @sprintf("min = %.2f", my) |> ls
	hline!([my], color=:red, linewidth=2)
	annotate!([(0.1, 0.95*my+0.05*maximum(y), Plots.text(annot_str, :red, :center, pointsize=7))])
	xlabel!(L"t")
	p_main = plot(p..., layout=(n_row, n_col))
	display(p_main)
	savefig(joinpath(fig_path,"plot_kernel_hat_sign.pdf"))
end


"""
Approximate kernel properties on non-encrypted data, for different degree and τ
"""
function plot_kernel_hat_properties()
	ls(s::String) = LaTeXString(s)
	# kernel parameters
	degree_lst = (4, 8, 14, 60, 100)#, 100)
	τ = 5:100
	n_row = length(degree_lst); n_col = 4
	p = Array{Plots.Plot{backend() |> typeof}, 2}(undef, n_col, n_row)
	xlab_str = L"\sqrt{\tau}"
	# plots
	for (row_idx, degree) in enumerate(degree_lst)
		# Col 1: 0.5*(length of interval where K(u) well behaving)
		col_idx = 1
		y1 = map(x -> _get_limits(degree, x)[2], τ)
		title_str = @sprintf("\$ l_{1}(d, \\tau),~ d=%d\$", degree) |> ls
		p[col_idx, row_idx] = plot(sqrt.(τ), [y1, sqrt.(τ)],
						 legend=false,#:right,
						 legendfontsize=5,
						 linewidth=2,
						 linestyle=[:solid :dash],
						 label=["" L"45^\circ"],
						 color=[:red :blue],
						 top_margin=0mm,
						 bottom_margin=1mm,
						 left_margin=-1mm,
						 right_margin=-1mm,
						 xtickfontsize=4,
						 ytickfontsize=3,
						 guidefontsize=7,
						 title=title_str,
						 titlefontsize=8)
		annotate!([(3, maximum(y1)-1,
		 Plots.text(L"45^\circ", :blue, :center, pointsize=7))])
		xlabel!(xlab_str)

		# Col 2: 0.5*(length of interval where K(u) nonzero)
		col_idx = 2
		y2 = map(x -> _get_zeropoints(degree, x)[2], τ)
		title_str = @sprintf("\$l_{2}(d, \\tau),~ d=%d\$", degree) |> ls
		# title_str =
		# 	if row_idx == 1
		# 		L"\mu(u\in [\pm\sqrt{\tau}] : \tilde{K}_{d,\tau}(u) > 0)/2"
		# 	else
		# 	 	@sprintf("\$ d=%d,\$", degree) |> ls
		#  end
		p[col_idx, row_idx] = plot(sqrt.(τ), [y2, sqrt.(τ)],
									legend=false,#:right,
									legendfontsize=5,
									linewidth=2,
									linestyle=[:solid :dash],
									label=["" L"45^\circ"],
									color=[:red :blue],
									 top_margin=0mm,
									 bottom_margin=1mm,
									 left_margin=-1mm,
									 right_margin=-1mm,
									 xtickfontsize=4,
									 ytickfontsize=3,
									 guidefontsize=7,
									 title=title_str,
									 titlefontsize=8)
		annotate!([(3, maximum(y1)-1,
					 Plots.text(L"45^\circ", :blue, :center, pointsize=7))])
		xlabel!(xlab_str)
		# ylabel!(L"0.5*\mu(u\in [-\tau, \tau] : \tilde{K}_{d,\tau}(u)\neq 0)")

		# Col 3: ratio
		col_idx = 3
		y3 = y1 ./ y2
		title_str = @sprintf("\$l_{1}(d, \\tau)/l_{2}(d, \\tau),~ d=%d\$", degree) |> ls
		p[col_idx, row_idx] = plot(sqrt.(τ), y3,
						 legend=false,
						 legendfontsize=5,
						 linewidth=2,
						 linestyle=:solid,
						 label=[""],
						 color=[:red :blue],
						 top_margin=0mm,
						 bottom_margin=1mm,
						 left_margin=-1mm,
						 right_margin=-1mm,
						 xtickfontsize=4,
						 ytickfontsize=3,
						 guidefontsize=7,
						 title=title_str,
						 titlefontsize=8)
		xlabel!(xlab_str)

		# Col 4: difference
		col_idx = 4
		y4 = y1 .- y2
		title_str = @sprintf("\$l_{1}(d, \\tau)-l_{2}(d, \\tau),~ d=%d\$", degree) |> ls
		p[col_idx, row_idx] = plot(sqrt.(τ), y4,
						 legend=false,
						 legendfontsize=5,
						 linewidth=2,
						 linestyle=:solid,
						 label=[""],
						 color=[:red :blue],
						 top_margin=0mm,
						 bottom_margin=1mm,
						 left_margin=-1mm,
						 right_margin=-1mm,
						 xtickfontsize=4,
						 ytickfontsize=3,
						 guidefontsize=7,
						 title=title_str,
						 titlefontsize=8)
		xlabel!(xlab_str)
	end
	title_str =  "\$\\textrm{Kernel Approximation Properties}\$" |> ls
	p_main = plot(p..., layout=(n_row, n_col), plot_title=title_str)
	display(p_main)
	savefig(joinpath(fig_path, "plot_kernel_hat_properties.pdf"))
end
