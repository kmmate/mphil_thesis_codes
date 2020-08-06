#!/usr/local/bin/julia

#=
Plots
- integrals I₁, I₂
- MISE for binding bandwidth with statistically optimal d, τ
- MISE for binding bandwidth with encryption- optimal d, τ

=#

const fig_path = "./figures"

##########################################################
#   Plots
##########################################################

"""
Plot the integrals
h₁int = ∫_{-√τ}^{√τ} Ktilde_{d,τ}(z)^2 dz
h₂int = ∫_{-√τ}^{√τ} Ktilde_{d,τ}(z)*z^2 dz
as a function of `d` and `τ`.
"""
function plot_h_integrals()
	ls(s::String) = LaTeXString(s)
	degree_lst = 4:2:40
	τ_lst1 = 1:800
	τ_lst2 = 1:1000
	n_row = 1; n_col = 2
	xlab_str = L"\tau"
	ylab_str = L"d"
	p = Array{Plots.Plot{backend() |> typeof}, 2}(undef, n_col, n_row)

	# plot options
	xtickfontsize = 8
	ytickfontsize = 8
	legendfontsize = 9
	guidefontsize = 11
	titlefontsize = 11
	GR.setcolormap(46)
	mycgrad = cgrad([RGB(r...) for r=eachrow(GR.colormap())])

	# h1int
	title_str = "\$(d, \\tau)\\mapsto I_1(d, \\tau)= \\int_{-\\sqrt{\\tau}}^{\\sqrt{\\tau}}\\tilde{K}_{d,\\tau}(z)^2 dz\$" |> ls
	p[1,1] = contour(τ_lst1, degree_lst, (x, y) -> h1(y, x),
					seriescolor=mycgrad,
	 				fill=true,
					top_margin=4mm,
					bottom_margin=2mm,
					left_margin=2mm,
					right_margin=2mm,
					xtickfontsize=xtickfontsize,
					ytickfontsize=ytickfontsize,
					guidefontsize=guidefontsize,
					title=title_str,
					titlefontsize=titlefontsize)
	xlabel!(xlab_str)
	ylabel!(ylab_str)
	title_str = "\$(d, \\tau)\\mapsto I_2(d, \\tau)=\\int_{-\\sqrt{\\tau}}^{\\sqrt{\\tau}}z^2\\tilde{K}_{d,\\tau}(z) dz\$" |> ls
	# note: h2(d, τ, 3.0)=_h2_integral(d, τ)
	p[2,1] = contour(τ_lst2, degree_lst, (x, y) -> h2(y, x, 3.0),
					seriescolor=mycgrad,
					fill=true,
					top_margin=4mm,
					bottom_margin=2mm,
					left_margin=2mm,
					right_margin=2mm,
					xtickfontsize=xtickfontsize,
					ytickfontsize=ytickfontsize,
					guidefontsize=guidefontsize,
					title=title_str,
					titlefontsize=titlefontsize)
	xlabel!(xlab_str)
	ylabel!(ylab_str)
	# main plot
	title_str =  "\$\\textrm{MISE Integrals}\$" |> ls
	p_main = plot(p..., layout=(n_row, n_col))
	display(p_main)
	savefig(joinpath(fig_path,"plot_h_integrals.pdf"))
end

function plot_optimal_parameters_bind()
	ls(s::String) = LaTeXString(s)
	n_plots = 6
	p = Array{Plots.Plot{backend() |> typeof}, 1}(undef, n_plots)
	# distribution parameters
	βf = 3.0
	low = 0.0
	# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#NOTE: if high_dict, low, βf, dpower, τpower changes, LaTeX code mustb e chenged !!!!!!! #NOTE
	# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	high_dict = Dict("maintext"=> 1.5,
					"appendix_1"=>0.01,
					"appendix_2"=>1.0,
					"appendix_3"=>1.2,
					"appendix_4"=>1.7,
					"appendix_5"=>10,
					"appendix_6"=>10000)
	# optimal parameters
	dpower = 4/5
	τpower = 2/5
	optimald(n::Integer) = round2(n^dpower)
	optimalτ(n::Integer) = round(n^τpower)
	# optimalτ(n::Integer) = round((high-low)^2 * n^τpower)

	# plot options
	xtickfontsize = 7
	ytickfontsize = 7
	legendfontsize = 9
	guidefontsize = 11
	titlefontsize = 11
	GR.setcolormap(46)
	mycgrad = cgrad([RGB(r...) for r=eachrow(GR.colormap())])

	for (name, high) in high_dict
		# ------ d_n
		τ_lst = 2:100
		n_lst = 10:2:200
		mise_bindd(n::Integer, τ::Real) = _mise_bound_bind(optimald(n), τ, n, βf, low, high)

		# -- optimal d_n
		xlab_str = L"n"
		# title_str = @sprintf("\$S=[%.0f, %.2f],~\\beta_f=%.0f;~~~n\\mapsto d_n=\\lfloor n^{%d/5}\\rceil_2\$", low, high, βf, dpower * 5) |> ls
		title_str = @sprintf("\$n\\mapsto d_n=\\lfloor n^{%d/5}\\rceil_2\$",dpower * 5) |> ls
		p[1] = plot(n_lst,
						optimald.(n_lst),
						legend=false,
						color=:blue,
		 				linewidth=2,
						top_margin=1mm,
						bottom_margin=2mm,
						left_margin=2mm,
						right_margin=2mm,
						xtickfontsize=xtickfontsize,
						ytickfontsize=ytickfontsize,
						guidefontsize=guidefontsize,
						title=title_str,
						titlefontsize=titlefontsize)
		xlabel!(xlab_str)
		# -- MISE
		xlab_str = L"\tau"
		ylab_str = L"n"
		title_str = "\$(n, \\tau)\\mapsto M_{d_n,\\tau}\\left(h_{\\textrm{bind}}^*(\\tau)\\right)\$" |> ls
		p[2] = contour(τ_lst, n_lst,
						(x, y) -> mise_bindd(y, x),
						seriescolor=mycgrad,
		 				fill=true,
						top_margin=3mm,
						bottom_margin=2mm,
						left_margin=2mm,
						right_margin=2mm,
						xtickfontsize=xtickfontsize,
						ytickfontsize=ytickfontsize,
						guidefontsize=guidefontsize,
						title=title_str,
						titlefontsize=titlefontsize)
		xlabel!(xlab_str)
		ylabel!(ylab_str)


		# ------ τ_n
		d_lst = 2:2:50
		n_lst = 10:1000
		mise_bindτ(d::Integer, n::Integer) = _mise_bound_bind(d, optimalτ(n), n, βf, low, high)

		# -- optimal τ_n
		xlab_str = L"n"
		title_str = @sprintf("\$n\\mapsto \\tau_n=\\lfloor n^{%d/5}\\rceil\$", τpower * 5) |> ls
		p[3] = plot(n_lst,
						optimalτ.(n_lst),
						legend=false,
						color=:blue,
		 				linewidth=2,
						top_margin=1mm,
						bottom_margin=2mm,
						left_margin=2mm,
						right_margin=2mm,
						xtickfontsize=xtickfontsize,
						ytickfontsize=ytickfontsize,
						guidefontsize=guidefontsize,
						title=title_str,
						titlefontsize=titlefontsize)
		xlabel!(xlab_str)
		# -- MISE τ_n
		xlab_str = L"d"
		ylab_str = L"n"
		title_str = "\$(n, d)\\mapsto M_{d,\\tau_n}\\left(h_{\\textrm{bind}}^*(\\tau_n)\\right)\$" |> ls
		p[4] = contour(d_lst, n_lst,
						(x, y) -> mise_bindτ(x, y),
						seriescolor=mycgrad,
		 				fill=true,
						top_margin=3mm,
						bottom_margin=2mm,
						left_margin=2mm,
						right_margin=2mm,
						xtickfontsize=xtickfontsize,
						ytickfontsize=ytickfontsize,
						guidefontsize=guidefontsize,
						title=title_str,
						titlefontsize=titlefontsize)
		xlabel!(xlab_str)
		ylabel!(ylab_str)

		# ------ d_n, τ_n
		mise_binddτ(n::Integer) = _mise_bound_bind(optimald(n), optimalτ(n), n, βf, low, high)
		bandwidthdτ(n::Integer) = _hoptimal_bind(optimalτ(n), low, high)
		# -- bandwidth d_n, τ_n
		xlab_str = L"n"
		title_str = "\$n\\mapsto h_{\\textrm{bind}}(\\tau_n)\$" |> ls
		p[5] = plot(n_lst,
						bandwidthdτ.(n_lst),
						legend=false,
						color=:blue,
		 				linewidth=2,
						top_margin=1mm,
						bottom_margin=2mm,
						left_margin=2mm,
						right_margin=2mm,
						xtickfontsize=xtickfontsize,
						ytickfontsize=ytickfontsize,
						guidefontsize=guidefontsize,
						title=title_str,
						titlefontsize=titlefontsize)
		xlabel!(xlab_str)
		# -- MISE d_n, τ_n
		n_lst = 10:1:100
		xlab_str = L"n"
		title_str = "\$n\\mapsto M_{d_n,\\tau_n}\\left(h_{\\textrm{bind}}^*(\\tau_n)\\right)\$" |> ls
		y = mise_binddτ.(n_lst)
		# fit transformed linear model
		x = [ones(length(n_lst), 1) log.(collect(n_lst))]
		coeffhat = x \ (log.(y))
		println("coeffhat = ", coeffhat)
		mise_hat(n::Integer) = exp(coeffhat[1, 1])*n^coeffhat[2, 1]
		y_fitted = mise_hat.(n_lst)
		p[6] = plot(n_lst,
					[y y_fitted],
					legend=(0.55, 0.71),
					linestyle=[:solid :dash],
					label=["" @sprintf("\$O\\left(n^{%.2f}\\right)\\textrm{ fitted}\$", coeffhat[2,1])],
					color=[:red :blue],
					linewidth=2,
					top_margin=3mm,
					bottom_margin=2mm,
					left_margin=2mm,
					right_margin=2mm,
					xtickfontsize=xtickfontsize,
					ytickfontsize=ytickfontsize,
					guidefontsize=guidefontsize,
					legendfontsize=legendfontsize,
					title=title_str,
					titlefontsize=titlefontsize)
		xlabel!(xlab_str)

		# main plot
		p_main = plot(p..., layout=@layout([a b;c d; e f]))
		display(p_main)
		savefig(joinpath(fig_path, @sprintf("plot_optimal_parameters_bind_%s.pdf", name)))
	end
end



function plot_optimal_parameters_bind_he()
	badexample = false
	ls(s::String) = LaTeXString(s)
	n_plots = 6
	p = Array{Plots.Plot{backend() |> typeof}, 1}(undef, n_plots)
	# distribution parameters
	βf = 3.0
	low = 0.0
	# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#NOTE: if high_dict, low, βf, dpower, τpower changes, LaTeX code mustb e chenged !!!!!!! #NOTE
	# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	high_dict = Dict("maintext"=> 100,
					"appendix_1"=>0.1,
					"appendix_2"=>1,
					"appendix_3"=>10,
					"appendix_4"=>100,
					"appendix_5"=>1000,
					"appendix_6"=>10000000) #NOTE: changed
	# optimal parameters
	(τpower, dvalue, figname_append) = if badexample
		(2/5, 100, "_bad_example")
	else
		(1/5, 14, "")
	end
	optimald(n::Integer) = dvalue
	optimalτ(n::Integer) = round(n^τpower)
	# optimalτ(n::Integer) = round((high-low)^2 * n^τpower)

	# plot options
	xtickfontsize = 7
	ytickfontsize = 7
	legendfontsize = 9
	guidefontsize = 11
	titlefontsize = 11
	GR.setcolormap(46)
	mycgrad = cgrad([RGB(r...) for r=eachrow(GR.colormap())])

	for (name, high) in high_dict
		# ------ d_n
		τ_lst = 2:1000
		n_lst = 10:2:1000
		mise_bindd(n::Integer, τ::Real) = _mise_bound_bind(optimald(n), τ, n, βf, low, high)

		# -- optimal d_n
		xlab_str = L"n"
		# title_str = @sprintf("\$S=[%.0f, %.2f],~\\beta_f=%.0f;~~~n\\mapsto d_n=\\lfloor n^{%d/5}\\rceil_2\$", low, high, βf, dpower * 5) |> ls
		title_str = @sprintf("\$n\\mapsto d_n=%d\$", dvalue) |> ls
		p[1] = plot(n_lst,
						optimald.(n_lst),
						legend=false,
						color=:blue,
		 				linewidth=2,
						top_margin=1mm,
						bottom_margin=2mm,
						left_margin=2mm,
						right_margin=2mm,
						xtickfontsize=xtickfontsize,
						ytickfontsize=ytickfontsize,
						guidefontsize=guidefontsize,
						title=title_str,
						titlefontsize=titlefontsize)
		xlabel!(xlab_str)
		# -- MISE
		xlab_str = L"\tau"
		ylab_str = L"n"
		title_str = "\$(n, \\tau)\\mapsto M_{d_n,\\tau}\\left(h_{\\textrm{bind}}^*(\\tau)\\right)\$" |> ls
		p[2] = contour(τ_lst, n_lst,
						(x, y) -> mise_bindd(y, x),
						seriescolor=mycgrad,
		 				fill=true,
						top_margin=3mm,
						bottom_margin=2mm,
						left_margin=2mm,
						right_margin=2mm,
						xtickfontsize=xtickfontsize,
						ytickfontsize=ytickfontsize,
						guidefontsize=guidefontsize,
						title=title_str,
						titlefontsize=titlefontsize)
		xlabel!(xlab_str)
		ylabel!(ylab_str)


		# ------ τ_n
		d_lst = 2:2:100
		n_lst = 10:1000
		mise_bindτ(d::Integer, n::Integer) = _mise_bound_bind(d, optimalτ(n), n, βf, low, high)

		# -- optimal τ_n
		xlab_str = L"n"
		title_str = @sprintf("\$n\\mapsto \\tau_n=\\lfloor n^{%d/5}\\rceil\$", τpower * 5) |> ls
		p[3] = plot(n_lst,
						optimalτ.(n_lst),
						legend=false,
						color=:blue,
		 				linewidth=2,
						top_margin=1mm,
						bottom_margin=2mm,
						left_margin=2mm,
						right_margin=2mm,
						xtickfontsize=xtickfontsize,
						ytickfontsize=ytickfontsize,
						guidefontsize=guidefontsize,
						title=title_str,
						titlefontsize=titlefontsize)
		xlabel!(xlab_str)
		# -- MISE τ_n
		xlab_str = L"d"
		ylab_str = L"n"
		title_str = "\$(n, d)\\mapsto M_{d,\\tau_n}\\left(h_{\\textrm{bind}}^*(\\tau_n)\\right)\$" |> ls
		p[4] = contour(d_lst, n_lst,
						(x, y) -> mise_bindτ(x, y),
						seriescolor=mycgrad,
		 				fill=true,
						top_margin=3mm,
						bottom_margin=2mm,
						left_margin=2mm,
						right_margin=2mm,
						xtickfontsize=xtickfontsize,
						ytickfontsize=ytickfontsize,
						guidefontsize=guidefontsize,
						title=title_str,
						titlefontsize=titlefontsize)
		xlabel!(xlab_str)
		ylabel!(ylab_str)

		# ------ d_n, τ_n
		mise_binddτ(n::Integer) = _mise_bound_bind(optimald(n), optimalτ(n), n, βf, low, high)
		bandwidthdτ(n::Integer) = _hoptimal_bind(optimalτ(n), low, high)
		# -- bandwidth d_n, τ_n
		xlab_str = L"n"
		title_str = "\$n\\mapsto h_{\\textrm{bind}}(\\tau_n)\$" |> ls
		p[5] = plot(n_lst,
						bandwidthdτ.(n_lst),
						legend=false,
						color=:blue,
		 				linewidth=2,
						top_margin=1mm,
						bottom_margin=2mm,
						left_margin=2mm,
						right_margin=2mm,
						xtickfontsize=xtickfontsize,
						ytickfontsize=ytickfontsize,
						guidefontsize=guidefontsize,
						title=title_str,
						titlefontsize=titlefontsize)
		xlabel!(xlab_str)
		# -- MISE d_n, τ_n
		n_lst = 10:1:1000
		xlab_str = L"n"
		title_str = "\$n\\mapsto M_{d_n,\\tau_n}\\left(h_{\\textrm{bind}}^*(\\tau_n)\\right)\$" |> ls
		y = mise_binddτ.(n_lst)
		# fit transformed linear model
		x = [ones(length(n_lst), 1) log.(collect(n_lst))]
		coeffhat = x \ (log.(y))
		println("coeffhat = ", coeffhat)
		#d_n=n, τ_n=n^(2/5) ⟹ coeffhat=[c, -1.73]. When (2/5)↑ ⟹ -1.73↓ but c↑; when (2/5)↓ ⟹ -1.73↑ and c↑
		# so that n^(2/5) appears best. TODO: But HOW CAN τ_n AFFECT c??
		mise_hat(n::Integer) = exp(coeffhat[1, 1])*n^coeffhat[2, 1] # optimalτ(n) / optimald(n)
		y_fitted = mise_hat.(n_lst)
		p[6] = plot(n_lst,
					[y y_fitted],
					legend=(0.55, 0.71),
					linestyle=[:solid :dash],
					label=["" @sprintf("\$O\\left(n^{%.2f}\\right)\\textrm{ fitted}\$", coeffhat[2,1])],
					color=[:red :blue],
					linewidth=2,
					top_margin=3mm,
					bottom_margin=2mm,
					left_margin=2mm,
					right_margin=2mm,
					xtickfontsize=xtickfontsize,
					ytickfontsize=ytickfontsize,
					guidefontsize=guidefontsize,
					legendfontsize=legendfontsize,
					title=title_str,
					titlefontsize=titlefontsize)
		xlabel!(xlab_str)

		# main plot
		p_main = plot(p..., layout=@layout([a b;c d; e f]))
		display(p_main)
		savefig(joinpath(fig_path, @sprintf("plot_optimal_parameters_bind_he%s_%s.pdf", figname_append, name)))
	end
end
