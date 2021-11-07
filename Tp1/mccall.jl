print("McCall, J. J. 'Economics of Information and Job Search' The Quarterly Journal of Economics, 1970, vol. 84, issue 1, 113-126\n")
print("Loading codes… ")

using Distributions, LinearAlgebra, StatsBase, PlotlyJS
export @L_str, AbstractLayout, AbstractTrace, Cycler, GenericTrace, Inset, Layout, Plot, PlotConfig, PlotlyBase, PlotlyFrame, PlotlyJS, Shape, Spec, Subplots, Template, add_hline!, add_hrect!, add_layout_image!, add_recession_bands!, add_shape!, add_trace!, add_vline!, add_vrect!, addtraces, addtraces!, attr, bar, barpolar, box, candlestick, carpet, choropleth, choroplethmapbox, circle, colors, cone, contour, contourcarpet, dataset, deletetraces, deletetraces!, densitymapbox, download_image, extendtraces, extendtraces!, fork, frame, funnel, funnelarea, heatmap, heatmapgl, histogram, histogram2d, histogram2dcontour, hline, icicle, image, indicator, isosurface, json, line, list_datasets, make_subplots, mesh3d, mgrid, movetraces, movetraces!, ohlc, parcats, parcoords, path, pie, plot, pointcloud, prependtraces, prependtraces!, purge!, react, react!, rect, redraw, redraw!, relayout, relayout!, restyle, restyle!, sankey, savefig, savejson, scatter, scatter3d, scattercarpet, scattergeo, scattergl, scattermapbox, scatterpolar, scatterpolargl, scatterternary, splom, stem, streamtube, sunburst, surface, table, templates, to_image, treemap, update, update!, update_annotations!, update_geos!, update_images!, update_mapboxes!, update_polars!, update_scenes!, update_shapes!, update_ternaries!, update_xaxes!, update_yaxes!, violin, vline, volume, waterfall

mutable struct McCall
	β::Float64
	γ::Float64
	θ::Float64
	b::Float64


	wgrid::Vector{Float64}
	pw::Vector{Float64}

	w_star::Float64
	v::Vector{Float64}
end

function McCall(;
	β = 0.96,
	γ = 0,
	b = 1,
	μw = 1,
	θ = 0,	# Inicializando
	σw = 0.05,
	wmin = 0.5,
	wmax = 2,
	Nw = 50)

	wgrid = range(wmin, wmax, length=Nw)

	w_star = first(wgrid)

	d = Normal(μw, σw)

	pw = [pdf(d, wv) for wv in wgrid]

	pw = pw / sum(pw)

	v = zeros(Nw)

	return McCall(β, γ,θ,b,wgrid, pw, w_star, v)
end

function u(c, mc::McCall)
	γ = mc.γ

	if γ == 1
		return log(c)
	else
		return c^(1-γ) / (1-γ)
	end
end

function R(w, mc::McCall)
	## Valor de aceptar una oferta w: R(w) = u(w) + β R(w)
	β = mc.β
	return u(w, mc) / (1-β)
end

function E_v(mc::McCall)
	## Valor esperado de la función de valor integrando sobre la oferta de mañana
	Ev = 0.0
	if mc.θ>0
		for jwp in eachindex(mc.wgrid)
			Ev +=mc.pw[jwp] * exp(-mc.θ*mc.v[jwp])
		end
		Ev = (-1/ mc.θ)*log(Ev)
		return Ev
	else
		for jwp in eachindex(mc.wgrid)
			Ev += mc.pw[jwp] * mc.v[jwp]
		end
		return Ev
	end
	
end

function update_v(ac, re, EV)
	## Actualizar la función de valor con max(aceptar, rechazar) si EV es falso o usando la forma cerrada con el extreme value si EV es verdadero
	if EV
		χ = 1
		# Probabilidad de aceptar
		prob = exp(ac/χ)/(exp(ac/χ)+exp(re/χ))
		return prob * ac + (1-prob) * re
	else
		return max(ac, re)
	end
end

function vf_iter!(new_v, mc::McCall, flag = 0; EV::Bool=true)
	## Una iteración de la ecuación de Bellman

	# El valor de rechazar la oferta es independiente del estado de hoy
	rechazar = u(mc.b, mc) + mc.β * E_v(mc)
	for (jw, wv) in enumerate(mc.wgrid)
		# El valor de aceptar la oferta sí depende de la oferta de hoy
		aceptar = R(wv, mc)

		# Para una oferta w, v(w) es lo mejor entre aceptar y rechazar
		new_v[jw] = update_v(aceptar, rechazar, EV)

		# El salario de reserva es la primera vez que aceptar es mejor que rechazar
		if flag == 0 && aceptar >= rechazar
			mc.w_star = wv
			flag = 1
		end
	end
end

function vfi!(mc::McCall; maxiter = 2000, tol = 1e-8, verbose::Bool=true)
	dist, iter = 1+tol, 0
	new_v = similar(mc.v)
	#EV = mc.θ > 0 ? false : true
	EV = false
	while dist > tol && iter < maxiter
		iter += 1
		vf_iter!(new_v, mc,EV=EV)
		dist = norm(mc.v - new_v)
		mc.v .= new_v
	end
	if iter == maxiter
		verbose && print("Stopped after ")
	else
		verbose && print("Finished in ")
	end
	verbose && print("$iter iterations.\nDist = $dist\n")

end


function simul(mc::McCall, flag = 0; maxiter = 2000, verbose::Bool=true)
	t = 0
	PESOS = Weights(mc.pw)
	while flag == 0 && t < maxiter
		t += 1
		wt = sample(mc.wgrid, PESOS)
		verbose && print("Salario en el período $t: $wt. ")
		verbose && sleep(0.1)
		wt >= mc.w_star ? flag = 1 : verbose && println("Sigo buscando")
	end
	
	(verbose && flag == 1) && println("Oferta aceptada en $t períodos")
	
	return t
end

function dist_T(mc::McCall, K = 100)
	Tvec = Vector{Int64}(undef, K)
	for jt in eachindex(Tvec)
		Tvec[jt] = simul(mc, verbose=false)
	end
	Tvec
end

print("✓\nConstructor mc = McCall(; β = 0.96, θ = 0,γ = 0, b = 1, μw = 1, σw = 0.05, wmin = 0.5, wmax = 2, Nw = 50\n")
print("Main loop vfi!(mc)\n")