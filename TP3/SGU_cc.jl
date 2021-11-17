# Este archivo tiene hechas modificaciones para que el modelo tenga shocks de noticias y shocks a la 
# productividad.

using Optim, Interpolations, Printf, LinearAlgebra, PlotlyJS, ColorSchemes, Distributions
using QuantEcon: tauchen

print("Economía pequeña y abierta con rigideces de salario. Inspirado por Schmitt-Grohé y Uribe (2016)\n")
print("Cargando códigos…")



abstract type SOE
end

struct SOEwr_cc <: SOE
	pars::Dict{Symbol, Float64}

	agrid::Vector{Float64}
	zgrid::Vector{Float64}
	ξgrid::Vector{Float64}
	Pξ::Matrix{Float64}

	v::Array{Float64, 4}
	gc::Array{Float64, 4}
	ga::Array{Float64, 4}
	
	pN::Array{Float64, 3}
	w::Array{Float64, 3}
	Ap::Array{Float64, 3}
	Y::Array{Float64, 3}
end

function SOEwr_cc(; β = 0.97, γ = 2, r = 0.02, ϖN = 0.55, η = 1/0.83-1, α = 0.67, wbar = 0.8, ρξ = 0.945, σξ = 0.025, Na = 40, Nξ = 21, amin = -1, amax = 4)

	ϖT = 1 - ϖN
	Nz=Nξ
	pars = Dict(:β=>β, :γ=>γ, :r=>r, :ϖN=>ϖN, :ϖT=>ϖT, :η=>η, :α=>α, :wbar=>wbar)

	agrid = cdf.(Beta(2,1), range(0,1, length=Na))
	agrid = amin .+ (amax-amin)*agrid

	ξchain= tauchen(Nξ, ρξ, σξ, 0, 3)
	zgrid = exp.(ξchain.state_values)
 	ξgrid = zgrid
	Pξ = ξchain.p

	v = ones(Na, Na, Nz, Nξ)
	gc = ones(Na, Na, Nz, Nξ)
	ga = ones(Na, Na, Nz, Nξ)

	pN = ones(Na, Nz, Nξ)
	w = ones(Na, Nz, Nξ)
	Ap = [av for av in agrid, zv in zgrid, ξv in ξgrid]
	Y  = [exp(zv) for av in agrid, zv in zgrid, ξv in ξgrid]

	return SOEwr_cc(pars, agrid, zgrid, ξgrid, Pξ, v, gc, ga, pN, w, Ap, Y)
end


price_index(pN, sw::SOE) = price_index(pN, 1, sw)
function price_index(pN, pT, sw::SOE)
	ϖN, ϖT, η = (sw.pars[sym] for sym in (:ϖN, :ϖT, :η))
	return (ϖN^(1/(1+η)) * pN^(η/(1+η)) + ϖT^(1/(1+η)) * pT^(η/(1+η)))^((1+η)/η)
end

function utility(c, sw::SOE)
	γ = sw.pars[:γ]
	cmin = 1e-3
	if c < cmin
		return utility(cmin,sw) + (c-cmin) * (cmin)^-γ
	else
		γ == 1 && return log(c)
		return c^(1-γ)/(1-γ)
	end
end

function expect_v(apv, Apv, zpv, pξ, itp_v, sw::SOE)
	Ev = 0.0
	for (jξp,ξpv) in enumerate(sw.ξgrid)
		prob = pξ[jξp]
		Ev += prob * itp_v(apv, Apv, zpv, ξpv)
	end
	return Ev
end


budget_constraint(apv, av, yv, r, pCv) = ( yv + av - apv/(1+r) ) / pCv

function eval_value(apv, av, yv, Apv, zpv, pξ, pCv, itp_v, sw::SOE)
	c = budget_constraint(apv, av, yv, sw.pars[:r], pCv)
	u = utility(c, sw)

	Ev = expect_v(apv, Apv, zpv, pξ, itp_v, sw)

	return u + sw.pars[:β] * Ev
end

function optim_value(av, yv, Apv, zpv, pξ, pCv, itp_v, sw::SOE)

	obj_f(x) = -eval_value(x, av, yv, Apv, zpv, pξ, pCv, itp_v, sw)
	amin, amax = extrema(sw.agrid)

	res = Optim.optimize(obj_f, amin, amax)

	apv = res.minimizer
	v  = -res.minimum

	c = budget_constraint(apv, av, yv, sw.pars[:r], pCv)

	return v, apv, c
end

function vf_iter!(new_v, sw::SOE)
	itp_v = interpolate((sw.agrid, sw.agrid, sw.zgrid, sw.ξgrid), sw.v, Gridded(Linear()))

	for jA in eachindex(sw.agrid), jz in eachindex(sw.zgrid), jξ in eachindex(sw.ξgrid)
		pNv = sw.pN[jA, jz, jξ]
		pCv = price_index(pNv, sw)

		Apv = sw.Ap[jA, jz, jξ]
		yv = sw.Y[jA, jz, jξ]

        zpv = sw.zgrid[jξ]
		pξ = sw.Pξ[jξ, :]

		for (ja, av) in enumerate(sw.agrid)
			v, apv, c = optim_value(av, yv, Apv, zpv, pξ, pCv, itp_v, sw)

			new_v[ja, jA, jz, jξ] = v
			sw.ga[ja, jA, jz, jξ] = apv
			sw.gc[ja, jA, jz, jξ] = c
		end
	end
end

function vfi!(sw::SOE; tol=1e-4, maxiter = 2000)
	iter, dist = 0, 1+tol
	new_v = similar(sw.v)

	upd_η = 1
	while iter < maxiter && dist > tol
		iter += 1
		
		vf_iter!(new_v, sw)

		dist = norm(new_v - sw.v) / (1+norm(sw.v))

		norm_v = norm(sw.v)

		print("Iteration $iter: dist = $(round(dist, sigdigits=3)) at ‖v‖ = $(round(norm_v, sigdigits=3))\n")

		sw.v .= sw.v + upd_η * (new_v - sw.v)
	end
	return dist
end

function labor_demand(zv, cT, w, sw::SOE)
	α, ϖN, ϖT, η = (sw.pars[sym] for sym in (:α, :ϖN, :ϖT, :η))

	hN = (α/w * ϖN / ϖT)^(1/(1+α*η)) * cT^(1+η)
	hT = (zv*α/w)^(1/(1-α))

	return (h = hN+hT, hN = hN, hT = hT)
end

function find_w(zv, cT, wbar, sw::SOE)
	dem = labor_demand(zv, cT, wbar, sw) 
	hN, hT, H = dem.hN, dem.hT, dem.h

	if H < 1
		wopt = wbar
	else
		f(w) = (labor_demand(zv, cT, w, sw).h - 1)^2
		res = Optim.optimize(f, wbar, max(2*wbar, 3))
		wopt = res.minimizer
		dem = labor_demand(zv, cT, wbar, sw) 
		hN, hT = dem.hN, dem.hT
	end

	return hN, hT, wopt
end

function diff_pN(pNv, pcC, zv, sw::SOE)
	α, ϖN, ϖT, η, wbar = (sw.pars[sym] for sym in (:α, :ϖN, :ϖT, :η, :wbar))

	pCv = price_index(pNv, sw)
	C = pcC / pCv

	cT = C * ϖT * (pCv)^η      

	hN, hT, wopt = find_w(zv, cT, wbar, sw)

	yN = hN^α
	yT = zv * hT^α

	pN_new = ϖN / ϖT * (cT/yN)^(1+η)

	output = pN_new * yN + yT

	return (F = (pN_new-pNv)^2, y = output, w = wopt)
end


function iter_pN!(new_p, sw::SOE; upd_η = 1)
	minp = 0.9 * minimum(sw.pN)
	maxp = 1.1 * maximum(sw.pN)
	for (jA, Av) in enumerate(sw.agrid), (jz, zv) in enumerate(sw.zgrid), (jξ, ξv) in enumerate(sw.ξgrid)
		pNg = sw.pN[jA, jz, jξ]
		pcC = sw.gc[jA, jA, jz, jξ] * price_index(pNg, sw)
		obj_f(x) = diff_pN(x, pcC, zv, sw).F

		res = Optim.optimize(obj_f, minp, maxp)

		p = sw.pN[jA, jz, jξ] * (1-upd_η) + res.minimizer * upd_η
		new_p[jA, jz, jξ] = p
		
		eq = diff_pN(p, pcC, zv, sw)
		sw.Y[jA, jz, jξ] = eq.y
		sw.w[jA, jz, jξ] = eq.w
	end
end

function iter_LoM!(sw::SOE; upd_η = 1)
	for jA in eachindex(sw.agrid), jz in eachindex(sw.zgrid), jξ in eachindex(sw.ξgrid)
		sw.Ap[jA, jz, jξ] = (1-upd_η) * sw.Ap[jA, jz, jξ] + upd_η * sw.ga[jA, jA, jz, jξ]
	end
end

function update_eqm!(new_p, sw::SOE; upd_η = 1)
	iter_pN!(new_p, sw)
	iter_LoM!(sw)
	dist = norm(new_p - sw.pN)
	sw.pN .= sw.pN + upd_η * (new_p - sw.pN)
	return dist
end

function comp_eqm!(sw::SOE; tol = 1e-6, maxiter = 2000)
	iter, dist = 0, 1+tol
	new_p = similar(sw.pN)
	tol_vfi = 1e-2

	while dist > tol && iter < maxiter
		iter += 1
		print("Outer Iteration $iter (tol = $(round(tol_vfi, sigdigits=3)))\n")
		dist_v = vfi!(sw, tol = tol_vfi)

		norm_p = norm(sw.pN)
		dist_p = update_eqm!(new_p, sw) / max(1,norm_p)

		dist = max(dist_p, dist_v)

		print("After $iter iterations, dist = $(round(dist_p, sigdigits=3)) at ‖pN‖ = $(round(norm_p, sigdigits=3))\n\n")

		tol_vfi = max(tol, tol_vfi * 0.9)
	end
end

function iter_simul!(tt, path, itp_gc, itp_ga, itp_w, itp_Y, itp_pN, At, zt, ξt, sw::SOE)
	ϖN, ϖT, η = (sw.pars[sym] for sym in (:ϖN, :ϖT, :η))
	w  = itp_w(At, zt, ξt)
	Y  = itp_Y(At, zt, ξt)
	pN = itp_pN(At, zt, ξt)
	pC = price_index(pN, sw)

	C = itp_gc(At, At, zt, ξt)

	cT = C * ϖT * (pC)^η
	cN = C * ϖN * (pC/pN)^η

	CA = Y - pC * C

	path[:CA][tt] = CA
	path[:pN][tt] = pN
	path[:w][tt]  = w
	path[:Y][tt]  = Y
	path[:C][tt]  = C
	path[:A][tt]  = At
	path[:z][tt]  = zt
	path[:ξ][tt]  = ξt

	A_new = itp_ga(At, At, zt, ξt)

	amin, amax = extrema(sw.agrid)
	A_new = max(amin, min(amax, A_new))

	ρξ, σξ = 0.945, 0.025
	ϵ_new = rand(Normal(0,1))
	ξ_new = exp(ρξ * log(ξt) + σξ * ϵ_new)

	ξmin, ξmax = extrema(sw.ξgrid)
	ξ_new = max(ξmin, min(ξmax, ξ_new))
	z_new = ξt

	return A_new, z_new, ξ_new
end

function simul(sw::SOE; T = 100)
	path = Dict(key => zeros(T) for key in [:w, :Y, :CA, :C, :pN, :A, :z, :ξ])
	itp_gc = interpolate((sw.agrid, sw.agrid, sw.zgrid, sw.ξgrid), sw.gc, Gridded(Linear()))
    itp_ga = interpolate((sw.agrid, sw.agrid, sw.zgrid, sw.ξgrid), sw.ga, Gridded(Linear()))

	itp_w = interpolate((sw.agrid, sw.zgrid, sw.ξgrid), sw.w, Gridded(Linear()))
    itp_Y = interpolate((sw.agrid, sw.zgrid, sw.ξgrid), sw.Y, Gridded(Linear()))
    itp_pN = interpolate((sw.agrid, sw.zgrid, sw.ξgrid), sw.pN, Gridded(Linear()))

	A0 = 0.0
	z0 = 1.0
	ξ0 = 1.0
	for tt in 1:T
		A0, z0, ξ0 = iter_simul!(tt, path, itp_gc, itp_ga, itp_w, itp_Y, itp_pN, A0, z0, ξ0, sw)
	end
	return path
end

print(" ✓\n")
print("Constructor sw = SOEwr_cc(; β = 0.97, γ = 2, r = 0.02, ϖN = 0.55, η = 1/0.83-1, α = 0.67, wbar = 0.8, ρz = 0.945, σz = 0.025, Na = 40, Nz = 21, amin = -0.5, amax = 10)\n")
print("Loop: comp_eqm!(sw; tol = 1e-3, maxiter = 2000)")