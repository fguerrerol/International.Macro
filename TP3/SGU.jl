####################################
### Schmitt-Grohé y Uribe (2016) ###
####################################

# Economía pequeña y abierta con rigideces de salario
# PROBLEMA DESCENTRALIZADO

# Listamos los paquetes
#using Pkg
#Pkg.add("ColorSchemes")
using Optim, Interpolations, Printf, LinearAlgebra, PlotlyJS, ColorSchemes, Distributions
using QuantEcon: tauchen

print("Economía pequeña y abierta con rigideces de salario. Inspirado por Schmitt-Grohé y Uribe (2016)\n")
print("Cargando códigos…")

# Declaro un tipo abstracto
abstract type SOE
end

# A la estructura, la declaramos como un subtipo del tipo abrascto SOE
# Ejemplo: Float64 <: Real --> eso te da true

# SOEwr: SOE con rigideces de salario
struct SOEwr <: SOE
	
	# Diccionario de parámetros
	pars::Dict{Symbol, Float64}

	# Grilla para el ahorro
	agrid::Vector{Float64}
	# Grilla para la productividad
	zgrid::Vector{Float64}
	# Matriz de transición para la productividad (depende de su valor pasado)
	Pz::Matrix{Float64}

	# Función de valor: depende de (a, A, z)
	v::Array{Float64, 3}
	# Función de política del consumo: depende de (a, A, z)
	gc::Array{Float64, 3}
	# Función de política del ahorro: depende de (a, A, z)
	ga::Array{Float64, 3}
	
	# Precio de los no transables: depende de (A, z)
	pN::Array{Float64, 2}
	# Salario: depende de (A, z)
	w::Array{Float64, 2}
	# Ley de movimiento del ahorro (A grande - a chica)
	Ap::Array{Float64, 2}
	# Producto: depende de (A, z)
	Y::Array{Float64, 2}
end

# Definimos el constructor
function SOEwr(; β = 0.97, γ = 2, r = 0.02, ϖN = 0.55, η = 1/0.83-1, α = 0.67, 
	wbar = 0.8, ρz = 0.945, σz = 0.025, Na = 40, Nz = 21, amin = -1, amax = 4)

	ϖT = 1 - ϖN

	pars = Dict(:β=>β, :γ=>γ, :r=>r, :ϖN=>ϖN, :ϖT=>ϖT, :η=>η, :α=>α, :wbar=>wbar)

	agrid = cdf.(Beta(2,1), range(0,1, length=Na))
	agrid = amin .+ (amax-amin)*agrid

	# Ahora el proceso de Markov es sobre la productividad y no sobre el producto
	zchain = tauchen(Nz, ρz, σz, 0, 3)
	zgrid = exp.(zchain.state_values)
	Pz = zchain.p
	
	# Inicializo las funciones de valor y las de política
	v = ones(Na, Na, Nz)
	gc = ones(Na, Na, Nz)
	ga = ones(Na, Na, Nz)

	pN = ones(Na, Nz)
	w = ones(Na, Nz)
	Ap = [av for av in agrid, zv in zgrid]
	
	# Producto:
	Y  = [exp(zv) for av in agrid, zv in zgrid]

	return SOEwr(pars, agrid, zgrid, Pz, v, gc, ga, pN, w, Ap, Y)
end

# Esta función te devuelve el precio de la canasta de bienes (compuesta por bien transable y no transable)
# El precio de los transables y de los no transables es un dato
# Ver slide 31 de 3_tnt-rigideces
price_index(pN, sw::SOE) = price_index(pN, 1, sw)
function price_index(pN, pT, sw::SOE)
	ϖN, ϖT, η = (sw.pars[sym] for sym in (:ϖN, :ϖT, :η))

	return (ϖN^(1/(1+η)) * pN^(η/(1+η)) + ϖT^(1/(1+η)) * pT^(η/(1+η)))^((1+η)/η)
end

# Función de utilidad
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

# Función que calcula el valor esperado de Vp
function expect_v(apv, Apv, pz, itp_v, sw::SOE)
	Ev = 0.0
	for (jzp, zpv) in enumerate(sw.zgrid)
		prob = pz[jzp]
		Ev += prob * itp_v(apv, Apv, zpv)
	end
	return Ev
end

# Función de la restricción presupuestaria
budget_constraint(apv, av, yv, r, pCv) = ( yv + av - apv/(1+r) ) / pCv

# Función que calcula la función de Bellman
function eval_value(apv, av, yv, Apv, pz, pCv, itp_v, sw::SOE)
	c = budget_constraint(apv, av, yv, sw.pars[:r], pCv)
	u = utility(c, sw)

	Ev = expect_v(apv, Apv, pz, itp_v, sw)

	return u + sw.pars[:β] * Ev
end

# Función que maximiza a la Ecuación de Bellman
function optim_value(av, yv, Apv, pz, pCv, itp_v, sw::SOE)

	# La función objetivo es la ecuación de bellman
	obj_f(x) = -eval_value(x, av, yv, Apv, pz, pCv, itp_v, sw)
	
	# defino los valores minimos y máximos posibles para el ahorro (según la grilla)
	amin, amax = extrema(sw.agrid)

	# Óptimo
	res = Optim.optimize(obj_f, amin, amax)

	# Nivel de ahorro óptimo
	apv = res.minimizer
	# Función de valor
	v  = -res.minimum
	# Consumo óptimo
	c = budget_constraint(apv, av, yv, sw.pars[:r], pCv)

	return v, apv, c
end


function vf_iter!(new_v, sw::SOE)
	itp_v = interpolate((sw.agrid, sw.agrid, sw.zgrid), sw.v, Gridded(Linear()))

	for jA in eachindex(sw.agrid), jz in eachindex(sw.zgrid)
		# Precio del no transable
		pNv = sw.pN[jA, jz]
		
		# Precio de al canasta de consumo
		pCv = price_index(pNv, sw)

		# Ahorro agregado 
		Apv = sw.Ap[jA, jz]

		# Producto
		yv = sw.Y[jA, jz]

		# Probailidad (transición a partir del nivel de productividad de hoy)
		pz = sw.Pz[jz, :]

		for (ja, av) in enumerate(sw.agrid)
			v, apv, c = optim_value(av, yv, Apv, pz, pCv, itp_v, sw)

			new_v[ja, jA, jz] = v
			sw.ga[ja, jA, jz] = apv
			sw.gc[ja, jA, jz] = c
		end
	end
end

# Iteración de la función de valor
function vfi!(sw::SOE; tol=1e-3, maxiter = 2000)
	# Inicializo la cantidad de iteraciones y la distancia
	iter, dist = 0, 1+tol

	# Inicializo un vector para ir guardando la función de valor
	new_v = similar(sw.v)

	upd_η = 1
	while iter < maxiter && dist > tol
		# Actualizo la cantidad de iteraciones
		iter += 1
		
		# Actualizo la función de valor, el ahorro y el consumo
		vf_iter!(new_v, sw)

		# Calculo la distancia y la norma
		dist = norm(new_v - sw.v) / (1+norm(sw.v))

		norm_v = norm(sw.v)

		print("Iteration $iter: dist = $(round(dist, sigdigits=3)) at ‖v‖ = $(round(norm_v, sigdigits=3))\n")

		sw.v .= sw.v + upd_η * (new_v - sw.v)
	end
	return dist
end

# Demanda de trabajo
function labor_demand(zv, cT, w, sw::SOE)
	α, ϖN, ϖT, η = (sw.pars[sym] for sym in (:α, :ϖN, :ϖT, :η))
	
	# Demanda de trabajo en el sector no transable
	hN = (α/w * ϖN / ϖT)^(1/(1+α*η)) * cT^(1+η)
	
	# Demanda de trabajo en el sector transable
	hT = (zv*α/w)^(1/(1-α))

	return (h = hN+hT, hN = hN, hT = hT)
end

#######################################
# TENGO UNA PREGUNTA DE ESTA FUNCIÓN: #
#######################################
# Encuentra el salario óptimo y las demandas de trabajo
function find_w(zv, cT, wbar, sw::SOE)

	# Dado zv, cT y xbar, te devuelve las demandas de trabajo
	dem = labor_demand(zv, cT, wbar, sw) 
	
	# Trabajo en no transable, en transable y trabajo total
	hN, hT, H = dem.hN, dem.hT, dem.h

	# Si la demanda de trabajo total es menor a 1 (salario menor al de wbar pero eso no puede pasar):
	if H < 1
		wopt = wbar
	# Si la demanda de trabajo es mayor a la oferta (inelastica: 1), el salario es mayor:
	else
		# La producción de no transables:
		f(w) = (labor_demand(zv, cT, w, sw).h - 1)^2
		
		# Calculo el óptimo de la ecuación de Bellman
		res = Optim.optimize(f, wbar, max(2*wbar, 3))
		
		# Salario óptimo
		wopt = res.minimizer
		
		# Demanda de trabajo: NO DEBERÍA SER WOPT?
		dem = labor_demand(zv, cT, wbar, sw) 
		hN, hT = dem.hN, dem.hT
	end

	return hN, hT, wopt
end

# Diferencia entre el precio de no transables anterior y el pN que surge
function diff_pN(pNv, pcC, zv, sw::SOE)
	α, ϖN, ϖT, η, wbar = (sw.pars[sym] for sym in (:α, :ϖN, :ϖT, :η, :wbar))

	# Precio de la canasta de consumo
	pCv = price_index(pNv, sw)

	# El consumo es igual al ratio del gasto en consumo sobre el precio del consumo
	C = pcC / pCv

	# El consumo de transables es igual al consumo total por la fracción
	cT = C * ϖT * (pCv)^η      # cᵢ = ϖᵢ (pᵢ/p)^(-η) C

	# Para esos consumos, la demanda de trabajo y el salario óptimo
	hN, hT, wopt = find_w(zv, cT, wbar, sw)

	# Producción
	yN = hN^α
	yT = zv * hT^α

	# Precio de los no transables dado los valores de la economía
	pN_new = ϖN / ϖT * (cT/yN)^(1+η)

	# Valor del outpur en precios de transables
	output = pN_new * yN + yT

	return (F = (pN_new-pNv)^2, y = output, w = wopt)
end

# Iteración del precio de los no transables
function iter_pN!(new_p, sw::SOE; upd_η = 1)

	# Valores mínimo y máximo del precio
	minp = 0.9 * minimum(sw.pN)
	maxp = 1.1 * maximum(sw.pN)

	for (jA, Av) in enumerate(sw.agrid), (jz, zv) in enumerate(sw.zgrid)
		# guess inicial de precio no transable
		pNg = sw.pN[jA, jz]
		# Dado el pNg, eso genera un cT y un h, eso te devuelve un pn
		pcC = sw.gc[jA, jA, jz] * price_index(pNg, sw)
		obj_f(x) = diff_pN(x, pcC, zv, sw).F

		res = Optim.optimize(obj_f, minp, maxp)

		p = sw.pN[jA, jz] * (1-upd_η) + res.minimizer * upd_η
		new_p[jA, jz] = p
		
		eq = diff_pN(p, pcC, zv, sw)
		sw.Y[jA, jz] = eq.y
		sw.w[jA, jz] = eq.w
	end
end

# Iteración de la ley de movimiento del ahorro
function iter_LoM!(sw::SOE; upd_η = 1)
	for jA in eachindex(sw.agrid), jz in eachindex(sw.zgrid)
		sw.Ap[jA, jz] = (1-upd_η) * sw.Ap[jA, jz] + upd_η * sw.ga[jA, jA, jz]
	end
end

# Actualización del equilibrio
function update_eqm!(new_p, sw::SOE; upd_η = 1)
	iter_pN!(new_p, sw)
	iter_LoM!(sw)
	dist = norm(new_p - sw.pN)
	sw.pN .= sw.pN + upd_η * (new_p - sw.pN)
	return dist
end

# Computo del equilibrio
function comp_eqm!(sw::SOE; tol = 1e-3, maxiter = 2000)
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

# Iterador del simulador
function iter_simul!(tt, path, itp_gc, itp_ga, itp_w, itp_Y, itp_pN, At, zt, sw::SOE)
	# ϖN, ϖT = (sw.pars[sym] for sym in (:ϖN, :ϖT))
	ϖN, ϖT, η = (sw.pars[sym] for sym in (:ϖN, :ϖT, :η))

	w  = itp_w(At, zt)
	Y  = itp_Y(At, zt)
	pN = itp_pN(At, zt)
	pC = price_index(pN, sw)

	C = itp_gc(At, At, zt)

	# este renglon tira error: "type SOEwr has no field ϖT"
	# cT = C * sw.ϖT * (pC)^sw.η 
	# cT = C * (1 - sw.ϖN) * (pC)^sw.η 
	cT = C * ϖT * (pC)^η 
	# cN = C * sw.ϖN * (pC/pN)^sw.η
	cN = C * ϖN * (pC/pN)^η

	CA = Y - pC * C

	path[:CA][tt] = CA
	path[:pN][tt] = pN
	path[:w][tt]  = w
	path[:Y][tt]  = Y
	path[:C][tt]  = C
	path[:A][tt]  = At
	path[:z][tt]  = zt

	A_new = itp_ga(At, At, zt)

	amin, amax = extrema(sw.agrid)
	A_new = max(amin, min(amax, A_new))

	ρz, σz = 0.945, 0.025
	ϵ_new = rand(Normal(0,1))
	z_new = exp(ρz * log(zt) + σz * ϵ_new)

	zmin, zmax = extrema(sw.zgrid)
	z_new = max(zmin, min(zmax, z_new))	

	return A_new, z_new
end

# Simulador
function simul(sw::SOE; T = 100)
	path = Dict(key => zeros(T) for key in [:w, :Y, :CA, :C, :pN, :A, :z])
	itp_gc = interpolate((sw.agrid, sw.agrid, sw.zgrid), sw.gc, Gridded(Linear()))
    itp_ga = interpolate((sw.agrid, sw.agrid, sw.zgrid), sw.ga, Gridded(Linear()))

	itp_w = interpolate((sw.agrid, sw.zgrid), sw.w, Gridded(Linear()))
    itp_Y = interpolate((sw.agrid, sw.zgrid), sw.Y, Gridded(Linear()))
    itp_pN = interpolate((sw.agrid, sw.zgrid), sw.pN, Gridded(Linear()))

	A0 = 0.0
	z0 = 1.0
	for tt in 1:T
		A0, z0 = iter_simul!(tt, path, itp_gc, itp_ga, itp_w, itp_Y, itp_pN, A0, z0, sw)
	end
	return path
end

print(" ✓\n")
print("Constructor sw = SOEwr(; β = 0.97, γ = 2, r = 0.02, ϖN = 0.55, η = 1/0.83-1, α = 0.67, wbar = 0.8, ρz = 0.945, σz = 0.025, Na = 40, Nz = 21, amin = -0.5, amax = 10)\n")
print("Loop: comp_eqm!(sw; tol = 1e-3, maxiter = 2000) \n")