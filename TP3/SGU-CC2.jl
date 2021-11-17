####################################
### Schmitt-Grohé y Uribe (2016) ###
####################################

# Economía pequeña y abierta con rigideces de salario
# PROBLEMA DESCENTRALIZADO

# Incluimos el código de Schmitt-Grohé y Uribe (2016): Economía pequeña y abierta con rigideces de salario
include("SGU.jl")

#######################################################
#   Sobre la interpretación de la cuenta corriente    #
#######################################################

print("Economía pequeña y abierta con rigideces de salario. Interpretación de la cuenta corriente. \n")
print("Cargando códigos…")

# En este problema, se agrega un estado adicional, entonces hay que crear un nuevo subtipo de SOE

# SOEwr_cc: SOE con rigideces de salario + cuenta corriente (se agrega shock de noticias sobre la productividad)
struct SOEwr_cc <: SOE
	
	# Diccionario de parámetros
	pars::Dict{Symbol, Float64}

	# Grilla para el ahorro
	agrid::Vector{Float64}
	# Grilla para la productividad
	zgrid::Vector{Float64}

	# Nuevo: introducción de los ξ
	# Grilla para los ξ
	ξgrid::Vector{Float64}
	# Matriz de transición para los ξ (depende de su valor pasado)
	Pξ::Matrix{Float64}

	# Función de valor: depende de (a, A, z, ξ), la vp va a depender de depende de (ap, Ap, ξ, ξp)
	v::Array{Float64, 4}
	# Función de política del consumo: depende de (a, A, z, ξ) 
	gc::Array{Float64, 4}
	# Función de política del ahorro: depende de (a, A, z, ξ)
	ga::Array{Float64, 4}
	
	# Precio de los no transables: depende de (A, z, ξ)
	pN::Array{Float64, 3}
	# Salario: depende de (A, z, ξ)
	w::Array{Float64, 3}
	# Ley de movimiento del ahorro: depende de (A, z, ξ)
	Ap::Array{Float64, 3}
	# Producto: depende de (A, z, ξ)
	Y::Array{Float64, 3}
end

# Definimos el constructor
function SOEwr_cc(; β = 0.97, γ = 2, r = 0.02, ϖN = 0.55, η = 1/0.83-1, α = 0.67, 
	wbar = 0.8, ρz = 0.945, ρξ = 0.945, σz = 0.025, σξ = 0.025, Na = 20, Nz = 11, Nξ = 11, 
	amin = -1, amax = 4)
	
    ϖT = 1 - ϖN

	pars = Dict(:β=>β, :γ=>γ, :r=>r, :ϖN=>ϖN, :ϖT=>ϖT, :η=>η, :α=>α, :wbar=>wbar)

	agrid = cdf.(Beta(2,1), range(0,1, length=Na))
	agrid = amin .+ (amax-amin)*agrid

	# Ahora el proceso de Markov es sobre los ξ
	ξchain = tauchen(Nξ, ρξ, σξ, 0, 3)
	ξgrid = exp.(ξchain.state_values)
	Pξ = ξchain.p
	zgrid = ξgrid

	# Inicializo las funciones de valor y las de política
	v = ones(Na, Na, Nz, Nξ)
	gc = ones(Na, Na, Nz, Nξ)
	ga = ones(Na, Na, Nz, Nξ)

	pN = ones(Na, Nz, Nξ)
	w = ones(Na, Nz, Nξ)
	Ap = [av for av in agrid, zv in zgrid, ξv in ξgrid]
	
	# Producto:
	Y  = [exp(zv) for av in agrid, zv in zgrid, ξv in ξgrid]

	return SOEwr_cc(pars, agrid, zgrid, ξgrid, Pξ, v, gc, ga, pN, w, Ap, Y)
end

# la función price_index sigue igual
# la función utility sigue igual

# Definimos un nuevo método para la función expect_v
function expect_v(apv, Apv, pξ, ξv, itp_v, sw::SOEwr_cc)
	
    # Inicializo el valor esperado
    Ev = zeros(length(sw.ξgrid))
    
    # Hoy, zv y ξv es conocido. Además, ξv = zpv (consigna). Lo que no se es ξpv
	Threads.@threads for i = 1:length(sw.ξgrid)  
		
		# Probabilidad de pasar a tener otra shock de noticias (pξ es un vector dado que estoy en ξv)
		prob = pξ[i]
		Ev[i] = prob * itp_v(apv, Apv, ξv, sw.ξgrid[i])
	end

	return sum(Ev)/length(Ev)
end

# la función de la restricción presupuestaria es la misma

# Defino un nuevo método para la función eval_value: Función que define la Ecuación de Bellman
# Lo que cambia es el input de expect_v: pasa de ser pz a pξ y se agrega ξv
function eval_value(apv, av, yv, Apv, pξ, ξv, pCv, itp_v, sw::SOEwr_cc)
	c = budget_constraint(apv, av, yv, sw.pars[:r], pCv)
	u = utility(c, sw)
	Ev = expect_v(apv, Apv, pξ, ξv, itp_v, sw)
	return u + sw.pars[:β] * Ev
end

# Defino un nuevo método para optim_value: función que maximiza a la Ecuación de Bellman
# Lo que cambia es que pasa de ser pz a pξ y se agrega ξv
function optim_value(av, yv, Apv, pξ, ξv, pCv, itp_v, sw::SOEwr_cc)

	# La función objetivo es la ecuación de bellman
	obj_f(x) = -eval_value(x, av, yv, Apv, pξ, ξv, pCv, itp_v, sw)
	
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

# Defino un nuevo método para vf_iter!
# lo que cambia:
# (1) se agrega sw.ξgrid como knot al interpolador
# (2) se agrega el índice jξ
# (3) se cambia pz por pξ (probabilidad de transición dado que estas en ξv)
function vf_iter!(new_v, sw::SOEwr_cc)
	itp_v = interpolate((sw.agrid, sw.agrid, sw.zgrid, sw.ξgrid), sw.v, Gridded(Linear()))
	
	for jA in eachindex(sw.agrid), jz in eachindex(sw.zgrid), (jξ, ξv) in enumerate(sw.ξgrid)
		# Precio del no transable
		pNv = sw.pN[jA, jz, jξ]
		
		# Precio de al canasta de consumo 
        # --> notar que si bien ahora pNv depende tmb de ξv la función price_index no cambia
		pCv = price_index(pNv, sw)

		# Ahorro agregado 
		Apv = sw.Ap[jA, jz, jξ]

		# Producto
		yv = sw.Y[jA, jz, jξ]

		# Probailidad (transición a partir del nivel de ξv de hoy)
		pξ = sw.Pξ[jξ, :]

		for (ja, av) in enumerate(sw.agrid)
			v, apv, c = optim_value(av, yv, Apv, pξ, ξv, pCv, itp_v, sw)

			new_v[ja, jA, jz, jξ] = v
			sw.ga[ja, jA, jz, jξ] = apv
			sw.gc[ja, jA, jz, jξ] = c
		end
	end
end

# la función vfi! es la misma 
# la función labor_demand es la misma
# la función find_w es la misma
# la función diff_pN es la misma

# Defino un nuevo método para la función iter_pN!
# Cambio: se agrega como índice a jξ
function iter_pN!(new_p, sw::SOEwr_cc; upd_η = 1)

	# Valores mínimo y máximo del precio
	minp = 0.9 * minimum(sw.pN)
	maxp = 1.1 * maximum(sw.pN)

	for (jA, Av) in enumerate(sw.agrid), (jz, zv) in enumerate(sw.zgrid), (jξ, ξv) in enumerate(sw.ξgrid)
		# guess inicial de precio no transable --> no estoy segura que haya que agregar jξ
		pNg = sw.pN[jA, jz, jξ]
		# Dado el pNg, eso genera un cT y un h, eso te devuelve un pn
		pcC = sw.gc[jA, jA, jz, jξ] * price_index(pNg, sw)
        # Definimos la función objetivo: recordar el output de diff_pN: F = (pN_new-pNv)^2
		obj_f(x) = diff_pN(x, pcC, zv, sw).F

		res = Optim.optimize(obj_f, minp, maxp)

		p = sw.pN[jA, jz, jξ] * (1-upd_η) + res.minimizer * upd_η
		new_p[jA, jz, jξ] = p
		
		eq = diff_pN(p, pcC, zv, sw)
		sw.Y[jA, jz, jξ] = eq.y
		sw.w[jA, jz, jξ] = eq.w
	end
end

# Definimos un nuevo método para la función que hace la iteración de la ley de movimiento del ahorro
function iter_LoM!(sw::SOEwr_cc; upd_η = 1)
	for jA in eachindex(sw.agrid), jz in eachindex(sw.zgrid), jξ in eachindex(sw.ξgrid)
		sw.Ap[jA, jz, jξ] = (1-upd_η) * sw.Ap[jA, jz, jξ] + upd_η * sw.ga[jA, jA, jz, jξ]
	end
end

# la función update_eqm! es la misma 
# la función comp_eqm! es la misma 

# Nuevo método del iterador del simulador
function iter_simul!(tt, path, itp_gc, itp_ga, itp_w, itp_Y, itp_pN, At, zt, ξt, sw::SOEwr_cc)
	# ϖN, ϖT = (sw.pars[sym] for sym in (:ϖN, :ϖT))
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

	A_new = itp_ga(At, At, zt, ξt)

	amin, amax = extrema(sw.agrid)
	A_new = max(amin, min(amax, A_new))

	ρz, σz = 0.945, 0.025
	ϵ_new = rand(Normal(0,1))
	z_new = exp(ρz * log(zt) + σz * ϵ_new)

	zmin, zmax = extrema(sw.zgrid)
	z_new = max(zmin, min(zmax, z_new))	

	return A_new, z_new
end

# Nuevo método para el simulador
# El único cambio es que le agregas sw.ξgrid dentro de los knots del interpolador
function simul(sw::SOEwr_cc; T = 100)
	path = Dict(key => zeros(T) for key in [:w, :Y, :CA, :C, :pN, :A, :z])
	itp_gc = interpolate((sw.agrid, sw.agrid, sw.zgrid, sw.ξgrid), sw.gc, Gridded(Linear()))
    itp_ga = interpolate((sw.agrid, sw.agrid, sw.zgrid, sw.ξgrid), sw.ga, Gridded(Linear()))

	itp_w = interpolate((sw.agrid, sw.zgrid, sw.ξgrid), sw.w, Gridded(Linear()))
    itp_Y = interpolate((sw.agrid, sw.zgrid, sw.ξgrid), sw.Y, Gridded(Linear()))
    itp_pN = interpolate((sw.agrid, sw.zgrid, sw.ξgrid), sw.pN, Gridded(Linear()))

	A0 = 0.0
	z0 = 1.0
	for tt in 1:T
		A0, z0 = iter_simul!(tt, path, itp_gc, itp_ga, itp_w, itp_Y, itp_pN, A0, z0, sw)
	end
	return path
end



print(" ✓\n")
print("Constructor sw = SOEwr_cc(; β = 0.97, γ = 2, r = 0.02, ϖN = 0.55, η = 1/0.83-1, α = 0.67, 
        wbar = 0.8, ρz = 0.945, ρξ = 0.945, σz = 0.025, σξ = 0.025, Na = 40, Nz = 21, Nξ = 21, 
        amin = -1, amax = 4)\n")
print("Loop: comp_eqm!(sw; tol = 1e-3, maxiter = 2000)\n")
