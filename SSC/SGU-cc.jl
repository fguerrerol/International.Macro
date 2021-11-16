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
	# Matriz de transición para la productividad (depende de su valor pasado)
	Pz::Matrix{Float64}

	# Nuevo: introducción de los ξ
	# Grilla para los ξ
	ξgrid::Vector{Float64}
	# Matriz de transición para los ξ (depende de su valor pasado)
	Pξ::Matrix{Float64}

	# Función de valor: depende de (a, A, z, ξ)
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
	wbar = 0.8, ρz = 0.945, ρξ = 0.945, σz = 0.025, σξ = 0.025, Na = 40, Nz = 21, Nξ = 21, 
	amin = -1, amax = 4)

	ϖT = 1 - ϖN

	pars = Dict(:β=>β, :γ=>γ, :r=>r, :ϖN=>ϖN, :ϖT=>ϖT, :η=>η, :α=>α, :wbar=>wbar)

	agrid = cdf.(Beta(2,1), range(0,1, length=Na))
	agrid = amin .+ (amax-amin)*agrid

	# Ahora el proceso de Markov es sobre la productividad y no sobre el producto
	zchain = tauchen(Nz, ρz, σz, 0, 3)
	zgrid = exp.(zchain.state_values)
	Pz = zchain.p

	# Inicializamos el vector de ξ
	# ξgrid = ones(Nz)

	# Ahora el proceso de Markov es sobre los ξ
	ξchain = tauchen(Nξ, ρξ, σξ, 0, 3)
	ξgrid = exp.(ξchain.state_values)
	Pξ = ξchain.p
	
	# Inicializo las funciones de valor y las de política
	v = ones(Na, Na, Nz, Nξ)
	gc = ones(Na, Na, Nz, Nξ)
	ga = ones(Na, Na, Nz, Nξ)

	pN = ones(Na, Nz, Nξ)
	w = ones(Na, Nz, Nξ)
	Ap = [av for av in agrid, zv in zgrid, ξv in ξgrid]
	
	# Producto:
	Y  = [exp(zv) for av in agrid, zv in zgrid, ξv in ξgrid]

	return SOEwr_cc(pars, agrid, zgrid, Pz, ξgrid, Pξ, v, gc, ga, pN, w, Ap, Y)
end


# la función price_index sigue igual
# la función utility sigue igual

# Definimos un nuevo método para la función expect_v
function expect_v(apv, Apv, pz, itp_v, sw::SOEwr_cc)
	Ev = 0.0
	for (jzp, zpv) in enumerate(sw.zgrid)
		ξv = zpv 
		# La ubicación del ξv en el ξgrid
		jξ = findall(ξ -> ξ == ξv, sw.ξgrid)
		# Dado que estoy en ξv, la probabilidad de estar en ξpv (esto es un vector)
		pξ = sw.Pξ[jξ, :]
				
        # Loop para los valores futuros de ξ
		for (jξp, ξpv) in enumerate(sw.ξgrid)
            # Probabilidad de estar en ξpv dado que hoy estoy en ξv
			probξ = pξ[jξp]
            # Valor de ξpv --> pero esto ya sale del enumerate 
			# ξpv = sw.ξgrid[jξp]
			# Probabilidad de pasar a tener otra productividad (pz es un vector dado que estoy en zv)
			probz = pz[jzp]
			# Valor esperado de V
			prob = probξ * probz
			Ev += prob * itp_v(apv, Apv, zpv, ξpv)
		end
	end
	return Ev
end

# la función de la restricción presupuestaria es la misma
# la función eval_value es la misma
# la función optim_value es la misma


# Definimos un nuevo método para la función vfi_iter!
function vf_iter!(new_v, sw::SOEwr_cc)
	# Agregamos el sw.ξgrid como dimensión del itp_v
	itp_v = interpolate((sw.agrid, sw.agrid, sw.zgrid, sw.ξgrid), sw.v, Gridded(Linear()))


	# Agregamos una iteración sobre ξ
	#for jA in eachindex(sw.agrid), jz in eachindex(sw.zgrid)
	for jA in eachindex(sw.agrid), (jz, zpv) in enumerate(sw.zgrid), (jξ, ξpv) in enumerate(sw.ξgrid)
		
		# Precio del no transable
		pNv = sw.pN[jA, jz, jξ]
		
		# Precio de al canasta de consumo
		pCv = price_index(pNv, sw)

		# Ahorro agregado 
		Apv = sw.Ap[jA, jz, jξ]

		# Producto
		yv = sw.Y[jA, jz, jξ]

		# Probailidad (transición a partir del nivel de productividad de hoy)
		pz = sw.Pz[jz, :]

		for (ja, av) in enumerate(sw.agrid)
			v, apv, c = optim_value(av, yv, Apv, pz, pCv, itp_v, sw)

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

# Defino un nuevo metodo para la función iter_pN!
function iter_pN!(new_p, sw::SOE; upd_η = 1)

	# Valores mínimo y máximo del precio
	minp = 0.9 * minimum(sw.pN)
	maxp = 1.1 * maximum(sw.pN)

	for (jA, Av) in enumerate(sw.agrid), (jz, zv) in enumerate(sw.zgrid), (jξ, ξv) in enumerate(sw.ξgrid)
		# guess inicial de precio no transable
		pNg = sw.pN[jA, jz, jξ]
		# Dado el pNg, eso genera un cT y un h, eso te devuelve un pn
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

# Nuevo método para la iteración de la ley de movimiento del ahorro
function iter_LoM!(sw::SOE; upd_η = 1)
	for jA in eachindex(sw.agrid), jz in eachindex(sw.zgrid), jξ in eachindex(sw.ξgrid)
		sw.Ap[jA, jz, jξ] = (1-upd_η) * sw.Ap[jA, jz, jξ] + upd_η * sw.ga[jA, jA, jz, jξ]
	end
end

# la función update_eqm! es la misma 
# la función comp_eqm! es la misma

print(" ✓\n")
print("Constructor sw = SOEwr_cc(; β = 0.97, γ = 2, r = 0.02, ϖN = 0.55, η = 1/0.83-1, α = 0.67, 
        wbar = 0.8, ρz = 0.945, ρξ = 0.945, σz = 0.025, σξ = 0.025, Na = 40, Nz = 21, Nξ = 21, 
        amin = -1, amax = 4)\n")
print("Loop: comp_eqm!(sw; tol = 1e-3, maxiter = 2000)")