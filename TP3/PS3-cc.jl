#####################
### PROBLEM SET 3 ###
#####################

# Averiguamos donde esta el directorio
pwd()
# Incluimos el código de CC
include("Codigo-julia\\SGU-CC.jl")

# Inicializamos el modelo
sw = SOEwr_cc();

# Corremos el modelo
#comp_eqm!(sw; tol = 1e-3, maxiter = 2000)
comp_eqm!(sw; tol = 1e-3, maxiter = 100)

# Shock de productividad

# Gráfica del producto en función de ξ y del nivel de productividad
ps3_3a =   plot(contour(
            x=sw.ξgrid, # horizontal axis
            y=sw.zgrid, # vertical axis
            z=sw.Y,
            colorscale="inferno",
            contours_start=0,
            contours_size=0.05
            ),
            Layout(title="Curvas de nivel del producto (shock prod.)", xaxis_title="ξ", yaxis_title="Productividad"))

savefig(ps3_3a, "Output/ps3_3a.jpeg")

# Gráfica del ahorro en función de ξ y del nivel de productividad
ps3_3b =   plot(contour(
			x=sw.ξgrid, # horizontal axis
			y=sw.zgrid, # vertical axis
			z=sw.Ap,
			colorscale="inferno",
			contours_start=0,
			contours_size=0.05
			),
			Layout(title="Curvas de nivel del A' (shock prod.)", xaxis_title="ξ", yaxis_title="Productividad"))
				
savefig(ps3_3b, "Output/ps3_3b.jpeg")

# Cuenta corriente en función de ξ y del nivel de productividad
# La cuenta corriente es: sw.Y - 
# pcC = sw.gc * price_index

savefig(ps3_3c, "Output/ps3_3c.jpeg")

# Shock de noticias sobre la productividad
# CREO QUE HAY QUE VOLVER A INICIALIZAR EL MODELO CON OTRO NIVEL DE PARÁMETRO

# Gráfica del producto en función de ξ y del nivel de productividad
ps3_3d =   plot(contour(
            x=sw.ξgrid, # horizontal axis
            y=sw.zgrid, # vertical axis
            z=sw.Y,
            colorscale="inferno",
            contours_start=0,
            contours_size=0.05
            ),
            Layout(title="Curvas de nivel del producto (shock prod.)", xaxis_title="ξ", yaxis_title="Productividad"))

savefig(ps3_3d, "Output/ps3_3d.jpeg")

# Gráfica del ahorro en función de ξ y del nivel de productividad
ps3_3e =   plot(contour(
			x=sw.ξgrid, # horizontal axis
			y=sw.zgrid, # vertical axis
			z=sw.Ap,
			colorscale="inferno",
			contours_start=0,
			contours_size=0.05
			),
			Layout(title="Curvas de nivel del A' (shock prod.)", xaxis_title="ξ", yaxis_title="Productividad"))

savefig(ps3_3e, "Output/ps3_3e.jpeg")

# Cuenta corriente en función de ξ y del nivel de productividad