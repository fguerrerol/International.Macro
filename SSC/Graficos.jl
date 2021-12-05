plot([scatter(x = eachindex(Results["Contagios"][:,1,2]),y = Results["Contagios"][:,8,jy],
    name= "Capacidad de testeo = $(Results["Betas"][jy])") for jy in eachindex(Results["Betas"])],
    Layout(title ="Proporcion de contagiados dado un  Beta = 0.5 y distintos niveles de testeo", 
    yaxis_title ="Proporcion de contagiados",
    xaxis_title ="Periodos"))


    plot([scatter(x = eachindex(Ab["Restrictos"][:,1,2]),y = Ab["Restrictos"][:,8,jy],
    name= "Beta = $(Ab["Betas"][jy])") for jy in eachindex(Ab["Betas"])],
    Layout(title ="Proporcion de retringidos para distintos niveles de Beta", 
    yaxis_title ="Proporcion de restringidos",
    xaxis_title ="Periodos"))