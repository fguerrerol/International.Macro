


    ####  Genero los modelos el cual el primero es con un salario
    ####  Primer modelo con salario Default

    M_1 = SOEwr();

    ### Corro el modelo para obtrner los distintos valores.

    comp_eqm!(M_1);




    M_2 = SOEwr(wbar = 0.0001);

    comp_eqm!(M_1);

    #### Corro el segundo modelo el cual tinee salairo W  muy cercnao a 0,
    #### con 0 el programa tiraba error