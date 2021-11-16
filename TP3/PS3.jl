
#### Modelo normal
A = SOEwr();
comp_eqm!(A)

### Modelo con el valor más bajo de w 

p3_1a = plot(contour(x =A.agrid,
                     y =A.zgrid,
                     Z = A.Y,
                     colorscale="inferno",
                     contours_start = 0 
                     contours_size = 0.05),
                     Layout(tite ="Curvas de nivel de prodcuto (wbar =0)",
                     xaxis_title="Ahorro",
                     yaxis_title="Productividad"))


B = SOEwr(wbar = 0);
comp_eqm!(B)

p3_1b = plot(contour(x =B.agrid,
                     y =B.zgrid,
                     Z = B.Y,
                     colorscale="inferno",
                     contours_start = 0, 
                     contours_size = 0.05),
                     Layout(tite ="Curvas de nivel de prodcuto (wbar =0)",
                     xaxis_title="Ahorro",
                     yaxis_title="Productividad"))