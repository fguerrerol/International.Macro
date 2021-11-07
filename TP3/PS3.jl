
#### Modelo normal
A = SOEwr();
comp_eqm!(A)
plot_cons(A)
### Modelo con el valor más bajo de w 

B = SOEwr(wbar = 0);
comp_eqm!(B)
plot_cons(B)