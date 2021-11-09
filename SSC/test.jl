function Cadena(N)
    global A
    A = N;
end


function Accion()
    global A
    print("Hola soy nuevo")
    println()
    print(A)
end

function Erase()
    global A
    A=0;
end

function Envelope(N)
    Cadena(N)
    Accion()
end