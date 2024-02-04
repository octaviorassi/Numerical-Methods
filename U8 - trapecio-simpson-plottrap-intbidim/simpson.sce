// Simpson
clc
function res = simpson(f, a, b, n)
    // Definimos la funcion
    deff("y = f(x)", "y = " + f)
    
    // Creamos los limites de los subintervalos. Si queremos
    // n subintervalos, tendremos n + 1 puntos.
    h = (b - a) / n
    for i = 1 : n + 1
        x(i) = a + h * (i - 1);
    end 
    
    res = f(x(1));
    
    for i = 2 : n
        mod = modulo(i, 2)
        res = res + (4 - mod * 2) * f(x(i))
    end
    
    res = (h / 3) * (res + f(x(n + 1)))
end


//
//res = simpson("x ^ 3", 3, 12, 300);
//disp(res)
//
//res_cuad = simpson("x ^ 3", 3, 12, 2);
//disp(res_cuad)
//disp(res - res_cuad)


