// Metodo del trapecio para n subintervalos
// * f es un string que representa el cuerpo de la funcion
// * a, b son los extremos del intervalo
// * n es la cantidad de subintervalos
function res = trapecio(f, a, b, n)
    // Definimos la funcion
    deff("y = f(x)", "y = " + f)
    
    // Creamos los limites de los subintervalos. Si queremos
    // n subintervalos, tendremos n + 1 puntos.
    h = (b - a) / n
    for i = 1 : n + 1
        x(i) = a + h * (i - 1);
    end    
    
    // Aplicamos la regla del trapecio
    res = 1/2 * f(x(1));    
    
    for i = 2 : n
        res = res + f(x(i));
    end    
    
    res = res + 1/2 * f(x(n + 1));
    res = res * h;
endfunction

// Ejemplo
//res = trapecio("x ^ 3", 3, 12, 300);
//disp(res)
//res_real = (12^4 - 3^4)/4; // 5100 aprox
//disp(res - res_real)




