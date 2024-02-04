// Graficador de phi_n
// Dados n + 1 puntos x_i grafica la funcion
//    phi_n(x) = (x - x_0)(x - x_1) ... (x - x_n)
// y traza una linea sobre el error maximo en valor absoluto

// Devuelve el poly phi_n y el error maximo en abs calculado.

function [phi_n, error_max] = graficaPhi(x)
    // Construyo el polinomio
    phi_n = poly(x, "x", "roots")
    
    // Genero el rango de ploteo
    x_min = min(x)
    x_max = max(x)
    
    x_rng = x_min : 0.01 : x_max
    y_rng = horner(phi_n, x_rng)
    
    // Pido los ejes y grafico phi_n
    subplot(1, 2, 1)
    a = gca()
    a.x_location = "origin"
    a.y_location = "origin"
    
    plot(x_rng, y_rng)
    
    // Pido los ejes y grafico abs(phi_n) y la cota del error
    subplot(1, 2, 2)
    a = gca()
    a.x_location = "origin"
    a.y_location = "origin"
    
    y_abs = abs(y_rng)
    error_max = max(y_abs)
    
    plot(x_rng, abs(y_rng))
    plot([x_min, x_max], [error_max, error_max], "red")

    
    // Muestro en pantalla el maximo error
    printf("El error maximo en valor absoluto es: %f\n", error_max)
endfunction

raices =[1, 2, 3, 4]
pol = graficaPhi(raices)
disp(pol)
for i = 1 : 4
    disp(horner(pol, i))
end
