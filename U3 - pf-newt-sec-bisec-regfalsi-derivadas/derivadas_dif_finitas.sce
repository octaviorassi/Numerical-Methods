clc // limpia la consola
clear // borra el contenido de la memoria

// Definici´n de la funci´on
function y = f(x)
    y = x * x * x * x;
endfunction

// Calculo de la derivada utilizando diferencias finitas
function y = dfa(f, x, h)
    y = (f(x + h) - f(x)) ./ h;
endfunction

// Implementacion de la derivada por diferencias finitas de orden
// n (forward finite differences)
function y = dfa_ordern(f, x, h, n)
    numerador = 0;
    
    // Construyo el numerador acorde a la expresion para
    // la forward difference de orden n
    for i = (0 : n)
        signo = (-1) ^ (n - i);
        combin = nchoosek(n, i);
        eval = f(x + i * h);
        
        producto = signo * combin * eval;
        numerador = numerador + producto;        
    end
    
    // El denominador de la n-esima derivada es h ^ n
    denominador = h ^ n;    
    
    // Construimos el cociente incremental
    y = numerador / denominador;
endfunction

// EJEMPLO. Derivada segunda de f(x) = x ^ 3
//                              f''(x) = 6x
// funciona para ciertos h, si h es muy pequenio el error es muy grande.
disp(dfa_ordern(f, 10, 0.0001, 2));

