function y = f(x)
    y = %e ^ x;
endfunction

function y = f_alt(x)
    y = 1 / (%e ^ x);
endfunction

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

// Implementar taylor() que calcule el valor de un polinomio de Taylor 
// de grado n de una funcion f en un punto x

// OBSERVACION. Conviene definir el polinomio alrededor del mismo punto
// es decir, hacer a == x.
function y = taylor(f, n, a, x)
    
    // El polinomio de Taylor de orden n de la funcion f alrededor de
    // un punto a evaluado en un punto x es:
    
//      ———  n          
//      \                f^(k)(a)   *   (x - a) ^ k
//       |               ————————————————
//      /                           k!
//      ———  k = 0 
    
    suma = 0;
    for k = (0 : n)
        fact =  factorial(k);
        deriv_evaluada = dfa_ordern(f, a, 0.0001, k);
        x_term = (x - a) ^ k;
        
        taylor_term = deriv_evaluada * x_term / fact;
        
        suma = suma + taylor_term;        
    end
    
    y = suma;    
endfunction 

// EJ6 - A
// CONSULTAR. La utilidad de Taylor es poder aproximar funciones
// trascendentales como e^x con un polinomio, sin tener que evaluar
// la funcion en si. Aca cuando defino el Taylor, estoy usando 
// la ufncion dfa_ordern que aproxima la derivada haciendo, entre otras
// operaciones, evaluaciones de la funcion f.
// Esta bien?
disp(taylor(f, 10, -2, -2));
disp(taylor(f_alt, 10, 2, 2));
