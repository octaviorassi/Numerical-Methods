exec("C:\Users\octav\OneDrive - frsn.utn.edu.ar\LCC\METODOS NUMERICOS\u4\gauss_final.sce") // Para calcular la inversa con Gauss
clc
function y = sumatoria_seidel(A, b, n, x0, x1, i)
    b_term = b(i);
    coef = 1/A(i, i);
        
    sumatoria1= 0;
    for j = 1 : i - 1
            sumatoria1 = sumatoria1 + A(i, j) * x1(j);
    end 
    
    sumatoria2 = 0;
    for j = i + 1 : n
        sumatoria2 = sumatoria2 + A(i, j) * x0(j);
    end            
                        
    y = coef * (b_term - sumatoria1 - sumatoria2);
endfunction

function [x1, iter] = gauss_seidel(A, b, x0, tol, maxiter)
    n = size(A)(1);
    x1 = zeros(x0);

    // Primer iteracion, para que no falle el norm(x1 - x0)
    for i = 1 : n          
        x1(i) = sumatoria_seidel(A, b, n, x0, x1, i)                 
    end

    iter = 1;       
    while norm(x1 - x0) > tol && iter < maxiter
       x0 = x1;
       
       for i = 1 : n          
            x1(i) = sumatoria_seidel(A, b, n, x0, x1, i)                 
       end       
            
       iter = iter + 1;
    end     
endfunction

// P = I - N_inv * A
function P = iteracion_gauss_seidel(A)
    N = A
    n = size(A)(1)
    
    for i = 1 : n
        for j = i + 1 : n
            N(i, j) = 0
        end
    end
    
    N_inv = gauss_inversa(N)
    
    P = eye(A) - N_inv * A
endfunction

// OTRA MANERA ES USANDO LA CONDICION SUFICIENTE ( Y NO NECESARIA )
// A es diagonal dominante
//                      ==>
// Gauss-Seidel converge para cualquier x0 inicial 
function y = convergencia_gauss_seidel(A)
    // Definimos la matriz N que es el resultado de poner en 0
    // los valores por encima de la diagonal de A (estricto)
    M = iteracion_gauss_seidel(A)
   
    p = max(abs(spec(M)))
    
    y = p < 1
    
    if y then
        printf("Converge para cualquier x0 inicial.\n")
    else
        printf("Convergencia no asegurada.\n")
    end
endfunction

//A = [1 -1 0; -1 2 -1; 0 -1 1.1]
//b = [0 1 0]
//
//[x1, iter] = gauss_seidel(A, b', [51 0 -33]', 0.01, 1000)
//
//disp(b)
//disp(A * x1)
//disp(iter)
//
//[x2, iter] = gauss_seidel(A, b', x1, 0.01, 100)
//
//disp(b)
//disp(A * x2)
//disp(iter)

