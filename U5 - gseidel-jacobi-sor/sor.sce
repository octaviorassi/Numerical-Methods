function y = sor_iter(A, b, w, x0, x1, n, i)
    coef = w / A(i, i)    
    unomenosw = (1 - w) * x0(i)    
    b_term = b(i)
    
    sumatoria = 0;
    for j = 1 : i - 1
        sumatoria = sumatoria + A(i, j) * x1(j)
    end
    
    for j = i + 1 : n
        sumatoria = sumatoria + A(i, j) * x0(j)
    end
    
    y = unomenosw + coef * (b_term - sumatoria)    
endfunction


// OBSERVACIONES.
// * Con w = 1, tenemos el metódo de Gauss Seidel.
// * Con 0 < w < 1 es un metódo de subrelajacion.  Permiten obtener la
// convergencia para sistemas donde Gauss-Seidel diverge.
// * Con w > 1m es un metódo de sobrerelajacion (SOR). Aceleran la
// convergencia para sistemas donde Gauss-Seidel ya converge.
function [x1, iter] = sor(A, b, w, x0, tol, maxiter)
    n = size(A)(1)
    
    // Calculo la primer iteracion, para que no falle norm(x1 - x0)
    x1(1 : n) = 0
    for i = 1 : n  
        x1(i) = sor_iter(A, b, w, x0, x1, n, i)
    end
    
    iter = 1    
    while norm(x1 - x0) > tol && iter < maxiter
        x0 = x1;
        
        for i = 1 : n  
            x1(i) = sor_iter(A, b, w, x0, x1, n, i)
        end
        
        iter = iter + 1;       
    end
    
    if iter < maxiter
        printf("\n Se alcanzo la tolerancia deseada. \n")
    else
        printf("\n Se alcanzo el maximo de iteraciones (%i).\n", maxiter)
    end    
endfunction

function w = omega_optimo(A)
    n = size(A)(1);    
    // diag(vector) === la matriz diagonal que tiene a vector como diag.
    T_j = eye(n, n) - diag(1 ./ diag(A)) * A; 
    p = max(abs(spec(T_j)));
   
    w = 2 / (1 + sqrt(1 - p**2))
endfunction

function [x1, iter] = sor_optimo(A, b, x0, tol, maxiter)
    w = omega_optimo(A)
    [x1, iter] = sor(A, b, w, x0, tol, maxiter)
endfunction
