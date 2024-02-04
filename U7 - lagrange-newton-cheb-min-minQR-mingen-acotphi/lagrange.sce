// Lagrange
function pol = L_k(X, k)
    n = length(X)
    pol = poly([1], "x", "coeff");
    x_k = X(k + 1);
    
    for i = 1 : k 
        numerador = poly([-X(i) 1], "x", "coeff")
        denom = x_k - X(i) //  x(i) != x(j)  =>  x(k) - x(i) != 0
        
        pol = pol * numerador / denom
    end    
    
    for i = k + 2 : n
        numerador = poly([-X(i) 1], "x", "coeff")
        denom = x_k - X(i) //  x(i) != x(j)  =>  x(k) - x(i) != 0
        
        pol = pol * numerador / denom
    end
    
endfunction

// disp(L_k([0 1 2], 1)) DEBERIA DAR 2x - x^2

// Asumo X, Y vectores fila
function pol = lagrange(X, Y)
    n = length(X)
    m = length(Y)
    
    if n <> m then
        error("lagrange - X e Y son de distinto tamanio")
        abort
    end
    
    pol = poly([0], "x", "coeff")
    
    for k = 0 : n - 1
        pol = pol + L_k(X, k) * Y(k + 1)
    end    
endfunction

//pol_lag = lagrange([0 1 2], [-1 -1 7]);
//disp(pol_lag)
