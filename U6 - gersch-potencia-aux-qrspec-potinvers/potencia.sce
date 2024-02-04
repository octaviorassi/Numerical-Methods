// Realiza siempre el maximo de iteraciones posibles.
// Devuelve el maximo autovalor de una matriz. Es decir, su radio
// espectral.
clc
function [lambd, z] = metodo_potencia(A, z0, maxiter)
    n = size(A)(1)
    
    // Repetimos la iteracion tantas veces como sea posible
    // (n --> inf)
    for i = 1 : maxiter
        w = A * z0
        
        z = w / norm(w, 'inf')
        
        z0 = z
    end
    
    // Ahora tenemos w^(n) y z^(n), queremos hacer una iteracion mas
    // para w y llegar a w^(n + 1)
    w = A * z0 // == A * z

    // Buscamos una componente no nula de w. 
    // En particular, buscamos la maxima (abs) para disminuir el error.
    k = 1
    
    for i = 2 : n
        if (abs(w(i)) > abs(w(k)))
            k = 1 // Actualizo el k
        end
    end
    
    // Debe existir k tal que w(k) no nulo pues autovector no nulo
    
    // lambda^(n + 1) = w^(n + 1)_k / z^(n)_k
    lambd = w(k) / z0(k)
endfunction

// ESTO NO FUNCIONA!!
function [autovals, autovects] = metodo_deflacion(A, z0, maxiter)
    n = size(A)(1);
    autovals = zeros(n, 1);
    autovects = zeros(n, n);

    for i = 1 : n
        lambd, z = metodo_potencia(A, z0, maxiter);
        autovals(i) = lambd;
        autovects(:,i) = z;
        A = A - lambd * (z * z') * (1 / norm(z, 'inf') ** 2);
    end
endfunction


//C = [6 4 4 1; 4 6 1 4; 4 1 6 4; 1 4 4 6]
//B = [12 1 3 4; 1 -3 1 5; 3 1 6 -2; 4 5 -2 -1]
//autovals, autovects = metodo_deflacion(B, [1 0 1 0]', 50)
//disp(autovals)
//disp(autovects)
//disp(spec(B))

A = [6 4 4 1; 4 6 1 4; 4 1 6 4; 1 4 4 6]
[lambd, z] = metodo_potencia(A, [1 1 1 1]', 100)

printf("Valor de lambd: %3.3f\n", lambd)
disp(z)
disp(lambd * z)
disp(A * z)
disp(A * z - lambd * z)

