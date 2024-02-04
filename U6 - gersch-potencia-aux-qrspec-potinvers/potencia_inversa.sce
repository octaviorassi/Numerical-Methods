// Metodo de la potencia inversa
// Necesitamos Gauss para resolver sistemas
// MEJORA ---> USAR GAUSS-SEIDEL

// Auxiliar. Funcion de remonte para triangulares superiores.
function b = remontesup(A, b)
    [n, m] = size(A)
    
    b(n) = b(n) / A(n, n)
    for i = n - 1 : -1 : 1
        suma = b(i)
        
        for j = i + 1 : n 
            suma = suma - A(i, j) * b(j)
        end
        
        b(i) = suma / A(i, i)        
    end   
endfunction

// Gauss version final
// Resuelve el sistema aplicando pivoteo parcial y calculando ademas
// las matrices L, U y P.
// Admite matrices B de n x m

// gauss_slup: A, B ---> [S, L, U, P]
// donde:
//  * A es la matriz cuadrada de coeficientes de n x n
//  * B es la matriz del lado derecho, de n x m
//  * S es una matriz de n x m donde la j-esima columna es la solucion
// al sistema resultante de considerar A como matriz de coeficientes y
// la j-esima columna de B como termino independiente
//  * L es la matriz triangular inferior
//  * U es la matriz triangular superior
//  * P es la matriz de permutaciones

//  Las matrices anteriores verifican la relacion
//                         P * A = L * U

// Observacion. Si A es una matriz invertible y B es eye(A), entonces
// S es la inversa de A.

function [S, L, U, P] = gauss_slup(A, B)
    [n, mA] = size(A) 
    [nB, m] = size(B)
    
    // Chequeos de dimensiones y compatibilidad del producto
    if n <> mA then
        error('gauss - La matriz A debe ser cuadrada');
        abort;
    elseif mA <> nB then
        error('gauss - dimensiones incompatibles entre A y B');
        abort;
    end;
    
    // a es la matriz ampliada del sistema.
    a = [A B]   
    L = eye(A)
    P = eye(A)
    
    for i = 1 : n - 1
        // Pivoteo. Busqueda de la fila con magnitud maxima en el pivot.
        fila_piv = i;
        for k = i + 1 : n
           if abs(a(k, i)) > abs(a(fila_piv, i))
              fila_piv = k; 
           end
        end
        
        // Swaps en las tres matrices: a, L, P
        if fila_piv <> i
            // Swap en a.
            tmp = a(i, i : n + m);
            a(i, i : n + m) = a(fila_piv, i : n + m);
            a(fila_piv, i : n + m) = tmp;
            
            // Swap en L.
            tmp = L(i, 1 : i - 1);
            L(i, 1 : i - 1) = L(fila_piv, 1 : i - 1);
            L(fila_piv, 1 : i - 1) = tmp;
            
            // Swap en P.
            tmp = P(i, 1 : n);
            P(i, 1 : n) = P (fila_piv, 1 : n);
            P (fila_piv, 1 : n) = tmp;        
        end
        
        // Eliminacion   
        for e = i + 1 : n       // e es la fila que quiero eliminar
            mult = a(e, i) / a(i, i)
            L(e, i) = mult
                            
            for j = i + 1 : n + m    // j es la columna en 'e'
                a(e, j) =  a(e, j) - mult * a(i, j)
            end

            a(e, i) = 0
        end
    end

    
    // Sustitucion regresiva para cada sistema
    r = 1 : n
    U = a(r, r)  

    for j = 1 : m
        b = a(r, j + n)
        S(j, r) = remontesup(U, b)
    end
    
    S = S'  
       
endfunction

// Realiza siempre el maximo de iteraciones posibles.
// Devuelve el minimo autovalor de una matriz.
clc
function [lambd, z] = potencia_inversa(A, z0, maxiter)
    // Realmente deberia usar una aproximacion a un autovalor
    // mu, y usar A - mu * I como matriz en los sistemas.
    // Creo que como esta escrito aca, uso mu = 0
    
    // Obs. Tambien puede hacerse analogo al potencia normal
    // pero en vez de resolver SEL, calcular inversas
    
    // En particular, 
    // z ^ (k + 1) =       (A - mu * I)^(-1) * z^(k)
    //                ------------------------------------
    //                  norm((A - mu * I)^(-1) * z^(k))    
    
    n = size(A)(1)
    z0 = z0'
    
    // Repetimos la iteracion tantas veces como sea posible
    // (n --> inf)
    for i = 1 : maxiter
        w = gauss_slup(A, z0)
       
        z = w / norm(w, 'inf')
        
        z0 = z
    end
    
    // Ahora tenemos w^(n) y z^(n), queremos hacer una iteracion mas
    // para w y llegar a w^(n + 1)
    w = gauss_slup(A, z0)

    // Buscamos una componente no nula de w. 
    // En particular, buscamos la maxima (abs) para disminuir el error.
    k = 1
    
    for i = 2 : n
        if (abs(w(i)) > abs(w(k)))
            k = i  // Actualizo el k
        end
    end
    
    // Debe existir k tal que w(k) no nulo pues autovector no nulo
    
    // lambda^(n + 1) = z^(n)_k / w^(n + 1)_k
    lambd = z0(k) / w(k)
endfunction

B = [12 1 3 4; 1 -3 1 5; 3 1 6 -2; 4 5 -2 -1]

[lambd, z] = potencia_inversa(B, [1 1 1 1], 200)
disp(lambd)
disp(spec(B))


