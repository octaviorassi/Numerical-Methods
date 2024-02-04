// Las funciones son phi1 = 1, phi2 = x, phi3 = x^3..
// X vector de preimagenes
// Y vector de imagenes
// n grado del polinomio por el que aproximaremos a minimos cuadrados
function [p, err] = minimos_cuadrados_poli(X, Y, n)
    // m es la cantidad de nodos
    m = length(X)
    
    // Creo una matriz A que es una columna de unos.
    // Obs. que A es la matriz de las funciones phi_i. En nuestra
    // aproximacion polinomial, definimos a phi_i como
    //                phi_{i}(x) = x ** (i - 1)
    // para i en {1, 2, ..., n + 1}
    // (como n es el grado del polinomio, necesitamos n + 1 funciones
    // pues la phi_1 se corresponde al termino independiente.
    A = ones(m, 1)
    
    // Ahora agrego n columnas a la derecha de A que resultan de 
    // evaluar a X en phi_j
    for j = 1 : n
        A = [A (X')**(j)];
    end    
    
    // Ahora por el teorema, los coeficientes del polinomio de
    // minimos cuadrados son la solucion a este sistema.    
    X = gauss_slup(A' * A, A' * Y');
    
    // Armamos el polinomio
    p = poly(X, "x", "coeff");
    
    // Calculamos el error
    err = norm(A*X - Y', 2);                
endfunction















// ***************** ****************** *****************
// ***************** AUXILIARES - GAUSS *****************
// ***************** ****************** *****************
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
            tmp = a(i, i : n + m); // De la diagonal hacia la der
            a(i, i : n + m) = a(fila_piv, i : n + m);
            a(fila_piv, i : n + m) = tmp;
            
            // Swap en L.
            tmp = L(i, 1 : i - 1); // De la diagonal hacia la izq
            L(i, 1 : i - 1) = L(fila_piv, 1 : i - 1);
            L(fila_piv, 1 : i - 1) = tmp;
            
            // Swap en P.
            tmp = P(i, 1 : n);     // Todo
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

