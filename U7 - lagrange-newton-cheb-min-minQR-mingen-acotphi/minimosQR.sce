// Minimo cuadrados utilizando factorizacion QR
// n es el grado del polinomio resultante
function [poli, err] = minimosQR(x, y, n, plotear)
    // Definimos la matriz A de funciones phi evaluadas en x
    m = length(x)
    A = ones(m, 1)
    
    for j = 1 : n
        A = [A ((x') ** j)]
    end
    
    // Descomponemos A = QR
    coefs = factQR_sol(A, y')
    
    // Definimos el polinomio de minimos cuadrados
    poli = poly(coefs, "x", "coeff") 
    
    // Calculamos el error
    // << REVISAR SI ESTO ESTA BIEN >>
     err = norm(A*x - y', 2);  
    
    if plotear == 1
       a = gca();
       a.x_location = "origin";
       a.y_location = "origin";
       
       x_min = min(x)
//       x_min = x_min - abs(x_min)/2
       
       x_max = max(x)
//       x_max = x_max + abs(x_max)/2
              
       x_rang = x_min : 0.01 : x_max
       y_rang = horner(poli, x_rang)
       
       scatter(x, y)
       plot(x_rang, y_rang)        
    end  
endfunction












// ***************** ***************** *****************
// *****************   AUXILIARES QR   *****************
// ***************** ***************** *****************

// Importamos la factorizacion QR
function c = col(A, k)
    c = A(:,k)
endfunction

// Factorizacion QR
// Q = (q1 | q2 | ... | qn)
// R es tal que A = QR

// Teorema. Asi definida, tenemos que
//                  Q' = inv(Q)
// Luego
//                 Ax = QRx = b
//                     <==>
//              inv(Q)QRx = inv(Q)b
//                     <==>
//                   Rx = inv(Q)b
//                     <==>
//                   Rx = Q'b
// Donde R es triangular superior, luego podemos resolverlo por
// sustitucion regresiva.
function [Q, R] = factQR(A)
    // Necesitamos calcular la base ortonormal por GS para Q
    // luego R deriva de esta

    [m, n] = size(A)
    Q = zeros(A)
    R = zeros(n, n)
    
    // Defino el primer vector columna de la base, q1 = a1/v1
    v_1 = norm(col(A, 1));
    Q(:, 1) = col(A, 1) / v_1
    R(1, 1) = v_1;
    
    for k = 2 : n
        acum = zeros(1 : m)';
        
        for i = 1 : k - 1
            // Calculo la sumatoria
            q_i = col(Q, i)
            a_k = col(A, k)
            
            // Calculo de paso el coeficiente ik de R
            r_ik = a_k' * q_i;
            acum = acum + r_ik .* q_i;
            
            R(i, k) = r_ik;
        end    
        
        // Calculo v_k
        v_k = norm(a_k - acum)
            
        // Defino la columna k de Q
        q_k = (a_k - acum) / v_k;
        Q(:, k) = q_k;     
        
        // Y colocamos en R el valor de la diagonal
        R(k, k) = v_k;
    end
endfunction

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

function res = factQR_sol(A, b)
    [Q, R] = factQR(A);
    
    b_prima = Q' * b;
    res = remontesup(R, b_prima);
endfunction



