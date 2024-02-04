// MinimosQR para phis cualquiera
// * x e y son preimagenes e imagenes
// * funcs es un array de funciones phi [phi_1, phi_2, ..., phi_n]
// que formaran las n columnas de la matriz A.
// * Devuelve un string funresult que es el string representante de 
// la funcion de minimos cuadrados de n parametros o coeficientes.
clc
function [funresult, err, coefs] = minimosQR_gen(x, y, funcs)
    m = length(x)
    n = size(funcs)(2)

    A = []
    
    // Generamos la matriz A de funciones phi
    for i = 1 : n
        // Creo la funcion f_i
        deff("y = f(x)", "y =" + funcs(i))
        
        // Agrego la columna a la matriz A
        A = [A map(x, f)']
    end
    
    // Descomponemos A = QR
    coefs = factQR_sol(A, y')
    disp(coefs)
    
    // Armamos el string de funcion resultante
    funresult = ""
    
    for i = 1 : n - 1
        funresult = funresult + string(coefs(i)) + "*" + funcs(i) + " + "
    end
    
    funresult = funresult + string(coefs(n)) + "*" + funcs(n)
    
    // Devolvemos el error
    err = norm(A * coefs - y', 2);
endfunction


// Funcion auxiliar map. Toma un vector y una funcion (no string)
// y devuelve el vector resultante de aplicarle la funcion.
//   (la utilidad es que f no necesita ser escrita con ./ .*)
function y = map(x, f)    
    // La defino como igual a x asi preserva el tipo
    m = length(x)
    
    // Hago esto para preservar el tipo de x
    y = x    
    for i = 1 : m
        y(i) = f(x(i))
    end
endfunction


// *************** *************** *************** **************
// *************** Importamos la factorizacion QR ***************
// *************** *************** *************** **************
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

// *************** *************** *************** ***************
// *************** ******   TESTEO EN EJ7    ***** ***************
// *************** *************** *************** ***************

//preimg = [0 0.15 0.31 0.5 0.6 0.75]
//img = [1 1.004 1.31 1.117 1.223 1.422]

//preimg = [4 4.2 4.5 4.7 5.1 5.5 5.9 6.3 6.8 7.1]
//img = [102.56 113.18 130.11 142.05 167.53 195.14 224.87 256.73 299.5 326.72]
//
//funstring = minimosQR_gen(preimg, img, ["1", "cos(x)", "sin(x)", "cos(2*x)", "sin(2*x)"])
//deff("y = f(x)", "y = " + funstring)
//
//a = gca()
//a.x_location = "origin"
//a.y_location = "origin"
//
////x_rng = 0 : 0.01 : 1
//x_rng = 3 : 0.01 : 8
//y_rng = f(x_rng)
//
//plot(x_rng, y_rng)
//scatter(preimg, img)



