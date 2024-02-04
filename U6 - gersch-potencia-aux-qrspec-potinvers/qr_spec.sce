// Factorizacion QR
// Q = (q1 | q2 | ... | qn)
// R es tal que A = QR
clc

function c = col(A, k)
    c = A(:,k)
endfunction

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

function [autovals] = qr_spec(A, maxiter)
    n = size(A)(1);
    autovals = zeros(n, 1);

    for i = 1 : maxiter
        [Q, R] = factQR(A);
        A = R * Q;
    end

    // The eigenvalues are the diagonal elements of the final matrix
    autovals = diag(A);
endfunction

B = [12 1 3 4; 1 -3 1 5; 3 1 6 -2; 4 5 -2 -1]
eigvals = qr_spec(B, 150)
disp(eigvals)
disp(spec(B))


