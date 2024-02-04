function [L, U, P] = factLU(A)
    [n, m] = size(A)
    
    if n <> m then
        error("A no es cuadrada. Error.");
        abort
    end
      
    
    L = eye(A)
    U = A
    P = eye(A)    
    for i = 1 : n - 1
         // Pivoteo.
         // Busqueda de la fila_piv con magnitud maxima en el pivot.
        fila_piv = i;
        for k = i + 1 : n
           if abs(U(k, i)) > abs(U(fila_piv, i))
              fila_piv = k; 
           end
        end
        
        // Swap en U.
        tmp = U(i, i : n);
        U(i, i : n) = U(fila_piv, i : n);
        U(fila_piv, i : n) = tmp;
        
        // Swap en L.
        tmp = L(i, 1 : i - 1);
        L(i, 1 : i - 1) = L(fila_piv, 1 : i - 1);
        L(fila_piv, 1 : i - 1) = tmp;
        
        // Swap en P.
        tmp = P(i, 1 : n);
        P(i, 1 : n) = P (fila_piv, 1 : n);
        P (fila_piv, 1 : n) = tmp;
        
        // Eliminacion Gaussiana
        for j = i + 1 : n
            L(j, i) = U(j, i) / U(i, i);
            U(j, i : n) = U(j, i : n) - L(j, i) .* U(i, i : n)  
        end       
    end
endfunction

//A = [2 1 1 0; 4 3 3 1; 8 7 9 5; 6 7 9 8]
//[L, U, P] = factLU(A);
//
//disp(L);
//disp(U);
//disp(P)
//disp(L * U)
//disp(P * A)

