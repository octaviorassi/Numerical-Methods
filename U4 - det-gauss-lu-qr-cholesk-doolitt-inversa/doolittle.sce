// Auxiliar. Funciones de remonte para triangulares.
clc
function b = remonteinf(A, b)
    [n, m] = size(A)
    b(1) = b(1) / A(1, 1);
    
    for i = 2 : n
        suma = b(i)
        
        for j = 1 : i - 1
            suma = suma - A(i, j) * b(j)
        end
        
        b(i) = suma / A(i, i)         
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

// OBSERVACION. La matriz 'A' debe ser invertible.
function [L, U] = doolittle(A)
    [n, m] = size(A)
    
    if n <> m then
        error("doolittle - A no es cuadrada\n")
        abort
    end
    
    L = eye(A)
    U = zeros(A)
    U(1, 1 : n) = A(1, 1 : n)
    
    for i = 2 : n        
        // Primero calculamos la proxima columna de L
        //  Dentro de este bucle, iteramos sobre las j-esimas filas
        // de la columna FIJA i - 1 de L (la notacion se da vuelta)
        col = i - 1;
        for j = i : n 
            acum = 0;
            for k = 1 : col - 1
               acum = acum + L(j, k) * U(k, col) 
            end
            
            L(j, col) = (A(j, col) - acum) / U(col, col);
        end
        
        // Ahora podemos calcular la i-esima fila de U
        for j = i : n
            acum = 0;
            for k = 1 : i - 1
               acum = acum + L(i, k) * U(k, j)
            end
            
            U(i, j) = A(i, j) - acum;
        end
    end    
endfunction

// Tenemos que (LU)x = b
// Luego su llamamos Ux = y
// Tenemos que Ux es un vector de n x 1. Resolvemos el sistema
//                          L * y = b
// donde y es un vector cualquiera.

// Esto nos da que valores tiene que verificar y. Como Ux es un vector
// tal que L * (Ux) = b, en particular Ux debe ser igual a ese y
// Es decir, nos resta calcular para que vector x se verifica que
//                          U x = y
// Y ese vector es la solucion de LUx = b

// Para calcular cada sistema, aprovechamos que son triangulares y
// utilizamos los remontes para inferiores y superiores.
function res = resuelve_sist_lu(L, U, b)
    y = remonteinf(L, b)
    res = remontesup(U, y)
endfunction

A = [1 2 3 4; 1 4 9 16; 1 8 27 64; 1 16 81 256] 
C = [5 6 6 8; 2 2 2 8; 6 6 2 8; 2 3 6 7]       
b = [2 10 44 190]'

[L, U] = doolittle(A)
//disp(L)
//disp(U)
//disp(L * U)
//disp(A)
x = resuelve_sist_lu(L, U, b)

printf("La solucion del sistema es \n")
disp(x')




