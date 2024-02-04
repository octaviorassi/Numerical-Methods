//exec("C:\Users\octav\OneDrive - frsn.utn.edu.ar\LCC\METODOS NUMERICOS\u4\gauss_sistemas.sce")
//clc
function res = determinante(A)
    [n, m] = size(A)
    
    if n <> m then
        error("No es una matriz cuadrada");
        abort;
    end
    
    // Eliminacion Gaussiana
    for i = 1 : n - 1               // i es la fila en la que estoy parado 
        for e = i + 1 : n           // e es la fila que voy a eliminar en
            mult = - A(e, i) / A(i, i)
            for j = i + 1 : n       // j es la columna en 'e' que estoy corrigiendo
               A(e, j) =  A(e, j) + mult * A(i, j)
            end
            
            // No hace falta para encontrar la solucion.
            // Ponemos en 0 los lugares que eliminamos, no hacemos
            // la cuenta para no perder precision.
            // a(e, i) = 0
        end
    end
    
    res = prod(diag(A))
    
endfunction

//A = [1 2 3; 3 -2 1; 4 2 -1]
//res = determinante(A);
//disp(res)
