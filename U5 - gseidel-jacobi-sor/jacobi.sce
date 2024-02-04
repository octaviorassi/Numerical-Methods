clc
function y = sumatoria_jacobi(A, b, x0, i)
    b_term = b(i);
    coef = 1/A(i, i);
        
    // Tambien puede separarse en dos sumatorias hasta i y desp
    // de i
    sumatoria = 0;
    for j = 1 : n
        if i <> j 
            sumatoria = sumatoria + A(i, j) * x0(j);
        end
    end              
                        
    y = coef * (b_term - sumatoria);  
endfunction

function [x1, cont] = jacobi(A, b, x0, tol, iter_max)
    n = size(A)(1);
    x1 = zeros(x0)    
    
    // Calculo la primer iteracion asi no falla el norm()
    for i = 1 : n
        x1(i) = sumatoria_jacobi(A, b, x0, i)                  
    end
    
    cont = 1;
       
    while norm(x1 - x0) > tol && cont < iter_max
       x0 = x1;
       
       // Piso todos los x^(cont) con x^(cont + 1)
       for i = 1 : n          
            x1(i) = sumatoria_jacobi(A, b, x0, i)                        
       end       
            
        cont = cont + 1;
    end        
endfunction

// Convergencia de Jacobi.
// Aplicando el corolario, convergencia_jacobi devuelve True si
// converge para cualquier x0 inicial, o False si existe un x0
// tal que no converge
// (EL COROLARIO ES UN SI Y SOLO SI)
function y = convergencia_jacobi(A)
    // M == I - inv(N) * A
    // Como N = diag(A), inv(N) son los reciprocos de la diagonal
    
    // En el teorema 3 puede verse que la matriz M es efectivamente
    // I - inv(N) * A como la hemos construido.
    n = size(A)(1)
    M = zeros(A)
    for i = 1 : n
       for j = 1 : n
          if i <> j
             M(i, j) = - A(i, j) / A(i, i) 
          end 
       end 
    end
    
    // Obs podria sacar spec(M) con factorizacion QR
    p = max(abs(spec(M)))  
    
    y = p < 1
    
    if y then
        printf("Converge para cualquier x0 inicial.\n")
    else
        printf("Convergencia no asegurada.\n")
    end
endfunction

// Aplica la condicion suficiente, que es mas debil que la
// condicion necesaria y suficiente
// Es el teorema 3
// Devuelve 1 si el teorema asegura la convergencia para todo x0
// o devuelve 0 si no puede asegurarlo.
function esDominante = convergencia_jacobi_condsuf(A)
    esDominante = 1
    i = 1
    n = size(A)(1)
    
    while esDominante && i <= n
        sumatoria = 0
        j = 1
        
        for j = 1 : i - 1
            sumatoria = sumatoria + abs(A(i, j))
        end
        for j = i + 1 : n
            sumatoria = sumatoria + abs(A(i, j))
        end
        
        esDominante = abs(A(i, i)) > sumatoria        
        i = i + 1
    end   
    
    if esDominante then
       printf("Convergencia asegurarada por condicion suficiente\n") 
    else
       i = i - 1
       printf("No puede asegurarse la convergencia. Falla la fila %i\n", i - 1)
       disp(A(i, :)) 
       printf("Sumatoria: %3.3f\nValor de la diagonal: %3.3f\n", sumatoria, abs(A(i, i)))
    end
endfunction

//
//
//
A = [1 -1 0; -1 2 -1; 0 -1 1.1]
//b = [0 1 0]
//
//[x1, iter] = jacobi(A, b', [51 0 -33]', 0.01, 1000)
//
//disp(b)
//disp(A * x1)
//disp(iter)
//
//[x2, iter] = jacobi(A, b', x1, 0.01, 100)
//
//disp(b)
//disp(A * x2)
//disp(iter)
//
//convergencia_jacobi(A)
disp(convergencia_jacobi_condsuf(A))
