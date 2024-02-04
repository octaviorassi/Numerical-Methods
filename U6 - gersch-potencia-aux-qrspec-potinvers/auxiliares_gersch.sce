function r = radio_fila(A, i)    
   r = 0
           
   for j = 1 :i - 1
       r = r + abs(A(i, j))
   end    
    
   for j = i + 1 : n
       r = r + abs(A(i, j))
   end
endfunction

function radios = radios_gersch(A)
    n = size(A)(1)
     
    for i = 1 : n
        radios(i) = radio_fila(A, i)
    end    
endfunction

function plot_autovalores(A)
    av = spec(A)
    n = length(av)
    y(1 : n) = 0
    
    plot(av, 0, 'X')
endfunction

function plot_diagonales(A)
    centros = diag(A)
    
    n = length(centros)
    
    y(1 : n) = 0
    
    plot(centros, y, 'o')    
endfunction

function p = pol_caracteristico(A)
   lambdaI = eye(A)
    
end
    
    
endfunction
