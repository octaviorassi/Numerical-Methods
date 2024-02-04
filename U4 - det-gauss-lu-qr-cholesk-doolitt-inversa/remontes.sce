// Resolver sistema triangular superior/inferior

// A es n x n
// b es n x 1
// NO ACEPTA CEROS EN LA DIAGONAL

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

//A = [1 2 3; 0 4 5; 0 0 6]
//b = [7 8 9]'
//D = [1 0 0; 2 3 0; 4 5 6]
//
//
//y = remontesup(A, b);
//disp(y)
//
//y2 = remonteinf(D, b)
//disp(y)
