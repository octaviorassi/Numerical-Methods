// Devuelve las raices del polinomio de chebyshev de orden n
clc
function r = chebyshev(n)
    coef = %pi / (2 * n);
    for k = 0 : n - 1
        k_odd = 2 * k + 1;
        
        r(k + 1) = cos(coef * k_odd);
    end
endfunction

function r_ab = chebyshev_ab(n, a, b)
    r = chebyshev(n);
    b_mas_a = b + a;
    b_menos_a = b - a;
    
    for k = 1 : n
        r_ab(k) = (b_mas_a + r(k) * b_menos_a) / 2
    end
endfunction

//disp(chebyshev(4))


