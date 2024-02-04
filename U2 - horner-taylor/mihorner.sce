// mihorner recibe:
// * p: un polinomio dado como vector de coeficientes. 
//     ejemplo. 1 + 2x^2 + 5x^3 ==== [1 2 5]
// * x0: el punto donde queremos evaluar al polinomio

// mihorner retorna el vector res = [a b] donde
//          a == p(x0)
//          b == p'(x0)
// e imprime en pantalla estos valores.
function res = mihorner(p, x0)    
    n = length(p) - 1; 
    
    res(n + 1) = p(n + 1);         
    res(n) = p(n) + x0 * res(n + 1)
    resderiv(n) = res(n + 1);     
    
    for i = n - 1 : -1 : 1 
        res(i) = p(i) + x0 * res(i + 1);
        resderiv(i) = res(i + 1) + x0 * resderiv(i + 1);
    end
    
    printf("p(%.f) = %.f \n", x0, res(1));
    printf("p_deriv(%.f) = %.f \n", x0, resderiv(1));      
endfunction

//mihorner([3 -5 4 2], 4);
//eval = mihorner([1 2], 3)
//disp(eval)
