exec("C:\Users\octav\OneDrive - frsn.utn.edu.ar\LCC\METODOS NUMERICOS\u3\metodos_u3.sce")

// EJ7

// 0 = 1 + x^2 − y^2 + e^x cos y
// 0 = 2xy + e^x sin y
// [x0, y0] = [-1, 4]
function res=sistema_ej7(X)
    res = [1 + X(1) ^ 2 - X(2) ^ 2 + (%e ^ X(1)) * cos(X(2));
           2 * X(1) * X(2) + (%e ^ X(1)) * sin(X(2))]
endfunction

//newtonSistema(sistema_ej7, [-1, 4]', %eps, 10)


// ej8
//0 = x2 + xy3 − 9
//0 = 3x2y − 4 − y3
// a) (1.2, 2.5) b) (−2, 2.5) c) (−1.2, −2.5) d) (2, −2.5).

function res=sistema_ej8(X)
    res = [X(1) ^ 2 + X(1) * (X(2) ^ 3) - 9;
           3 * (X(1) ^ 2) * X(2) - 4 - (X(2) ^ 3)]
endfunction

// en la primer eval. estoy mostrando los resultados y que al evaluarlo
// da casi 0 
//disp(sistema_ej8(newtonSistema(sistema_ej8, [-1.2, 2.5]', %eps, 50)'))
//newtonSistema(sistema_ej8, [-2, 2.5]', %eps, 50)
//newtonSistema(sistema_ej8, [-1.2, -2.5]', %eps, 50)
//newtonSistema(sistema_ej8, [2, -2.5]', %eps, 50)



// ej9
// modelamos el sistema como
function res=sistema_ej9(K)
    res = [K(1) * exp(K(2)) + K(3) - 10;
           K(1) * exp(2 * K(2)) + 2 * K(3) - 12;
           K(1) * exp(3 * K(2)) + 3 * K(3) - 15]
endfunction


K = newtonSistema(sistema_ej9, [2, 2, 2]', %eps, 500)
disp(sistema_ej9(K))
// a)) k_1 = 8.7712864
//     k_2 = 0.2696954
//     k_3 = -1.372283
//  (para d = 1ft = 12in)


// convertimos la ecuacion en un problema de punto fijo
//       _____________________________________________
// __  /  _____________________500___________________ = r
//   \/       pi * (k(1)* exp(k(2) * r) + k(3) * r)

func_pf = "sqrt(500 / (%pi * (K(1) * exp(K(2) * x) + K(3) * x)))"
r_calc = puntofijo(func_pf, 0.2, 500, %eps);
disp(r_calc)

// luego debe tener radio al menos r_calc.

