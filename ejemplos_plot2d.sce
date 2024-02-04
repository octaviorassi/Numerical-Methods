exec("C:\Users\octav\OneDrive - frsn.utn.edu.ar\LCC\METODOS NUMERICOS\u3\metodos_u3.sce");


x = -5:0.001:5;

// funciones a - b - c - d - e
y1 = cos(x) .* cosh(x) + 1;
y2 = 2 .* sin(x) - x .^ 2;
y3 = (%e ^ (-x)) - x .^ 4;
y4 = log(x) - x + 1;
y5 = ((x .^ 2) / 4) - sin(x);
y6 = 2 ^ (x - 1);
y7 = x;

y8 = ((%e ^ x) - x) / 2;
y9 = (%e ^ x) / 3
// raices aproximadas graficamente
// y1) -1.85 y 1.85
// y2) 0 y 1,4
// y3) -1.4 y 0.8
// y4) -0.28 y 1 (descarto -0.28 por log de algo negativo)
// y5) 0 y 1.95

for i = 1:6
    subplot(2,3,i);
    a = gca();
    a.x_location = "origin";
    a.y_location = "origin";    
end

subplot(2,3,1);
a = gca();
//plot(x, y8);
plot(x, y7);
plot(x, y9);

subplot(2,3,2);
a = gca();
plot(x, y2);

subplot(2,3,3);
a = gca();
plot(x, y3);

subplot(2,3,4);
a = gca();
plot(x, y4);

subplot(2,3,5);
a = gca();
plot(x, y5);

subplot(2,3,6);
a = gca();
plot(x, y6);
plot(x, y7);


// 2) De ser posible, hallar todas las ra´ıces de las funciones del ´ıtem anterior utilizando el m´etodos de
//la bisecci´on y la secante, ambos con una precisi´on de 10−6
// Comparar la cantidad de iteraciones obtenidas en cada metodo.
// function res=bisect(fun, a, b, tol, maxiter)
// function res=secante(fun, x0, x1, tol)

disp("f1, raiz 1")
disp(secante("cos(x) .* cosh(x) + 1", -2, -1.7, 10 ^ (-6)));
disp(bisect("cos(x) .* cosh(x) + 1", -2, -1, 10 ^ (-6), 1000));

// f1 raiz 2
disp("f1 raiz 2")
disp(secante("cos(x) .* cosh(x) + 1", 2, 1.7, 10 ^ (-6)));
disp(bisect("cos(x) .* cosh(x) + 1", 1, 3, 10 ^ (-6), 1000));

// f2 raiz 1
disp("f2 raiz 1")
disp(secante("2 .* sin(x) - x .^ 2", -0.5, 0.5, 10 ^ (-6)));
disp(bisect("2 .* sin(x) - x .^ 2", -0.2, 0.2, 10 ^ (-6), 1000));

// f2 raiz 2
disp("f2 raiz 2")
disp(secante("2 .* sin(x) - x .^ 2", 1.2, 1.6, 10 ^ (-6)));
disp(bisect("2 .* sin(x) - x .^ 2", -0.5, 1.5, 10 ^ (-6), 1000));

// f3 raiz 1
disp("f3 raiz 1")
disp(secante("(%e ^ (-x)) - x .^ 4", -1.6, -1.2, 10 ^ (-6)));

// f3 raiz 2
disp("f3 raiz 2")
disp(secante("(%e ^ (-x)) - x .^ 4", 0.7, 0.9, 10 ^ (-6)));

// f4 raiz 2 (la raiz negativa no tiene sentido creo)
disp("f4 raiz 2")
disp(secante("log(x) - x + 1", 0.5, 1, 10 ^ (-6)));

// f5 raiz 1
disp("f5 raiz 1")
disp(bisect("((x .^ 2) / 4) - sin(x)", -1, 1, 10 ^ (-6), 100));

// aca usando los mismos a y b con la sec me aproxima la otra raiz
disp(secante("((x .^ 2) / 4) - sin(x)", -1, 1, 10 ^ (-6))); 

// aca eligiendo a y b para que la secante pase mas o menos por el 0
// hago que me encuentre la raiz en 0
disp(secante("((x .^ 2) / 4) - sin(x)", -2, 0.5, 10 ^ (-6))); 

// f5 raiz 2
disp("f5 raiz 2")
disp(bisect("((x .^ 2) / 4) - sin(x)", 1, 3, 10 ^ (-6), 100));
disp(secante("((x .^ 2) / 4) - sin(x)", 1, 3, 10 ^ (-6))); 
 





