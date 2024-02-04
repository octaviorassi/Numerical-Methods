// IMPORTANTE. Todos los metodos asumen que las funciones son pasadas
// como strings. Por ejemplo,
// >> f = "3 * x ^ 2"
// >> newton(f, x0, tol, maxiter)

// DOCUMENTACION deff. deff simplemente se reemplaza de la sig manera:

// deff(functionHeadline, functionBody)
//
//                     <===> 
//
// function functionHeadline
//      functionBody
// endfunction

function y = puntofijo(g, x0, maxiter, tol)
    deff("y=g(x)","y="+g)   // La iteracion de punto fijo
    
    cont = 0;        
    while abs(g(x0) - x0) > tol && cont < maxiter
        x0 = g(x0)
        cont = cont + 1;
    end
    
    if cont >= maxiter then
        printf("puntofijo finalizo por alcanzar el maximo de iteraciones\n");
    else
        printf("puntofijo finalizo por alcanzar la tolerancia requerida\n");
    end
    
    y = x0;   
endfunction

function y = puntofijo_alt(g, x0, maxiter, tol)
   deff("y=g(x)","y="+g)   // La iteracion de punto fijo
   x1 = g(x0)
   
   cont = 1;
   while abs(g(x1) - x1) > tol && cont < maxiter
        x0 = x1;
        x1 = g(x1);
        cont = cont + 1;
    end
    
    if cont >= maxiter then
        printf("puntofijo finalizo por alcanzar el maximo de iteraciones");
    else
        printf("puntofijo finalizo por alcanzar la tolerancia requerida");
    end
endfunction


// si z es tal que f(z) = 0 y x0 esta 'cerca' de z y verifica que
// f'(x0) != 0 ents 
function salida=newton(fun, x0, tol, maxiter)
    deff("y=f(x)","y="+fun);

    i = 0;
    x1 = x0 - f(x0)/numderivative(f, x0)
    
    while abs(x1-x0) > tol && i < maxiter
        i = i + 1;
        x0 = x1;
        x1 = x0 - f(x0)/numderivative(f, x0)
    end;
    
    if abs(x1 - x0) > tol then
        disp("Tolerancia");
    else 
        disp("Iteraciones");
    end    
    salida = x1;
endfunction

function salida=newton_pf(fun, x0, tol, maxiter)
    deff("y=f(x)","y="+fun);
    deff("y=g(x)0","y=" + "x -" + "f(x)/numderivative(f, x)")
    
    salida = puntofijo(g, x0, tol, maxiter)    
endfunction

// REPRESENTA EL SISTEMA
// x^2 - 2y
function res=sistema(X)
    res = [X(1)^2 - 2 * X(2) - 2; X(1) * X(2) - 2]
endfunction

function salida=newtonSistema(f, x0, tol, maxiter)   
    i = 0;
    
    J = numderivative(f, x0);
    J = inv(J);
    xn = x0 - J* f(x0);
    
    while norm(xn-x0) > tol && i < maxiter
        
        i = i + 1;
        x0 = xn;
        
        J = numderivative(f, x0);
        J = inv(J);
        xn = x0 - J* f(x0);
                
    end;
    
    if norm(xn - x0) > tol then
        disp("Tolerancia");
    else 
        disp("Iteraciones");
        disp("Cantidad de operaciones: ", i);
    end    
    
    disp(xn);
    salida = xn;
endfunction


function res=bisect(fun, a, b, tol, maxiter)
    deff("y=f(x)","y="+fun);  
      
    c = (a + b) / 2;
    i = 0;
    
    // es mejor comparar f(c) ?
    while abs(f(c)) > tol && i < maxiter
//  while abs(b - c) > tol && i < maxiter
        i = i + 1;
        
        if (f(b) * f(c) < 0) then
            a = c;
        else
            b = c;
        end;
        
        c = (a + b) / 2;        
    end;
    
    res = c;
endfunction

function res=secante(fun, x0, x1, tol)
    deff("y=f(x)","y="+fun);   
    
    x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
        
    while abs(f(x2)) > tol
        x0 = x1; // x
        x1 = x2;     
           
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));                
    end
    
    res = x2;        
endfunction

function res=regulafalsi(fun, a, b, tol, maxit)
    deff("y=f(x)","y="+fun);
    //    Metodo de la secante
    c = b - f(b) * (b - a) / (f(b) - f(a));    
    i = 0;
    while abs(f(c)) > tol && i < maxit
        i = i + 1;
        // Metodo de la biseccion
        if f(a) * f(b) < 0 then
            b = c;
        else 
            a = c;
        end
        
        c = b - f(b) * (b - a) / (f(b) - f(a)); 
    end
    
    res = c;
endfunction



