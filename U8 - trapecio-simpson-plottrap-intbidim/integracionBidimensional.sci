// * fun esta dada como string
// * c y d como funciones
// a y b son constantes, d y c son funciones de una variable.
// n y m cantidad de particiones por intervalo.
function integral = integralBidimensionalTrapecio (fun, a, b, c, d, n ,m)
    h = (b - a) / n;

    deff("y = f(x,y)", "y = " + fun);
    deff("y = c(x)", "y =" + c);
    deff("y = d(x)", "y =" + d);  
   
    integral = 0
    xj = a;
  
    // Primer iteracion con el 1/2
    integral = integral  + 1/2 * trapecio(construirFuncion(fun, xj), c(xj), d(xj), m); 
    xj = xj + h;

    for j = 1 : n - 1
        integral = integral  + trapecio(construirFuncion(fun, xj), c(xj), d(xj), m); 
        xj = xj + h;
    end

    // Ultima iteracion con el 1/2
    integral = integral  + 1/2 * trapecio(construirFuncion(fun, b), c(b), d(b), m); 

    // Todo por h
    integral = h * integral;
endfunction


// * fun es la funcion en x e y dada como string
// * a y b son constantes, d y c son funciones de una variable.
// * c y d es una funcion en una variable dada como string tambien
// * n y m cantidad de particiones por intervalo.
function integral = integralBidimensionalSimpson (fun, a, b, c, d, n, m)
    h = (b - a) / n;
  
    // Calculo los limites de los intervalos
    for i = 1 : n + 1
        x(i) = a + h * (i - 1);
    end 
  
    deff("y = f(x,y)", "y = " + fun);  
    deff("y = c(x)", "y =" + c);
    deff("y = d(x)", "y =" + d);

    integral = simpson(construirFuncion(fun, a), c(a), d(a), m)
  
    for j = 2 : n
        mod = modulo(j, 2)
        integral = integral  + (4 - mod * 2) * simpson(construirFuncion(fun, x(j)), c(x(j)), d(x(j)), m);
    end

    integral = (h / 3) * (integral + simpson(construirFuncion(fun, b), c(b), d(b), m));
endfunction

// construirFuncion recibe la funcion f(x, y)  como string
// y la transforma en otra funcion h(x) dada como string tambien
// tal que
//                    h(x) === f(xj, x)
// o sea, le cambia el nombre de la variable 'y' a 'x', y evalua
// a la variable 'x' original en 'xj'.
function fxj = construirFuncion(fun, xj)
  fxj = strsubst(fun, 'x', string(xj)); // Evaluo x en xj
  fxj = strsubst(fxj, 'y', "x");        // Intercambio y con x
endfunction

















// ************** **************** **************** **************
// ************** METODOS DEL TRAPECIO Y DE SIMPSON **************
// ************** **************** **************** **************


// Metodo del trapecio para n subintervalos
// * f es un string que representa el cuerpo de la funcion
// * a, b son los extremos del intervalo
// * n es la cantidad de subintervalos
function res = trapecio(f, a, b, n)
    // Definimos la funcion
    deff("y = f(x)", "y = " + f)
    
    // Creamos los limites de los subintervalos. Si queremos
    // n subintervalos, tendremos n + 1 puntos.
    h = (b - a) / n
    for i = 1 : n + 1
        x(i) = a + h * (i - 1);
    end    
    
    // Aplicamos la regla del trapecio
    res = 1/2 * f(x(1));    
    
    for i = 2 : n
        res = res + f(x(i));
    end    
    
    res = res + 1/2 * f(x(n + 1));
    res = res * h;
endfunction

function res = simpson(f, a, b, n)
    // Definimos la funcion
    deff("y = f(x)", "y = " + f)
    
    // Creamos los limites de los subintervalos. Si queremos
    // n subintervalos, tendremos n + 1 puntos.
    h = (b - a) / n
    for i = 1 : n + 1
        x(i) = a + h * (i - 1);
    end 
    
    res = f(x(1));
    
    for i = 2 : n
        mod = modulo(i, 2)
        res = res + (4 - mod * 2) * f(x(i))
    end
    
    res = (h / 3) * (res + f(x(n + 1)))
end

// ************** **************** **************** *************
// Ejemplo ejercicio 5 de la practica

f = "sin(x + y)"
a = 0;
b = 2;
c = "0";
d = "1";
n = 50;

// 1.63
I = integralBidimensionalSimpson(f, a, b, c, d, n, n)
disp(I)


