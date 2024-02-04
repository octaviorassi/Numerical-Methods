// Plot trapecios.
function plot_trapecios(f, a, b, n)
    deff("y = f(x)", "y = " + f)
    
    // Creamos los limites de los subintervalos. Si queremos
    // n subintervalos, tendremos n + 1 puntos.
    h = (b - a) / n
    for i = 1 : n + 1
        x(i) = a + h * (i - 1);
    end 
    
    y = f(x);
    
    // Ploteamos
    axes = gca();
    axes.x_location = "origin";
    axes.y_location = "origin"; 
        
    // Esta va a ser la verdadera funcion,
    // la ploteamos en azul.
    x_preciso = a : 0.01 : b;
    y_preciso = f(x_preciso)    
    plot(x_preciso, y_preciso, 'linewidth', 2, 'color', 'blue')
    
    // Ploteamos tambien los trapecios
    // en naranja
    plot(x, y, 'linewidth', 2, 'color', 'orange')
    for i = 1 : n + 1
    // Imprimo la linea vertical de (x(i), 0) a (x(i), y(i))
        plot([x(i), x(i)], [0, y(i)], 'r--')
    end
    scatter(x, y)
endfunction

//plot_trapecios("tan(x)", 0, %pi/2 - 0.1, 3);
