function circ(r, x, y)
    w = 2 * r
    h = 2 * r
    
    xarc(x - r, y + r, w, h, 0, 360 * 64)
    gce().thickness = 2
endfunction

// SI LA MATRIZ ES DEFINIDA POSITIVA ENTS LOS AUTOVALS SON POS
// Y VALE LA INVERSA, SI SON AUTOVALS POS ENTS ES DEF POSITIV.

// gersch plottea:
// * Los circulos de Gersch
// * Sus radios
// * Los centros de los circulos con 'o'
// * Los autovalores de A con 'x'
function gersch(A)
    [n, m] = size(A)
    
    centros = diag(A)
    radios = sum(abs(A), 'c') - abs(centros)    
    
    // (min_x, min_y) es la esquina inferior izquierda
    // (max_x, max_y) es la esquina superior derecha    
    min_x = min(centros - radios)
    max_x = max(centros + radios)
    min_xr = round(min_x - 1)
    max_xr = round(max_x + 1)
    
    min_y = min(-radios)
    max_y = max(radios)    
    min_yr = round(min_y - 1)
    max_yr = round(max_y + 1)
    
    rectangulo = [min_xr, min_yr, max_xr, max_yr]
    
    // Acotamos el radio espectral
    cota = max(abs(min_x), abs(max_x))    
    printf("p(A) <= %3.3f\n", cota)   
    
    a = gca()
    a.x_location = "origin"
//  a.y_location = "origin"
    plot2d(real(spec(A)), imag(spec(A)), -1, "031", "", rectangulo);
    xgrid();
    
    // Imprimo los circculos
    for i = 1 : n
        circ(radios(i), centros(i), 0)
    end
    
    // Imprimo los centros
    y(1 : n) = 0    
    plot(centros, y, 'O')
    
    // Imprimo los autovalores
    av = spec(A)    
    plot(av, 0, 'X')

    // Imprimimos los radios
    for i = 1 : n
        // Definimos los vectores para graficar los radios
        x_line = [centros(i), centros(i)];
        y_line = [0, radios(i)];
    
        plot(x_line, y_line, 'LineWidth', 2);
    
        // Imprimimos el largo de los radios al lado de cada uno
        xstring(centros(i), radios(i)/2, sprintf('Radio = %.2f', radios(i)));
    end
    
    // Imprimo la linea vertical que acota al radio espectral
    x_line = [cota, cota]
    y_line = [min_yr, max_yr]
    
    plot(x_line, y_line, 'LineWidth', 1.5)
    
    if real(av) == av then
        printf("Todos los autovalores son reales.\n")
    end

endfunction

// Muestra en pantallas las cotas para los autovalores
function cotas_av(A)
    n = size(A)(1)
    
    for i = 1 : n
        acum = 0
        for j = 1 : i - 1
            acum = acum + abs(A(i, j)) 
        end
        
        for j = i + 1 : n
            acum = acum + abs(A(i, j))
        end
        
        printf("|z - %3.3f| <= %3.3f\n", A(i,i), acum);             
    end
endfunction

A = [1 0 0; -1 0 1; -1 -1 2]
B = [1 0 0; -0.1 0 0.1; -0.1 -0.1 2]
C = [1 0 0; -0.25 0 0.25; -0.25 -0.25 2]
D = [4 -1 0; -1 4 -1; -1 -1 4]
E = [3 2 1; 2 3 0; 1 0 3]
F = [4.75 2.25 -0.25; 2.25 4.75 1.25; -0.25 1.25 4.75]

M = [1 2; 1 -1]

gersch(M)
gersch(M')
