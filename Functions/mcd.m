function resultado = mcd(numeros)
%MCD Calcula el máximo común divisor de un conjunto de enteros
%
%   resultado = MCD(numeros)
%
%   Esta función devuelve el máximo común divisor (MCD) de un vector de
%   números enteros positivos, usando la función incorporada gcd de MATLAB.
%
%   Entradas:
%       - numeros: vector fila o columna de enteros positivos
%
%   Salidas:
%       - resultado: máximo común divisor de todos los elementos del vector
%

    % Comprobamos si el vector está vacío
    if isempty(numeros)
        error('El vector de entrada está vacío.');
    end

    % Inicializamos el resultado con el primer número
    resultado = numeros(1);
    
    % Iteramos para calcular el MCD acumulado con cada elemento
    for k = 2:length(numeros)
        resultado = gcd(resultado, numeros(k)); % gcd es la función de MATLAB para dos números
    end
end
