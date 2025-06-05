function funcEval = evaluateFunction(func, var, range)
%% evaluateFunction Evalúa una función simbólica sobre un conjunto de valores numéricos.
%
% Esta función toma una expresión simbólica definida en una variable simbólica
% y devuelve su evaluación numérica en cada punto de un vector dado.
%
% Sintaxis:
%   funcEval = evaluateFunction(func, var, range)
%
% Entradas:
%   func  : Expresión simbólica a evaluar (por ejemplo, sin(x), x^2 + cos(x), etc.)
%   var   : Variable simbólica en la que está definida la función (por ejemplo, x)
%   range : Vector de valores numéricos (1D) en los que se desea evaluar la función
%
% Salida:
%   funcEval : Vector de valores numéricos resultantes de evaluar func(var) en cada punto de range
%              Tiene el mismo tamaño y forma que range
%
% Nota:
%   - Es útil para graficar o comparar funciones simbólicas con datos numéricos.
%   - Internamente usa subs() para sustituir la variable y double() para obtener el valor numérico.
%
% Ejemplo:
%   syms x
%   f = x^2 * sin(x);
%   t = linspace(0, pi, 100);
%   y = evaluateFunction(f, x, t);

    % Inicializar vector de salida con ceros del mismo tamaño que el rango
    funcEval = zeros(size(range));

    % Evaluar la función simbólica en cada valor del rango
    for k = 1:length(range)
        funcEval(k) = double(subs(func, var, range(k)));
    end
end
