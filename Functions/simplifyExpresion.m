function simplifiedExpresion = simplifyExpresion(expresion, numLoops, numSteps)
%% simplifyExpresion Simplificación adaptativa de expresiones simbólicas.
%
% Esta función aplica un proceso iterativo de simplificación simbólica con ajustes
% dinámicos en el número de pasos de la función simplify(), intentando mejorar
% la forma cerrada o legibilidad de una expresión matemática.
%
% Sintaxis:
%   simplifiedExpresion = simplifyExpresion(expresion, numLoops, numSteps)
%
% Entradas:
%   expresion : Expresión simbólica que se desea simplificar.
%   numLoops  : Número máximo de iteraciones (intentos) de simplificación.
%   numSteps  : Número inicial de pasos que se pasarán a simplify() en la primera iteración.
%               Este número se reducirá progresivamente en cada ciclo.
%
% Salida:
%   simplifiedExpresion : Expresión simbólica resultante tras aplicar simplificaciones adaptativas.
%
% Estrategia:
%   - Se utiliza simplify() con la opción 'preferReal' y sin restricciones analíticas.
%   - En cada iteración, el número de pasos se reduce (división exponencial).
%   - El proceso se detiene si la expresión no cambia respecto a la iteración anterior.
%
% Ejemplo:
%   syms x
%   expr = sin(x)^2 + cos(x)^2;
%   simpler = simplifyExpresion(expr, 4, 100);

    % Inicializar la expresión simplificada con la original
    simplifiedExpresion = expresion;

    % Iterar hasta numLoops veces con pasos decrecientes
    for j = 1:numLoops
        % Reducir progresivamente el número de pasos de simplify, mínimo 5
        stepsUsed = max(5, round(numSteps / (2^(j - 1))));

        % Intentar simplificación preferentemente real e ignorando restricciones analíticas
        newExpresion = simplify(simplifiedExpresion, ...
                                'Criterion', 'preferReal', ...
                                'Steps', stepsUsed, ...
                                'IgnoreAnalyticConstraints', true);

        % Si no hay cambio, finalizar el proceso anticipadamente
        if isequal(newExpresion, simplifiedExpresion)
            break;
        end

        % Actualizar la expresión para la siguiente iteración
        simplifiedExpresion = newExpresion;
    end
end
