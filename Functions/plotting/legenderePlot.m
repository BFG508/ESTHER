function legenderePlot(N, funcHandle, domain)
%% legenderePlot Representa los polinomios de Legendre P₀ a P_N aplicados a una función.
%
% Sintaxis:
%   legenderePlot(N, funcHandle)
%   legenderePlot(N, funcHandle, domain)
%
% Entradas:
%   N          : Grado máximo del polinomio de Legendre a representar
%   funcHandle : Función anónima, por ejemplo @(x) sin(x)
%   domain     : (Opcional) Dominio de evaluación [xmin, xmax] (por defecto [-1, 1])
%
% Ejemplo:
%   legenderePlot(4, @(x) sin(x))
%   legenderePlot(7, @(i) sin(i), [-pi/2, pi/2])
%   legenderePlot(12, @(x) cos(x), [0, pi])

    % --- Dominio por defecto ---
    if nargin < 3 || isempty(domain)
        domain = linspace(-1, 1, 1000);
    end

    % --- Detección automática del nombre de variable ---
    funcStrFull = func2str(funcHandle);                    
    varDetected = regexp(funcStrFull, '@\((\w)\)', 'tokens', 'once');
    varName = varDetected{1};                              
    funcLabel = ['$' varName '$'];                         

    % --- Preparación del string para etiquetas ---
    funcStr = regexprep(funcStrFull, '@\(.*?\)', '');      
    funcStr = strrep(funcStr, '^', '^{');
    funcStr = strrep(funcStr, '.*', '\cdot');
    funcStr = regexprep(funcStr, '([0-9\.]+)\{', '$1}{');
    funcLaTeX = ['$' funcStr '$'];                         

    % --- Evaluación de la función base ---
    x_vals = linspace(min(domain), max(domain), 1000);
    fx_vals = funcHandle(x_vals);

    % --- Cálculo de los polinomios de Legendre ---
    P_vals = zeros(N + 1, numel(x_vals));
    syms t
    for n = 0:N
        Pn = (1 / (2^n * factorial(n))) * diff((t^2 - 1)^n, t, n);
        P_vals(n + 1, :) = double(subs(Pn, t, fx_vals));
    end

    % --- Gráfico ---
    cmap = lines(N + 1);
    figure; hold on;
    for n = 0:N
        plot(x_vals, P_vals(n + 1, :), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('$P_{%d}(%s)$', n, funcStr), ...
             'Color', cmap(n + 1, :));
    end

    xlabel(funcLabel, 'Interpreter', 'latex');
    ylabel(sprintf('$P_n(%s)$', funcStr), 'Interpreter', 'latex');
    title(['Polinomios de Legendre evaluados sobre ', funcLaTeX], 'Interpreter', 'latex');

    % --- Leyenda con 5 elementos por fila ---
    legendHandle = legend('Interpreter', 'latex', 'Location', 'northoutside');
    set(legendHandle, 'Orientation', 'horizontal', ...
                      'NumColumns', 5, ...
                      'FontSize', 11);

    grid on;
    xlim([min(x_vals), max(x_vals)]);
    hold off;
end
