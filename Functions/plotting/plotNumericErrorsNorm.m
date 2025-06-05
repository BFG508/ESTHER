function plotNumericErrorsNorm(func1, func2, funcNumeric, range, legend1, legend2, title_text, xlabel_text, ylabel_text, color1, color2)
%% plotNumericErrorsNorm Grafica la norma del error entre dos funciones y una solución numérica de referencia.
%
% Esta función representa la norma euclidiana de los errores cometidos por dos funciones
% respecto a un mismo conjunto de datos numéricos de referencia.
%
% Sintaxis:
%   plotNumericErrorsNorm(func1, func2, funcNumeric, range)
%   plotNumericErrorsNorm(..., legend1, legend2)
%   plotNumericErrorsNorm(..., title_text, xlabel_text, ylabel_text)
%   plotNumericErrorsNorm(..., color1, color2)
%
% Entradas obligatorias:
%   func1        : Primer conjunto de valores exactos (vector, Nx3 o 3xN)
%   func2        : Segundo conjunto de valores exactos
%   funcNumeric  : Solución numérica de referencia (mismo tamaño que func1 y func2)
%   range        : Vector de puntos donde se evalúan los errores
%
% Entradas opcionales:
%   legend1      : Etiqueta para el error de func1 [por defecto: '']
%   legend2      : Etiqueta para el error de func2 [por defecto: '']
%   title_text   : Título del gráfico (LaTeX) [por defecto: '']
%   xlabel_text  : Etiqueta del eje X (LaTeX) [por defecto: '']
%   ylabel_text  : Etiqueta del eje Y (LaTeX) [por defecto: '']
%   color1       : Color de la curva de func1 [por defecto: 'b']
%   color2       : Color de la curva de func2 [por defecto: 'r']

    % Valores por defecto
    if nargin < 5 || isempty(legend1), legend1 = ''; end
    if nargin < 6 || isempty(legend2), legend2 = ''; end
    if nargin < 7 || isempty(title_text), title_text = ''; end
    if nargin < 8 || isempty(xlabel_text), xlabel_text = ''; end
    if nargin < 9 || isempty(ylabel_text), ylabel_text = ''; end
    if nargin < 10 || isempty(color1), color1 = 'b'; end
    if nargin < 11 || isempty(color2), color2 = 'r'; end

    % Calcular errores
    funcError1 = funcNumeric - func1;
    funcError2 = funcNumeric - func2;

    % Norma del error para func1
    if size(funcError1, 2) == 3
        normError1 = vecnorm(funcError1, 2, 2);
    elseif size(funcError1, 1) == 3
        normError1 = vecnorm(funcError1, 2, 1);
    else
        normError1 = abs(funcError1);
    end

    % Norma del error para func2
    if size(funcError2, 2) == 3
        normError2 = vecnorm(funcError2, 2, 2);
    elseif size(funcError2, 1) == 3
        normError2 = vecnorm(funcError2, 2, 1);
    else
        normError2 = abs(funcError2);
    end

    % Crear figura y graficar
    figure;
    hold on;
    plot(range, normError1, 'Color', color1, 'LineWidth', 1.5);
    plot(range, normError2, 'Color', color2, 'LineWidth', 1.5);
    hold off;

    % Estética
    grid on;
    box off;

    % Leyenda y etiquetas
    legendHandle = legend({legend1, legend2}, 'Interpreter', 'latex', 'Location', 'northoutside');
    set(legendHandle, 'Orientation', 'horizontal');

    title(title_text, 'Interpreter', 'latex');
    xlabel(xlabel_text, 'Interpreter', 'latex');
    ylabel(ylabel_text, 'Interpreter', 'latex');

    xlim([min(range), max(range)]);
end
