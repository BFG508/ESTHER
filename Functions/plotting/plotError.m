function plotError(var, func, func_approx, range, title_text, xlabel_text, ylabel_text, color)
%% plotError Representa el error entre una función exacta y su aproximación simbólica.
%
% Esta función evalúa y representa el error entre dos funciones simbólicas sobre un rango de valores.
%
% Sintaxis:
%   plotError(var, func, func_approx, range)
%   plotError(..., title_text)
%   plotError(..., title_text, xlabel_text, ylabel_text)
%   plotError(..., title_text, xlabel_text, ylabel_text, color)
%
% Entradas obligatorias:
%   var          : Variable simbólica (por ejemplo, x)
%   func         : Función simbólica exacta
%   func_approx  : Función simbólica aproximada
%   range        : Vector de valores numéricos donde se evalúa el error
%
% Entradas opcionales:
%   title_text   : Título del gráfico (formato LaTeX) [por defecto: '']
%   xlabel_text  : Etiqueta del eje X (LaTeX) [por defecto: '']
%   ylabel_text  : Etiqueta del eje Y (LaTeX) [por defecto: '']
%   color        : Color de la curva del error [por defecto: 'b']

    % Asignar valores por defecto si no se proporcionan
    if nargin < 5 || isempty(title_text)
        title_text = '';
    end
    if nargin < 6 || isempty(xlabel_text)
        xlabel_text = '';
    end
    if nargin < 7 || isempty(ylabel_text)
        ylabel_text = '';
    end
    if nargin < 8 || isempty(color)
        color = 'b';
    end

    % Calcular y evaluar el error
    funcError = func - func_approx;
    errorVal = evaluateFunction(funcError, var, range);

    % Crear figura y graficar
    figure;
    plot(range, errorVal, 'Color', color, 'LineWidth', 1.5);
    grid on;
    box off;

    % Etiquetas y título en formato LaTeX
    title(title_text, 'Interpreter', 'latex');
    xlabel(xlabel_text, 'Interpreter', 'latex');
    ylabel(ylabel_text, 'Interpreter', 'latex');

    % Ajuste de eje X
    xlim([min(range), max(range)]);
end
