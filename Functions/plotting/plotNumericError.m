function plotNumericError(func, funcNumeric, range, title_text, xlabel_text, ylabel_text, color)
%% plotNumericError Representa el error entre una función exacta y su aproximación numérica.
%
% Esta función calcula el error punto a punto entre un vector de valores exactos
% y uno aproximado, y lo representa gráficamente sobre un rango dado.
%
% Sintaxis:
%   plotNumericError(func, funcNumeric, range)
%   plotNumericError(..., title_text)
%   plotNumericError(..., title_text, xlabel_text, ylabel_text)
%   plotNumericError(..., title_text, xlabel_text, ylabel_text, color)
%
% Entradas obligatorias:
%   func         : Vector de valores exactos
%   funcNumeric  : Vector de valores aproximados
%   range        : Vector de valores donde se ha evaluado la función
%
% Entradas opcionales:
%   title_text   : Título del gráfico (formato LaTeX) [por defecto: '']
%   xlabel_text  : Etiqueta del eje X (formato LaTeX) [por defecto: '']
%   ylabel_text  : Etiqueta del eje Y (formato LaTeX) [por defecto: '']
%   color        : Color o estilo de línea [por defecto: 'b']

    % Asignar valores por defecto si no se proporcionan
    if nargin < 4 || isempty(title_text)
        title_text = '';
    end
    if nargin < 5 || isempty(xlabel_text)
        xlabel_text = '';
    end
    if nargin < 6 || isempty(ylabel_text)
        ylabel_text = '';
    end
    if nargin < 7 || isempty(color)
        color = 'b';
    end

    % Calcular el error punto a punto
    funcError = funcNumeric - func;

    % Crear figura y graficar
    figure;
    plot(range, funcError, 'Color', color, 'LineWidth', 1.5);
    grid on;
    box off;

    % Etiquetas y título en LaTeX
    title(title_text, 'Interpreter', 'latex');
    xlabel(xlabel_text, 'Interpreter', 'latex');
    ylabel(ylabel_text, 'Interpreter', 'latex');

    % Ajustar límites del eje X
    xlim([min(range), max(range)]);
end
