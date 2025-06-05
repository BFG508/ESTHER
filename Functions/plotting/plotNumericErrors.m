function plotNumericErrors(func1, func2, funcNumeric, range, legend1, legend2, title_text, xlabel_text, ylabel_text, color1, color2)
%% plotNumericErrors Representa dos errores numéricos con respecto a una referencia común.
%
% Esta función calcula y grafica los errores entre un conjunto de datos numéricos
% y dos referencias diferentes, en el mismo gráfico.
%
% Sintaxis:
%   plotNumericErrors(func1, func2, funcNumeric, range)
%   plotNumericErrors(..., legend1, legend2)
%   plotNumericErrors(..., title_text, xlabel_text, ylabel_text)
%   plotNumericErrors(..., color1, color2)
%
% Entradas obligatorias:
%   func1        : Primer conjunto de valores exactos (vector)
%   func2        : Segundo conjunto de valores exactos (vector)
%   funcNumeric  : Valores numéricos comunes (vector)
%   range        : Vector con los puntos de evaluación
%
% Entradas opcionales:
%   legend1      : Leyenda para func1 [por defecto: '']
%   legend2      : Leyenda para func2 [por defecto: '']
%   title_text   : Título del gráfico [por defecto: '']
%   xlabel_text  : Etiqueta del eje X [por defecto: '']
%   ylabel_text  : Etiqueta del eje Y [por defecto: '']
%   color1       : Color de la primera curva [por defecto: 'b']
%   color2       : Color de la segunda curva [por defecto: 'r']

    % Valores por defecto
    if nargin < 5 || isempty(legend1), legend1 = ''; end
    if nargin < 6 || isempty(legend2), legend2 = ''; end
    if nargin < 7 || isempty(title_text), title_text = ''; end
    if nargin < 8 || isempty(xlabel_text), xlabel_text = ''; end
    if nargin < 9 || isempty(ylabel_text), ylabel_text = ''; end
    if nargin < 10 || isempty(color1), color1 = 'b'; end
    if nargin < 11 || isempty(color2), color2 = 'r'; end

    % Calcular errores punto a punto
    funcError1 = funcNumeric - func1;
    funcError2 = funcNumeric - func2;

    % Crear figura
    figure;
    hold on;
    plot(range, funcError1, 'Color', color1, 'LineWidth', 1.5);
    plot(range, funcError2, 'Color', color2, 'LineWidth', 1.5);
    hold off;

    % Ajustes visuales
    grid on;
    box off;

    % Leyenda
    legendHandle = legend({legend1, legend2}, ...
                          'Interpreter', 'latex', ...
                          'Location', 'northoutside');
    set(legendHandle, 'Orientation', 'horizontal');

    % Etiquetas y título
    title(title_text, 'Interpreter', 'latex');
    xlabel(xlabel_text, 'Interpreter', 'latex');
    ylabel(ylabel_text, 'Interpreter', 'latex');

    % Ajustar eje X
    xlim([min(range), max(range)]);
end
