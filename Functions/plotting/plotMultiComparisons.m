function plotMultiComparisons(yFunctions, xFunction, titleText, labels, legendsCell, colorCell)
%% plotMultiComparisons Representa múltiples funciones respecto a una variable común.
%
% Esta función genera una gráfica 2D donde se representan varias curvas
% en función de una variable independiente común (por ejemplo, tiempo).
%
% Sintaxis:
%   plotMultiComparisons(yFunctions, xFunction)
%   plotMultiComparisons(..., titleText)
%   plotMultiComparisons(..., titleText, labels)
%   plotMultiComparisons(..., titleText, labels, legendsCell)
%   plotMultiComparisons(..., titleText, labels, legendsCell, colorCell)
%
% Entradas obligatorias:
%   yFunctions   : Celda con vectores {y1, y2, ..., yn} para cada función dependiente
%   xFunction    : Vector con la variable independiente común (por ejemplo, tiempo)
%
% Entradas opcionales:
%   titleText    : Título del gráfico (LaTeX) [por defecto: '']
%   labels       : Celda con etiquetas {ylabel, xlabel} (LaTeX) [por defecto: {'', ''}]
%   legendsCell  : Celda con etiquetas de leyenda [por defecto: {'f_1', ..., 'f_n'}]
%   colorCell    : Celda con colores para las curvas [por defecto: MATLAB default]

    % Número de curvas a representar
    nCurves = numel(yFunctions);

    % Valores por defecto
    if nargin < 3 || isempty(titleText)
        titleText = '';
    end
    if nargin < 4 || isempty(labels)
        labels = {'', ''};
    end
    if nargin < 5 || isempty(legendsCell)
        legendsCell = arrayfun(@(k) sprintf('$f_{%d}$', k), 1:nCurves, 'UniformOutput', false);
    end
    if nargin < 6 || isempty(colorCell)
        defaultColors = get(gca, 'ColorOrder');
        colorCell = arrayfun(@(k) defaultColors(mod(k-1, size(defaultColors, 1))+1, :), ...
                             1:nCurves, 'UniformOutput', false);
    end

    % Crear y configurar la figura
    figure;
    hold on;
    grid on;

    % Dibujar cada función
    for i = 1:nCurves
        plot(xFunction, yFunctions{i}, 'Color', colorCell{i}, 'LineWidth', 1.5);
    end

    % Añadir etiquetas y título
    title(titleText, 'Interpreter', 'latex');
    ylabel(labels{1}, 'Interpreter', 'latex');
    xlabel(labels{2}, 'Interpreter', 'latex');

    % Añadir leyenda
    legendHandle = legend(legendsCell, ...
                          'Interpreter', 'latex', ...
                          'Location', 'northoutside');
    set(legendHandle, 'Orientation', 'horizontal');

    % Ajustes estéticos
    set(gca, 'FontSize', 12);
    xlim([min(xFunction), max(xFunction)]);
end
