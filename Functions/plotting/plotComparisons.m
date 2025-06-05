function plotComparisons(dataCell, labelCell, titleStr, legendCell, colorCell)
%% plotComparisons Genera una figura comparativa 2D de dos curvas con formato personalizable en LaTeX.
%
% Esta función permite comparar gráficamente dos conjuntos de datos en función
% de un eje común, con opciones de personalización para etiquetas, colores, título y leyenda.
%
% Sintaxis:
%   plotComparisons(dataCell)
%   plotComparisons(dataCell, labelCell)
%   plotComparisons(dataCell, labelCell, titleStr)
%   plotComparisons(dataCell, labelCell, titleStr, legendCell)
%   plotComparisons(dataCell, labelCell, titleStr, legendCell, colorCell)
%
% Entradas:
%   dataCell     : (Obligatoria) Celda de 1x3 con vectores numéricos {x, y1, y2}
%                   - x  : Vector del eje X
%                   - y1 : Primer conjunto de datos
%                   - y2 : Segundo conjunto de datos
%
%   labelCell    : (Opcional) Celda con etiquetas en LaTeX para los ejes {xlabel, ylabel1, ylabel2}
%                   - Por defecto: {'', '', ''}
%
%   titleStr     : (Opcional) Cadena para el título de la figura
%                   - Por defecto: ''
%
%   legendCell   : (Opcional) Celda con dos cadenas en LaTeX para la leyenda de las curvas
%                   - Por defecto: {'', ''}
%
%   colorCell    : (Opcional) Celda con los estilos de línea o color para cada curva
%                   - Por defecto: {'b', 'r'}
%
% Ejemplo de uso:
%   plotComparisons({t, x1, x2}, ...
%                   {'$t$', '$x_1$', '$x_2$'}, ...
%                   'Comparación de señales', ...
%                   {'$x_1(t)$', '$x_2(t)$'}, ...
%                   {'k--', 'r-'});

    % Valores por defecto para entradas opcionales
    if nargin < 2 || isempty(labelCell)
        labelCell = {'', '', ''};
    end
    if nargin < 3 || isempty(titleStr)
        titleStr = '';
    end
    if nargin < 4 || isempty(legendCell)
        legendCell = {'', ''};
    end
    if nargin < 5 || isempty(colorCell)
        colorCell = {'b', 'r'};
    end

    % Extraer datos de entrada
    x  = dataCell{1};
    y1 = dataCell{2};
    y2 = dataCell{3};

    % Etiquetas de ejes
    xlabelStr   = labelCell{1};
    ylabelArray = {labelCell{2}, labelCell{3}};

    % Crear figura y trazar curvas
    figure;
    hold on;
    plot(x, y1, 'Color', colorCell{1}, 'LineWidth', 1.5);
    plot(x, y2, 'Color', colorCell{2}, 'LineWidth', 1.5);
    hold off;

    % Etiquetar eje X
    xlabel(xlabelStr, 'Interpreter', 'latex');

    % Etiqueta del eje Y combinada como texto rotado
    ylabel(''); % se elimina etiqueta automática
    text(-0.13, 0.5, ylabelArray, ...
         'Units', 'normalized', ...
         'Rotation', 90, ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...
         'Interpreter', 'latex', ...
         'FontSize', 12);

    % Título de la figura
    title(titleStr, 'Interpreter', 'latex');

    % Añadir leyenda
    legendHandle = legend(legendCell, ...
                          'Interpreter', 'latex', ...
                          'Location', 'northoutside');
    set(legendHandle, 'Orientation', 'horizontal');

    % Activar rejilla
    grid on;
end
