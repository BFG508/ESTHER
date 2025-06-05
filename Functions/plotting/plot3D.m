function plot3D(dataCell, viewArray, titleStr, labelCell, lineColor)
%% plot3D Representa una trayectoria 3D con formato personalizable y etiquetas en LaTeX.
%
% Esta función genera una gráfica tridimensional de una trayectoria espacial o trayectoria 
% general en el espacio 3D, a partir de los vectores {x, y, z}. Permite personalizar el color 
% de la curva, el ángulo de vista, el título de la figura y las etiquetas de los ejes.
%
% Sintaxis:
%   plot3D(dataCell)
%   plot3D(dataCell, viewArray)
%   plot3D(dataCell, viewArray, titleStr)
%   plot3D(dataCell, viewArray, titleStr, labelCell)
%   plot3D(dataCell, viewArray, titleStr, labelCell, lineColor)
%
% Entradas:
%   dataCell   : (Obligatorio) Celda de tamaño 1x3 que contiene los vectores {x, y, z}
%                 - Cada uno debe ser un vector numérico del mismo tamaño.
%
%   viewArray  : (Opcional) Vector [azimuth, elevation] que define la vista 3D de la cámara.
%                 - Por defecto: 3 (equivalente a view(3), vista isométrica).
%
%   titleStr   : (Opcional) Título de la figura (en formato LaTeX).
%                 - Por defecto: '' (sin título).
%
%   labelCell  : (Opcional) Celda con las etiquetas en LaTeX para los ejes {xlabel, ylabel, zlabel}.
%                 - Por defecto: {'', '', ''}.
%
%   lineColor  : (Opcional) Color o estilo de línea para la curva (compatible con plot3).
%                 - Por defecto: 'b' (azul).
%
% Ejemplo de uso:
%   plot3D({x, y, z}, [-152 27], 'Órbita 3D', ...
%          {'$x$ (km)', '$y$ (km)', '$z$ (km)'}, 'r--');

    % Asignar valores por defecto si no se proporcionan
    if nargin < 2 || isempty(viewArray), viewArray = 3; end
    if nargin < 3 || isempty(titleStr), titleStr = ''; end
    if nargin < 4 || isempty(labelCell), labelCell = {'', '', ''}; end
    if nargin < 5 || isempty(lineColor), lineColor = 'b'; end

    % Extraer componentes de posición desde la celda
    x = dataCell{1};
    y = dataCell{2};
    z = dataCell{3};

    % Extraer etiquetas para los ejes
    xlabelStr = labelCell{1};
    ylabelStr = labelCell{2};
    zlabelStr = labelCell{3};

    % Crear figura nueva
    figure;
    plot3(x, y, z, 'Color', lineColor, 'LineWidth', 1.5);

    % Ajustar ángulo de vista 3D
    view(viewArray);

    % Añadir título con formato LaTeX
    title(titleStr, 'Interpreter', 'latex');

    % Etiquetas de ejes en LaTeX
    xlabel(xlabelStr, 'Interpreter', 'latex');
    ylabel(ylabelStr, 'Interpreter', 'latex');
    zlabel(zlabelStr, 'Interpreter', 'latex');

    % Activar rejilla
    grid on;
end
