function plotNumericErrorNorm(func, funcNumeric, range, title_text, xlabel_text, ylabel_text, color)
%% plotNumericErrorNorm Representa la norma del error entre dos conjuntos de datos vectoriales o escalares.
%
% Esta función calcula y representa gráficamente la norma del error entre dos conjuntos
% de datos (exactos y aproximados), que pueden ser escalares, vectores o matrices Nx3/3xN.
%
% Sintaxis:
%   plotNumericErrorNorm(func, funcNumeric, range)
%   plotNumericErrorNorm(..., title_text)
%   plotNumericErrorNorm(..., title_text, xlabel_text, ylabel_text)
%   plotNumericErrorNorm(..., title_text, xlabel_text, ylabel_text, color)
%
% Entradas obligatorias:
%   func         : Datos exactos (vector, matriz Nx3 o 3xN)
%   funcNumeric  : Datos aproximados (mismo tamaño que func)
%   range        : Vector de puntos donde se evalúa el error
%
% Entradas opcionales:
%   title_text   : Título del gráfico (LaTeX) [por defecto: '']
%   xlabel_text  : Etiqueta del eje X (LaTeX) [por defecto: '']
%   ylabel_text  : Etiqueta del eje Y (LaTeX) [por defecto: '']
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

    % Calcular norma del error según orientación
    if size(funcError, 2) == 3       % Caso Nx3 (vectores fila)
        normError = vecnorm(funcError, 2, 2);
    elseif size(funcError, 1) == 3   % Caso 3xN (vectores columna)
        normError = vecnorm(funcError, 2, 1);
    else                             % Caso escalar o vector simple
        normError = abs(funcError);
    end

    % Crear figura y graficar
    figure;
    plot(range, normError, 'Color', color, 'LineWidth', 1.5);
    grid on;
    box off;

    % Añadir etiquetas y título con formato LaTeX
    title(title_text, 'Interpreter', 'latex');
    xlabel(xlabel_text, 'Interpreter', 'latex');
    ylabel(ylabel_text, 'Interpreter', 'latex');

    % Ajustar límites del eje X
    xlim([min(range), max(range)]);
end
