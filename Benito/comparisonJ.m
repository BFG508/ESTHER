%% Sección de configuración inicial
close all; clear; clc;  % Limpia el entorno de trabajo

% Constantes fundamentales
mu = 3.986004418e14;  % Parámetro gravitacional terrestre [m³/s²]
J2 = 1.082635854e-3;  % Segundo coeficiente armónico zonal
J3 = -2.532435346e-6; % Tercer coeficiente armónico zonal
R  = 6378.137e3;      % Radio terrestre ecuatorial [m]

% Ángulos orbitales iniciales
Omega   = deg2rad(0); % Ascensión recta del nodo ascendente
omega   = deg2rad(0); % Argumento del perigeo
theta_0 = deg2rad(0); % Anomalía verdadera inicial

% Parámetros orbitales variables
a_values   = R + [400e3, 700e3, 1400e3];                 % Semiejes mayores [m]
e_values   = [0, 0.01, 0.1];                             % Excentricidades [-]
INC_values = deg2rad([0, 30, 45, asind(sqrt(4/5)), 90]); % Inclinaciones [rad]

% Preparación de combinaciones
num_a         = numel(a_values);         % Número de semiejes mayores
num_e         = numel(e_values);         % Número de excentricidades
num_inc       = numel(INC_values);       % Número de inclinaciones
combinaciones = num_a * num_e * num_inc; % Total de combinaciones

% Inicialización de estructuras de datos
t_values = cell(num_a, 1);          % Cell array para vectores temporales
posECI   = zeros(3, combinaciones); % Matriz para posiciones ECI
velECI   = zeros(3, combinaciones); % Matriz para velocidades ECI
n_values = zeros(1, num_a);         % Vector de movimientos medios
T_values = zeros(1, num_a);         % Vector de periodos orbitales

% Cálculo de parámetros temporales
for ia = 1:num_a
    n_values(ia) = sqrt(mu / a_values(ia)^3); % Movimiento medio [rad/s]
    T_values(ia) = 2*pi/n_values(ia);         % Periodo orbital [s]
    t_values{ia} = 0:0.1:2*T_values(ia)';     % Vector temporal para cada semieje mayor
end

%% Generación de condiciones iniciales
idx = 1;  % Contador global de combinaciones
for ia = 1:num_a
    a = a_values(ia);
    for ie = 1:num_e
        e = e_values(ie);
        for iinc = 1:num_inc
            INC = INC_values(iinc);
            
            % Conversión de elementos orbitales a ECI
            [pos, vel] = OE2ECI(mu, a, e, INC, Omega, omega, theta_0);
            
            % Almacenamiento de posición y velocidad
            posECI(:, idx) = pos;
            velECI(:, idx) = vel;
            idx = idx + 1;
        end
    end
end

% Extracción de componentes individuales
x0I  = posECI(1, :); % Componente x inicial
y0I  = posECI(2, :); % Componente y inicial
z0I  = posECI(3, :); % Componente z inicial
dx0I = velECI(1, :); % Velocidad x inicial
dy0I = velECI(2, :); % Velocidad y inicial
dz0I = velECI(3, :); % Velocidad z inicial

%% Simulación parcial de trayectorias analíticas
% soluciones_J2   = cell(num_a, num_e, num_inc);
% soluciones_J3   = cell(num_a, num_e, num_inc);
soluciones_J2J3 = cell(num_a, num_e, num_inc);

for ia = 1:num_a
    a = a_values(ia);
    n = n_values(ia);
    current_t = t_values{ia};
    
    for ie = 1:num_e
        e = e_values(ie);
        for iinc = 1:num_inc
            INC = INC_values(iinc);
            
            % Vector inicial nulo para solución general
            X0 = [0; 0; 0; 0; 0; 0];
            
            % Cálculo de la solución general analítica
            [xJ2J3, yJ2J3, zJ2J3, xdotJ2J3, ydotJ2J3, zdotJ2J3] = computeGeneralSolution(...
                J2, J3, R, a, e, INC, Omega, omega, theta_0, X0, n, current_t);

            % [xJ2, yJ2, zJ2, xdotJ2, ydotJ2, zdotJ2] = computeGeneralSolution(...
            %     J2, 0, R, a, e, INC, Omega, omega, theta_0, X0, n, current_t);
            % 
            % [xJ3, yJ3, zJ3, xdotJ3, ydotJ3, zdotJ3] = computeGeneralSolution(...
            %     0, J3, R, a, e, INC, Omega, omega, theta_0, X0, n, current_t);
            
            % Almacenamiento de la solución analítica
            soluciones_J2J3{ia, ie, iinc} = struct(...
                't', current_t, ...
                'a', a, 'e', e, 'inc', INC, ...
                'xJ', xJ2J3, 'yJ', yJ2J3, 'zJ', zJ2J3,...
                'xdotJ', xdotJ2J3, 'ydotJ', ydotJ2J3, 'zdotJ', zdotJ2J3);

            % soluciones_J2{ia, ie, iinc} = struct(...
            %     't', current_t, ...
            %     'a', a, 'e', e, 'inc', INC, ...
            %     'xJ', xJ2, 'yJ', yJ2, 'zJ', zJ2,...
            %     'xdotJ', xdotJ2, 'ydotJ', ydotJ2, 'zdotJ', zdotJ2);
            % 
            % soluciones_J3{ia, ie, iinc} = struct(...
            %     't', current_t, ...
            %     'a', a, 'e', e, 'inc', INC, ...
            %     'xJ', xJ3, 'yJ', yJ3, 'zJ', zJ3,...
            %     'xdotJ', xdotJ3, 'ydotJ', ydotJ3, 'zdotJ', zdotJ3);
        end
    end
end

%% Resolución de ODEs
% Configuración del solver de ODEs
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-5); % Tolerancias numéricas

% Inicialización de almacenamiento de soluciones numéricas
dXdt_kepler      = cell(1, combinaciones); % Soluciones keplerianas
dXdt_cowell_J2J3 = cell(1, combinaciones); % Soluciones con J2/J3
dXdt_hodei       = cell(1, combinaciones); % Soluciones con variables polares-nodales

% Simulación de dinámica orbital para cada combinación
idx = 1;  % Índice global para combinaciones
for ia = 1:num_a
    current_t = t_values{ia};
    
    for ie = 1:num_e
        for iinc = 1:num_inc
            % Configuración del estado inicial para ODE
            X0I = [posECI(:,idx); velECI(:,idx)];  % Vector de estado [pos; vel]
            
            % Dinámica kepleriana (sin perturbaciones)
            [~, sol_kepler] = ode113(@(t, X) odeKepler(X, mu), current_t, X0I, opts);
            dXdt_kepler{idx} = sol_kepler;
            
            % Dinámica con perturbaciones J2/J3 (Cowell)
            dXdt_cowell_J2J3{idx} = cowellZonals(...
                X0I, current_t, mu, R, [J2; J3; 0; 0; 0], opts);

            % Dinámica con variables polares-nodales
            dXdt_hodei{idx} = hodeiMotionCartesian(current_t, X0I, mu, R, 0, J2, opts);
            
            idx = idx + 1;
        end
    end
end

%% Organización final en estructura con componentes separados
resultados = cell(num_a, num_e, num_inc);
idx = 1;

for ia = 1:num_a
    for ie = 1:num_e
        for iinc = 1:num_inc

            % Extraer datos de las soluciones numéricas
            analytical_data  = soluciones_J2J3{ia, ie, iinc};
            kepler_data      = dXdt_kepler{idx};
            cowell_data      = dXdt_cowell_J2J3{idx};
            hodei_data       = dXdt_hodei{idx};
            
            % Almacenar cada componente por separado
            resultados{ia, ie, iinc} = struct(...
                't',   t_values{ia},...
                'a',   a_values(ia),...
                'e',   e_values(ie),...
                'inc', INC_values(iinc), ...
                'analytical_x',  kepler_data(:, 1) + analytical_data.xJ', ...
                'analytical_y',  kepler_data(:, 2) + analytical_data.yJ', ...
                'analytical_z',  kepler_data(:, 3) + analytical_data.zJ', ...
                'analytical_dx', kepler_data(:, 4) + analytical_data.xdotJ', ...
                'analytical_dy', kepler_data(:, 5) + analytical_data.ydotJ', ...
                'analytical_dz', kepler_data(:, 6) + analytical_data.zdotJ', ...
                'cowell_x',  cowell_data(:, 1), ...
                'cowell_y',  cowell_data(:, 2), ...
                'cowell_z',  cowell_data(:, 3), ...
                'cowell_dx', cowell_data(:, 4), ...
                'cowell_dy', cowell_data(:, 5), ...
                'cowell_dz', cowell_data(:, 6),  ...
                'hodei_x',  hodei_data(:, 1), ...
                'hodei_y',  hodei_data(:, 2), ...
                'hodei_z',  hodei_data(:, 3), ...
                'hodei_dx', hodei_data(:, 4), ...
                'hodei_dy', hodei_data(:, 5), ...
                'hodei_dz', hodei_data(:, 6) ...
            );
            
            idx = idx + 1;
        end
    end
end

%% Generación de gráficos individuales
if ~exist('comparaciones', 'dir')
    mkdir('comparaciones'); 
end
set(0, 'DefaultFigureVisible', 'off');

plasma = [0.20 0.02 0.40;
          0.32 0.02 0.60;
          0.56 0.07 0.65;
          0.82 0.29 0.48;
          0.95 0.52 0.27];

componentes = {
    'x'  '$x(t)$'          'km'
    'y'  '$y(t)$'          'km' 
    'z'  '$z(t)$'          'km'
    'dx' '$\dot{x}(t)$'    'm/s'
    'dy' '$\dot{y}(t)$'    'm/s'
    'dz' '$\dot{z}(t)$'    'm/s'
};

for ia = 1:num_a
    for ie = 1:num_e
        for k = 1:size(componentes, 1)
            var = componentes{k, 1};
            leyendas = cell(1, num_inc);
            figure; hold on;
            for iinc = 1:num_inc
                datos = resultados{ia, ie, iinc};
                error_comp = datos.(['cowell_' var]) - datos.(['analytical_' var]);
                if k < 4
                    plot(datos.t, error_comp / 1e3, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
                else
                    plot(datos.t, error_comp, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
                end
                leyendas{iinc} = sprintf('i = %.2f$^\\circ$', rad2deg(datos.inc));
            end

            % Formatear parámetros para el nombre del archivo
            a_clean = sprintf('%.0f', datos.a/1e3);
            e_str = sprintf('%.3f', datos.e);
            e_parts = strsplit(e_str, '.');
            e_clean = e_parts{2};

            % Título y etiquetas
            subtitulo = sprintf('a = %.0f km, e = %.3f', datos.a/1e3, datos.e);
            subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
            title(['Error de ' componentes{k,2}], 'Interpreter', 'latex');
            xlabel('$t$ (s)', 'Interpreter', 'latex');
            ylabel(['Error de ' componentes{k,2} ' (' componentes{k,3} ')'], 'Interpreter', 'latex');
            legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
            xlim([min(datos.t), max(datos.t)]);
            grid on; 
            box off;

            % Guardar figura
            filenamePNG = sprintf('%s_a%s_e%s.png', var, a_clean, e_clean);
            saveas(gcf, fullfile('comparaciones', filenamePNG));
            % filenameSVG = sprintf('%s_a%s_e%s.svg', var, a_clean, e_clean);
            % saveas(gcf, fullfile('comparaciones', filenameSVG), 'svg');
            close(gcf);
        end
    end
end

%% Generación de gráficos con múltiples variables
for ia = 1:num_a
    for ie = 1:num_e
        leyendas = cell(1, num_inc);

        % Error de posición
        figure; hold on;
        for iinc = 1:num_inc
            datos      = resultados{ia, ie, iinc};
            
            cowell_r     = [datos.cowell_x, datos.cowell_y, datos.cowell_z];
            analytical_r = [datos.analytical_x, datos.analytical_y, datos.analytical_z];
            error_pos    = vecnorm(cowell_r - analytical_r, 2, 2);

            plot(datos.t, error_pos/1e3, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
            leyendas{iinc} = sprintf('i = %.2f$^\\circ$', rad2deg(datos.inc));
        end
        hold off;
        
        % Título y etiquetas
        subtititulo = sprintf('a = %.0f km, e = %.3f', datos.a/1e3, datos.e);
        subtitle(subtititulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||r(t)||$', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||r(t)||$ (km)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos.t), max(datos.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos.a/1e3);
        e_str = sprintf('%.3f', datos.e);
        e_parts = strsplit(e_str, '.');
        e_clean = e_parts{2};
        filenamePNG = sprintf('r_a%s_e%s.png', a_clean, e_clean);
        saveas(gcf, fullfile('comparaciones', filenamePNG));
        % filenameSVG = sprintf('r_a%s_e%s.svg', a_clean, e_clean);
        % saveas(gcf, fullfile('comparaciones', filenameSVG), 'svg');
        close(gcf);

        % Error de velocidad
        figure; hold on;
        for iinc = 1:num_inc
            datos      = resultados{ia, ie, iinc};
            
            cowell_v     = [datos.cowell_dx, datos.cowell_dy, datos.cowell_dz];
            analytical_v = [datos.analytical_dx, datos.analytical_dy, datos.analytical_dz];
            error_vel    = vecnorm(cowell_v - analytical_v, 2, 2);

            plot(datos.t, error_vel, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
            leyendas{iinc} = sprintf('i = %.2f$^\\circ$', rad2deg(datos.inc));
        end
        hold off;
        
        % Título y etiquetas
        subtititulo = sprintf('a = %.0f km, e = %.3f', datos.a/1e3, datos.e);
        subtitle(subtititulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||v(t)||$', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||v(t)||$ (m/s)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos.t), max(datos.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos.a/1e3);
        e_str = sprintf('%.3f', datos.e);
        e_parts = strsplit(e_str, '.');
        e_clean = e_parts{2};
        filenamePNG = sprintf('v_a%s_e%s.png', a_clean, e_clean);
        saveas(gcf, fullfile('comparaciones', filenamePNG));
        % filenameSVG = sprintf('v_a%s_e%s.svg', a_clean, e_clean);
        % saveas(gcf, fullfile('comparaciones', filenameSVG), 'svg');
        close(gcf);
    end
end

%% Comparaciones con Kechichian

plasma_complement = [
    0.80 0.98 0.60;
    0.68 0.98 0.40;
    0.44 0.93 0.35;
    0.18 0.71 0.52;
    0.05 0.48 0.73
];

%% Comparaciones con Hodei

for ia = 1:num_a
    for ie = 1:num_e
        for k = 1:size(componentes, 1)
            var = componentes{k, 1};
            leyendas = cell(1, 2*num_inc);
            figure; hold on;
            for iinc = 1:num_inc
                datos = resultados{ia, ie, iinc};
                % Error analytical - cowell
                error_B = datos.(['analytical_' var]) - datos.(['cowell_' var]);
                % Error hodei - cowell
                error_H = datos.(['hodei_' var]) - datos.(['cowell_' var]);
                % Escalado para posición
                if k < 4
                    plot(datos.t, error_B / 1e3, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
                    plot(datos.t, error_H / 1e3, 'Color', plasma_complement(iinc, :), 'LineWidth', 1.5);
                else
                    plot(datos.t, error_B, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
                    plot(datos.t, error_H, 'Color', plasma_complement(iinc, :), 'LineWidth', 1.5);
                end
                leyendas{2*iinc-1} = sprintf('B %.2f$^\\circ$', rad2deg(datos.inc));
                leyendas{2*iinc}   = sprintf('H %.2f$^\\circ$', rad2deg(datos.inc));
            end

            % Formatear parámetros para el nombre del archivo
            a_clean = sprintf('%.0f', datos.a/1e3);
            e_str = sprintf('%.3f', datos.e);
            e_parts = strsplit(e_str, '.');
            e_clean = e_parts{2};

            % Título y etiquetas
            subtitulo = sprintf('a = %.0f km, e = %.3f', datos.a/1e3, datos.e);
            subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
            title(['Error de ' componentes{k,2}], 'Interpreter', 'latex');
            xlabel('$t$ (s)', 'Interpreter', 'latex');
            ylabel(['Error de ' componentes{k,2} ' (' componentes{k,3} ')'], 'Interpreter', 'latex');
            legend(leyendas, 'Interpreter', 'latex', 'Location', 'eastoutside', 'Orientation', 'vertical');
            xlim([min(datos.t), max(datos.t)]);
            grid on; box off;

            % Guardar figura
            filenamePNG = sprintf('hodei_%s_a%s_e%s.png', var, a_clean, e_clean);
            saveas(gcf, fullfile('comparaciones', filenamePNG));
            % filenameSVG = sprintf('hodei_%s_a%s_e%s.svg', var, a_clean, e_clean);
            % saveas(gcf, fullfile('comparaciones', filenameSVG), 'svg');
            close(gcf);
        end
    end
end
