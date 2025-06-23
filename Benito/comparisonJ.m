%% Sección de configuración inicial
close all; clear; clc;  % Limpia el entorno de trabajo

% Constantes fundamentales
mu = 3.986004418e14;  % Parámetro gravitacional terrestre  [m³/s²]
J2 = 1.082635854e-3;  % Segundo coeficiente armónico zonal [-]
J3 = -2.532435346e-6; % Tercer coeficiente armónico zonal  [-]
R  = 6378.137e3;      % Radio terrestre ecuatorial         [m]

% Ángulos orbitales iniciales
Omega   = deg2rad(0); % Ascensión recta del nodo ascendente
omega   = deg2rad(0); % Argumento del perigeo
theta_0 = deg2rad(0); % Anomalía verdadera inicial

% Parámetros orbitales variables
a_values   = R + [400e3, 700e3, 1400e3];                 % Semiejes mayores [m]
e_values   = [0, 0.01, 0.1];                             % Excentricidades  [-]
INC_values = deg2rad([0, 30, 45, asind(sqrt(4/5)), 90]); % Inclinaciones    [rad]

% Preparación de combinaciones
num_a         = numel(a_values);         % Número de semiejes mayores
num_e         = numel(e_values);         % Número de excentricidades
num_inc       = numel(INC_values);       % Número de inclinaciones
combinations  = num_a * num_e * num_inc; % Total de combinaciones

% Inicialización de estructuras de datos
t_values = cell(num_a, 1);         % Cell array para vectores temporales
posECI   = zeros(3, combinations); % Matriz para posiciones ECI
velECI   = zeros(3, combinations); % Matriz para velocidades ECI
n_values = zeros(1, num_a);        % Vector de movimientos medios
T_values = zeros(1, num_a);        % Vector de periodos orbitales

% Cálculo de parámetros temporales
for ia = 1:num_a
    n_values(ia) = sqrt(mu / a_values(ia)^3); % Movimiento medio [rad/s]
    T_values(ia) = 2*pi/n_values(ia);         % Periodo orbital [s]
    t_values{ia} = 0:0.1:2*T_values(ia)';     % Vector temporal para cada semieje mayor
end

%% Representación de casos de estudios
if ~exist('figuras', 'dir')
    mkdir('figuras'); 
end

set(0, 'DefaultFigureVisible', 'on');

a_km = a_values / 1e3;
inc_deg = rad2deg(INC_values);

[a_grid, inc_grid] = meshgrid(a_km, inc_deg);

casos1 = reshape(1:num_a*num_inc, num_inc, num_a);
casos2 = reshape((num_a*num_inc+1):(2*num_a*num_inc), num_inc, num_a);
casos3 = reshape((2*num_a*num_inc+1):(3*num_a*num_inc), num_inc, num_a);

casos = {casos1, casos2, casos3};

figure
sgtitle('\bf{Matriz de casos de estudios}', 'Interpreter', 'latex')

for i = 1:3
    subplot(1,3,i)
    hold on
    
    scatter(a_grid(:), inc_grid(:), 40, 'filled', 'MarkerFaceColor', [0.56 0.07 0.65])
    
    for j = 1:numel(a_grid)
        offset_x = range(a_km) * 0.05;
        text(a_grid(j) + offset_x, inc_grid(j), num2str(casos{i}(j)), ...
             'FontSize', 14, 'Interpreter', 'latex')
    end
    
    hold off
    title(['Excentricidad ', num2str(e_values(i))], 'FontSize', 12, 'Interpreter', 'latex')
    xlabel('Semieje mayor (km)', 'FontSize', 12, 'Interpreter', 'latex')
    if i == 1
        ylabel('Inclinaci\''on ($^\circ$)', 'FontSize', 12, 'Interpreter', 'latex')
    end

    xlim([min(a_km), max(a_km)])
    ylim([min(inc_deg)-5, max(inc_deg)+5])
    
    set(gca,'XTick', a_km)
    set(gca,'YTick', inc_deg)
    
    grid on
    box off
end

set(gcf, 'WindowState', 'fullscreen');
pause(1); 
filenameSVG = 'casos_de_estudio';
print(gcf, fullfile('figuras', filenameSVG), '-dsvg', '-vector')
close(gcf);

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
soluciones_J2            = cell(num_a, num_e, num_inc);
soluciones_J3            = cell(num_a, num_e, num_inc);
soluciones_J2J3          = cell(num_a, num_e, num_inc);
soluciones_kechichian_J2 = cell(num_a, num_e, num_inc);
soluciones_kechichian_J3 = cell(num_a, num_e, num_inc);

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

            [xJ2, yJ2, zJ2, xdotJ2, ydotJ2, zdotJ2] = computeGeneralSolution(...
                J2, 0, R, a, e, INC, Omega, omega, theta_0, X0, n, current_t);

            [xJ3, yJ3, zJ3, xdotJ3, ydotJ3, zdotJ3] = computeGeneralSolution(...
                0, J3, R, a, e, INC, Omega, omega, theta_0, X0, n, current_t);

            [xKJ2, yKJ2, zKJ2, xdotKJ2, ydotKJ2, zdotKJ2] = computeGeneralSolutionKechichian(...
                J2, 0, R, a, e, INC, Omega, omega, theta_0, X0, n, current_t);

            [xKJ3, yKJ3, zKJ3, xdotKJ3, ydotKJ3, zdotKJ3] = computeGeneralSolutionKechichian(...
                0, J3, R, a, e, INC, Omega, omega, theta_0, X0, n, current_t);

            % Almacenamiento de la solución analítica
            soluciones_J2J3{ia, ie, iinc} = struct(...
                't', current_t, ...
                'a', a, 'e', e, 'inc', INC, ...
                'xJ', xJ2J3, 'yJ', yJ2J3, 'zJ', zJ2J3,...
                'xdotJ', xdotJ2J3, 'ydotJ', ydotJ2J3, 'zdotJ', zdotJ2J3);

            soluciones_J2{ia, ie, iinc} = struct(...
                't', current_t, ...
                'a', a, 'e', e, 'inc', INC, ...
                'xJ', xJ2, 'yJ', yJ2, 'zJ', zJ2,...
                'xdotJ', xdotJ2, 'ydotJ', ydotJ2, 'zdotJ', zdotJ2);

            soluciones_J3{ia, ie, iinc} = struct(...
                't', current_t, ...
                'a', a, 'e', e, 'inc', INC, ...
                'xJ', xJ3, 'yJ', yJ3, 'zJ', zJ3,...
                'xdotJ', xdotJ3, 'ydotJ', ydotJ3, 'zdotJ', zdotJ3);

            soluciones_kechichian_J2{ia, ie, iinc} = struct(...
                't', current_t, ...
                'a', a, 'e', e, 'inc', INC, ...
                'xJ', xKJ2, 'yJ', yKJ2, 'zJ', zKJ2,...
                'xdotJ', xdotKJ2, 'ydotJ', ydotKJ2, 'zdotJ', zdotKJ2);
            
            soluciones_kechichian_J3{ia, ie, iinc} = struct(...
                't', current_t, ...
                'a', a, 'e', e, 'inc', INC, ...
                'xJ', xKJ3, 'yJ', yKJ3, 'zJ', zKJ3,...
                'xdotJ', xdotKJ3, 'ydotJ', ydotKJ3, 'zdotJ', zdotKJ3);
        end
    end
end

%% Resolución de ODEs
% Configuración del solver de ODEs
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-5); % Tolerancias numéricas

% Inicialización de almacenamiento de soluciones numéricas
dXdt_kepler      = cell(1, combinations); % Soluciones keplerianas
dXdt_cowell_J2   = cell(1, combinations); % Soluciones con J2
dXdt_cowell_J3   = cell(1, combinations); % Soluciones con J3
dXdt_cowell_J2J3 = cell(1, combinations); % Soluciones con J2/J3
dXdt_hodei       = cell(1, combinations); % Soluciones con variables polares-nodales

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

            dXdt_cowell_J2{idx} = cowellZonals(...
                X0I, current_t, mu, R, [J2; 0; 0; 0; 0], opts);

            dXdt_cowell_J3{idx} = cowellZonals(...
                X0I, current_t, mu, R, [0; J3; 0; 0; 0], opts);

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
            analytical_data_J2   = soluciones_J2{ia, ie, iinc};
            analytical_data_J3   = soluciones_J3{ia, ie, iinc};
            analytical_data_J2J3 = soluciones_J2J3{ia, ie, iinc};
            kechichian_data_J2   = soluciones_kechichian_J2{ia, ie, iinc};
            kechichian_data_J3   = soluciones_kechichian_J3{ia, ie, iinc};
            kepler_data          = dXdt_kepler{idx};
            cowell_data_J2       = dXdt_cowell_J2{idx};
            cowell_data_J3       = dXdt_cowell_J3{idx};
            cowell_data_J2J3     = dXdt_cowell_J2J3{idx};
            hodei_data           = dXdt_hodei{idx};
            
            % Almacenar cada componente por separado
            resultados{ia, ie, iinc} = struct(...
                't',   t_values{ia},...
                'a',   a_values(ia),...
                'e',   e_values(ie),...
                'inc', INC_values(iinc), ...
                'analytical_x_J2',  kepler_data(:, 1) + analytical_data_J2.xJ', ...
                'analytical_y_J2',  kepler_data(:, 2) + analytical_data_J2.yJ', ...
                'analytical_z_J2',  kepler_data(:, 3) + analytical_data_J2.zJ', ...
                'analytical_dx_J2', kepler_data(:, 4) + analytical_data_J2.xdotJ', ...
                'analytical_dy_J2', kepler_data(:, 5) + analytical_data_J2.ydotJ', ...
                'analytical_dz_J2', kepler_data(:, 6) + analytical_data_J2.zdotJ', ...
                'analytical_x_J3',  kepler_data(:, 1) + analytical_data_J3.xJ', ...
                'analytical_y_J3',  kepler_data(:, 2) + analytical_data_J3.yJ', ...
                'analytical_z_J3',  kepler_data(:, 3) + analytical_data_J3.zJ', ...
                'analytical_dx_J3', kepler_data(:, 4) + analytical_data_J3.xdotJ', ...
                'analytical_dy_J3', kepler_data(:, 5) + analytical_data_J3.ydotJ', ...
                'analytical_dz_J3', kepler_data(:, 6) + analytical_data_J3.zdotJ', ...
                'analytical_x_J2J3',  kepler_data(:, 1) + analytical_data_J2J3.xJ', ...
                'analytical_y_J2J3',  kepler_data(:, 2) + analytical_data_J2J3.yJ', ...
                'analytical_z_J2J3',  kepler_data(:, 3) + analytical_data_J2J3.zJ', ...
                'analytical_dx_J2J3', kepler_data(:, 4) + analytical_data_J2J3.xdotJ', ...
                'analytical_dy_J2J3', kepler_data(:, 5) + analytical_data_J2J3.ydotJ', ...
                'analytical_dz_J2J3', kepler_data(:, 6) + analytical_data_J2J3.zdotJ', ...
                'kechichian_x_J2',  kepler_data(:, 1) + kechichian_data_J2.xJ, ...
                'kechichian_y_J2',  kepler_data(:, 2) + kechichian_data_J2.yJ, ...
                'kechichian_z_J2',  kepler_data(:, 3) + kechichian_data_J2.zJ, ...
                'kechichian_dx_J2', kepler_data(:, 4) + kechichian_data_J2.xdotJ, ...
                'kechichian_dy_J2', kepler_data(:, 5) + kechichian_data_J2.ydotJ, ...
                'kechichian_dz_J2', kepler_data(:, 6) + kechichian_data_J2.zdotJ, ...
                'kechichian_x_J3',  kepler_data(:, 1) + kechichian_data_J3.xJ, ...
                'kechichian_y_J3',  kepler_data(:, 2) + kechichian_data_J3.yJ, ...
                'kechichian_z_J3',  kepler_data(:, 3) + kechichian_data_J3.zJ, ...
                'kechichian_dx_J3', kepler_data(:, 4) + kechichian_data_J3.xdotJ, ...
                'kechichian_dy_J3', kepler_data(:, 5) + kechichian_data_J3.ydotJ, ...
                'kechichian_dz_J3', kepler_data(:, 6) + kechichian_data_J3.zdotJ, ...
                'cowell_x_J2',  cowell_data_J2(:, 1), ...
                'cowell_y_J2',  cowell_data_J2(:, 2), ...
                'cowell_z_J2',  cowell_data_J2(:, 3), ...
                'cowell_dx_J2', cowell_data_J2(:, 4), ...
                'cowell_dy_J2', cowell_data_J2(:, 5), ...
                'cowell_dz_J2', cowell_data_J2(:, 6), ...
                'cowell_x_J3',  cowell_data_J3(:, 1), ...
                'cowell_y_J3',  cowell_data_J3(:, 2), ...
                'cowell_z_J3',  cowell_data_J3(:, 3), ...
                'cowell_dx_J3', cowell_data_J3(:, 4), ...
                'cowell_dy_J3', cowell_data_J3(:, 5), ...
                'cowell_dz_J3', cowell_data_J3(:, 6),  ...
                'cowell_x_J2J3',  cowell_data_J2J3(:, 1), ...
                'cowell_y_J2J3',  cowell_data_J2J3(:, 2), ...
                'cowell_z_J2J3',  cowell_data_J2J3(:, 3), ...
                'cowell_dx_J2J3', cowell_data_J2J3(:, 4), ...
                'cowell_dy_J2J3', cowell_data_J2J3(:, 5), ...
                'cowell_dz_J2J3', cowell_data_J2J3(:, 6),  ...
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

%% Conversión a elementos orbitales
for ia = 1:num_a
    for ie = 1:num_e
        for iinc = 1:num_inc
            datos = resultados{ia, ie, iinc};
            N = length(datos.t);
            
            % Preparar matrices de posición y velocidad (Nx3)
            % Método analítico J2
            r_analytical_J2 = [datos.analytical_x_J2, datos.analytical_y_J2, datos.analytical_z_J2];
            v_analytical_J2 = [datos.analytical_dx_J2, datos.analytical_dy_J2, datos.analytical_dz_J2];
            
            % Método analítico J2+J3
            r_analytical_J2J3 = [datos.analytical_x_J2J3, datos.analytical_y_J2J3, datos.analytical_z_J2J3];
            v_analytical_J2J3 = [datos.analytical_dx_J2J3, datos.analytical_dy_J2J3, datos.analytical_dz_J2J3];
            
            % Cowell J2
            r_cowell_J2 = [datos.cowell_x_J2, datos.cowell_y_J2, datos.cowell_z_J2];
            v_cowell_J2 = [datos.cowell_dx_J2, datos.cowell_dy_J2, datos.cowell_dz_J2];
            
            % Cowell J2+J3
            r_cowell_J2J3 = [datos.cowell_x_J2J3, datos.cowell_y_J2J3, datos.cowell_z_J2J3];
            v_cowell_J2J3 = [datos.cowell_dx_J2J3, datos.cowell_dy_J2J3, datos.cowell_dz_J2J3];
            
            % Calcular elementos orbitales (llamada vectorizada)
            [a_analytical_J2, e_analytical_J2, i_analytical_J2, ...
                RAAN_analytical_J2, omega_analytical_J2, theta_analytical_J2, ...
                AOLon_analytical_J2, AOLat_analytical_J2] = ECI2OE(mu, r_analytical_J2, v_analytical_J2);
            [a_analytical_J2J3, e_analytical_J2J3, i_analytical_J2J3, ...
                RAAN_analytical_J2J3, omega_analytical_J2J3, theta_analytical_J2J3, ...
                AOLon_analytical_J2J3, AOLat_analytical_J2J3] = ECI2OE(mu, r_analytical_J2J3, v_analytical_J2J3);
            [a_cowell_J2, e_cowell_J2, i_cowell_J2, ...
                RAAN_cowell_J2, omega_cowell_J2, theta_cowell_J2, ...
                AOLon_cowell_J2, AOLat_cowell_J2] = ECI2OE(mu, r_cowell_J2, v_cowell_J2);
            [a_cowell_J2J3, e_cowell_J2J3, i_cowell_J2J3, ...
                RAAN_cowell_J2J3, omega_cowell_J2J3, theta_cowell_J2J3, ...
                AOLon_cowell_J2J3, AOLat_cowell_J2J3] = ECI2OE(mu, r_cowell_J2J3, v_cowell_J2J3);
            
            % Almacenar resultados (estructura vectorizada)
            % Analítico J2
            resultados{ia, ie, iinc}.oe_analytical_J2 = struct(...
                'a', a_analytical_J2, 'e', e_analytical_J2, 'i', i_analytical_J2, ...
                'Omega', RAAN_analytical_J2, 'omega', omega_analytical_J2, 'theta', theta_analytical_J2, ...
                'AOLon', AOLon_analytical_J2, 'AOLat', AOLat_analytical_J2);
            
            % Analítico J2+J3
            resultados{ia, ie, iinc}.oe_analytical_J2J3 = struct(...
                'a', a_analytical_J2J3, 'e', e_analytical_J2J3, 'i', i_analytical_J2J3, ...
                'Omega', RAAN_analytical_J2J3, 'omega', omega_analytical_J2J3, 'theta', theta_analytical_J2J3, ...
                'AOLon', AOLon_analytical_J2J3, 'AOLat', AOLat_analytical_J2J3);
            
            % Cowell J2
            resultados{ia, ie, iinc}.oe_cowell_J2 = struct(...
                'a', a_cowell_J2, 'e', e_cowell_J2, 'i', i_cowell_J2, ...
                'Omega', RAAN_cowell_J2, 'omega', omega_cowell_J2, 'theta', theta_cowell_J2, ...
                'AOLon', AOLon_cowell_J2, 'AOLat', AOLat_cowell_J2);
            
            % Cowell J2+J3
            resultados{ia, ie, iinc}.oe_cowell_J2J3 = struct(...
                'a', a_cowell_J2J3, 'e', e_cowell_J2J3, 'i', i_cowell_J2J3, ...
                'Omega', RAAN_cowell_J2J3, 'omega', omega_cowell_J2J3, 'theta', theta_cowell_J2J3, ...
                'AOLon', AOLon_cowell_J2J3, 'AOLat', AOLat_cowell_J2J3);
        end
    end
end

%% Generación de gráficos individuales
set(0, 'DefaultFigureVisible', 'off');

plasma = [0.20 0.02 0.40;
          0.32 0.02 0.60;
          0.56 0.07 0.65;
          0.82 0.29 0.48;
          0.95 0.52 0.27];

componentes = {
    'x'   '$x(t)$'        'km'
    'y'   '$y(t)$'        'km' 
    'z'   '$z(t)$'        'km'
    'dx'  '$\dot{x}(t)$'  'km/s'
    'dy'  '$\dot{y}(t)$'  'km/s'
    'dz'  '$\dot{z}(t)$'  'km/s'
};

if ~exist('figuras/eciJ2J3vsJ2J3', 'dir')
    mkdir('figuras', 'eciJ2J3vsJ2J3');
end

for ia = 1:num_a
    for ie = 1:num_e
        for k = 1:size(componentes, 1)
            var = componentes{k, 1};
            leyendas = cell(1, num_inc);
            figure; hold on;
            for iinc = 1:num_inc
                datos = resultados{ia, ie, iinc};
                error_comp = datos.(['cowell_' var '_J2J3']) - datos.(['analytical_' var '_J2J3']);
                plot(datos.t, error_comp / 1e3, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
                leyendas{iinc} = sprintf('$i$ = %.2f$^\\circ$', rad2deg(datos.inc));
            end

            % Formatear parámetros para el nombre del archivo
            a_clean = sprintf('%.0f', datos.a/1e3);
            e_str = sprintf('%.3f', datos.e);
            e_parts = strsplit(e_str, '.');
            e_clean = e_parts{2};

            % Título y etiquetas
            subtitulo = sprintf('$a$ = %.0f km, $e$ = %.3f, $J_3 \\neq 0$', datos.a/1e3, datos.e);
            subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
            title(['Error de ' componentes{k,2}], 'Interpreter', 'latex');
            xlabel('$t$ (s)', 'Interpreter', 'latex');
            ylabel(['Error de ' componentes{k,2} ' (' componentes{k,3} ')'], 'Interpreter', 'latex');
            legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
            xlim([min(datos.t), max(datos.t)]);
            grid on; 
            box off;

            % Guardar figura
            filenamePNG = sprintf('%s_a%s_e%s_J2J3vsJ2J3.png', var, a_clean, e_clean);
            saveas(gcf, fullfile('figuras/eciJ2J3vsJ2J3', filenamePNG));
            filenameSVG = sprintf('%s_a%s_e%s_J2J3vsJ2J3.svg', var, a_clean, e_clean);
            saveas(gcf, fullfile('figuras/eciJ2J3vsJ2J3', filenameSVG), 'svg');
            close(gcf);
        end
    end
end

if ~exist('figuras/eciJ2vsJ2J3', 'dir')
    mkdir('figuras', 'eciJ2vsJ2J3');
end

for ia = 1:num_a
    for ie = 1:num_e
        for k = 1:size(componentes, 1)
            var = componentes{k, 1};
            leyendas = cell(1, num_inc);
            figure; hold on;
            for iinc = 1:num_inc
                datos = resultados{ia, ie, iinc};
                error_comp = datos.(['cowell_' var '_J2J3']) - datos.(['analytical_' var '_J2']);
                plot(datos.t, error_comp / 1e3, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
                leyendas{iinc} = sprintf('$i$ = %.2f$^\\circ$', rad2deg(datos.inc));
            end

            % Formatear parámetros para el nombre del archivo
            a_clean = sprintf('%.0f', datos.a/1e3);
            e_str = sprintf('%.3f', datos.e);
            e_parts = strsplit(e_str, '.');
            e_clean = e_parts{2};

            % Título y etiquetas
            subtitulo = sprintf('$a$ = %.0f km, $e$ = %.3f, $J_3 \\neq 0$', datos.a/1e3, datos.e);
            subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
            title(['Error de ' componentes{k,2} ' del modelo truncado'], 'Interpreter', 'latex');
            xlabel('$t$ (s)', 'Interpreter', 'latex');
            ylabel(['Error de ' componentes{k,2} ' (' componentes{k,3} ')'], 'Interpreter', 'latex');
            legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
            xlim([min(datos.t), max(datos.t)]);
            grid on; 
            box off;

            % Guardar figura
            filenamePNG = sprintf('%s_a%s_e%s_J2vsJ2J3.png', var, a_clean, e_clean);
            saveas(gcf, fullfile('figuras/eciJ2vsJ2J3', filenamePNG));
            filenameSVG = sprintf('%s_a%s_e%s_J2vsJ2J3.svg', var, a_clean, e_clean);
            saveas(gcf, fullfile('figuras/eciJ2vsJ2J3', filenameSVG), 'svg');
            close(gcf);
        end
    end
end

if ~exist('figuras/eciJ2vsJ2', 'dir')
    mkdir('figuras', 'eciJ2vsJ2');
end

for ia = 1:num_a
    for ie = 1:num_e
        for k = 1:size(componentes, 1)
            var = componentes{k, 1};
            leyendas = cell(1, num_inc);
            figure; hold on;
            for iinc = 1:num_inc
                datos = resultados{ia, ie, iinc};
                error_comp = datos.(['cowell_' var '_J2']) - datos.(['analytical_' var '_J2']);
                plot(datos.t, error_comp / 1e3, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
                leyendas{iinc} = sprintf('$i$ = %.2f$^\\circ$', rad2deg(datos.inc));
            end

            % Formatear parámetros para el nombre del archivo
            a_clean = sprintf('%.0f', datos.a/1e3);
            e_str = sprintf('%.3f', datos.e);
            e_parts = strsplit(e_str, '.');
            e_clean = e_parts{2};

            % Título y etiquetas
            subtitulo = sprintf('$a$ = %.0f km, $e$ = %.3f, $J_3 = 0$', datos.a/1e3, datos.e);
            subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
            title(['Error de ' componentes{k,2} ' del modelo simplificado'], 'Interpreter', 'latex');
            xlabel('$t$ (s)', 'Interpreter', 'latex');
            ylabel(['Error de ' componentes{k,2} ' (' componentes{k,3} ')'], 'Interpreter', 'latex');
            legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
            xlim([min(datos.t), max(datos.t)]);
            grid on; 
            box off;

            % Guardar figura
            filenamePNG = sprintf('%s_a%s_e%s_J2vsJ2.png', var, a_clean, e_clean);
            saveas(gcf, fullfile('figuras/eciJ2vsJ2', filenamePNG));
            filenameSVG = sprintf('%s_a%s_e%s_J2vsJ2.svg', var, a_clean, e_clean);
            saveas(gcf, fullfile('figuras/eciJ2vsJ2', filenameSVG), 'svg');
            close(gcf);
        end
    end
end

%% Gráficos de diferencia en elementos orbitales (MODIFICAR ESTA PARTE PARA RAAN, AOP Y TA)

set(0, 'DefaultFigureVisible', 'off');

plasma = [0.20 0.02 0.40;
          0.32 0.02 0.60;
          0.56 0.07 0.65;
          0.82 0.29 0.48;
          0.95 0.52 0.27];

elementos_orbitales = {
    'a'      '$a(t)$'       'km'        'SMA'
    'e'      '$e(t)$'       '-'         'ECC'
    'i'      '$i(t)$'       '$^\circ$'  'INC'
    'omega'  '$\Omega(t)$'  '$^\circ$'  'RAAN'
    'omega'  '$\omega(t)$'  '$^\circ$'  'AOP'
    'theta'  '$\theta(t)$'  '$^\circ$'  'TA' 
    'AOLon'  '$\varpi(t)$'  '$^\circ$'  'AOLon'
    'AOLat'  '$u(t)$'       '$^\circ$'  'AOLat' 
};

if ~exist('figuras/oeJ2J3vsJ2J3', 'dir')
    mkdir('figuras', 'oeJ2J3vsJ2J3');
end

for ia = 1:num_a
    for ie = 1:num_e
        for k = 1:size(elementos_orbitales, 1)
            elem = elementos_orbitales{k, 1};
            elem_name = elementos_orbitales{k, 4};
            leyendas = cell(1, num_inc);
            figure; hold on;
            for iinc = 1:num_inc
                datos = resultados{ia, ie, iinc};
                if strcmp(elem, 'Omega') || strcmp(elem, 'omega') || strcmp(elem, 'theta') || strcmp(elem, 'AOLon') || strcmp(elem, 'AOLat')
                    diff_data = wrapTo2Pi(datos.oe_analytical_J2J3.(elem)) - wrapTo2Pi(datos.oe_cowell_J2J3.(elem));
                    diff_data = wrapTo360(rad2deg(diff_data));
                elseif strcmp(elem, 'a')
                    diff_data = (datos.oe_analytical_J2J3.(elem) - datos.oe_cowell_J2J3.(elem) ) / 1e3;
                else
                    diff_data = datos.oe_analytical_J2J3.(elem) - datos.oe_cowell_J2J3.(elem);
                end
                
                plot(datos.t, diff_data, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
                leyendas{iinc} = sprintf('$i$ = %.2f$^\\circ$', rad2deg(datos.inc));
            end
            hold off;
            
            % Formatear parámetros para nombre
            a_clean = sprintf('%.0f', datos.a/1e3);
            e_str = sprintf('%.3f', datos.e);
            e_parts = strsplit(e_str, '.');
            e_clean = e_parts{2};
            
            % Configuración del gráfico
            subtitulo = sprintf('$a$ = %.0f km, $e$ = %.3f, $J_3 \\neq 0$', datos.a/1e3, datos.e);
            subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
            title(['Error de ' elementos_orbitales{k,2}], 'Interpreter', 'latex');
            xlabel('$t$ (s)', 'Interpreter', 'latex');
            ylabel(['Error de ' elementos_orbitales{k,2} ' (' elementos_orbitales{k,3} ')'], 'Interpreter', 'latex');
            legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
            xlim([min(datos.t), max(datos.t)]);
            grid on; 
            box off;
            
            % Guardar figura en directorio 'oe'
            filenamePNG = sprintf('%s_a%s_e%s_J2J3vsJ2J3.png', elem_name, a_clean, e_clean);
            saveas(gcf, fullfile('figuras/oeJ2J3vsJ2J3', filenamePNG));
            filenameSVG = sprintf('%s_a%s_e%s_J2J3vsJ2J3.svg', elem_name, a_clean, e_clean);
            saveas(gcf, fullfile('figuras/oeJ2J3vsJ2J3', filenameSVG), 'svg');
            close(gcf);
        end
    end
end

if ~exist('figuras/oeJ2vsJ2J3', 'dir')
    mkdir('figuras', 'oeJ2vsJ2J3');
end

for ia = 1:num_a
    for ie = 1:num_e
        for k = 1:size(elementos_orbitales, 1)
            elem = elementos_orbitales{k, 1};
            elem_name = elementos_orbitales{k, 4};
            leyendas = cell(1, num_inc);
            figure; hold on;
            for iinc = 1:num_inc
                datos = resultados{ia, ie, iinc};
                if strcmp(elem, 'Omega') || strcmp(elem, 'omega') || strcmp(elem, 'theta') || strcmp(elem, 'AOLon') || strcmp(elem, 'AOLat')
                    diff_data = wrapTo2Pi(datos.oe_analytical_J2.(elem)) - wrapTo2Pi(datos.oe_cowell_J2J3.(elem));
                    diff_data = wrapTo360(rad2deg(diff_data));
                elseif strcmp(elem, 'a')
                    diff_data = (datos.oe_analytical_J2.(elem) - datos.oe_cowell_J2J3.(elem) ) / 1e3;
                else
                    diff_data = datos.oe_analytical_J2.(elem) - datos.oe_cowell_J2J3.(elem);
                end
                
                plot(datos.t, diff_data, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
                leyendas{iinc} = sprintf('$i$ = %.2f$^\\circ$', rad2deg(datos.inc));
            end
            hold off;
            
            % Formatear parámetros para nombre
            a_clean = sprintf('%.0f', datos.a/1e3);
            e_str = sprintf('%.3f', datos.e);
            e_parts = strsplit(e_str, '.');
            e_clean = e_parts{2};
            
            % Configuración del gráfico
            subtitulo = sprintf('$a$ = %.0f km, $e$ = %.3f, $J_3 \\neq 0$', datos.a/1e3, datos.e);
            subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
            title(['Error de ' elementos_orbitales{k,2} ' del modelo truncado'], 'Interpreter', 'latex');
            xlabel('$t$ (s)', 'Interpreter', 'latex');
            ylabel(['Error de ' elementos_orbitales{k,2} ' (' elementos_orbitales{k,3} ')'], 'Interpreter', 'latex');
            legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
            xlim([min(datos.t), max(datos.t)]);
            grid on; 
            box off;
            
            % Guardar figura en directorio 'oe'
            filenamePNG = sprintf('%s_a%s_e%s_J2vsJ2J3.png', elem_name, a_clean, e_clean);
            saveas(gcf, fullfile('figuras/oeJ2vsJ2J3', filenamePNG));
            filenameSVG = sprintf('%s_a%s_e%s_J2vsJ2J3.svg', elem_name, a_clean, e_clean);
            saveas(gcf, fullfile('figuras/oeJ2vsJ2J3', filenameSVG), 'svg');
            close(gcf);
        end
    end
end

if ~exist('figuras/oeJ2vsJ2', 'dir')
    mkdir('figuras', 'oeJ2vsJ2');
end

for ia = 1:num_a
    for ie = 1:num_e
        for k = 1:size(elementos_orbitales, 1)
            elem = elementos_orbitales{k, 1};
            elem_name = elementos_orbitales{k, 4};
            leyendas = cell(1, num_inc);
            figure; hold on;
            for iinc = 1:num_inc
                datos = resultados{ia, ie, iinc};
                if strcmp(elem, 'Omega') || strcmp(elem, 'omega') || strcmp(elem, 'theta') || strcmp(elem, 'AOLon') || strcmp(elem, 'AOLat')
                    diff_data = wrapTo2Pi(datos.oe_analytical_J2.(elem)) - wrapTo2Pi(datos.oe_cowell_J2.(elem));
                    diff_data = wrapTo360(rad2deg(diff_data));
                elseif strcmp(elem, 'a')
                    diff_data = (datos.oe_analytical_J2.(elem) - datos.oe_cowell_J2.(elem) ) / 1e3;
                else
                    diff_data = datos.oe_analytical_J2.(elem) - datos.oe_cowell_J2.(elem);
                end
                
                plot(datos.t, diff_data, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
                leyendas{iinc} = sprintf('$i$ = %.2f$^\\circ$', rad2deg(datos.inc));
            end
            hold off;
            
            % Formatear parámetros para nombre
            a_clean = sprintf('%.0f', datos.a/1e3);
            e_str = sprintf('%.3f', datos.e);
            e_parts = strsplit(e_str, '.');
            e_clean = e_parts{2};
            
            % Configuración del gráfico
            subtitulo = sprintf('$a$ = %.0f km, $e$ = %.3f, $J_3$ = 0', datos.a/1e3, datos.e);
            subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
            title(['Error de ' elementos_orbitales{k,2}], 'Interpreter', 'latex');
            xlabel('$t$ (s)', 'Interpreter', 'latex');
            ylabel(['Error de ' elementos_orbitales{k,2} ' (' elementos_orbitales{k,3} ')'], 'Interpreter', 'latex');
            legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
            xlim([min(datos.t), max(datos.t)]);
            grid on; 
            box off;
            
            % Guardar figura en directorio 'oe'
            filenamePNG = sprintf('%s_a%s_e%s_J2vsJ2.png', elem_name, a_clean, e_clean);
            saveas(gcf, fullfile('figuras/oeJ2vsJ2', filenamePNG));
            filenameSVG = sprintf('%s_a%s_e%s_J2vsJ2.svg', elem_name, a_clean, e_clean);
            saveas(gcf, fullfile('figuras/oeJ2vsJ2', filenameSVG), 'svg');
            close(gcf);
        end
    end
end

%% Generación de gráficos con múltiples variables
if ~exist('figuras/modulosJ2J3vsJ2J3', 'dir')
    mkdir('figuras', 'modulosJ2J3vsJ2J3');
end

for ia = 1:num_a
    for ie = 1:num_e
        leyendas = cell(1, num_inc);

        % Error de posición
        figure; hold on;
        for iinc = 1:num_inc
            datos      = resultados{ia, ie, iinc};
            
            cowell_r     = [datos.cowell_x_J2J3, datos.cowell_y_J2J3, datos.cowell_z_J2J3];
            analytical_r = [datos.analytical_x_J2J3, datos.analytical_y_J2J3, datos.analytical_z_J2J3];
            error_pos    = vecnorm(cowell_r - analytical_r, 2, 2);

            plot(datos.t, error_pos/1e3, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
            leyendas{iinc} = sprintf('$i$ = %.2f$^\\circ$', rad2deg(datos.inc));
        end
        hold off;
        
        % Título y etiquetas
        subtititulo = sprintf('$a$ = %.0f km, $e$ = %.3f, $J_3 \\neq 0$', datos.a/1e3, datos.e);
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
        filenamePNG = sprintf('r_a%s_e%s_J2J3vsJ2J3.png', a_clean, e_clean);
        saveas(gcf, fullfile('figuras/modulosJ2J3vsJ2J3', filenamePNG));
        filenameSVG = sprintf('r_a%s_e%s_J2J3vsJ2J3.svg', a_clean, e_clean);
        saveas(gcf, fullfile('figuras/modulosJ2J3vsJ2J3', filenameSVG), 'svg');
        close(gcf);

        % Error de velocidad
        figure; hold on;
        for iinc = 1:num_inc
            datos      = resultados{ia, ie, iinc};
            
            cowell_v     = [datos.cowell_dx_J2J3, datos.cowell_dy_J2J3, datos.cowell_dz_J2J3];
            analytical_v = [datos.analytical_dx_J2J3, datos.analytical_dy_J2J3, datos.analytical_dz_J2J3];
            error_vel    = vecnorm(cowell_v - analytical_v, 2, 2);

            plot(datos.t, error_vel/1e3, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
            leyendas{iinc} = sprintf('$i$ = %.2f$^\\circ$', rad2deg(datos.inc));
        end
        hold off;
        
        % Título y etiquetas
        subtititulo = sprintf('$a$ = %.0f km, $e$ = %.3f, $J_3 \\neq 0$', datos.a/1e3, datos.e);
        subtitle(subtititulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||v(t)||$', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||v(t)||$ (km/s)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos.t), max(datos.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos.a/1e3);
        e_str = sprintf('%.3f', datos.e);
        e_parts = strsplit(e_str, '.');
        e_clean = e_parts{2};
        filenamePNG = sprintf('v_a%s_e%s_J2J3vsJ2J3.png', a_clean, e_clean);
        saveas(gcf, fullfile('figuras/modulosJ2J3vsJ2J3', filenamePNG));
        filenameSVG = sprintf('v_a%s_e%s_J2J3vsJ2J3.svg', a_clean, e_clean);
        saveas(gcf, fullfile('figuras/modulosJ2J3vsJ2J3', filenameSVG), 'svg');
        close(gcf);
    end
end

if ~exist('figuras/modulosJ2vsJ2J3', 'dir')
    mkdir('figuras', 'modulosJ2vsJ2J3');
end

for ia = 1:num_a
    for ie = 1:num_e
        leyendas = cell(1, num_inc);

        % Error de posición
        figure; hold on;
        for iinc = 1:num_inc
            datos      = resultados{ia, ie, iinc};
            
            cowell_r     = [datos.cowell_x_J2J3, datos.cowell_y_J2J3, datos.cowell_z_J2J3];
            analytical_r = [datos.analytical_x_J2, datos.analytical_y_J2, datos.analytical_z_J2];
            error_pos    = vecnorm(cowell_r - analytical_r, 2, 2);

            plot(datos.t, error_pos/1e3, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
            leyendas{iinc} = sprintf('$i$ = %.2f$^\\circ$', rad2deg(datos.inc));
        end
        hold off;
        
        % Título y etiquetas
        subtititulo = sprintf('$a$ = %.0f km, $e$ = %.3f, $J_3 \\neq 0$', datos.a/1e3, datos.e);
        subtitle(subtititulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||r(t)||$ del modelo truncado', 'Interpreter', 'latex');
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
        filenamePNG = sprintf('r_a%s_e%s_J2vsJ2J3.png', a_clean, e_clean);
        saveas(gcf, fullfile('figuras/modulosJ2vsJ2J3', filenamePNG));
        filenameSVG = sprintf('r_a%s_e%s_J2vsJ2J3.svg', a_clean, e_clean);
        saveas(gcf, fullfile('figuras/modulosJ2vsJ2J3', filenameSVG), 'svg');
        close(gcf);

        % Error de velocidad
        figure; hold on;
        for iinc = 1:num_inc
            datos      = resultados{ia, ie, iinc};
            
            cowell_v     = [datos.cowell_dx_J2J3, datos.cowell_dy_J2J3, datos.cowell_dz_J2J3];
            analytical_v = [datos.analytical_dx_J2, datos.analytical_dy_J2, datos.analytical_dz_J2];
            error_vel    = vecnorm(cowell_v - analytical_v, 2, 2);

            plot(datos.t, error_vel/1e3, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
            leyendas{iinc} = sprintf('$i$ = %.2f$^\\circ$', rad2deg(datos.inc));
        end
        hold off;
        
        % Título y etiquetas
        subtititulo = sprintf('$a$ = %.0f km, $e$ = %.3f, $J_3 \\neq 0$', datos.a/1e3, datos.e);
        subtitle(subtititulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||v(t)||$ del modelo truncado', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||v(t)||$ (km/s)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos.t), max(datos.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos.a/1e3);
        e_str = sprintf('%.3f', datos.e);
        e_parts = strsplit(e_str, '.');
        e_clean = e_parts{2};
        filenamePNG = sprintf('v_a%s_e%s_J2vsJ2J3.png', a_clean, e_clean);
        saveas(gcf, fullfile('figuras/modulosJ2vsJ2J3', filenamePNG));
        filenameSVG = sprintf('v_a%s_e%s_J2vsJ2J3.svg', a_clean, e_clean);
        saveas(gcf, fullfile('figuras/modulosJ2vsJ2J3', filenameSVG), 'svg');
        close(gcf);
    end
end

if ~exist('figuras/modulosJ2vsJ2', 'dir')
    mkdir('figuras', 'modulosJ2vsJ2');
end

for ia = 1:num_a
    for ie = 1:num_e
        leyendas = cell(1, num_inc);

        % Error de posición
        figure; hold on;
        for iinc = 1:num_inc
            datos      = resultados{ia, ie, iinc};
            
            cowell_r     = [datos.cowell_x_J2, datos.cowell_y_J2, datos.cowell_z_J2];
            analytical_r = [datos.analytical_x_J2, datos.analytical_y_J2, datos.analytical_z_J2];
            error_pos    = vecnorm(cowell_r - analytical_r, 2, 2);

            plot(datos.t, error_pos/1e3, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
            leyendas{iinc} = sprintf('$i$ = %.2f$^\\circ$', rad2deg(datos.inc));
        end
        hold off;
        
        % Título y etiquetas
        subtititulo = sprintf('$a$ = %.0f km, $e$ = %.3f, $J_3 = 0$', datos.a/1e3, datos.e);
        subtitle(subtititulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||r(t)||$ del modelo simplificado', 'Interpreter', 'latex');
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
        filenamePNG = sprintf('r_a%s_e%s_J2vsJ2.png', a_clean, e_clean);
        saveas(gcf, fullfile('figuras/modulosJ2vsJ2', filenamePNG));
        filenameSVG = sprintf('r_a%s_e%s_J2vsJ2.svg', a_clean, e_clean);
        saveas(gcf, fullfile('figuras/modulosJ2vsJ2', filenameSVG), 'svg');
        close(gcf);

        % Error de velocidad
        figure; hold on;
        for iinc = 1:num_inc
            datos      = resultados{ia, ie, iinc};
            
            cowell_v     = [datos.cowell_dx_J2, datos.cowell_dy_J2, datos.cowell_dz_J2];
            analytical_v = [datos.analytical_dx_J2, datos.analytical_dy_J2, datos.analytical_dz_J2];
            error_vel    = vecnorm(cowell_v - analytical_v, 2, 2);

            plot(datos.t, error_vel/1e3, 'Color', plasma(iinc, :), 'LineWidth', 1.5);
            leyendas{iinc} = sprintf('$i$ = %.2f$^\\circ$', rad2deg(datos.inc));
        end
        hold off;
        
        % Título y etiquetas
        subtititulo = sprintf('$a$ = %.0f km, $e$ = %.3f, $J_3 = 0$', datos.a/1e3, datos.e);
        subtitle(subtititulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||v(t)||$ del modelo simplificado', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||v(t)||$ (km/s)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos.t), max(datos.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos.a/1e3);
        e_str = sprintf('%.3f', datos.e);
        e_parts = strsplit(e_str, '.');
        e_clean = e_parts{2};
        filenamePNG = sprintf('v_a%s_e%s_J2vsJ2.png', a_clean, e_clean);
        saveas(gcf, fullfile('figuras/modulosJ2vsJ2', filenamePNG));
        filenameSVG = sprintf('v_a%s_e%s_J2vsJ2.svg', a_clean, e_clean);
        saveas(gcf, fullfile('figuras/modulosJ2vsJ2', filenameSVG), 'svg');
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

% Obtener paleta de colores
color_analytical_min = plasma(2, :);
color_comparison_min = plasma_complement(2, :);
color_analytical_max = plasma(3, :);
color_comparison_max = plasma_complement(3, :);

if ~exist('figuras/kechichianJ2vsJ2', 'dir')
    mkdir('figuras', 'kechichianJ2vsJ2');
end

for ia = 1:num_a
    for iinc = 1:num_inc
        figure; hold on;
        
        % Primer valor de excentricidad (mínimo)
        ie_min = 1;
        datos_min = resultados{ia, ie_min, iinc};
        
        % Cálculo de errores para excentricidad mínima
        cowell_r_min     = [    datos_min.cowell_x_J2,     datos_min.cowell_y_J2,     datos_min.cowell_z_J2];
        analytical_r_min = [datos_min.analytical_x_J2, datos_min.analytical_y_J2, datos_min.analytical_z_J2];
        kechichian_r_min = [datos_min.kechichian_x_J2, datos_min.kechichian_y_J2, datos_min.kechichian_z_J2];

        cowell_v_min     = [    datos_min.cowell_dx_J2,     datos_min.cowell_dy_J2,     datos_min.cowell_dz_J2];
        analytical_v_min = [datos_min.analytical_dx_J2, datos_min.analytical_dy_J2, datos_min.analytical_dz_J2];
        kechichian_v_min = [datos_min.kechichian_dx_J2, datos_min.kechichian_dy_J2, datos_min.kechichian_dz_J2];
        
        error_r_analytical_min = vecnorm(cowell_r_min - analytical_r_min, 2, 2);
        error_r_kechichian_min = vecnorm(cowell_r_min - kechichian_r_min, 2, 2);

        error_v_analytical_min = vecnorm(cowell_v_min - analytical_v_min, 2, 2);
        error_v_kechichian_min = vecnorm(cowell_v_min - kechichian_v_min, 2, 2);
        
        % Último valor de excentricidad (máximo)
        ie_max = num_e;
        datos_max = resultados{ia, ie_max, iinc};
        
        % Cálculo de errores para excentricidad máxima
        cowell_r_max     = [    datos_max.cowell_x_J2,     datos_max.cowell_y_J2,    datos_max.cowell_z_J2];
        analytical_r_max = [datos_max.analytical_x_J2, datos_max.analytical_y_J2, datos_max.analytical_z_J2];
        kechichian_r_max = [datos_max.kechichian_x_J2, datos_max.kechichian_y_J2, datos_max.kechichian_z_J2];

        cowell_v_max     = [    datos_max.cowell_dx_J2,     datos_max.cowell_dy_J2,    datos_max.cowell_dz_J2];
        analytical_v_max = [datos_max.analytical_dx_J2, datos_max.analytical_dy_J2, datos_max.analytical_dz_J2];
        kechichian_v_max = [datos_max.kechichian_dx_J2, datos_max.kechichian_dy_J2, datos_max.kechichian_dz_J2];
        
        error_r_analytical_max = vecnorm(cowell_r_max - analytical_r_max, 2, 2);
        error_r_kechichian_max = vecnorm(cowell_r_max - kechichian_r_max, 2, 2);

        error_v_analytical_max = vecnorm(cowell_v_max - analytical_v_max, 2, 2);
        error_v_kechichian_max = vecnorm(cowell_v_max - kechichian_v_max, 2, 2);
        
        % Configurar leyenda
        leyendas = {
            sprintf(    'Benito ($e$ = %.3f)', datos_min.e)
            sprintf('Kechichian ($e$ = %.3f)', datos_min.e)
            sprintf(    'Benito ($e$ = %.3f)', datos_max.e)
            sprintf('Kechichian ($e$ = %.3f)', datos_max.e)
        };

        % Graficar (líneas discontinuas para primera excentricidad)
        plot(datos_min.t, error_r_analytical_min/1e3, '--', 'Color', color_analytical_min, 'LineWidth', 1.5);
        plot(datos_min.t, error_r_kechichian_min/1e3, '--', 'Color', color_comparison_min, 'LineWidth', 1.5);
        
        % Graficar (líneas continuas para última excentricidad)
        plot(datos_max.t, error_r_analytical_max/1e3, '-', 'Color', color_analytical_max, 'LineWidth', 1.5);
        plot(datos_max.t, error_r_kechichian_max/1e3, '-', 'Color', color_comparison_max, 'LineWidth', 1.5);
        
        hold off;
        
        % Configuración del gráfico
        subtitulo = sprintf('$a$ = %.0f km, $i$ = %.2f$^\\circ$, $J_3 = 0$', datos_min.a/1e3, rad2deg(datos_min.inc));
        subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||r(t)||$', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||r(t)||$ (km)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos_min.t), max(datos_min.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos_min.a/1e3);
        inc_deg = rad2deg(datos_min.inc);
        inc_clean = sprintf('%.0f', inc_deg);
        filenamePNG = sprintf('r_a%s_i%s_kechichianJ2vsJ2.png', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/kechichianJ2vsJ2', filenamePNG));        
        filenameSVG = sprintf('r_a%s_i%s_kechichianJ2vsJ2.svg', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/kechichianJ2vsJ2', filenameSVG), 'svg');
        close(gcf);


        figure; hold on;
        % Graficar (líneas discontinuas para primera excentricidad)
        plot(datos_min.t, error_v_analytical_min/1e3, '--', 'Color', color_analytical_min, 'LineWidth', 1.5);
        plot(datos_min.t, error_v_kechichian_min/1e3, '--', 'Color', color_comparison_min, 'LineWidth', 1.5);
        
        % Graficar (líneas continuas para última excentricidad)
        plot(datos_max.t, error_v_analytical_max/1e3, '-', 'Color', color_analytical_max, 'LineWidth', 1.5);
        plot(datos_max.t, error_v_kechichian_max/1e3, '-', 'Color', color_comparison_max, 'LineWidth', 1.5);
        
        hold off;
        
        % Configuración del gráfico
        subtitulo = sprintf('$a$ = %.0f km, $i$ = %.2f$^\\circ$, $J_3 = 0$', datos_min.a/1e3, rad2deg(datos_min.inc));
        subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||v(t)||$', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||v(t)||$ (km/s)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos_min.t), max(datos_min.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos_min.a/1e3);
        inc_deg = rad2deg(datos_min.inc);
        inc_clean = sprintf('%.0f', inc_deg);
        filenamePNG = sprintf('v_a%s_i%s_kechichianJ2vsJ2.png', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/kechichianJ2vsJ2', filenamePNG));        
        filenameSVG = sprintf('v_a%s_i%s_kechichianJ2vsJ2.svg', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/kechichianJ2vsJ2', filenameSVG), 'svg');
        close(gcf);
    end
end

if ~exist('figuras/kechichianJ3vsJ3', 'dir')
    mkdir('figuras', 'kechichianJ3vsJ3');
end

for ia = 1:num_a
    for iinc = 1:num_inc
        figure; hold on;
        
        % Primer valor de excentricidad (mínimo)
        ie_min = 1;
        datos_min = resultados{ia, ie_min, iinc};
        
        % Cálculo de errores para excentricidad mínima
        cowell_r_min     = [    datos_min.cowell_x_J3,     datos_min.cowell_y_J3,     datos_min.cowell_z_J3];
        analytical_r_min = [datos_min.analytical_x_J3, datos_min.analytical_y_J3, datos_min.analytical_z_J3];
        kechichian_r_min = [datos_min.kechichian_x_J3, datos_min.kechichian_y_J3, datos_min.kechichian_z_J3];

        cowell_v_min     = [    datos_min.cowell_dx_J3,     datos_min.cowell_dy_J3,     datos_min.cowell_dz_J3];
        analytical_v_min = [datos_min.analytical_dx_J3, datos_min.analytical_dy_J3, datos_min.analytical_dz_J3];
        kechichian_v_min = [datos_min.kechichian_dx_J3, datos_min.kechichian_dy_J3, datos_min.kechichian_dz_J3];
        
        error_r_analytical_min = vecnorm(cowell_r_min - analytical_r_min, 2, 2);
        error_r_kechichian_min = vecnorm(cowell_r_min - kechichian_r_min, 2, 2);

        error_v_analytical_min = vecnorm(cowell_v_min - analytical_v_min, 2, 2);
        error_v_kechichian_min = vecnorm(cowell_v_min - kechichian_v_min, 2, 2);
        
        % Último valor de excentricidad (máximo)
        ie_max = num_e;
        datos_max = resultados{ia, ie_max, iinc};
        
        % Cálculo de errores para excentricidad máxima
        cowell_r_max     = [    datos_max.cowell_x_J3,     datos_max.cowell_y_J3,     datos_max.cowell_z_J3];
        analytical_r_max = [datos_max.analytical_x_J3, datos_max.analytical_y_J3, datos_max.analytical_z_J3];
        kechichian_r_max = [datos_max.kechichian_x_J3, datos_max.kechichian_y_J3, datos_max.kechichian_z_J3];

        cowell_v_max     = [    datos_max.cowell_dx_J3,     datos_max.cowell_dy_J3,     datos_max.cowell_dz_J3];
        analytical_v_max = [datos_max.analytical_dx_J3, datos_max.analytical_dy_J3, datos_max.analytical_dz_J3];
        kechichian_v_max = [datos_max.kechichian_dx_J3, datos_max.kechichian_dy_J3, datos_max.kechichian_dz_J3];
        
        error_r_analytical_max = vecnorm(cowell_r_max - analytical_r_max, 2, 2);
        error_r_kechichian_max = vecnorm(cowell_r_max - kechichian_r_max, 2, 2);

        error_v_analytical_max = vecnorm(cowell_v_max - analytical_v_max, 2, 2);
        error_v_kechichian_max = vecnorm(cowell_v_max - kechichian_v_max, 2, 2);
        
        % Configurar leyenda
        leyendas = {
            sprintf(    'Benito ($e$ = %.3f)', datos_min.e)
            sprintf('Kechichian ($e$ = %.3f)', datos_min.e)
            sprintf(    'Benito ($e$ = %.3f)', datos_max.e)
            sprintf('Kechichian ($e$ = %.3f)', datos_max.e)
        };

        % Graficar (líneas discontinuas para primera excentricidad)
        plot(datos_min.t, error_r_analytical_min/1e3, '--', 'Color', color_analytical_min, 'LineWidth', 1.5);
        plot(datos_min.t, error_r_kechichian_min/1e3, '--', 'Color', color_comparison_min, 'LineWidth', 1.5);
        
        % Graficar (líneas continuas para última excentricidad)
        plot(datos_max.t, error_r_analytical_max/1e3, '-', 'Color', color_analytical_max, 'LineWidth', 1.5);
        plot(datos_max.t, error_r_kechichian_max/1e3, '-', 'Color', color_comparison_max, 'LineWidth', 1.5);
        
        hold off;
        
        % Configuración del gráfico
        subtitulo = sprintf('$a$ = %.0f km, $i$ = %.2f$^\\circ$, $J_2 = 0$', datos_min.a/1e3, rad2deg(datos_min.inc));
        subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||r(t)||$', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||r(t)||$ (km)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos_min.t), max(datos_min.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos_min.a/1e3);
        inc_deg = rad2deg(datos_min.inc);
        inc_clean = sprintf('%.0f', inc_deg);
        filenamePNG = sprintf('r_a%s_i%s_kechichianJ3vsJ3.png', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/kechichianJ3vsJ3', filenamePNG));        
        filenameSVG = sprintf('r_a%s_i%s_kechichianJ3vsJ3.svg', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/kechichianJ3vsJ3', filenameSVG), 'svg');
        close(gcf);


        figure; hold on;
        % Graficar (líneas discontinuas para primera excentricidad)
        plot(datos_min.t, error_v_analytical_min/1e3, '--', 'Color', color_analytical_min, 'LineWidth', 1.5);
        plot(datos_min.t, error_v_kechichian_min/1e3, '--', 'Color', color_comparison_min, 'LineWidth', 1.5);
        
        % Graficar (líneas continuas para última excentricidad)
        plot(datos_max.t, error_v_analytical_max/1e3, '-', 'Color', color_analytical_max, 'LineWidth', 1.5);
        plot(datos_max.t, error_v_kechichian_max/1e3, '-', 'Color', color_comparison_max, 'LineWidth', 1.5);
        
        hold off;
        
        % Configuración del gráfico
        subtitulo = sprintf('$a$ = %.0f km, $i$ = %.2f$^\\circ$, $J_2 = 0$', datos_min.a/1e3, rad2deg(datos_min.inc));
        subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||v(t)||$', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||v(t)||$ (km/s)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos_min.t), max(datos_min.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos_min.a/1e3);
        inc_deg = rad2deg(datos_min.inc);
        inc_clean = sprintf('%.0f', inc_deg);
        filenamePNG = sprintf('v_a%s_i%s_kechichianJ3vsJ3.png', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/kechichianJ3vsJ3', filenamePNG));        
        filenameSVG = sprintf('v_a%s_i%s_kechichianJ3vsJ3.svg', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/kechichianJ3vsJ3', filenameSVG), 'svg');
        close(gcf);
    end
end

%% Comparaciones con Hodei

if ~exist('figuras/hodeiJ2vsJ2', 'dir')
    mkdir('figuras', 'hodeiJ2vsJ2');
end

for ia = 1:num_a
    for iinc = 1:num_inc
        figure; hold on;
        
        % Primer valor de excentricidad (mínimo)
        ie_min = 1;
        datos_min = resultados{ia, ie_min, iinc};
        
        % Cálculo de errores para excentricidad mínima
        cowell_r_min     = [    datos_min.cowell_x_J2,     datos_min.cowell_y_J2,     datos_min.cowell_z_J2];
        analytical_r_min = [datos_min.analytical_x_J2, datos_min.analytical_y_J2, datos_min.analytical_z_J2];
        hodei_r_min      = [        datos_min.hodei_x,         datos_min.hodei_y,         datos_min.hodei_z];

        cowell_v_min     = [    datos_min.cowell_dx_J2,     datos_min.cowell_dy_J2,     datos_min.cowell_dz_J2];
        analytical_v_min = [datos_min.analytical_dx_J2, datos_min.analytical_dy_J2, datos_min.analytical_dz_J2];
        hodei_v_min      = [        datos_min.hodei_dx,         datos_min.hodei_dy,         datos_min.hodei_dz];
        
        error_r_analytical_min = vecnorm(cowell_r_min - analytical_r_min, 2, 2);
        error_r_hodei_min      = vecnorm(cowell_r_min -      hodei_r_min, 2, 2);

        error_v_analytical_min = vecnorm(cowell_v_min - analytical_v_min, 2, 2);
        error_v_hodei_min      = vecnorm(cowell_v_min -      hodei_v_min, 2, 2);
        
        % Último valor de excentricidad (máximo)
        ie_max = num_e;
        datos_max = resultados{ia, ie_max, iinc};
        
        % Cálculo de errores para excentricidad máxima
        cowell_r_max     = [    datos_max.cowell_x_J2,     datos_max.cowell_y_J2,     datos_max.cowell_z_J2];
        analytical_r_max = [datos_max.analytical_x_J2, datos_max.analytical_y_J2, datos_max.analytical_z_J2];
        hodei_r_max      = [        datos_max.hodei_x,         datos_max.hodei_y,         datos_max.hodei_z];

        cowell_v_max     = [    datos_max.cowell_dx_J2,     datos_max.cowell_dy_J2,     datos_max.cowell_dz_J2];
        analytical_v_max = [datos_max.analytical_dx_J2, datos_max.analytical_dy_J2, datos_max.analytical_dz_J2];
        hodei_v_max      = [        datos_max.hodei_dx,         datos_max.hodei_dy,         datos_max.hodei_dz];
        
        error_r_analytical_max = vecnorm(cowell_r_max - analytical_r_max, 2, 2);
        error_r_hodei_max      = vecnorm(cowell_r_max -      hodei_r_max, 2, 2);

        error_v_analytical_max = vecnorm(cowell_v_max - analytical_v_max, 2, 2);
        error_v_hodei_max      = vecnorm(cowell_v_max -      hodei_v_max, 2, 2);
        
        % Graficar (líneas discontinuas para primera excentricidad)
        plot(datos_min.t, error_r_analytical_min/1e3, '--', 'Color', color_analytical_min, 'LineWidth', 1.5);
        plot(datos_min.t,      error_r_hodei_min/1e3, '--', 'Color', color_comparison_min, 'LineWidth', 1.5);
        
        % Graficar (líneas continuas para última excentricidad)
        plot(datos_max.t, error_r_analytical_max/1e3, '-', 'Color', color_analytical_max, 'LineWidth', 1.5);
        plot(datos_max.t,      error_r_hodei_max/1e3, '-', 'Color', color_comparison_max, 'LineWidth', 1.5);
        
        hold off;
        
        % Configurar leyenda
        leyendas = {
            sprintf('Benito ($e$ = %.3f)', datos_min.e)
            sprintf( 'Hodei ($e$ = %.3f)', datos_min.e)
            sprintf('Benito ($e$ = %.3f)', datos_max.e)
            sprintf( 'Hodei ($e$ = %.3f)', datos_max.e)
        };
        
        % Configuración del gráfico
        subtitulo = sprintf('$a$ = %.0f km, $i$ = %.2f$^\\circ$, $J_3 = 0$', datos_min.a/1e3, rad2deg(datos_min.inc));
        subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||r(t)||$', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||r(t)||$ (km)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos_min.t), max(datos_min.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos_min.a/1e3);
        inc_deg = rad2deg(datos_min.inc);
        inc_clean = sprintf('%.0f', inc_deg);
        filenamePNG = sprintf('r_a%s_i%s_hodeiJ2vsJ2.png', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/hodeiJ2vsJ2', filenamePNG));        
        filenameSVG = sprintf('r_a%s_i%s_hodeiJ2vsJ2.svg', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/hodeiJ2vsJ2', filenameSVG), 'svg');
        close(gcf);


        figure; hold on;
        % Graficar (líneas discontinuas para primera excentricidad)
        plot(datos_min.t, error_v_analytical_min/1e3, '--', 'Color', color_analytical_min, 'LineWidth', 1.5);
        plot(datos_min.t,      error_v_hodei_min/1e3, '--', 'Color', color_comparison_min, 'LineWidth', 1.5);
        
        % Graficar (líneas continuas para última excentricidad)
        plot(datos_max.t, error_v_analytical_max/1e3, '-', 'Color', color_analytical_max, 'LineWidth', 1.5);
        plot(datos_max.t,      error_v_hodei_max/1e3, '-', 'Color', color_comparison_max, 'LineWidth', 1.5);
        
        hold off;
        
        % Configuración del gráfico
        subtitulo = sprintf('$a$ = %.0f km, $i$ = %.2f$^\\circ$, $J_3 = 0$', datos_min.a/1e3, rad2deg(datos_min.inc));
        subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||v(t)||$', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||v(t)||$ (km/s)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos_min.t), max(datos_min.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos_min.a/1e3);
        inc_deg = rad2deg(datos_min.inc);
        inc_clean = sprintf('%.0f', inc_deg);
        filenamePNG = sprintf('v_a%s_i%s_hodeiJ2vsJ2.png', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/hodeiJ2vsJ2', filenamePNG));        
        filenameSVG = sprintf('v_a%s_i%s_hodeiJ2vsJ2.svg', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/hodeiJ2vsJ2', filenameSVG), 'svg');
        close(gcf);
    end
end

if ~exist('figuras/hodeiJ2J3vsJ2J3', 'dir')
    mkdir('figuras', 'hodeiJ2J3vsJ2J3');
end

for ia = 1:num_a
    for iinc = 1:num_inc
        figure; hold on;
        
        % Primer valor de excentricidad (mínimo)
        ie_min = 1;
        datos_min = resultados{ia, ie_min, iinc};
        
        % Cálculo de errores para excentricidad mínima
        cowell_r_min     = [    datos_min.cowell_x_J2J3,     datos_min.cowell_y_J2J3,     datos_min.cowell_z_J2J3];
        analytical_r_min = [datos_min.analytical_x_J2J3, datos_min.analytical_y_J2J3, datos_min.analytical_z_J2J3];
        hodei_r_min      = [          datos_min.hodei_x,           datos_min.hodei_y,           datos_min.hodei_z];

        cowell_v_min     = [    datos_min.cowell_dx_J2J3,     datos_min.cowell_dy_J2J3,     datos_min.cowell_dz_J2J3];
        analytical_v_min = [datos_min.analytical_dx_J2J3, datos_min.analytical_dy_J2J3, datos_min.analytical_dz_J2J3];
        hodei_v_min      = [          datos_min.hodei_dx,           datos_min.hodei_dy,           datos_min.hodei_dz];
        
        error_r_analytical_min = vecnorm(cowell_r_min - analytical_r_min, 2, 2);
        error_r_hodei_min      = vecnorm(cowell_r_min -      hodei_r_min, 2, 2);

        error_v_analytical_min = vecnorm(cowell_v_min - analytical_v_min, 2, 2);
        error_v_hodei_min      = vecnorm(cowell_v_min -      hodei_v_min, 2, 2);
        
        % Último valor de excentricidad (máximo)
        ie_max = num_e;
        datos_max = resultados{ia, ie_max, iinc};
        
        % Cálculo de errores para excentricidad máxima
        cowell_r_max     = [    datos_max.cowell_x_J2J3,     datos_max.cowell_y_J2J3,     datos_max.cowell_z_J2J3];
        analytical_r_max = [datos_max.analytical_x_J2J3, datos_max.analytical_y_J2J3, datos_max.analytical_z_J2J3];
        hodei_r_max      = [          datos_max.hodei_x,           datos_max.hodei_y,           datos_max.hodei_z];

        cowell_v_max     = [    datos_max.cowell_dx_J2J3,     datos_max.cowell_dy_J2J3,     datos_max.cowell_dz_J2J3];
        analytical_v_max = [datos_max.analytical_dx_J2J3, datos_max.analytical_dy_J2J3, datos_max.analytical_dz_J2J3];
        hodei_v_max      = [          datos_max.hodei_dx,           datos_max.hodei_dy,           datos_max.hodei_dz];
        
        error_r_analytical_max = vecnorm(cowell_r_max - analytical_r_max, 2, 2);
        error_r_hodei_max      = vecnorm(cowell_r_max -      hodei_r_max, 2, 2);

        error_v_analytical_max = vecnorm(cowell_v_max - analytical_v_max, 2, 2);
        error_v_hodei_max      = vecnorm(cowell_v_max -      hodei_v_max, 2, 2);
        
        % Graficar (líneas discontinuas para primera excentricidad)
        plot(datos_min.t, error_r_analytical_min/1e3, '--', 'Color', color_analytical_min, 'LineWidth', 1.5);
        plot(datos_min.t,      error_r_hodei_min/1e3, '--', 'Color', color_comparison_min, 'LineWidth', 1.5);
        
        % Graficar (líneas continuas para última excentricidad)
        plot(datos_max.t, error_r_analytical_max/1e3, '-', 'Color', color_analytical_max, 'LineWidth', 1.5);
        plot(datos_max.t,      error_r_hodei_max/1e3, '-', 'Color', color_comparison_max, 'LineWidth', 1.5);
        
        hold off;
        
        % Configurar leyenda
        leyendas = {
            sprintf('Benito ($e$ = %.3f)', datos_min.e)
            sprintf(' Hodei ($e$ = %.3f)', datos_min.e)
            sprintf('Benito ($e$ = %.3f)', datos_max.e)
            sprintf(' Hodei ($e$ = %.3f)', datos_max.e)
        };
        
        % Configuración del gráfico
        subtitulo = sprintf('$a$ = %.0f km, $i$ = %.2f$^\\circ$, $J_3 \\neq 0$', datos_min.a/1e3, rad2deg(datos_min.inc));
        subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||r(t)||$', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||r(t)||$ (km)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos_min.t), max(datos_min.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos_min.a/1e3);
        inc_deg = rad2deg(datos_min.inc);
        inc_clean = sprintf('%.0f', inc_deg);
        filenamePNG = sprintf('r_a%s_i%s_hodeiJ2J3vsJ2J3.png', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/hodeiJ2J3vsJ2J3', filenamePNG));      
        filenameSVG = sprintf('r_a%s_i%s_hodeiJ2J3vsJ2J3.svg', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/hodeiJ2J3vsJ2J3', filenameSVG), 'svg');
        close(gcf);


        figure; hold on;
        % Graficar (líneas discontinuas para primera excentricidad)
        plot(datos_min.t, error_v_analytical_min/1e3, '--', 'Color', color_analytical_min, 'LineWidth', 1.5);
        plot(datos_min.t,      error_v_hodei_min/1e3, '--', 'Color', color_comparison_min, 'LineWidth', 1.5);
        
        % Graficar (líneas continuas para última excentricidad)
        plot(datos_max.t, error_v_analytical_max/1e3, '-', 'Color', color_analytical_max, 'LineWidth', 1.5);
        plot(datos_max.t,      error_v_hodei_max/1e3, '-', 'Color', color_comparison_max, 'LineWidth', 1.5);
        
        hold off;
        
        % Configuración del gráfico
        subtitulo = sprintf('$a$ = %.0f km, $i$ = %.2f$^\\circ$, $J_3 = 0$', datos_min.a/1e3, rad2deg(datos_min.inc));
        subtitle(subtitulo, 'Interpreter', 'latex', 'FontSize', 10);
        title('Error de $||v(t)||$', 'Interpreter', 'latex');
        xlabel('$t$ (s)', 'Interpreter', 'latex');
        ylabel('Error de $||v(t)||$ (km/s)', 'Interpreter', 'latex');
        legend(leyendas, 'Interpreter', 'latex', 'Location', 'best');
        xlim([min(datos_min.t), max(datos_min.t)]);
        grid on; 
        box off;
        
        % Guardar figura
        a_clean = sprintf('%.0f', datos_min.a/1e3);
        inc_deg = rad2deg(datos_min.inc);
        inc_clean = sprintf('%.0f', inc_deg);
        filenamePNG = sprintf('v_a%s_i%s_hodeiJ2J3vsJ2J3.png', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/hodeiJ2J3vsJ2J3', filenamePNG));        
        filenameSVG = sprintf('v_a%s_i%s_hodeiJ2J3vsJ2J3.svg', a_clean, inc_clean);
        saveas(gcf, fullfile('figuras/hodeiJ2J3vsJ2J3', filenameSVG), 'svg');
        close(gcf);
    end
end