function [r, theta, nu, R, Theta, N] = cartesian2PolarNodal(x, y, z, dx, dy, dz)
%% cartesian2PolarNodal Conversión de coordenadas cartesianas a variables polares-nodales.
%
% Esta función transforma la posición y velocidad en coordenadas cartesianas
% a un conjunto de variables polares-nodales. Estas variables son útiles en 
% dinámica orbital para describir el estado orbital de un cuerpo en un 
% sistema no inercial ligado al plano orbital.
%
% Entradas:
%   x, y, z     : Componentes de posición en el sistema cartesiano [m]
%   dx, dy, dz  : Componentes de velocidad en el sistema cartesiano [m/s]
%
% Salidas:
%   r      : Radio orbital (norma del vector de posición) [m]
%   theta  : Argumento de latitud (posición angular en el plano orbital) [rad]
%   nu     : Ascensión recta del nodo ascendente (RAAN) [rad]
%   R      : Velocidad radial (proyección de la velocidad en la dirección radial) [m/s]
%   Theta  : Módulo del momento angular orbital [m²/s]
%   N      : Componente z del momento angular (Theta * cos(i)) [m²/s]
%
% Nota:
%   La conversión asume un sistema inercial con origen en el cuerpo central
%   y el eje Z como referencia para el plano ecuatorial.

    % --- Paso 1: Vectores de posición y velocidad ---
    r_vec = [x; y; z]; % Vector de posición [m]
    v_vec = [dx; dy; dz]; % Vector de velocidad [m/s]
    
    % --- Paso 2: Módulo del vector de posición ---
    r = norm(r_vec); % Radio orbital [m]
    
    % --- Paso 3: Momento angular ---
    h_vec = cross(r_vec, v_vec); % Vector del momento angular específico [m²/s]
    Theta = norm(h_vec); % Módulo del momento angular [m²/s]
    h_hat = h_vec / Theta;
    N = h_vec(3); % Componente sobre el eje Z [m²/s]
    
    % --- Paso 4: Línea de nodos ---
    n_vec = cross([0; 0; 1], h_vec);
    n_norm = norm(n_vec);
    
    if n_norm < norm(h_vec) * eps
        nu = 0; % Caso especial: órbitas ecuatoriales
    else
        nu = wrapTo2Pi(atan2(n_vec(2), n_vec(1)));
    end
    
    % --- Paso 5: Base del plano orbital usando el nodo ascendente ---
    if n_norm > eps
        n_hat = n_vec / n_norm; % Dirección del nodo ascendente
        m_hat = cross(h_hat, n_hat); % Perpendicular en el plano orbital
    else
        % Caso ecuatorial
        n_hat = [1; 0; 0];
        m_hat = [0; 1; 0];
    end
    
    % --- Paso 6: Argumento de latitud ---
    r_n = dot(r_vec, n_hat);
    r_m = dot(r_vec, m_hat);
    theta = wrapTo2Pi(atan2(r_m, r_n));
    
    % --- Paso 7: Velocidad radial ---
    R = dot(r_vec, v_vec) / r;
end
