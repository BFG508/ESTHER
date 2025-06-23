function [SMA, ECC, INC, RAAN, AOP, TA, AOLon, AOLat] = ECI2OE(mu, r, v)
% ECI2OE Convierte vectores de posición y velocidad ECI en elementos orbitales clásicos.
%
% Esta función toma los vectores de posición y velocidad en el sistema Earth-Centered Inertial (ECI)
% y calcula los elementos orbitales keplerianos clásicos, incluyendo parámetros derivados.
%
% ENTRADAS:
%   mu  : Parámetro gravitacional [m³/s²]
%   r   : Vector de posición en marco ECI [m] (3x1 o Nx3)
%   v   : Vector de velocidad en marco ECI [m/s] (mismas dimensiones que r)
%
% SALIDAS:
%   SMA    : Semieje mayor [m]
%   ECC    : Excentricidad (adimensional)
%   INC    : Inclinación [rad]
%   RAAN   : Ascensión Recta del Nodo Ascendente [rad]
%   AOP    : Argumento del Periapsis [rad]
%   TA     : Anomalía Verdadera [rad]
%   AOLat  : Argumento de Latitud [rad] (TA + AOP)
%   AOLon  : Longitud del Periapsis [rad] (RAAN + AOP)

    % Ajustar dimensiones a Nx3 para manejar múltiples entradas
    if size(r, 1) == 3 && size(r, 2) ~= 3
        r = r';
    end
    if size(v, 1) == 3 && size(v, 2) ~= 3
        v = v';
    end
    N = size(r, 1);  % Número de conjuntos de datos
    
    % 1. Calcular momento angular específico (h = r × v)
    h = cross(r, v, 2);
    h_norm = sqrt(sum(h.^2, 2));  % Módulo del momento angular
    
    % 2. Calcular vector excentricidad (e = (v × h)/μ - r/|r|)
    r_norm = sqrt(sum(r.^2, 2));
    e_vec = cross(v, h, 2) ./ mu - r ./ r_norm;
    ECC = sqrt(sum(e_vec.^2, 2));  % Excentricidad
    
    % 3. Calcular semieje mayor mediante energía específica
    v_norm = sqrt(sum(v.^2, 2));
    E = v_norm.^2/2 - mu./r_norm;  % Energía específica
    SMA = -mu./(2*E);              % Semieje mayor
    
    % 4. Calcular inclinación (i = acos(h_z / |h|))
    arg_inc = h(:,3) ./ h_norm;
    INC = acos(arg_inc);
    
    % 5. Calcular RAAN (Ω = acos(h_x / -h_y))
    RAAN = atan2(h(:,1), -h(:,2));
    
    % 6. Calcular argumento del periapsis (ω = acos(n·e / |n||e|))
    arg1_aop = e_vec(:, 3) ./ sin(INC);
    arg2_aop = e_vec(:, 1) .* cos(RAAN) + e_vec(:, 2) .* sin(RAAN);
    AOP = atan2(arg1_aop, arg2_aop);
    
    % 7. Calcular anomalía verdadera (ν = acos(e·r / |e||r|))
    TA = zeros(N, 1);
    valid_nu = (ECC > 1e-10);
    if any(valid_nu)
        e_dot_r = dot(e_vec(valid_nu,:), r(valid_nu,:), 2);
        arg_ta = e_dot_r ./ (ECC(valid_nu) .* r_norm(valid_nu));
        arg_ta = min(max(arg_ta, -1), 1);  % Clamping
        TA(valid_nu) = acos(arg_ta);
    end

    % 8. Calcular argumento de latitud y longitud
    AOLat =   TA + AOP;
    AOLon = RAAN + AOP;
end
