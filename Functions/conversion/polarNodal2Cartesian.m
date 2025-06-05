function [x, y, z, dx, dy, dz] = polarNodal2Cartesian(r, theta, nu, R, Theta, N)
%% polarNodal2Cartesian Conversión de variables polares-nodales a coordenadas cartesianas.
%
% Esta función transforma un conjunto de variables polares-nodales orbitales
% (r, θ, ν, R, Θ, N) a coordenadas cartesianas en un marco inercial.
% El algoritmo es compatible con múltiples instantes, operando sobre vectores columna (Nx1).
%
% Sintaxis:
%   [x, y, z, dx, dy, dz] = polarNodal2Cartesian(r, theta, nu, R, Theta, N)
%
% Entradas:
%   r      : Radio orbital (distancia al cuerpo central) [m]                   - (Nx1)
%   theta  : Argumento de latitud (posición angular en el plano orbital) [rad] - (Nx1)
%   nu     : Ascensión recta del nodo ascendente (RAAN) [rad]                  - (Nx1)
%   R      : Velocidad radial [m/s]                                            - (Nx1)
%   Theta  : Momento angular orbital total [m²/s]                              - (Nx1)
%   N      : Componente del momento angular sobre el eje Z [m²/s]              - (Nx1)
%
% Salidas:
%   x, y, z    : Posición en coordenadas cartesianas [m]         - (Nx1)
%   dx, dy, dz : Velocidad en coordenadas cartesianas [m/s]      - (Nx1)
%
% Notas:
%   - Si N/Theta excede el dominio de acos (es decir, |cos(i)| > 1), se lanza un error.
%   - La conversión asume que el eje Z es la dirección del momento angular del cuerpo central.
%   - La rotación se realiza aplicando primero una rotación en Z(nu) y luego una rotación en X(i).
%

    % --- Validación inicial ---
    assert(all(r > 0), 'El radio orbital debe ser positivo.');
    assert(all(Theta > 0), 'El momento angular total debe ser positivo.');

    % --- Cálculo de inclinación ---
    cosi = N ./ Theta;
    i = acos(cosi);

    % --- Direcciones unitarias en el plano orbital ---
    er_x = cos(theta);
    er_y = sin(theta);
    et_x = -sin(theta);
    et_y = cos(theta);

    % --- Posición y velocidad en marco orbital ---
    r_orb_x = r .* er_x;
    r_orb_y = r .* er_y;

    v_orb_x = R .* er_x + (Theta ./ r) .* et_x;
    v_orb_y = R .* er_y + (Theta ./ r) .* et_y;

    % --- Parámetros de rotación ---
    cos_nu = cos(nu);
    sin_nu = sin(nu);
    cos_i  = cos(i);
    sin_i  = sin(i);

    % --- Rotación Z(nu) seguida de X(i) para posiciÃ³n ---
    x = cos_nu .* r_orb_x - sin_nu .* cos_i .* r_orb_y;
    y = sin_nu .* r_orb_x + cos_nu .* cos_i .* r_orb_y;
    z = sin_i .* r_orb_y;

    % --- Rotación Z(nu) seguida de X(i) para velocidad ---
    dx = cos_nu .* v_orb_x - sin_nu .* cos_i .* v_orb_y;
    dy = sin_nu .* v_orb_x + cos_nu .* cos_i .* v_orb_y;
    dz = sin_i .* v_orb_y;
end
