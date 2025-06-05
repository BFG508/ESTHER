function [dXdt, dXdt_Kepler, dXdt_J2J3_I] = odeMotion(t, X0, X0I, mu, J2, J3, R, a, e, n, INC, Omega, omega, theta_0, opts)
%% odeMotion Combina movimiento kepleriano y movimiento relativo perturbado por J2–J3.
%
% Esta función calcula la evolución temporal del estado orbital en coordenadas inerciales,
% combinando el movimiento kepleriano del centro de masas y el movimiento relativo linealizado
% del satélite respecto al cuerpo central, incluyendo perturbaciones zonales J2 y J3.
%
% Sintaxis:
%   [dXdt, dXdt_Kepler, dXdt_J2J3_I] = odeMotion(t, X0, X0I, mu, J2, J3, R, a, e, n, INC, Omega, omega, theta_0)
%   [dXdt, dXdt_Kepler, dXdt_J2J3_I] = odeMotion(..., opts)
%
% Entradas:
%   t        : Vector de tiempos de evaluación [Nx1] (s)
%   X0       : Estado inicial relativo en sistema orbital (Euler-Hill) [6x1]
%   X0I      : Estado inicial inercial kepleriano del centro de masas [6x1]
%   mu       : Parámetro gravitacional del cuerpo central [km³/s² o m³/s²]
%   J2, J3   : Coeficientes zonales del potencial gravitacional
%   R        : Radio ecuatorial del cuerpo central [km o m]
%   a        : Semieje mayor de la órbita de referencia [km o m]
%   e        : Excentricidad de la órbita de referencia
%   n        : Velocidad angular media (mean motion) [rad/s]
%   INC      : Inclinación orbital [rad]
%   Omega    : Ascensión recta del nodo ascendente (RAAN) [rad]
%   omega    : Argumento del perigeo [rad]
%   theta_0  : Anomalía verdadera inicial [rad]
%   opts     : (Opcional) Estructura de opciones para el integrador ODE (por defecto: vacío = estándar de MATLAB)
%
% Salidas:
%   dXdt         : Movimiento total en sistema inercial [Nx6] (posición + velocidad)
%   dXdt_Kepler  : Movimiento kepleriano puro en sistema inercial [Nx6]
%   dXdt_J2J3_I  : Movimiento perturbado por J2/J3 en sistema inercial [Nx6]
%
% Notas:
%   - Internamente se integran por separado las trayectorias: kepleriana e inercialización del
%     movimiento relativo perturbado (modelo CW elíptico con J2 y J3).
%   - La transformación se realiza desde el sistema orbital a inercial mediante 3 rotaciones sucesivas.
%

    % --- Si no se especifica 'opts', usar configuración por defecto ---
    if nargin < 15 || isempty(opts)
        opts = []; 
    end

    % --- Movimiento kepleriano del centro de masas ---
    [~, dXdt_Kepler] = ode113(@(tOde, X) odeKepler(X, mu), t, X0I, opts);

    % --- Movimiento relativo perturbado por J2 y J3 (modelo CW elíptico) ---
    [~, dXdt_J2J3] = ode113(@(tOde, X) ...
        odeClohessyWiltshireECCJ2J3(tOde, X, J2, J3, R, a, e, n, INC, omega, theta_0), ...
        t, X0, opts);

    % --- Cálculo del ángulo de latitud corregido ---
    delta = 2 * (theta_0 - e * (1 + sqrt(1 - e^2)) * sin(theta_0));
    sd2 = sin(delta / 2);
    cd2 = cos(delta / 2);
    theta_c = omega + theta_0 - 2 * e * sd2;

    % --- Corrección angular u(t) dependiente del tiempo ---
    t = t(:);
    snt = sin(n * t);
    cnt = cos(n * t);
    u = theta_c + n * t + 2 * e * (sd2 * cnt + cd2 * snt);

    % --- Transformaciones a sistema inercial ---
    Rz_RAAN = [cos(Omega), -sin(Omega), 0;
               sin(Omega),  cos(Omega), 0;
                        0,           0, 1];

    Rx_INC = [1,        0,         0;
              0, cos(INC), -sin(INC);
              0, sin(INC),  cos(INC)];

    dXdt_J2J3_I = zeros(length(t), 6);

    for j = 1:length(t)
        Rz_u = [cos(u(j)), -sin(u(j)), 0;
                sin(u(j)),  cos(u(j)), 0;
                        0,          0, 1];

        % Rotación total: sistema Euler-Hill → sistema inercial
        R_EUtoI = Rz_RAAN * Rx_INC * Rz_u;

        % Extraer posición y velocidad relativas perturbadas
        rJ2J3 = dXdt_J2J3(j, 1:3)';
        vJ2J3 = dXdt_J2J3(j, 4:6)';

        % Transformar a marco inercial
        rJ2J3_I = R_EUtoI * rJ2J3;
        vJ2J3_I = R_EUtoI * vJ2J3;

        % Almacenar resultado transformado
        dXdt_J2J3_I(j, :) = [rJ2J3_I; vJ2J3_I]';
    end

    % --- Composición del movimiento total ---
    dXdt = dXdt_Kepler + dXdt_J2J3_I;
end
