function dXdt = odeClohessyWiltshireECCJ2J3(tOde, X, J2, J3, R, a, e, n, INC, omega, theta_0)
%% odeClohessyWiltshireECCJ2J3 Movimiento relativo en órbita elíptica con efectos J2 y J3.
%
% Esta función modela el movimiento relativo linealizado entre dos satélites
% en una órbita elíptica utilizando un modelo extendido de Clohessy-Wiltshire,
% que incluye las perturbaciones gravitatorias zonales J2 y J3 del cuerpo central.
%
% Sintaxis:
%   dXdt = odeClohessyWiltshireECCJ2J3(t, X, J2, J3, R, a, e, n, INC, omega, theta_0)
%
% Entradas:
%   tOde     : Tiempo actual de integración [s]
%   X        : Vector de estado relativo [x; y; z; dx; dy; dz] (posición y velocidad) [6x1]
%   J2, J3   : Coeficientes zonales del potencial gravitacional (sin unidades)
%   R        : Radio ecuatorial del cuerpo central (e.g. Tierra) [m]
%   a        : Semieje mayor de la órbita de referencia [m]
%   e        : Excentricidad de la órbita de referencia
%   n        : Movimiento medio orbital (mean motion) [rad/s]
%   INC      : Inclinación orbital de la órbita de referencia [rad]
%   omega    : Argumento del perigeo [rad]
%   theta_0  : Anomalía verdadera inicial [rad]
%
% Salida:
%   dXdt     : Derivadas del vector de estado [dx; dy; dz; ddx; ddy; ddz] [6x1]
%
% Ejemplo de uso con un integrador numérico:
%   [t, Y] = ode45(@(t, X) odeClohessyWiltshireECCJ2J3(t, X, J2, J3, R, a, e, n, INC, omega, theta_0), tspan, X0);
%
% Referencias:
%   - Clohessy-Wiltshire modificado para órbitas elípticas.
%   - Modelo de perturbaciones zonales basado en desarrollo armónico del campo gravitacional.
%   - Corrección de fase y oscilación para incluir efectos de la excentricidad.
%

    %% --- Términos orbitales con corrección por excentricidad ---
    r4 =  (a * (1 - e^2)^4)   / (4 * e * cos(theta_0) + 1);   % Normalización para J2
    r5 = -(a^2 * (1 - e^2)^5) / (5 * e * cos(theta_0) + 1);   % Normalización para J3

    %% --- Corrección angular para órbitas elípticas ---
    delta = 2 * (theta_0 - e * (1 + sqrt(1 - e^2)) * sin(theta_0));
    sd2 = sin(delta/2);
    cd2 = cos(delta/2);

    theta_c = omega + theta_0 - 2 * e * sd2;  % Ángulo corregido base

    %% --- Oscilaciones instantáneas en el ángulo u corregido ---
    snt = sin(n * tOde);
    cnt = cos(n * tOde);
    u = theta_c + n * tOde + 2 * e * (cd2 * snt + sd2 * cnt);  % Argumento de latitud corregido

    %% --- Aceleraciones perturbativas zonales ---
    % Radial (x)
    fr = (3 * J2 * R^2 * n^2 * (3 * sin(INC)^2 * sin(u)^2 - 1)) / (2 * r4) ...
       + (2 * J3 * R^3 * n^2 * sin(INC) * sin(u) * ...
         (5 * sin(INC)^2 * sin(u)^2 - 3)) / r5;

    % Transversal (y)
    ftheta = - (3 * J2 * R^2 * n^2 * sin(2*u) * sin(INC)^2) / (2 * r4) ...
             + (3 * J3 * R^3 * n^2 * sin(INC) * cos(u) * ...
               (5 * sin(INC)^2 * (cos(u)^2 - 1) + 1)) / (2 * r5);

    % Normal (z)
    fh = - (3 * J2 * R^2 * n^2 * cos(INC) * sin(INC) * sin(u)) / r4 ...
         + (3 * J3 * R^3 * n^2 * cos(INC) * ...
           (5 * sin(u)^2 * (cos(INC)^2 - 1) + 1)) / (2 * r5);

    %% --- Estado relativo actual ---
    x  = X(1); y  = X(2); z  = X(3);
    dx = X(4); dy = X(5); dz = X(6);

    %% --- Dinámica relativa perturbada (modelo CW extendido) ---
    ddx = fr     + 3 * n^2 * x + 2 * n * dy;
    ddy = ftheta - 2 * n * dx;
    ddz = fh     - n^2 * z;

    %% --- Salida: derivadas del estado relativo ---
    dXdt = [dx; dy; dz; ddx; ddy; ddz];
end
