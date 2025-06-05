function dXdt = odeClohessyWiltshireJ2J3(tOde, X, J2, J3, R, a, n, INC, omega, theta_0)
%% odeClohessyWiltshireJ2J3 Movimiento relativo perturbado con J2 y J3 (modelo CW modificado).
%
% Esta función define el sistema de ecuaciones diferenciales del movimiento relativo
% linealizado de un satélite respecto a una órbita de referencia circular usando
% el modelo de Clohessy-Wiltshire (Hill), incluyendo perturbaciones zonales J2 y J3.
%
% Sintaxis:
%   dXdt = odeClohessyWiltshireJ2J3(t, X, J2, J3, R, a, n, INC, omega, theta_0)
%
% Entradas:
%   tOde     : Tiempo actual de integración [s]
%   X        : Vector de estado relativo en el sistema Hill [6x1]:
%                - [x; y; z; dx; dy; dz], posición y velocidad relativas en m y m/s
%   J2, J3   : Coeficientes armónicos zonales del potencial gravitacional (sin unidades)
%   R        : Radio ecuatorial del cuerpo central (e.g., Tierra) [m]
%   a        : Semieje mayor de la órbita de referencia [m]
%   n        : Velocidad angular media de la órbita de referencia [rad/s]
%   INC      : Inclinación orbital de la órbita de referencia [rad]
%   omega    : Argumento del perigeo de la órbita de referencia [rad]
%   theta_0  : Anomalía verdadera inicial de la órbita de referencia [rad]
%
% Salida:
%   dXdt     : Derivada del estado relativo (6x1), en formato [dx; dy; dz; ddx; ddy; ddz]
% 
% Ejemplo de uso con un integrador numérico:
%   [t, Y] = ode45(@(t, X) odeClohessyWiltshireJ2J3(t, X, J2, J3, R, a, n, INC, omega, theta_0), tspan, X0);
% 
% Referencias:
%   - Modelo de Hill-Clohessy-Wiltshire para movimiento relativo
%   - Incorporación de efectos de J2 y J3 proyectados en el sistema local Euler-Hill
%   - Válido para órbitas cuasi-circulares y perturbaciones moderadas

    %% --- Términos perturbativos precomputados ---
    J2_term = J2 * R^2 * n^2 / a;
    J3_term = J3 * R^3 * n^2 / (2 * a^2);

    %% Ángulo argumental instantáneo: u = ω + θ
    u = omega + theta_0 + n * tOde;

    %% --- Aceleraciones perturbativas debido a J2 ---
    fr_J2     =   J2_term * (3/2) * (3*sin(u)^2*sin(INC)^2 - 1);     % Radial
    ftheta_J2 = - J2_term * (3/2) * sin(2*u) * sin(INC)^2;           % Transversal (along-track)
    fh_J2     = - J2_term * 3 * cos(INC) * sin(u) * sin(INC);        % Normal (out-of-plane)

    %% --- Aceleraciones perturbativas debido a J3 ---
    fr_J3     =   J3_term * 4 * sin(INC) * sin(u) * ...
                       (5 * sin(u)^2 * sin(INC)^2 - 3);              % Radial
    ftheta_J3 = - J3_term * 3 * cos(u) * sin(INC) * ...
                       (5 * sin(u)^2 * sin(INC)^2 - 1);              % Transversal
    fh_J3     = - J3_term * 3 * cos(INC) * ...
                       (5 * sin(u)^2 * sin(INC)^2 - 1);              % Normal

    %% --- Aceleración total (suma de J2 + J3) ---
    fr     = fr_J2     + fr_J3;
    ftheta = ftheta_J2 + ftheta_J3;
    fh     = fh_J2     + fh_J3;

    %% --- Estado relativo ---
    x  = X(1);  y  = X(2);  z  = X(3);
    dx = X(4);  dy = X(5);  dz = X(6);

    %% --- Dinámica relativa con perturbaciones ---
    ddx = fr     + 3*n^2*x + 2*n*dy;  % Eje radial-local
    ddy = ftheta - 2*n*dx;            % Eje transversal-local
    ddz = fh     - n^2*z;             % Eje normal-local

    %% --- Derivadas del estado ---
    dXdt = [dx; dy; dz; ddx; ddy; ddz];
end
