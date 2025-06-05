function dX = odeHodei(X, mu, Re, alpha, J2)
%% odeHodei Dinámica orbital con efecto J2 y empuje bajo en dirección radial.
%
% Esta función define un sistema de ecuaciones diferenciales en coordenadas orbitales no inerciales,
% incorporando la perturbación del coeficiente zonal J2 y un empuje continuo de baja intensidad
% aplicado en la dirección radial. El modelo es útil para estudiar la evolución orbital de satélites
% o sondas espaciales con propulsión continua y perturbaciones gravitacionales.
%
% Sintaxis:
%   dX = odeHodei(X, mu, Re, alpha, J2)
%
% Entradas:
%   X     : Vector de estado (6x1), con las siguientes componentes:
%           X(1) = r      → Radio orbital (m)
%           X(2) = θ      → Argumento de latitud (rad)
%           X(3) = ν      → Ascensión recta del nodo ascendente (rad)
%           X(4) = R      → Velocidad radial (m/s)
%           X(5) = Θ      → Módulo del momento angular orbital total (m²/s)
%           X(6) = N      → Proyección del momento angular sobre el eje Z (m²/s)
%
%   mu    : Parámetro gravitacional estándar (μ = G·M) del cuerpo central [m³/s²]
%   Re    : Radio ecuatorial del cuerpo central [m]
%   alpha : Aceleración continua de bajo empuje en la dirección radial [m/s²]
%           - Representa una propulsión constante en la dirección del radio
%   J2    : Coeficiente de achatamiento zonal del cuerpo central
%
% Salida:
%   dX    : Derivada del estado (6x1), con componentes:
%           [dr/dt; dθ/dt; dν/dt; dR/dt; dΘ/dt; dN/dt]
% 
% Ejemplo de uso con un integrador numérico:
%   [t, Y] = ode45(@(t,X) odeHodei(X, mu, Re, alpha, J2), tspan, X0);

    % === Estado extraído del vector X ===
    r     = X(1);   % Radio orbital
    theta = X(2);   % Argumento de latitud (en plano orbital)
    nu    = X(3);   % Ascensión recta del nodo ascendente (RAAN)
    R     = X(4);   % Velocidad radial
    Theta = X(5);   % Momento angular orbital total
    N     = X(6);   % Proyección del momento angular sobre el eje Z

    % === Cálculos auxiliares ===
    p     = Theta^2 / mu;                 % Semi-latus rectum (parámetro orbital)
    sigma = 0.5 * J2 * (Re / p)^2;        % Factor de perturbación J2
    k     = p/r - 1;                      % Factor de forma orbital
    c     = N / Theta;                    % cos(inclinación)
    s     = sqrt(1 - c^2);                % sin(inclinación)

    % === Derivadas del estado ===
    drdt     = R;
    dthetadt = Theta/r^2 * (1 + 3*sigma*(1 + k)*c^2*(1 - cos(2*theta)));
    dnudt    = Theta/r^2 * (-3*sigma*(1 + k)*c*(1 - cos(2*theta)));
    dRdt     = Theta^2/r^3 * ...
               (k/(1 + k) - 3*sigma*(1 + k)*(1 - 1.5*s^2*(1 - cos(2*theta)))) + alpha;
    dThetadt = Theta/r^2 * (-3*sigma*(1 + k)*s^2*(1 - cos(2*theta)));
    dNdt     = 0;  % N es constante (conservación del momento angular proyectado)

    % === Salida del sistema dinámico ===
    dX = [drdt; dthetadt; dnudt; dRdt; dThetadt; dNdt];
end
