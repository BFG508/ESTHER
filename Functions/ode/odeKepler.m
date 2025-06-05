function dXdt = odeKepler(X, mu)
%% odeKepler Sistema dinámico del problema de Kepler (dos cuerpos).
%
% Esta función define las ecuaciones diferenciales de primer orden para
% el movimiento de un cuerpo bajo la atracción gravitatoria de otro,
% según la ley de la gravitación universal de Newton.
%
% Sintaxis:
%   dXdt = odeKepler(X, mu)
%
% Entradas:
%   X   : Vector de estado de dimensión 6x1 [r; v], donde:
%           - r = posición tridimensional [x; y; z]
%           - v = velocidad tridimensional [vx; vy; vz]
%
%   mu  : Parámetro gravitacional estándar (μ = G·M), donde:
%           - G es la constante de gravitación universal
%           - M es la masa del cuerpo central (por ejemplo, un planeta)
%
% Salida:
%   dXdt : Derivada temporal del estado [dr/dt; dv/dt], también de tamaño 6x1
%
% Ejemplo de uso con un integrador numérico:
%   [t, Y] = ode45(@(t, X) odeKepler(X, mu), tspan, X0);

    % Descomposición del estado
    r = X(1:3);  % Posición [x; y; z]
    v = X(4:6);  % Velocidad [vx; vy; vz]

    % Dinámica orbital
    dr = v;                            % dr/dt = v
    dv = -mu / norm(r)^3 * r;          % dv/dt = aceleración gravitacional

    % Ensamblaje del vector de derivadas
    dXdt = [dr; dv];
end
