function dXdt = odeCowell(X, a)
%% odeCowell Ecuaciones dinámicas orbitales mediante el método de Cowell.
%
% Esta función implementa el método de Cowell para integrar directamente
% la ecuación del movimiento de un cuerpo bajo cualquier aceleración total dada.
% Es especialmente útil para simulaciones donde se desea incluir perturbaciones
% no conservativas o efectos adicionales.
%
% Sintaxis:
%   dXdt = odeCowell(X, a)
%
% Entradas:
%   X : Vector de estado (6x1) con:
%         - X(1:3) = [x; y; z]   → posición en coordenadas cartesianas (m)
%         - X(4:6) = [vx; vy; vz] → velocidad en coordenadas cartesianas (m/s)
%
%   a : Vector de aceleración total (3x1) aplicada al cuerpo, en m/s²
%         - Puede incluir componentes gravitacionales y perturbativas
%
% Salida:
%   dXdt : Derivada del estado (6x1), compuesta por:
%           - dr/dt = v
%           - dv/dt = a
%
% Ejemplo de uso con un integrador numérico:
%   [t, Y] = ode45(@(t, X) odeCowell(X, a_total(t,X)), tspan, X0);

    % Extraer velocidad desde el estado
    v = X(4:6);  % [vx; vy; vz]

    % Devolver el sistema de ecuaciones diferenciales
    dXdt = [v; a];
end
