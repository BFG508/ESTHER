function dXdt = cowellZonals(X0, tSpan, mu, R, zonals, opts)
%% cowellZonals Integración orbital con método de Cowell incluyendo zonales J2–J6.
%
% Esta función simula la propagación orbital de un cuerpo bajo la influencia de la 
% atracción gravitacional central (ley de Newton) y las perturbaciones armónicas 
% zonales J2 a J6, utilizando el método directo de Cowell.
%
% Sintaxis:
%   dXdt = cowellZonals(X0, tSpan, mu, R, zonals)
%   dXdt = cowellZonals(X0, tSpan, mu, R, zonals, opts)
%
% Entradas:
%   X0     : Vector de estado inicial [6x1] → [x; y; z; vx; vy; vz]
%   tSpan  : Intervalo de integración [t0, tf] (s)
%   mu     : Parámetro gravitacional estándar (G·M) del cuerpo central [m³/s² o km³/s²]
%   R      : Radio ecuatorial del cuerpo central [m o km]
%   zonals : Vector de coeficientes zonales [J2, J3, J4, J5, J6]
%   opts   : (Opcional) Estructura de opciones para el integrador ODE (ej. creada con odeset)
%            - Si no se proporciona, se utilizarán las opciones por defecto de MATLAB
%
% Salida:
%   dXdt : Matriz [Nx6] con estados (posición y velocidad) a lo largo de la integración

    % Si no se especifica 'opts', usar configuración por defecto
    if nargin < 6 || isempty(opts)
        opts = [];
    end

    % Aceleración central (Kepler)
    aKepler = @(X) -mu * X(1:3) / norm(X(1:3))^3;

    % Aceleración por perturbaciones zonales (J2 a J6)
    aZonals = @(X) accelerationJ(X(1:3), mu, R, zonals(1), zonals(2), zonals(3), zonals(4), zonals(5));

    % Dinámica combinada total
    fCowell = @(t, X) odeCowell(X, aKepler(X(1:3)) + aZonals(X(1:3)));

    % Integración temporal con ode113
    [~, dXdt] = ode113(fCowell, tSpan, X0, opts);
end
