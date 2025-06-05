function dXdt = hodeiMotionCartesian(t, X0, mu, Re, alpha, J2, opts)
%% hodeiMotionCartesian Convierte el movimiento orbital desde coordenadas cartesianas a nodales-polares,
% simula su dinámica bajo bajo empuje radial constante y perturbación J2, 
% y devuelve el resultado nuevamente en coordenadas cartesianas.
%
% Esta función utiliza una transformación a variables polares-nodales para 
% calcular la dinámica mediante el modelo hodeiMotion, que incluye la 
% perturbación gravitacional J2 y un empuje radial constante. Posteriormente, 
% los resultados se reconvierten a coordenadas cartesianas.
%
% INPUTS:
%   t       - Tiempo actual de evaluación [s]
%   X0      - Vector de estado en coordenadas cartesianas:
%             [x; y; z; dx; dy; dz]
%   mu      - Parámetro gravitacional estándar del cuerpo central [m^3/s^2]
%   Re      - Radio ecuatorial del cuerpo central [m]
%   alpha   - Aceleración constante de bajo empuje en dirección radial [m/s^2]
%   J2      - Coeficiente zonal de segundo orden del potencial gravitatorio
%   opts    - (Opcional) Estructura con opciones del integrador ODE
%
% OUTPUT:
%   dXdt    - Derivada temporal del estado en coordenadas cartesianas:
%             [dx; dy; dz; ddx; ddy; ddz]

    % Extraer componentes del estado en cartesianas
    x0  = X0(1);
    y0  = X0(2);
    z0  = X0(3);
    dx0 = X0(4);
    dy0 = X0(5);
    dz0 = X0(6);

    % Transformar a coordenadas polares-nodales
    [r0, theta0, nu0, R0, Theta0, N0] = cartesian2PolarNodal(x0, y0, z0, dx0, dy0, dz0);
    X0_PN = [r0, theta0, nu0, R0, Theta0, N0];

    % Calcular dinámica en coordenadas polares-nodales
    if nargin < 7
        dXdt_PN = hodeiMotion(t, X0_PN, mu, Re, alpha, J2);
    else
        dXdt_PN = hodeiMotion(t, X0_PN, mu, Re, alpha, J2, opts);
    end

    % Convertir derivadas a coordenadas cartesianas
    [x, y, z, dx, dy, dz] = polarNodal2Cartesian(dXdt_PN(:, 1), dXdt_PN(:, 2), dXdt_PN(:, 3), ...
                                                 dXdt_PN(:, 4), dXdt_PN(:, 5), dXdt_PN(:, 6));

    % Ensamblar resultado
    dXdt = [x, y, z, dx, dy, dz];
end
