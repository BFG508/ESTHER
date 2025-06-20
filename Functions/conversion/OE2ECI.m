function [posECI, velECI] = OE2ECI(mu, SMA, ECC, INC, RAAN, AOP, TA)
% OE2ECI Convierte elementos orbitales clásicos en vectores de posición y velocidad en el marco ECI.
%
%   Esta función toma los elementos orbitales clásicos (COEs) de una nave 
%   y los convierte en vectores de posición y velocidad en el sistema de coordenadas 
%   Earth-Centered Inertial (ECI).
%
%   ENTRADAS:
%       - SMA  : Semieje mayor [m]
%       - ECC  : Excentricidad (adimensional)
%       - INC  : Inclinación [rad]
%       - RAAN : Ascensión Recta del Nodo Ascendente [rad]
%       - AOP  : Argumento del Periápside [rad]
%       - TA   : Anomalía Verdadera [rad]
%
%   SALIDAS:
%       - posECI : Vector de posición en el marco ECI [m]
%       - velECI : Vector de velocidad en el marco ECI [m/s]

% Calcular el semi-latus rectum (parámetro de la sección cónica)
semilatusRectum = SMA * (1 - ECC^2);

% Calcular el momento angular específico
angularMomentum = sqrt(mu * semilatusRectum);

% Calcular el vector de posición en el marco perifocal
perifocalPos = (semilatusRectum / (1 + ECC * cos(TA))) * [cos(TA); sin(TA); 0];

% Calcular el vector de velocidad en el marco perifocal
perifocalVel = (mu / angularMomentum) * [ECC*sin(TA); 1 + ECC*cos(TA); 0];

% Calcular la matriz de transformación del marco perifocal al ECI
dcmPQW2ECI = getDcmPQW2ECI(INC, RAAN, AOP);

% Transformar posición y velocidad del marco perifocal al ECI
posECI = dcmPQW2ECI' * perifocalPos;
velECI = dcmPQW2ECI' * perifocalVel;

end
