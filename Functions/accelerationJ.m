function acceleration = accelerationJ(position, mu, R, J2, J3, J4, J5, J6)
%% accelerationJ Aceleración perturbativa por armónicos zonales J2–J6 en un campo gravitacional.
%
% Esta función calcula la aceleración perturbativa total sobre un satélite debido a los
% coeficientes armónicos zonales J2 hasta J6 del potencial gravitacional, aplicando las
% fórmulas de Cunningham y Escobal.
%
% Sintaxis:
%   acceleration = accelerationJ(position, mu, R, J2, J3, J4, J5, J6)
%
% Entradas:
%   position : Vector de posición [x; y; z] en coordenadas cartesianas [m o km]
%   mu       : Parámetro gravitacional del cuerpo central (G·M) [m³/s² o km³/s²]
%   R        : Radio ecuatorial del cuerpo central [m o km]
%   J2–J6    : Coeficientes zonales del desarrollo armónico del potencial (sin unidades)
%
% Salida:
%   acceleration : Vector [ax; ay; az] con la aceleración total debido a J2–J6 [m/s² o km/s²]
%
% Referencias:
%   - Cunningham, G.F. (1970). "On the Computation of the Spherical Harmonic Terms Needed
%     During the Numerical Integration of the Orbital Motion of an Artificial Satellite"
%   - Escobal, P.R. (1976). "Methods of Orbit Determination"

    %% --- Descomposición de posición ---
    x = position(1);
    y = position(2);
    z = position(3);
    r = norm(position);  % Magnitud del vector posición

    %% --- Contribución de J2 (Cunningham) ---
    J2factor = - (3/2) * J2 * mu * R^2 / r^5;
    aI_J2 = J2factor * x * (1 - 5*(z/r)^2);
    aJ_J2 = J2factor * y * (1 - 5*(z/r)^2);
    aK_J2 = J2factor * z * (3 - 5*(z/r)^2);

    %% --- Contribución de J3 (Cunningham) ---
    J3factor = - (5/2) * J3 * mu * R^3 / r^7;
    aI_J3 = J3factor * x * z * (3 - 7*(z/r)^2);
    aJ_J3 = J3factor * y * z * (3 - 7*(z/r)^2);
    aK_J3 = J3factor * (6*z^2 - 7*(z^2/r)^2 - (3/5)*r^2);

    %% --- Contribución de J4 (Cunningham) ---
    J4factor = (15/8) * J4 * mu * R^4 / r^7;
    aI_J4 = J4factor * x * (1 - 14*(z/r)^2 + 21*(z/r)^4);
    aJ_J4 = J4factor * y * (1 - 14*(z/r)^2 + 21*(z/r)^4);
    aK_J4 = J4factor * z * (5 - (70/3)*(z/r)^2 + 21*(z/r)^4);

    %% --- Contribución de J5 (Escobal) ---
    J5factor = (3/8) * J5 * mu * R^5 / r^9;
    aI_J5 = J5factor * x * z * (35 - 210*(z/r)^2 + 231*(z/r)^4);
    aJ_J5 = J5factor * y * z * (35 - 210*(z/r)^2 + 231*(z/r)^4);
    aK_J5 = J5factor * (z^2 * (105 - 315*(z/r)^2 + 231*(z/r)^4) - 5*r^2);

    %% --- Contribución de J6 (Escobal) ---
    J6factor = - (1/16) * J6 * mu * R^6 / r^9;
    aI_J6 = J6factor * x * (35 - 945*(z/r)^2 + 3465*(z/r)^4 - 3003*(z/r)^6);
    aJ_J6 = J6factor * y * (35 - 945*(z/r)^2 + 3465*(z/r)^4 - 3003*(z/r)^6);
    aK_J6 = J6factor * z * (245 - 2205*(z/r)^2 + 4851*(z/r)^4 - 3003*(z/r)^6);

    %% --- Aceleración total perturbada ---
    acceleration = [aI_J2 + aI_J3 + aI_J4 + aI_J5 + aI_J6;
                    aJ_J2 + aJ_J3 + aJ_J4 + aJ_J5 + aJ_J6;
                    aK_J2 + aK_J3 + aK_J4 + aK_J5 + aK_J6];
end
