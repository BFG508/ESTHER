function dcmPQW2I = getDcmPQW2ECI(inclination, RAAN, argumentOfPeriapsis)
% getDcmPQW2I Calcula la Matriz de Cosenos Directores (DCM) del marco Perifocal (PQW) al marco Inercial (I).
%
%   Esta función calcula la matriz de transformación que convierte un vector del 
%   marco de coordenadas perifocal (PQW) al marco inercial (ECI) usando una secuencia de rotaciones pasivas.
%
%   ENTRADAS:
%       - inclination          : Inclinación orbital [rad]
%       - RAAN                 : Ascensión Recta del Nodo Ascendente (Ω) [rad]
%       - argumentOfPeriapsis  : Argumento del Periapside (ω) [rad]
%
%   SALIDA:
%       - dcmPQW2I : Matriz 3×3 de cosenos directores (DCM) que transforma un vector de PQW a ECI.
%
%   SECUENCIA DE ROTACIONES:
%       - 1ª rotación: Rotación pasiva alrededor del eje Z por el Argumento del Periapside (ω).
%       - 2ª rotación: Rotación pasiva alrededor del eje X por la Inclinación (i).
%       - 3ª rotación: Rotación pasiva alrededor del eje Z por la RAAN (Ω).

% Matriz de rotación alrededor del eje Z por la RAAN (Ω)
R_RAAN = [cos(RAAN), -sin(RAAN), 0;
          sin(RAAN),  cos(RAAN), 0;
                  0,          0, 1];

% Matriz de rotación alrededor del eje X por la inclinación (i)
R_INC = [1,                0,                 0;
         0, cos(inclination), -sin(inclination);
         0, sin(inclination),  cos(inclination)];

% Matriz de rotación alrededor del eje Z por el argumento del periapside (ω)
R_u = [cos(argumentOfPeriapsis), -sin(argumentOfPeriapsis), 0;
       sin(argumentOfPeriapsis),  cos(argumentOfPeriapsis), 0;
                              0,                         0, 1];

% Producto de las matrices de rotación: de ECI a PQW
dcmI2PQW = R_RAAN * R_INC * R_u;

% La transpuesta da la transformación de PQW a ECI
dcmPQW2I = dcmI2PQW';

end
