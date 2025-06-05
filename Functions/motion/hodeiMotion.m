function dXdt = hodeiMotion(tspan, X0, mu, Re, alpha, J2, optiset)
%% hodeiMotion Simula el movimiento orbital de un satélite bajo la acción combinada de
% un empuje radial constante y la perturbación zonal J2, utilizando un 
% modelo de primer orden basado en una transformación canónica.
%
% La integración se realiza en un conjunto de variables primadas que
% eliminan los efectos oscilatorios de corto período debidos a J2, y se 
% transforma posteriormente de nuevo a variables polares-nodales reales.
%
% INPUTS:
%   tspan   - Vector [t0 tf] o array de tiempos para la integración [s]
%   X0      - Vector de estado inicial en variables polares-nodales:
%             [r; theta; nu; R; Theta; N]
%   mu      - Parámetro gravitacional estándar del cuerpo central [m^3/s^2]
%   Re      - Radio ecuatorial del cuerpo central [m]
%   alpha   - Aceleración constante de bajo empuje en dirección radial [m/s^2]
%   J2      - Coeficiente zonal de segundo orden del potencial gravitatorio
%   optiset - (Opcional) Estructura con opciones del integrador ODE (ode113)
%
% OUTPUT:
%   dXdt    - Matriz Nx6 con la evolución temporal del estado en variables:
%             [r, theta, nu, R, Theta, N]

    % === 1. Estado inicial en variables polares-nodales ===
    r     = X0(1); 
    theta = X0(2); 
    nu    = X0(3);
    R     = X0(4); 
    Theta = X0(5); 
    N     = X0(6);

    % === 2. Cálculo de parámetros intermedios ===
    p     = Theta^2 / mu;
    sigma = 1/2 * J2 * (Re / p)^2;
    k     = p / r - 1;
    q     = p * R / Theta;
    c     = N / Theta;
    s     = sqrt(1 - c^2);

    % === 3. Transformación a variables primadas ===
    r_p     = r     + sigma * p * (1 - 3/2*s^2 - 1/2*s^2*cos(2*theta));
    theta_p = theta + sigma * ((1 - 6*c^2)*q + (1 - 2*c^2)*q*cos(2*theta) ...
                            - (1/4 + k - (7/4 + 3*k)*c^2)*sin(2*theta));
    nu_p    = nu    + sigma * c * ((3 + cos(2*theta))*q - (3/2 + 2*k)*sin(2*theta));
    R_p     = R     + sigma * Theta / p * (1 + k)^2 * s^2 * sin(2*theta);
    Theta_p = Theta - sigma * Theta * s^2 * ((3/2 + 2*k)*cos(2*theta) + q * sin(2*theta));
    N_p     = N     - 0;

    % === 4. Integración en espacio primado ===
    p_p         = Theta_p^2 / mu;
    sigma_p     = 1/2 * J2 * (Re / p_p)^2;
    c_p         = N_p / Theta_p;
    s_p         = sqrt(1 - c_p^2);
    Theta_tilde = Theta_p * sqrt(1 - sigma_p*(2 - 3*s_p^2));

    X0_p = [r_p; R_p; theta_p; nu_p];
    if nargin < 7
        [t, dXdt_p] = ode113(@(t, X) [X(2); ...
                                      Theta_tilde^2 / X(1)^3 - mu / X(1)^2 + alpha; ...
                                            Theta_p / X(1)^2 * (1 - sigma_p * (1 - 6*c_p^2)); ...
                                      - 3 * Theta_p / X(1)^2 * sigma_p * c_p], ...
                             tspan, X0_p);
    else
        [t, dXdt_p] = ode113(@(t, X) [X(2); ...
                                      Theta_tilde^2 / X(1)^3 - mu / X(1)^2 + alpha; ...
                                            Theta_p / X(1)^2 * (1 - sigma_p * (1 - 6*c_p^2)); ...
                                      - 3 * Theta_p / X(1)^2 * sigma_p * c_p], ...
                             tspan, X0_p, optiset);
    end
    r_p_vec     = dXdt_p(:,1);
    R_p_vec     = dXdt_p(:,2);
    theta_p_vec = dXdt_p(:,3);
    nu_p_vec    = dXdt_p(:,4);

    % === 5. Transformación inversa a variables originales ===
    dXdt = zeros(length(t), 6);
    for j = 1:length(t)
        r_pj     = r_p_vec(j);
        theta_pj = theta_p_vec(j);
        nu_pj    = nu_p_vec(j);
        R_pj     = R_p_vec(j);
        
        k_p = p_p / r_pj - 1;
        q_p = p_p * R_pj / Theta_p;
        
        r     = r_pj      - sigma_p * p_p * (1 - 3/2*s_p^2 - 1/2*s_p^2*cos(2*theta_pj));
        theta = theta_pj  - sigma_p * ((1 - 6*c_p^2)*q_p + (1 - 2*c_p^2)*q_p*cos(2*theta_pj) ...
                 - (1/4 + k_p - (7/4 + 3*k_p)*c_p^2)*sin(2*theta_pj));
        nu    = nu_pj     - sigma_p * c_p * ((3 + cos(2*theta_pj))*q_p - (3/2 + 2*k_p)*sin(2*theta_pj));
        R     = R_pj      - sigma_p * Theta_p/p_p * (1 + k_p)^2 * s_p^2 * sin(2*theta_pj);
        Theta = Theta_p + sigma_p * Theta_p * s_p^2 * ((3/2 + 2*k_p)*cos(2*theta_pj) + q_p * sin(2*theta_pj));
        N     = N_p     + 0;

        dXdt(j,:) = [r, theta, nu, R, Theta, N];
    end
end
