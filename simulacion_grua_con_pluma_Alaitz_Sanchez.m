
clear; clc; close all;

% 1. DEFINICIÓN DE PARÁMETROS DEL SISTEMA
M = 1000;      % Masa de la pluma (kg)
L = 20;        % Longitud de la pluma (m)
m = 50;       % Masa del gancho (kg)
l = 5;         % Longitud del cable (m)
H = 30;        % Altura inicial del CM (m)
g = 9.81;      % Aceleración de la gravedad (m/s^2)

ICM = (1/12) * M * L^2; % Momento de inercia

% Condiciones Iniciales
x0      = 0;            
y0      = H;            
theta0  = deg2rad(10);  
phi0    = deg2rad(50);  

dx0     = 0;            
dy0     = 0;            
dtheta0 = 0;            
dphi0   = 0.2;          % Impulso inicial para ver la oscilación suave

q0 = [x0; y0; theta0; phi0; dx0; dy0; dtheta0; dphi0];

% 2. RESOLUCIÓN NUMÉRICA CON EVENTO DE IMPACTO EN EL SUELO
% Definimos un tiempo máximo amplio, pero ode45 se detendrá cuando Y_CM llegue a un valor seguro
tspan = [0, 10]; 

% Opciones de ode45 configuradas con un Evento (detección de suelo Y = 0)
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'Events', @ground_event);

ode_sys = @(t, q) [
    q(5); ... 
    q(6); ... 
    q(7); ... 
    q(8); ... 
    [ (M+m),       0,       m*l*cos(q(4)); ...
        0,       (M+m),     m*l*sin(q(4)); ...
      cos(q(4)), sin(q(4)),      l         ] \ ...
    [ m*l*sin(q(4))*(q(8)^2); ...
     -m*l*cos(q(4))*(q(8)^2) - (M+m)*g; ...
     -g*sin(q(4)) ]; ...
    0 
];

[t, sol, te, ye, ie] = ode45(ode_sys, tspan, q0, options);

% Extracción de resultados
x     = sol(:, 1);  y_pos = sol(:, 2);
theta = sol(:, 3);  phi   = sol(:, 4);
dx    = sol(:, 5);  dy    = sol(:, 6);
dtheta= sol(:, 7);  dphi  = sol(:, 8);

xh = x + l * sin(phi);
yh = y_pos - l * cos(phi);

Px = (M+m)*dx + m*l.*dphi.*cos(phi); 
E  = 0.5*(M+m)*(dx.^2 + dy.^2) + m*l.*dphi.*(dx.*cos(phi) + dy.*sin(phi)) + ...
     0.5*m*l^2.*dphi.^2 + 0.5*ICM*dtheta.^2 + (M+m)*g.*y_pos - m*g*l.*cos(phi); 

% 3. ANIMACIÓN Y VISUALIZACIÓN
fig = figure('Name', 'Caída Libre con Parada Automática en el Suelo', ...
             'Position', [50, 50, 1200, 650], 'Color', 'white');

for k = 1:length(t)
    %Animación 
    sp1 = subplot(2, 2, [1, 3]);
    clf(fig);
    
    sp1 = subplot(2, 2, [1, 3]); hold(sp1, 'on'); grid(sp1, 'on'); axis(sp1, 'equal');
    axis(sp1, [-15, 25, -5, 35]);
    xlabel(sp1, 'Posición Horizontal X (m)', 'FontWeight', 'bold');
    ylabel(sp1, 'Posición Vertical Y (m)', 'FontWeight', 'bold');
    title(sp1, sprintf('Caída Libre (Impacto al suelo en t = %.2f s)', t(end)), 'FontSize', 12);
    
    x1 = x(k) - (L/2)*cos(theta(k)); y1 = y_pos(k) - (L/2)*sin(theta(k));
    x2 = x(k) + (L/2)*cos(theta(k)); y2 = y_pos(k) + (L/2)*sin(theta(k));
    
    plot(sp1, [0, 0], [0, H], 'k--', 'LineWidth', 1.5); 
    line(sp1, [-15, 25], [0, 0], 'Color', 'k', 'LineWidth', 2); 
    
    plot(sp1, x(1:k), y_pos(1:k), 'r:', 'LineWidth', 1.2); 
    plot(sp1, xh(1:k), yh(1:k), 'b--', 'LineWidth', 1);    
    
    plot(sp1, [x1, x2], [y1, y2], 'Color', [0.85 0.55 0.1], 'LineWidth', 6);
    plot(sp1, x(k), y_pos(k), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    plot(sp1, [x(k), xh(k)], [y_pos(k), yh(k)], 'k-', 'LineWidth', 2);
    plot(sp1, xh(k), yh(k), 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 7);
    
    legend(sp1, {'Eje Torre', 'Suelo', 'Trayectoria CM', 'Trayectoria Gancho', ...
                 'Pluma (M, L)', 'CM (L/2)', 'Gancho (m, l)'}, 'Location', 'northeast');

    % Ángulos 
    sp2 = subplot(2, 2, 2);
    plot(sp2, t(1:k), rad2deg(phi(1:k)), 'b-', 'LineWidth', 1.5); hold(sp2, 'on'); grid(sp2, 'on');
    plot(sp2, t(1:k), rad2deg(theta(1:k)), 'r-', 'LineWidth', 1.5);
    xlabel(sp2, 'Tiempo (s)'); ylabel(sp2, 'Ángulo (grados)');
    title(sp2, 'Evolución Ángulos \phi y \theta');
    legend(sp2, {'\phi (Gancho)', '\theta (Pluma)'}, 'Location', 'southeast');
    xlim(sp2, [0, t(end)]);

    %Conservación 
    sp3 = subplot(2, 2, 4);
    plot(sp3, t(1:k), Px(1:k), 'm-', 'LineWidth', 1.5); hold(sp3, 'on'); grid(sp3, 'on');
    plot(sp3, t(1:k), (E(1:k) - E(1))/1000, 'g-', 'LineWidth', 1.5);
    xlabel(sp3, 'Tiempo (s)'); ylabel(sp3, 'Magnitud');
    title(sp3, 'Verificación de Constantes');
    legend(sp3, {'Momento P_x', '\DeltaEnergía E'}, 'Location', 'southeast');
    xlim(sp3, [0, t(end)]);

    drawnow;
end

% FUNCIÓN DE EVENTO PARA DETECTAR EL SUELO
function [value, isterminal, direction] = ground_event(~, q)
    % El evento salta cuando el Centro de Masas de la pluma llega a Y = 0 
    y_cm = q(2);
    
    value = y_cm;      % Queremos que este valor llegue a 0
    isterminal = 1;    % 1 = Detener la integración inmediatamente
    direction = -1;    % -1 = Solo cuando va cayendo (dirección negativa)
end