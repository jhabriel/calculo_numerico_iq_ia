% EJ_13_EDO_SEGUNDO_ORDEN_REDUCCION
% Oscilador amortiguado (EDO de segundo orden):
%
%     x'' + 2 zeta w x' + w^2 x = 0,     x(0) = 1,  x'(0) = 0
%
% Reduccion de orden: u1 = x (posicion) y u2 = x' (velocidad). Entonces
%
%     u1' = u2
%     u2' = -2 zeta w u2 - w^2 u1
%
% que es el sistema lineal  u' = AA u  con  u = [u1; u2]  y
%
%     AA = [  0      1        ]
%          [ -w^2  -2 zeta w  ]
%
% Los tres metodos se aplican al vector u (con paso h):
%
%   Euler explicito :  u_{i+1} = u_i + h AA u_i
%   Euler implicito :  (I - h AA) u_{i+1} = u_i      (sistema lineal en cada paso)
%   Heun            :  predictor u* = u_i + h AA u_i,
%                      corrector u_{i+1} = u_i + h/2 [ AA u_i + AA u* ]
%
% La historia se guarda en matrices de (n+1) x 2: una fila por instante,
% la columna 1 es u1 = x y la columna 2 es u2 = x'.
%
% Requisito: haber ejecutado `iniciar` al abrir MATLAB/Octave (ver README).
% Ejecutar:  run('ejemplos/ej_13_edo_segundo_orden_reduccion.m')

clear; clc; close all;

w = 2 * pi;  zeta = 0.1;
AA = [ 0,    1;
      -w^2, -2*zeta*w ];
m = 2;                                  % numero de variables
u0 = [1; 0];                            % u1(0) = x(0) = 1,  u2(0) = x'(0) = 0
t0 = 0; tf = 5;
n = 500;
h = (tf - t0) / n;

wd = w * sqrt(1 - zeta^2);
u1ex = @(t) exp(-zeta*w*t) .* (cos(wd*t) + zeta*w/wd * sin(wd*t));   % x(t) exacta

t = t0 + (0:n)' * h;

%% Euler explicito
UE = zeros(n+1, m);  UE(1, :) = u0';
u = u0;
for i = 1:n
    u = u + h * AA * u;
    UE(i+1, :) = u';
end

%% Euler implicito: (I - h AA) u_{i+1} = u_i
UI = zeros(n+1, m);  UI(1, :) = u0';
u = u0;
for i = 1:n
    u = resolver_lu(eye(m) - h * AA, u);
    UI(i+1, :) = u';
end

%% Heun
UH = zeros(n+1, m);  UH(1, :) = u0';
u = u0;
for i = 1:n
    f1 = AA * u;
    up = u + h * f1;                    % predictor
    f2 = AA * up;
    u = u + h / 2 * (f1 + f2);          % corrector
    UH(i+1, :) = u';
end

%% Resultados
fprintf('Oscilador amortiguado, h = %.3f, error maximo |u1 - u1_exacta|:\n', h);
fprintf('   Euler explicito : %.3e  (la amplitud crece artificialmente)\n', max(abs(UE(:, 1) - u1ex(t))));
fprintf('   Euler implicito : %.3e  (amortigua de mas)\n', max(abs(UI(:, 1) - u1ex(t))));
fprintf('   Heun            : %.3e\n', max(abs(UH(:, 1) - u1ex(t))));

%% Grafica
figure;
plot(t, u1ex(t), 'k-', t, UE(:, 1), 'r--', t, UI(:, 1), 'g--', t, UH(:, 1), 'b-.', 'LineWidth', 1.2);
xlabel('t'); ylabel('u_1(t) = x(t)'); grid on;
legend('Exacta', 'E. explicito', 'E. implicito', 'Heun');
title('x'''' + 2\zeta\omega x'' + \omega^2 x = 0');
