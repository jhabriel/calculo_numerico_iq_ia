% EJ_14_SISTEMA_EDO_LINEAL
% Modelo compartimental (farmacocinetica de dos compartimentos con absorcion
% oral y eliminacion):
%
%   dA/dt  = -ka A                            (tracto gastrointestinal)
%   dC1/dt =  ka A - (k10 + k12) C1 + k21 C2  (sangre)
%   dC2/dt =  k12 C1 - k21 C2                 (tejido)
%   dE/dt  =  k10 C1                          (eliminado)
%
% La masa total A + C1 + C2 + E se conserva. Con y = [A; C1; C2; E] es un
% sistema lineal  y' = K y,  donde K es la matriz de coeficientes.
%
% Los tres metodos se aplican al vector y (con paso h):
%
%   Euler explicito :  y_{i+1} = y_i + h K y_i
%   Euler implicito :  (I - h K) y_{i+1} = y_i       (sistema lineal en cada paso)
%   Heun            :  predictor y* = y_i + h K y_i,
%                      corrector y_{i+1} = y_i + h/2 [ K y_i + K y* ]
%
% La historia se guarda en matrices de (n+1) x 4: una fila por instante y
% una columna por variable (A, C1, C2, E).
%
% Requisito: haber ejecutado `iniciar` al abrir MATLAB/Octave (ver README).
% Ejecutar:  run('ejemplos/ej_14_sistema_edo_lineal.m')

clear; clc; close all;

ka = 1.2; k10 = 0.25; k12 = 0.4; k21 = 0.2;      % 1/h
K = [-ka,        0,          0,    0;
      ka, -(k10 + k12),     k21,   0;
       0,       k12,       -k21,   0;
       0,       k10,         0,    0];
m = 4;                                           % numero de variables
y0 = [100; 0; 0; 0];                             % dosis de 100 mg en A
t0 = 0; tf = 24;                                 % h
n = 240;
h = (tf - t0) / n;

t = t0 + (0:n)' * h;

%% Euler explicito
ZE = zeros(n+1, m);  ZE(1, :) = y0';
y = y0;
for i = 1:n
    y = y + h * K * y;
    ZE(i+1, :) = y';
end

%% Euler implicito: (I - h K) y_{i+1} = y_i
ZI = zeros(n+1, m);  ZI(1, :) = y0';
y = y0;
for i = 1:n
    y = resolver_lu(eye(m) - h * K, y);
    ZI(i+1, :) = y';
end

%% Heun
ZH = zeros(n+1, m);  ZH(1, :) = y0';
y = y0;
for i = 1:n
    f1 = K * y;
    yp = y + h * f1;                    % predictor
    f2 = K * yp;
    y = y + h / 2 * (f1 + f2);          % corrector
    ZH(i+1, :) = y';
end

%% Referencia: solucion exacta y(t) = expm(K t) y0
Zex = zeros(n+1, m);
for i = 1:n+1
    Zex(i, :) = (expm(K * t(i)) * y0)';
end

%% Resultados
fprintf('Modelo compartimental, h = %.2f h, error maximo en C1:\n', h);
fprintf('   Euler explicito : %.3e\n', max(abs(ZE(:, 2) - Zex(:, 2))));
fprintf('   Euler implicito : %.3e\n', max(abs(ZI(:, 2) - Zex(:, 2))));
fprintf('   Heun            : %.3e\n', max(abs(ZH(:, 2) - Zex(:, 2))));
fprintf('   Masa total en t = %g h (Heun): %.6f mg (inicial 100)\n', tf, sum(ZH(end, :)));

[Cmax, imax] = max(ZH(:, 2));
fprintf('   Pico en sangre (Heun): %.2f mg a las %.1f h\n', Cmax, t(imax));

%% Grafica
figure;
plot(t, Zex(:, 2), 'k-', t, ZH(:, 2), 'b--', t, Zex(:, 3), 'k-', t, ZH(:, 3), 'r--', 'LineWidth', 1.2);
xlabel('t (h)'); ylabel('cantidad (mg)'); grid on;
legend('Exacta', 'Heun C_1 (sangre)', '', 'Heun C_2 (tejido)');
title('Modelo de dos compartimentos');
