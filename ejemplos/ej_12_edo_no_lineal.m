% EJ_12_EDO_NO_LINEAL
% Reaccion de segundo orden (EDO no lineal):   dC/dt = -k C^2,   C(0) = C0
% Solucion analitica:                          C(t) = C0 / (1 + k C0 t)
%
% Mismos tres metodos que en ej_11. La diferencia esta en Euler implicito: como
% f es no lineal en C, la ecuacion
%
%   C_{i+1} = C_i + h f(t_{i+1}, C_{i+1})
%
% no se puede despejar. En cada paso hay que hallar la raiz de
%
%   G(z) = z - C_i - h f(t_{i+1}, z) = 0,     G'(z) = 1 - h df/dC(t_{i+1}, z)
%
% con el metodo de Newton (funcion newton), usando la derivada analitica df/dC.
%
% Requisito: haber ejecutado `iniciar` al abrir MATLAB/Octave (ver README).
% Ejecutar:  run('ejemplos/ej_12_edo_no_lineal.m')

clear; clc; close all;

k = 0.5;             % L/(mol s)
C0 = 2;              % mol/L
t0 = 0; tf = 10;     % s
n = 20;              % numero de pasos
h = (tf - t0) / n;

f  = @(t, C) -k * C.^2;              % dC/dt = f(t, C)
df = @(t, C) -2 * k * C;             % derivada de f respecto de C
Cex = @(t) C0 ./ (1 + k * C0 * t);   % solucion analitica

t = t0 + (0:n)' * h;                 % instantes de tiempo (columna)
Ce = zeros(n+1, 1);                  % Euler explicito
Ci = zeros(n+1, 1);                  % Euler implicito
Ch = zeros(n+1, 1);                  % Heun
Ce(1) = C0; Ci(1) = C0; Ch(1) = C0;

for i = 1:n
    % Euler explicito
    Ce(i+1) = Ce(i) + h * f(t(i), Ce(i));

    % Euler implicito: raiz de G(z) = 0 con Newton
    G  = @(z) z - Ci(i) - h * f(t(i+1), z);
    dG = @(z) 1 - h * df(t(i+1), z);
    Ci(i+1) = newton(G, dG, Ci(i), 1e-12);

    % Heun (predictor-corrector)
    f1 = f(t(i), Ch(i));
    Cp = Ch(i) + h * f1;                         % predictor
    f2 = f(t(i+1), Cp);
    Ch(i+1) = Ch(i) + h / 2 * (f1 + f2);         % corrector
end

%% Resultados
fprintf('n = %d pasos (h = %.2f s)\n', n, h);
fprintf('  t     exacta   E.explicito  E.implicito   Heun\n');
idx = 1:4:n+1;
fprintf('%5.1f  %8.5f  %9.5f  %9.5f  %9.5f\n', [t(idx), Cex(t(idx)), Ce(idx), Ci(idx), Ch(idx)]');

fprintf('\nError maximo |C - C_exacta|:\n');
fprintf('  Euler explicito : %.3e\n', max(abs(Ce - Cex(t))));
fprintf('  Euler implicito : %.3e\n', max(abs(Ci - Cex(t))));
fprintf('  Heun            : %.3e\n', max(abs(Ch - Cex(t))));

%% Grafica
figure;
tt = linspace(t0, tf, 200);
plot(tt, Cex(tt), 'k-', t, Ce, 'ro--', t, Ci, 'gs--', t, Ch, 'b^--', 'LineWidth', 1.2);
xlabel('t (s)'); ylabel('C (mol/L)'); grid on;
legend('Exacta', 'E. explicito', 'E. implicito', 'Heun');
title('2A \rightarrow productos (2do orden), n = 20');
