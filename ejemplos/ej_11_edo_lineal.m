% EJ_11_EDO_LINEAL
% Reaccion de primer orden (EDO lineal):   dC/dt = -k C,   C(0) = C0
% Solucion analitica:                      C(t) = C0 exp(-k t)
%
% Se resuelve con tres metodos de un paso, con paso h y t_{i+1} = t_i + h:
%
%   Euler explicito :  C_{i+1} = C_i + h f(t_i, C_i)
%   Euler implicito :  C_{i+1} = C_i + h f(t_{i+1}, C_{i+1})
%   Heun            :  predictor  C* = C_i + h f(t_i, C_i)
%                      corrector  C_{i+1} = C_i + h/2 [ f(t_i, C_i) + f(t_{i+1}, C*) ]
%
% En Euler implicito la incognita C_{i+1} aparece a ambos lados. Como f es
% lineal se puede despejar:
%
%   C_{i+1} = C_i + h (-k C_{i+1})   =>   C_{i+1} = C_i / (1 + h k)
%
% La historia de cada metodo se guarda en un vector de n+1 componentes.
%
% Ejecutar:  run('ejemplos/ej_11_edo_lineal.m')

clear; clc; close all;

k = 0.5;             % 1/s
C0 = 2;              % mol/L
t0 = 0; tf = 10;     % s
n = 20;              % numero de pasos
h = (tf - t0) / n;

f = @(t, C) -k * C;                  % dC/dt = f(t, C)
Cex = @(t) C0 * exp(-k * t);         % solucion analitica

t = t0 + (0:n)' * h;                 % instantes de tiempo (columna)
Ce = zeros(n+1, 1);                  % Euler explicito
Ci = zeros(n+1, 1);                  % Euler implicito
Ch = zeros(n+1, 1);                  % Heun
Ce(1) = C0; Ci(1) = C0; Ch(1) = C0;

for i = 1:n
    % Euler explicito
    Ce(i+1) = Ce(i) + h * f(t(i), Ce(i));

    % Euler implicito (despeje para f lineal)
    Ci(i+1) = Ci(i) / (1 + h * k);

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
title('A \rightarrow productos (1er orden), n = 20');
