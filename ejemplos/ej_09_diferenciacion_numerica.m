% EJ_09_DIFERENCIACION_NUMERICA
% Primera y segunda derivada de f(x) = exp(x) sin(x) en x0 = 1 con diferencias
% finitas, y estudio del error al reducir el paso h.
%
%   f'(x)  = exp(x) (sin x + cos x)
%   f''(x) = 2 exp(x) cos x
%
% Notacion: h = paso, f_k = f(x0 + k h)  (f_0 = f(x0), f_1 = f(x0+h),
% f_{-1} = f(x0-h), etc.)
%
% PRIMERA DERIVADA
%   adelante  O(h)    f'(x0) ~ ( f_1 - f_0 ) / h
%   adelante  O(h^2)  f'(x0) ~ ( -3 f_0 + 4 f_1 - f_2 ) / (2h)
%   atras     O(h)    f'(x0) ~ ( f_0 - f_{-1} ) / h
%   atras     O(h^2)  f'(x0) ~ ( 3 f_0 - 4 f_{-1} + f_{-2} ) / (2h)
%   central   O(h^2)  f'(x0) ~ ( f_1 - f_{-1} ) / (2h)
%   central   O(h^4)  f'(x0) ~ ( -f_2 + 8 f_1 - 8 f_{-1} + f_{-2} ) / (12h)
%
% SEGUNDA DERIVADA
%   adelante  O(h)    f''(x0) ~ ( f_0 - 2 f_1 + f_2 ) / h^2
%   adelante  O(h^2)  f''(x0) ~ ( 2 f_0 - 5 f_1 + 4 f_2 - f_3 ) / h^2
%   atras     O(h)    f''(x0) ~ ( f_0 - 2 f_{-1} + f_{-2} ) / h^2
%   atras     O(h^2)  f''(x0) ~ ( 2 f_0 - 5 f_{-1} + 4 f_{-2} - f_{-3} ) / h^2
%   central   O(h^2)  f''(x0) ~ ( f_1 - 2 f_0 + f_{-1} ) / h^2
%   central   O(h^4)  f''(x0) ~ ( -f_2 + 16 f_1 - 30 f_0 + 16 f_{-1} - f_{-2} ) / (12 h^2)
%
% Ejecutar:  run('ejemplos/ej_09_diferenciacion_numerica.m')

clear; clc; close all;

f = @(x) exp(x) .* sin(x);
x0 = 1;
d1 = exp(x0) * (sin(x0) + cos(x0));        % derivada exacta
d2 = 2 * exp(x0) * cos(x0);                % segunda derivada exacta

h = 0.4 ./ 2.^(0:5);                       % pasos: 0.4, 0.2, 0.1, ...

% Valores de f en los puntos x0 + k h (h es un vector: se evalua para todos
% los h a la vez)
fm3 = f(x0 - 3*h);  fm2 = f(x0 - 2*h);  fm1 = f(x0 - h);
f0  = f(x0);
f1  = f(x0 + h);    f2  = f(x0 + 2*h);  f3  = f(x0 + 3*h);

%% Primera derivada
d1_ad1 = (f1 - f0) ./ h;                              % adelante O(h)
d1_ad2 = (-3*f0 + 4*f1 - f2) ./ (2*h);                % adelante O(h^2)
d1_at1 = (f0 - fm1) ./ h;                             % atras    O(h)
d1_at2 = (3*f0 - 4*fm1 + fm2) ./ (2*h);               % atras    O(h^2)
d1_ce2 = (f1 - fm1) ./ (2*h);                         % central  O(h^2)
d1_ce4 = (-f2 + 8*f1 - 8*fm1 + fm2) ./ (12*h);        % central  O(h^4)

%% Segunda derivada
d2_ad1 = (f0 - 2*f1 + f2) ./ h.^2;                    % adelante O(h)
d2_ad2 = (2*f0 - 5*f1 + 4*f2 - f3) ./ h.^2;           % adelante O(h^2)
d2_at1 = (f0 - 2*fm1 + fm2) ./ h.^2;                  % atras    O(h)
d2_at2 = (2*f0 - 5*fm1 + 4*fm2 - fm3) ./ h.^2;        % atras    O(h^2)
d2_ce2 = (f1 - 2*f0 + fm1) ./ h.^2;                   % central  O(h^2)
d2_ce4 = (-f2 + 16*f1 - 30*f0 + 16*fm1 - fm2) ./ (12*h.^2);   % central O(h^4)

%% Errores absolutos
E1 = abs([d1_ad1; d1_ad2; d1_at1; d1_at2; d1_ce2; d1_ce4] - d1);   % 6 x length(h)
E2 = abs([d2_ad1; d2_ad2; d2_at1; d2_at2; d2_ce2; d2_ce4] - d2);

nombres = {'adelante O(h)', 'adelante O(h^2)', 'atras O(h)', ...
           'atras O(h^2)', 'central O(h^2)', 'central O(h^4)'};
teorico = [1 2 1 2 2 4];

% Orden observado con los dos h mas pequenos: log2(e(h) / e(h/2))
% (con h grandes el error todavia no sigue su comportamiento asintotico)
orden = @(e) log2(e(end-1) / e(end));

%% Resultados
fprintf('f''(1)  exacta = %.8f\n', d1);
fprintf('f''''(1) exacta = %.8f\n\n', d2);

fprintf('Primera derivada, errores absolutos:\n');
fprintf('     h      ad O(h)   ad O(h^2)   at O(h)   at O(h^2)  cen O(h^2)  cen O(h^4)\n');
fprintf('%8.5f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', [h; E1]);

fprintf('\nSegunda derivada, errores absolutos:\n');
fprintf('     h      ad O(h)   ad O(h^2)   at O(h)   at O(h^2)  cen O(h^2)  cen O(h^4)\n');
fprintf('%8.5f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', [h; E2]);

fprintf('\nOrden observado (teorico entre parentesis):\n');
for j = 1:6
    fprintf('  %-16s  f'': %5.2f (%d)    f'''': %5.2f (%d)\n', nombres{j}, ...
            orden(E1(j, :)), teorico(j), orden(E2(j, :)), teorico(j));
end

%% Graficas
figure;
loglog(h, E1(1,:), 'o-', h, E1(2,:), 'o--', h, E1(3,:), 's-', h, E1(4,:), 's--', ...
       h, E1(5,:), '^-', h, E1(6,:), '^--', 'LineWidth', 1.3);
xlabel('h'); ylabel('error absoluto'); grid on;
legend(nombres, 'Location', 'southeast');
title('Primera derivada: error de las diferencias finitas');

figure;
loglog(h, E2(1,:), 'o-', h, E2(2,:), 'o--', h, E2(3,:), 's-', h, E2(4,:), 's--', ...
       h, E2(5,:), '^-', h, E2(6,:), '^--', 'LineWidth', 1.3);
xlabel('h'); ylabel('error absoluto'); grid on;
legend(nombres, 'Location', 'southeast');
title('Segunda derivada: error de las diferencias finitas');
