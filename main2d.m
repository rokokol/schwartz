clc; clear; clf; close all;

n = 50; m = 50;

x = linspace(0, 1, n); y = linspace(0, 1, m);
u_left  = zeros(m,1);       % u(0,y)
u_right = sin(pi * y)';     % u(1,y)
u_bottom = zeros(n,1);      % u(x,0)
u_top    = zeros(n,1);      % u(x,1)

% --- Координаты и правая часть (без изменений) ---
[Y, X] = meshgrid(y, x);
F = sin(pi * X) .* sin(pi * Y);
F_internal = F(2:end-1, 2:end-1);
[~, ~, U] = solve2d(x, y, u_left, u_bottom, u_right, u_top, F_internal);
U_real =  ((sinh(pi * X) ./ sinh(pi)) - (1 / (2 * pi^2)) * sin(pi * X)) .* sin(pi * Y);

figure(1); grid on
surf(X, Y, U)
title('Численное решение');
xlabel('x'); ylabel('y')

figure(2); grid on
surf(X, Y, U_real)
title('Аналитическое решение');
xlabel('x'); ylabel('y')

