clc; clear; clf; close all;

n = 5; m = 5;
n_internal_x = n - 2;
n_internal_y = m - 2;
n_internal = n_internal_x * n_internal_y;

x = linspace(0, 1, n); y = linspace(0, 1, m);
hx = x(2) - x(1);
hy = y(2) - y(1);

% --- Координаты и правая часть (без изменений) ---
[X_internal, Y_internal] = meshgrid(x(2:end-1), y(2:end-1));
[Y, X] = meshgrid(y, x);
F = (sin(pi * X_internal) .* sin(pi * Y_internal));

e = ones(n_internal, 1);
diagonals = [hy^2*e, hx^2*e, -2*(hx^2+hy^2)*e, hx^2*e, hy^2*e];
offsets = [-n_internal_y, -1, 0, 1, n_internal_y];
A = spdiags(diagonals, offsets, n_internal, n_internal);
% Убрать ложные связи между столбцами
b = hx^2*e;
b(n_internal_y:n_internal_y:n_internal) = 0;
A = spdiags(b, 1, A); % для +1 диагонали
A = spdiags(b, -1, A);% для -1 диагонали

u_internal = A \ (F(:) * hx^2*hy^2);
U_internal = reshape(u_internal, m - 2, n - 2)';
U = zeros(n, m);
U(2:end-1, 2:end-1) = U_internal;
U_real = -1/(2*pi^2) * sin(pi*X) .* sin(pi*Y);

figure(1); grid on
surf(X, Y, U)

figure(2); grid on
surf(X, Y, U_real)

