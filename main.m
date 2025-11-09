clc; clear; clf; close all;

%% init
n = 13; a = 0; b = 1 * pi; u_0 = 0; u_n = 0;
intersect_points = 3; accuracy = 1e-2;
x = linspace(a, b, n); f = sin(x)';

A = diag(ones(n - 3, 1), -1) + diag(ones(n - 3, 1), 1) - 2 * diag(ones(n - 2, 1));
u_real = A\f(2:end-1);

subplot(1, 3, 1); hold on; grid on;
plot(x, [u_0 u_real' u_n])
legend(sprintf("Прогонка на %i узлах", n))

%% calculations
subplot(1, 3, 2); hold on; grid on;
[u_approx, E] = schwartz_2d(n, a, b, intersect_points, f, u_0, u_n, accuracy, true);

subplot(1, 3, 3); hold on; grid on;
plot(E')
legend("left error", "right error", "MAE")

