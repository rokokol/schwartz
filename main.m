clc; clear; clf; close all;
% NOTE: Я воспринимаю это как УЧЕБУ и РАБОТУ

% NOTE: На параметры накладываются следующие ограничения при симметричном методе:
% n >= 5 && n % 2 == 1
% intersect_points_radius >= 2

% NOTE: Вопросы:
% - Получается, что intersect_points_radius=1 невозможна: на рисунке видно, что информации для корректировки точек будет слишком мало. При intersect_points_radius=2 имеем ситуацию, в которой доступен минимальный возможный "зазор" для оценки погрешности и, соответственно, коррективки краевых условий
% - Как можно было бы реализовать метод для четных n? Уместно ли это с точки зрения памяти и скорости? Кажется, что эффективнее при четном n просто отбрасывать лишнюю точку

%% init
n = 7; a = 0; b = 1 * pi; u_0 = 0; u_n = 0; h = (b - a) / n;
intersect_points_radius = 2; accuracy = 1e-2;
x = linspace(a, b, n); f = sin(x)';

A = diag(ones(n - 3, 1), -1) + diag(ones(n - 3, 1), 1) - 2 * diag(ones(n - 2, 1));
u_exact = A\(f(2:end-1) .* h^2);
u_exact = [u_0 u_exact' u_n];

%% calculations
subplot(1, 3, 2); hold on; grid on;
[U, E] = schwartz_2d_symmetrical(n, a, b, intersect_points_radius, f, u_0, u_n, accuracy, true);

mae = E(:, end);
mae_diff = max(abs(U(:, 2:end-1) - u_exact(2:end - 1)), [], 2);

[val_min, idx_min] = min(mae);
[val_max, idx_max] = max(mae);

subplot(1, 3, 1); hold on; grid on;
plot(x, u_exact, x, U(end, :))
legend(sprintf("Прогонка на %i узлах", n), sprintf("Шварц на %i узлах", n))

fprintf("%f => %f, min=%f (idx=%d), max=%f (idx=%d)\n", ...
  mae(1), mae(end), val_min, idx_min, val_max, idx_max);

[val_min, idx_min] = min(mae_diff);
[val_max, idx_max] = max(mae_diff);

fprintf("%f => %f, min=%f (idx=%d), max=%f (idx=%d)\n", ...
  mae_diff(1), mae_diff(end), val_min, idx_min, val_max, idx_max);

subplot(1, 3, 3); hold on; grid on;
plot(E)
plot(mae_diff)
plot(idx_min, val_min, "o")

legend("left error", "right error", "MAE", "real error", sprintf("Minimum real error at %i iteration", idx_min))

