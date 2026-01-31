solve1d(101, 0, pi, 0, 0, 1e-3)

function = solve1d(n, a, b, u_0, u_n, accuracy, target_ax)
% n >= 5 && n % 2 == 1
% intersect_points_radius >= 2

h = (b - a) / n;
intersect_points_radius = floor(n / 4); 
x = linspace(a, b, n); f = sin(x)';

A = diag(ones(n - 3, 1), -1) + diag(ones(n - 3, 1), 1) - 2 * diag(ones(n - 2, 1));
u_exact = A\(f(2:end-1) .* h^2);
% u_exact = seidel(A, f(2:end-1) .* h^2, zeros([n-2, 1]), 1e-6, 100000);
u_exact = [u_0 u_exact' u_n];

%% calculations
subplot(1, 3, 2); hold on; grid on;
tic
% [U, E] = seidel_schwartz_2d_symmetrical(n, a, b, intersect_points_radius, f, u_0, u_n, accuracy, true, 1);
[U, E] = schwartz_2d_symmetrical(n, a, b, intersect_points_radius, f, u_0, u_n, accuracy, true);
toc 

mae = E(:, end);
mae_diff = max(abs(U(:, 2:end-1) - u_exact(2:end - 1)), [], 2);

[val_min, idx_min] = min(mae);
[val_max, idx_max] = max(mae);

subplot(1, 3, 1); hold on; grid on;
plot(x, u_exact, x, U(end, :))
legend(sprintf("Прогонка на %i узлах", n), sprintf("Шварц на %i узлах", n))

subplot(1, 3, 3); hold on; grid on;
fprintf("%f => %f, min=%f (idx=%d), max=%f (idx=%d)\n", ...
    mae(1), mae(end), val_min, idx_min, val_max, idx_max);

[val_min, idx_min] = min(mae_diff);
[val_max, idx_max] = max(mae_diff);

fprintf("%f => %f, min=%f (idx=%d), max=%f (idx=%d)\n", ...
  mae_diff(1), mae_diff(end), val_min, idx_min, val_max, idx_max);

plot(E)
plot(mae_diff)
plot(idx_min, val_min, "o")

legend("left error", "right error", "MAE", "real error", sprintf("Minimum real error at %i iteration", idx_min))
  
