clc; clear; clf; close all;

%% init
n = 13; a = 0; b = 2 * pi; h = (b - a) / n; % n >= 13
intersect_points = 3;
norm_mod = h * intersect_points;
center = floor(n / 2);
left_end = center + intersect_points - 2; right_start = center - intersect_points + 2;

u_0 = 0; u_n = 0;
% accuracy = 0.0476;
accuracy = 1e-3;

x = linspace(a, b, n);
pivot = x(center + 1)
f = sin(x)' .* h ^ 2;
f = f(2:end-1);
f(1) = f(1) - u_0; f(end) = f(end) - u_n;

left_x = x(1:left_end + 1); right_x = x(right_start + 1:end);
left_un = 0; right_u0 = 0;
pleft_un = 0; pright_u0 = 0;
left_f = f(1:left_end); right_f = f(right_start:end);

A = diag(ones(n - 3, 1), -1) + diag(ones(n - 3, 1), 1) - 2 * diag(ones(n - 2, 1));
A_c = diag(ones(numel(left_f) - 1, 1), -1) + diag(ones(numel(left_f) - 1, 1), 1) - 2 * diag(ones(numel(left_f), 1));
u_real = A\f;

figure(1);
subplot(1, 3, 1); hold on; grid on;
plot(x, [u_0 u_real' u_n])
legend(sprintf("Прогонка на %i узлах", n))

%% calculations
subplot(1, 3, 2); hold on; grid on;
u_approx = zeros(n, 1);
error = Inf;
E = [];
i = 1;

% NOTE: из-за неправильного выбора опорной точки
%       даже при минимально возможном n=13 программа
%       сходится медленнее, чем при n=131.
%       Однако при n=1341 из-за слишком большой разницей
%       между n и intersect_points скорость сходимости падает

% TODO: починить выбор опорной точки и проверить гипотезу
%       о том, что чем больше разница между n и intersect_points
%       тем ниже скорость сходимости

while abs(error) > accuracy && i < 1000
    % while i < 16
    left_f(end) = left_f(end) + pleft_un - left_un;
    right_f(1) = right_f(1) + pright_u0 - right_u0;
    
    left_u = A_c\left_f; right_u = A_c\right_f;
    pleft_un = left_un; pright_u0 = right_u0;
    left_un = right_u(2 * intersect_points - 3); right_u0 = left_u(end - 2 * intersect_points + 4); % Корректировка
    plot([left_x x(left_end + 2)], [u_0 left_u' pleft_un]) %, '-', 'DisplayName', sprintf('left %d', i))
    plot([x(right_start) right_x], [pright_u0 right_u' u_n]) %, '-', 'DisplayName', sprintf('right %d', i))
    
    left_err = right_u(1) - left_un; right_err = left_u(end) - right_u0; %#ok<*NOPTS>
    error = max(abs(left_err) + abs(right_err)) / norm_mod
    E = [E [left_err / norm_mod, right_err / norm_mod, error]']; %#ok<AGROW>
    
    % modif = 1 + 1 / i^2; left_un = modif * left_un; right_u0 = modif * right_u0;
    plot(left_x(end), left_un, "o", right_x(1), right_u0, "o")
    
    i = i + 1;
end

i = i - 1
legend(sprintf("Пересечение: %i точек за %i итераций", intersect_points, i))

subplot(1, 3, 3); hold on; grid on;
plot(E')
legend("left error", "right error", "MAE")

