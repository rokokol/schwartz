% NOTE:
% Скорее всего что-то не так с методом Зейделя
function [U, E] = seidel_schwartz_2d_symmetrical(n, a, b, intersect_points_radius, f, u_0, u_n, accuracy, trace, iters)
%SCWARTZ_2D_SYMMETRICAL реализация симметричного двумерного метода Шварца
% n >= 5 && n % 2 == 1
% intersect_points_radius >= 2

if trace
  hold on; grid on;
end

h = (b - a) / n; norm_mod = h * intersect_points_radius;
center = floor(n / 2); x = linspace(a, b, n);
left_end = center + intersect_points_radius - 2; right_start = center - intersect_points_radius + 2;
left_piv = 2 * intersect_points_radius - 2; right_piv = - 2 * intersect_points_radius + 3;

f = f(2:end-1); f = f .* h ^ 2;
f(1) = f(1) - u_0; f(end) = f(end) - u_n;

left_x = x(1:left_end + 1); right_x = x(right_start + 1:end);
left_un = 0; right_u0 = 0;
pleft_un = 0; pright_u0 = 0;
left_f = f(1:left_end); right_f = f(right_start:end);
left_u = rand(size(left_f)); right_u = rand(size(right_f));
A_c = diag(ones(numel(left_f) - 1, 1), -1) + diag(ones(numel(left_f) - 1, 1), 1) - 2 * diag(ones(numel(left_f), 1));

%% calculations
error = Inf;
E = [];
U = [];

while abs(error) > accuracy && size(E, 1) < 1000
  left_f(end) = left_f(end) + pleft_un - left_un;
  right_f(1) = right_f(1) + pright_u0 - right_u0;
  
  left_u = seidel(A_c, left_f, left_u, accuracy, iters); right_u = seidel(A_c, right_f, right_u, accuracy, iters);
  pleft_un = left_un; pright_u0 = right_u0;
  left_un = right_u(left_piv); right_u0 = left_u(end + right_piv); % Корректировка
  
  left_err = right_u(1) - left_un; right_err = left_u(end) - right_u0;
  error = max(abs(left_err) + abs(right_err)) / norm_mod;
  E = [E; [left_err / norm_mod, right_err / norm_mod, error]];
  
  if trace
    plot([left_x x(left_end + 2)], [u_0 left_u' pleft_un])
    plot([x(right_start) right_x], [pright_u0 right_u' u_n])
    plot(right_x(left_piv), left_un, "o", left_x(end + right_piv), right_u0, "o")
  end
  
  u_approx = zeros([1 n]);
  u_approx(1:left_end + 2) = [u_0 left_u' left_un];
  u_approx(right_start:end) = u_approx(right_start:end) + [right_u0 right_u' u_n];
  u_approx(right_start:left_end + 2) = u_approx(right_start:left_end + 2) / 2;
  U = [U; u_approx]; %#ok<*AGROW>
end

if trace
  legend(sprintf("Пересечение: %i точек за %i итераций", intersect_points_radius, size(E, 1)))
end
end

