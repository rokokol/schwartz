function [u_approx, E] = schwartz_2d(n, a, b, f, u_0, u_n, accuracy, trace)
h = (b - a) / n; norm_mod = h * intersect_points;
center = floor(n / 2); x = linspace(a, b, n);
left_end = center + intersect_points - 2; right_start = center - intersect_points + 2;

f = f .* h ^ 2;
f = f(2:end-1);
f(1) = f(1) - u_0; f(end) = f(end) - u_n;

left_x = x(1:left_end + 1); right_x = x(right_start + 1:end);
left_un = 0; right_u0 = 0;
pleft_un = 0; pright_u0 = 0;
left_f = f(1:left_end); right_f = f(right_start:end);
A_c = diag(ones(numel(left_f) - 1, 1), -1) + diag(ones(numel(left_f) - 1, 1), 1) - 2 * diag(ones(numel(left_f), 1));

%% calculations
u_approx = zeros(n, 1);
error = Inf;
E = [];

while abs(error) > accuracy && size(E, 2) < 1000
  left_f(end) = left_f(end) + pleft_un - left_un;
  right_f(1) = right_f(1) + pright_u0 - right_u0;
  
  left_u = A_c\left_f; right_u = A_c\right_f;
  pleft_un = left_un; pright_u0 = right_u0;
  
  left_err = right_u(1) - left_un; right_err = left_u(end) - right_u0;
  error = max(abs(left_err) + abs(right_err)) / norm_mod;
  E = [E [left_err / norm_mod, right_err / norm_mod, error]']; %#ok<AGROW>
  
  if trace
    plot([left_x x(left_end + 2)], [u_0 left_u' pleft_un])
    plot([x(right_start) right_x], [pright_u0 right_u' u_n])
    
    
    plot(left_x(end), left_un, "o", right_x(1), right_u0, "o")
  end
end
end

