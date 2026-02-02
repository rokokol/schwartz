% TODO: мне нужно провести следующие эксперименты:
% - попробовать найти порядок ошибки \alpha при помощи n/h^\alpha = const
% - сделать динамическое вычисление числа точек пересечения?

function [U, e, ie] = schwartz_1d_symmetrical(x, intersect_points_radius, f, u_0, u_n, accuracy, trace, target_ax)
%SCWARTZ_2D_SYMMETRICAL реализация симметричного двумерного метода Шварца
% n >= 5
% intersect_points_radius >= 2
  n = numel(x); center = floor(n / 2); 
  left_end = center + intersect_points_radius - 2; right_start = center - intersect_points_radius + 2;
  left_piv = 2 * intersect_points_radius - 2; right_piv = - 2 * intersect_points_radius + 3;
  h = diff(x); 
  f = f(2:end-1);
  f(1) = f(1) - u_0; f(end) = f(end) - u_n;

  left_x = x(1:left_end + 1); right_x = x(right_start + 1:end);
  left_f = f(1:left_end); right_f = f(right_start:end);

  h_L = h(1:end-1);      
  h_R = h(2:end);
  sum_h = h_L + h_R;

  glob_main  = -2 ./ (h_L .* h_R);
  glob_lower =  2 ./ (h_L .* sum_h);
  glob_upper =  2 ./ (h_R .* sum_h);
  idx_L = 1:left_end; 

  m_L = glob_main(idx_L);
  l_L = glob_lower(idx_L); 
  u_L = glob_upper(idx_L);
  A_l = diag(m_L) + ...
        diag(l_L(2:end), -1) + ... 
        diag(u_L(1:end-1), 1);

  idx_R = right_start:length(f); 
  m_R = glob_main(idx_R);
  l_R = glob_lower(idx_R);
  u_R = glob_upper(idx_R);
  A_r = diag(m_R) + ...
        diag(l_R(2:end), -1) + ...
        diag(u_R(1:end-1), 1); 

  error = Inf;
  e = [];
  ie = [];
  U = [];
  
  if trace
      hold(target_ax, 'on'); 
  end

  iter = 0;
  base_left_f = left_f(end);
  base_right_f = right_f(1); 
  left_un = 0; right_u0 = 0;
  coeff_boundary_l = u_L(end); 
  coeff_boundary_r = l_R(1);
  while abs(error) > accuracy && iter < 100
      iter = iter + 1;
      
      left_f(end) = base_left_f; 
      right_f(1) = base_right_f;  

      left_f(end) = left_f(end) - coeff_boundary_l * left_un;
      right_f(1) = right_f(1) - coeff_boundary_r * right_u0;
      
      left_u = A_l\left_f; right_u = A_r\right_f;
      pleft_un = left_un; pright_u0 = right_u0; 
      left_un = right_u(left_piv); right_u0 = left_u(end + right_piv); 

      abs_diff = abs([left_u(end + right_piv:end); left_un] - [right_u0; right_u(1:left_piv)]);
      error = max(abs_diff);
      integ_error = dot(abs_diff, h(numel(left_u) + right_piv:numel(left_u) + 1));
      e = [e; error];
      ie = [ie; integ_error];
      
      if trace
          plot(target_ax, [left_x x(left_end + 2)], [u_0 left_u' pleft_un], 'Color', [0 0.4470 0.7410 0.3]); 
          plot(target_ax, [x(right_start) right_x], [pright_u0 right_u' u_n], 'Color', [0.8500 0.3250 0.0980 0.3]); 
          % Точки стыковки
          plot(target_ax, right_x(left_piv), left_un, "b.", 'MarkerSize', 5);
          plot(target_ax, left_x(end + right_piv), right_u0, "r.", 'MarkerSize', 5);
      end
      
      u_approx = zeros([1 n]);
      u_approx(1:left_end + 2) = [u_0 left_u' left_un];
      u_approx(right_start:end) = u_approx(right_start:end) + [right_u0 right_u' u_n];
      u_approx(right_start:left_end + 2) = u_approx(right_start:left_end + 2) / 2;
      U = [U; u_approx]; 
  end
end

