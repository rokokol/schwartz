% TODO: мне нужно провести следующие эксперименты:
% - добавить сгущение точек у границы
% - попробовать найти порядок ошибки \alpha при помощи n/h^\alpha = const
% - ввести мнимую точку для нечетных n
% - сделать динамическое вычисление числа точек пересечения?

function [U, E] = schwartz_1d_symmetrical(x, intersect_points_radius, f, u_0, u_n, accuracy, trace, target_ax)
%SCWARTZ_2D_SYMMETRICAL реализация симметричного двумерного метода Шварца
% n >= 5 && n % 2 == 0
% intersect_points_radius >= 2
  
  n = numel(x); center = floor(n / 2); 
  left_end = center + intersect_points_radius - 2; right_start = center - intersect_points_radius + 2;
  left_piv = 2 * intersect_points_radius - 2; right_piv = - 2 * intersect_points_radius + 3;
  h = diff(x); h = h(1:end-1);
  f = f(2:end-1); f = f .* (h .^ 2)';
  f(1) = f(1) - u_0; f(end) = f(end) - u_n;

  left_x = x(1:left_end + 1); right_x = x(right_start + 1:end);
  left_un = 0; right_u0 = 0;
  pleft_un = 0; pright_u0 = 0;
  left_f = f(1:left_end); right_f = f(right_start:end);
  A_c = diag(ones(numel(left_f) - 1, 1), -1) + diag(ones(numel(left_f) - 1, 1), 1) - 2 * diag(ones(numel(left_f), 1));

  error = Inf;
  E = [];
  U = [];
  
  if trace
      hold(target_ax, 'on'); % Используем target_ax вместо gcf
  end

  iter = 0;
  while abs(error) > accuracy && iter < 100
      iter = iter + 1;
      
      left_f(end) = left_f(end) + pleft_un - left_un;
      right_f(1) = right_f(1) + pright_u0 - right_u0;
      
      left_u = A_c\left_f; right_u = A_c\right_f;
      pleft_un = left_un; pright_u0 = right_u0;
      left_un = right_u(left_piv); right_u0 = left_u(end + right_piv); 
      
      left_err = right_u(1) - left_un; right_err = left_u(end) - right_u0;
      error = max(abs(left_err) + abs(right_err));
      E = [E; [left_err, right_err, error]];
      
      if trace
          % ВАЖНО: передаем первый аргумент (ось) в каждую команду plot
          plot(target_ax, [left_x x(left_end + 2)], [u_0 left_u' pleft_un], 'Color', [0 0.4470 0.7410 0.3]); % Полупрозрачный синий
          plot(target_ax, [x(right_start) right_x], [pright_u0 right_u' u_n], 'Color', [0.8500 0.3250 0.0980 0.3]); % Полупрозрачный красный
          % Точки стыковки
          plot(target_ax, right_x(left_piv), left_un, "b.", 'MarkerSize', 5);
          plot(target_ax, left_x(end + right_piv), right_u0, "r.", 'MarkerSize', 5);
          
          drawnow limitrate; % Обязательно, чтобы интерфейс обновлялся в реальном времени!
      end
      
      u_approx = zeros([1 n]);
      u_approx(1:left_end + 2) = [u_0 left_u' left_un];
      u_approx(right_start:end) = u_approx(right_start:end) + [right_u0 right_u' u_n];
      u_approx(right_start:left_end + 2) = u_approx(right_start:left_end + 2) / 2;
      U = [U; u_approx]; 
  end
end

