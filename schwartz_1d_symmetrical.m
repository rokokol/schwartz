% TODO: мне нужно провести следующие эксперименты:
% - попробовать найти порядок ошибки \alpha при помощи n/h^\alpha = const
% - сделать динамическое вычисление числа точек пересечения?

function [U, e] = schwartz_1d_symmetrical(x, bpoints, intersect_points_radius, f, u_0, u_n, accuracy, trace, target_ax)
%SCWARTZ_1D_SYMMETRICAL реализация мультипликативного одномерного метода Шварца
% n >= 5
% intersect_points_radius >= 2

  n = numel(x); 
  N_int = n - 2; % Количество внутренних узлов
  step = floor(N_int / bpoints); 
  
  h = diff(x); 
  f_int = f(2:end-1); % Правая часть только для внутренних узлов
  
  A = cell(1, bpoints);
  starts = zeros(1, bpoints);   
  ends_arr = zeros(1, bpoints); 

  h_L = h(1:end-1);      
  h_R = h(2:end);
  sum_h = h_L + h_R;
  
  glob_main  = -2 ./ (h_L .* h_R);
  glob_lower =  2 ./ (h_L .* sum_h);
  glob_upper =  2 ./ (h_R .* sum_h);

  for i = 1:bpoints
    base_s = (i - 1) * step + 1;
    if i == bpoints
        base_e = N_int; 
    else
        base_e = i * step;
    end
    
    % Важно случайно не вылететь за границы массива
    starts(i) = max(1, base_s - intersect_points_radius);
    ends_arr(i) = min(N_int, base_e + intersect_points_radius);
    
    idx_curr = starts(i) : ends_arr(i); 
    m_curr = glob_main(idx_curr);
    l_curr = glob_lower(idx_curr); 
    u_curr = glob_upper(idx_curr);
    A{i} = diag(m_curr) + ...
           diag(l_curr(2:end), -1) + ... 
           diag(u_curr(1:end-1), 1);
  end 

  u_curr = zeros(n, 1); 
  u_curr(1) = u_0;      
  u_curr(end) = u_n;    

  error = Inf;
  iter = 0;
  e = [];
  U = [];

  base_colors = [
      0.95, 0.05, 0.05; 
      0.40, 0.75, 0.00; 
      1.00, 0.80, 0.00; 
      0.10, 0.55, 0.75; 
      0.80, 0.20, 0.60; 
      0.20, 0.80, 0.70; 
      1.00, 0.40, 0.00  
  ];
  colors = repmat(base_colors, ceil(bpoints / size(base_colors, 1)), 1);
  
  if trace
      hold(target_ax, 'on'); 
  end

  while error > accuracy && iter < 100
      iter = iter + 1;
      
      % Аддитивный или мультипликативный проход по подобластям
      [u_next, max_diff] = schwartz_mult_step(u_curr, A, f_int, starts, ends_arr, glob_lower, glob_upper, bpoints);
      
      error = max_diff;
      e = [e; error];
      U = [U; u_next']; 
      
      u_curr = u_next;
      
      % Рисуем по подобластям
      if trace
          for i = 1:bpoints
              idx_plot = starts(i) : ends_arr(i) + 2; 
              c = [colors(i, :), 0.5]; 
              plot(target_ax, x(idx_plot), u_curr(idx_plot), 'Color', c, 'LineWidth', 1.5);
          end
      end
  end
end

