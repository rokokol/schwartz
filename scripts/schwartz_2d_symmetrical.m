function [U_hist, e, iter] = schwartz_2d_symmetrical(x, y, bpoints, intersect_points_radius, F, u_left, u_bottom, u_right, u_top, accuracy, is_mult, trace, target_ax)
  n = numel(x); 
  m = numel(y);
  N_int_x = n - 2; 
  N_int_y = m - 2;
  step = floor(N_int_x / bpoints); 
  
  [cx_L, cx_R, cy_B, cy_T, cx_main, cy_main] = get_2d_coeffs(x, y);
  F_int = F(2:end-1, 2:end-1); 
  
  A_cells = cell(1, bpoints);
  starts = zeros(1, bpoints);   
  ends_arr = zeros(1, bpoints); 

  % Decomposition along X-axis
  % TODO: add decomposition along longer axis
  for i = 1:bpoints
    base_s = (i - 1) * step + 1;
    if i == bpoints
        base_e = N_int_x; 
    else
        base_e = i * step;
    end
    
    starts(i) = max(1, base_s - intersect_points_radius);
    ends_arr(i) = min(N_int_x, base_e + intersect_points_radius);
    
    idx_x = starts(i) : ends_arr(i);
    A_cells{i} = build_2d_matrix(idx_x, cx_L, cx_R, cx_main, cy_B, cy_T, cy_main);
  end 

  U_curr = zeros(n, m); 
  U_curr(1, :) = u_left;
  U_curr(end, :) = u_right;
  U_curr(:, 1) = u_bottom;
  U_curr(:, end) = u_top;

  iter = 0;
  e = [];
  U_hist = {}; 
  
  if is_mult
      method = @schwartz_mult_step2d;
  else
      method = @schwartz_add_step2d;
  end
  
  if trace
      hold(target_ax, 'on'); 
  end

  flag = true;
  while flag
      iter = iter + 1;
      
      [U_next, max_diff] = method(U_curr, A_cells, F_int, starts, ends_arr, bpoints, cx_L, cx_R, cy_B, cy_T, N_int_y);
      
      error = max_diff;
      e = [e; error];
      U_hist{iter} = U_next; 
      U_curr = U_next;

      if trace
         cla(target_ax);
         surf(target_ax, x, y, U_curr');
         title(target_ax, sprintf('Итерация %d', iter));
         drawnow;
      end

      flag = error > accuracy && iter < 1000;
  end
end
