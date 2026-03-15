function [U_new, max_diff] = schwartz_mult_step2d(U_old, A_cells, F_int, starts, ends_arr, bpoints, cx_L, cx_R, cy_B, cy_T, n_internal_y)
    U_new = U_old;
    for i = 1:bpoints
        idx_x = starts(i) : ends_arr(i);
        len_x = length(idx_x);
        
        RHS = F_int(idx_x, :); % size len_x x n_internal_y
        
        left_u_val = U_new(starts(i), 2:end-1); % 1 x n_internal_y
        RHS(1, :) = RHS(1, :) - cx_L(idx_x(1)) * left_u_val;
        
        right_u_val = U_new(ends_arr(i) + 2, 2:end-1); % 1 x n_internal_y
        RHS(end, :) = RHS(end, :) - cx_R(idx_x(end)) * right_u_val;
        
        bottom_u_val = U_new(idx_x + 1, 1); % len_x x 1 (column vector)
        top_u_val = U_new(idx_x + 1, end);  % len_x x 1
        
        RHS(:, 1) = RHS(:, 1) - cy_B(1) * bottom_u_val;
        RHS(:, end) = RHS(:, end) - cy_T(end) * top_u_val;
        
        u_local_flat = A_cells{i} \ RHS(:);
        U_local = reshape(u_local_flat, len_x, n_internal_y); 
        
        U_new(idx_x + 1, 2:end-1) = U_local;
    end
    
    max_diff = max(max(abs(U_new(2:end-1, 2:end-1) - U_old(2:end-1, 2:end-1))));
end
