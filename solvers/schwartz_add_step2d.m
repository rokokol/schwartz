function [U_new, max_diff] = schwartz_add_step2d(U_old, A_cells, F_int, starts, ends_arr, bpoints, cx_L, cx_R, cy_B, cy_T, n_internal_y)
    U_new = U_old;
    U_new_sum = zeros(size(U_old));
    weights = zeros(size(U_old));
    
    U_locals = cell(bpoints, 1);
    
    parfor i = 1:bpoints
        idx_x = starts(i) : ends_arr(i);
        len_x = length(idx_x);
        RHS = F_int(idx_x, :); 
        
        left_u_val = U_old(starts(i), 2:end-1); 
        RHS(1, :) = RHS(1, :) - cx_L(idx_x(1)) * left_u_val;
        
        right_u_val = U_old(ends_arr(i) + 2, 2:end-1);
        RHS(end, :) = RHS(end, :) - cx_R(idx_x(end)) * right_u_val;
        
        bottom_u_val = U_old(idx_x + 1, 1); 
        top_u_val = U_old(idx_x + 1, end);  
        
        RHS(:, 1) = RHS(:, 1) - cy_B(1) * bottom_u_val;
        RHS(:, end) = RHS(:, end) - cy_T(end) * top_u_val;
        
        u_local_flat = A_cells{i} \ RHS(:);
        U_locals{i} = reshape(u_local_flat, len_x, n_internal_y);
    end
    
    for i = 1:bpoints
        idx_x = starts(i) : ends_arr(i);
        U_new_sum(idx_x + 1, 2:end-1) = U_new_sum(idx_x + 1, 2:end-1) + U_locals{i};
        weights(idx_x + 1, 2:end-1) = weights(idx_x + 1, 2:end-1) + 1;
    end
    
    mask = weights > 0;
    U_new(mask) = U_new_sum(mask) ./ weights(mask);
    
    max_diff = max(max(abs(U_new(2:end-1, 2:end-1) - U_old(2:end-1, 2:end-1))));
end
