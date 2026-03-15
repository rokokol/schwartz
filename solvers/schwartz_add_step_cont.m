function [u_new, max_diff] = schwartz_add_step(u_old, A_cells, f_int, starts, ends_arr, glob_lower, glob_upper, bpoints)
    u_new_sum = zeros(size(u_old));
    weights = zeros(size(u_old));
    
    for i = 1:bpoints
        idx_int = starts(i) : ends_arr(i); 
        f_local = f_int(idx_int);
        
        left_u_val = u_old(starts(i)); 
        f_local(1) = f_local(1) - glob_lower(idx_int(1)) * left_u_val;
        
        right_u_val = u_old(ends_arr(i) + 2);
        f_local(end) = f_local(end) - glob_upper(idx_int(end)) * right_u_val;
        
        u_local = A_cells{i} \ f_local;
        
        u_new_sum(idx_int + 1) = u_new_sum(idx_int + 1) + u_local;
        weights(idx_int + 1) = weights(idx_int + 1) + 1;
    end
    
    u_new = u_old;
    mask = weights > 0;
    u_new(mask) = u_new_sum(mask) ./ weights(mask);
    
    max_diff = max(abs(u_new(2:end-1) - u_old(2:end-1)));
end