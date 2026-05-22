function [u_new, max_diff] = schwartz_mult_step(u_old, A_cells, f_int, starts, ends_arr, glob_lower, glob_upper, bpoints)
    u_new = u_old; 
    
    for i = 1:bpoints
        idx_int = starts(i) : ends_arr(i); 
        f_local = f_int(idx_int);
        
        left_u_val = u_new(starts(i)); 
        f_local(1) = f_local(1) - glob_lower(idx_int(1)) * left_u_val;
        
        right_u_val = u_new(ends_arr(i) + 2);
        f_local(end) = f_local(end) - glob_upper(idx_int(end)) * right_u_val;
        
        u_local = A_cells{i} \ f_local;
        
        u_new(idx_int + 1) = u_local; 
    end
    
    max_diff = max(abs(u_new(2:end-1) - u_old(2:end-1)));
end
