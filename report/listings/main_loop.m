flag = true;
while flag
    iter = iter + 1;

    [u_next, max_diff] = method(u_curr, A, f_int, ...
                                starts, ends_arr,   ...
                                glob_lower, glob_upper, bpoints);
    e      = [e;      max_diff];
    U      = [U;      u_next' ];
    u_curr = u_next;

    flag = max_diff > accuracy && iter < 100;
end
