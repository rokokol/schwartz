function U = solve2d(x, y, u_left, u_bottom, u_right, u_top, F_int)
    n = numel(x); 
    m = numel(y);

    [cx_L, cx_R, cy_B, cy_T, cx_main, cy_main] = get_2d_coeffs(x, y);
    A = build_2d_matrix(1:(n-2), cx_L, cx_R, cx_main, cy_B, cy_T, cy_main);

    RHS = F_int;

    % Create robust boundary vectors regardless of row/col orientation
    u_l = reshape(u_left(2:end-1), 1, []);    % 1 x n_internal_y
    u_r = reshape(u_right(2:end-1), 1, []);   
    u_b = reshape(u_bottom(2:end-1), [], 1);  % n_internal_x x 1
    u_t = reshape(u_top(2:end-1), [], 1);     

    RHS(1, :)   = RHS(1, :)   - cx_L(1) * u_l;    
    RHS(end, :) = RHS(end, :) - cx_R(end) * u_r;   
    RHS(:, 1)   = RHS(:, 1)   - cy_B(1) * u_b;  
    RHS(:, end) = RHS(:, end) - cy_T(end) * u_t;    

    u_internal = A \ RHS(:);

    U_internal = reshape(u_internal, n-2, m-2);

    U = zeros(n, m);
    U(2:end-1, 2:end-1) = U_internal;

    U(:, 1)   = reshape(u_bottom, [], 1);  % y = 0
    U(:, end) = reshape(u_top, [], 1);     % y = 1
    U(1, :)   = reshape(u_left, 1, []);    % x = 0
    U(end, :) = reshape(u_right, 1, []);   % x = 1
end
