function [cx_L, cx_R, cy_B, cy_T, cx_main, cy_main] = get_2d_coeffs(x, y)
    hx = diff(x);
    hx_L = hx(1:end-1);
    hx_R = hx(2:end);
    sum_hx = hx_L + hx_R;
    cx_L = 2 ./ (hx_L .* sum_hx);
    cx_R = 2 ./ (hx_R .* sum_hx);
    cx_main = -2 ./ (hx_L .* hx_R);
    
    hy = diff(y);
    hy_B = hy(1:end-1);
    hy_T = hy(2:end);
    sum_hy = hy_B + hy_T;
    cy_B = 2 ./ (hy_B .* sum_hy);
    cy_T = 2 ./ (hy_T .* sum_hy);
    cy_main = -2 ./ (hy_B .* hy_T);
end
