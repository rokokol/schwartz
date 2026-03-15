function A = build_2d_matrix(x_idx, cx_L, cx_R, cx_main, cy_B, cy_T, cy_main)
    len_x = length(x_idx);
    len_y = length(cy_B);
    N = len_x * len_y;
    
    I = zeros(N * 5, 1);
    J = zeros(N * 5, 1);
    V = zeros(N * 5, 1);
    idx = 0;
    
    for j = 1:len_y
        for i = 1:len_x
            row = i + (j-1)*len_x;
            
            % Main diagonal
            idx = idx + 1;
            I(idx) = row; J(idx) = row;
            V(idx) = cx_main(x_idx(i)) + cy_main(j);
            
            % Left neighbor (X)
            if i > 1
                idx = idx + 1;
                I(idx) = row; J(idx) = row - 1;
                V(idx) = cx_L(x_idx(i));
            end
            
            % Right neighbor (X)
            if i < len_x
                idx = idx + 1;
                I(idx) = row; J(idx) = row + 1;
                V(idx) = cx_R(x_idx(i));
            end
            
            % Bottom neighbor (Y)
            if j > 1
                idx = idx + 1;
                I(idx) = row; J(idx) = row - len_x;
                V(idx) = cy_B(j);
            end
            
            % Top neighbor (Y)
            if j < len_y
                idx = idx + 1;
                I(idx) = row; J(idx) = row + len_x;
                V(idx) = cy_T(j);
            end
        end
    end
    
    A = sparse(I(1:idx), J(1:idx), V(1:idx), N, N);
end
