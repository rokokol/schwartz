function [X, Y, U] = solve2d(x, y, u_left, u_bottom, u_right, u_top, F)
n = numel(x); m = numel(y);
n_internal_x = n - 2;
n_internal_y = m - 2;
n_internal = n_internal_x * n_internal_y;
hx = x(2) - x(1);
hy = y(2) - y(1);

[Y, X] = meshgrid(y, x);

RHS = hx^2 * hy^2 * F;
RHS(:,1)   = RHS(:,1)   - hy^2 * u_left(2:end-1);    % i = 1 (по x) сосед слева
RHS(:,end) = RHS(:,end) - hy^2 * u_right(2:end-1);   % i = n-2 (по x) сосед справа

RHS(1,:)   = RHS(1,:)   - hx^2 * u_bottom(2:end-1)'; % j = 1 (по y) сосед снизу
RHS(end,:) = RHS(end,:) - hx^2 * u_top(2:end-1)';    % j = m-2 (по y) сосед сверху

e = ones(n_internal, 1);
diagonals = [hy^2*e, hx^2*e, -2*(hx^2+hy^2)*e, hx^2*e, hy^2*e];
offsets = [-n_internal_y, -1, 0, 1, n_internal_y];
A = spdiags(diagonals, offsets, n_internal, n_internal);
% Убрать ложные связи между столбцами
b = hx^2*e;
b(n_internal_y:n_internal_y:n_internal) = 0;
A = spdiags(b, 1, A); 
A = spdiags(b, -1, A);

u_internal = A \ RHS(:);
U_internal = reshape(u_internal, m - 2, n - 2)';
U = zeros(n, m);
U(2:end-1, 2:end-1) = U_internal;
U(:, 1)   = u_bottom;  % y = 0
U(:, end) = u_top;     % y = 1
U(1, :)   = u_left;    % x = 0
U(end, :) = u_right;   % x = 1
end

