function [x, k, err_hist] = seidel(A, b, x0, tol, max_iter)
% GAUSS_SEIDEL Решает систему Ax=b методом Зейделя.
%   Входные данные:
%   A          - Матрица коэффициентов N x N
%   b          - Вектор правой части N x 1
%   x0         - Начальное приближение (N x 1)
%   tol        - Допуск сходимости (относительная норма)
%   max_iter   - Максимальное число итераций
%
%   Выходные данные:
%   x          - Найденное решение
%   k          - Число выполненных итераций
%   err_hist   - История относительной ошибки

N = length(b);
x = x0(:); 
err_hist = [];

if any(abs(diag(A)) < eps)
  error('Ошибка: Диагональный элемент равен нулю. Деление на ноль невозможно.');
end

for k = 1:max_iter
  x_old = x;   
  for i = 1:N
    sigma = 0;
    
    for j = 1:N
      if i ~= j
        sigma = sigma + A(i, j) * x(j);
      end
    end
    x(i) = (b(i) - sigma) / A(i, i);
  end
  
  error_current = norm(x - x_old, inf) / norm(x, inf);
  err_hist(k) = error_current;
  
  if error_current < tol
    break;
  end
end

if k == max_iter && error_current >= tol
  warning('Метод не сошелся за максимальное число итераций.');
end
end
