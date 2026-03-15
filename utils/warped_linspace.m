function res = warped_linspace(a, b, n, k)
    x = linspace(-1, 1, n);
    p = exp(k);
    warped_x = sign(x) .* (abs(x) .^ p);
    
    res = a + (b - a) * (warped_x + 1) / 2;
end

