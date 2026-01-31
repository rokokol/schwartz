function SchwartzApp
fig = uifigure('Name', 'Метод Шварца: Исследование', 'Position', [100 100 1200 600]);

gl = uigridlayout(fig, [1, 2]);
gl.ColumnWidth = {300, '1x'};

panel = uipanel(gl, 'Title', 'Параметры задачи');

pLayout = uigridlayout(panel, [11, 1]);
pLayout.RowHeight = {22, 25, 22, 25, 22, 25, 22, 30, 22, 25, 50};
pLayout.RowSpacing = 5;

uilabel(pLayout, 'Text', 'Кол-во узлов (N, нечетное):');
nField = uieditfield(pLayout, 'numeric', 'Value', 101, 'Limits', [5, Inf]);

uilabel(pLayout, 'Text', 'Точность (Accuracy):');
accField = uieditfield(pLayout, 'numeric', 'Value', 1e-2, 'ValueDisplayFormat', '%1.1e');

uilabel(pLayout, 'Text', 'Радиус пересечения r(n)');
radField = uieditfield(pLayout, 'text', 'Value', 'floor(n / 4)');

uilabel(pLayout, 'Text', 'Интервал [a, b]:');

abGrid = uigridlayout(pLayout, [1, 2]);
abGrid.Padding = [0 0 0 0];
abGrid.ColumnSpacing = 10;

aField = uieditfield(abGrid, 'numeric', 'Value', 0);
aField.Placeholder = 'a';

bField = uieditfield(abGrid, 'numeric', 'Value', pi);
bField.Placeholder = 'b';

uilabel(pLayout, 'Text', 'Исходная функция f(x):');
funcField = uieditfield(pLayout, 'text', 'Value', 'sin(x)');

uilabel(pLayout, 'Text', 'Распределение узлов x:');
distSlider = uislider(pLayout, 'Limits', [-2 2]);

btn = uibutton(pLayout, 'Text', 'РАССЧИТАТЬ', ...
    'ButtonPushedFcn', @(btn, event) onRun());
btn.BackgroundColor = [0.30, 0.75, 0.93];
btn.FontWeight = 'bold';
btn.FontColor = 'white';

plotsGrid = uigridlayout(gl, [1, 3]);

ax1 = uiaxes(plotsGrid); title(ax1, 'Сравнение решений'); grid(ax1, 'on');
ax2 = uiaxes(plotsGrid); title(ax2, 'Ход итераций');      grid(ax2, 'on');
ax3 = uiaxes(plotsGrid); title(ax3, 'Ошибки');            grid(ax3, 'on');

    function onRun()
        n = round(nField.Value);

        if mod(n, 2) == 0
            n = n + 1;
            nField.Value = n;
        end

        acc = accField.Value;

        a = aField.Value;
        b = bField.Value;

        if a >= b
            uialert(fig, 'Начало отрезка a должно быть меньше конца b', 'Ошибка интервала');
            return;
        end

        try
            radius_fcn = str2func(['@(n) ' radField.Value]);
            rad = round(radius_fcn(n));
            if rad < 2, rad = 2; end
        catch
            uialert(fig, 'Радиус r(n) должен быть не меньше чем 2 и не больше чем n', 'Неправильный радиус');
            return;
        end

        x = warped_linspace(a, b, n, distSlider.Value);
        h = diff(x); h = h(1:end-1);

        try
            function_fcn = str2func(['@(x) ' funcField.Value]);
            f_vec = function_fcn(x)';
        catch
            uialert(fig, 'Твоя f(x) полное говно', 'Неправильная функция');
            return;
            % f_vec = (x ./ x)'
            % f_vex(1) = 1
        end

        u_0 = 0; u_n = 0;

        cla(ax1); cla(ax2); cla(ax3);

        A = diag(ones(n - 3, 1), -1) + diag(ones(n - 3, 1), 1) - 2 * diag(ones(n - 2, 1));
        rhs = f_vec(2:end-1) .* (h .^ 2)';
        u_exact_inner = A \ rhs;
        u_exact = [u_0 u_exact_inner' u_n];

        tic;
        [U, E] = schwartz_1d_symmetrical(x, rad, f_vec, u_0, u_n, acc, true, ax2);
        time_elapsed = toc;

        hold(ax1, 'on');
        plot(ax1, x, u_exact, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Точное');
        plot(ax1, x, U(end, :), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Шварц');
        legend(ax1);
        title(ax1, sprintf('N=%d, Radius=%d', n, rad));

        real_relative_err = max(abs(U(:, 2:end-1) - u_exact(2:end - 1) ./ U(:, 2:end-1)), [], 2);
        real_err = max(abs(U(:, 2:end-1) - u_exact(2:end - 1)), [], 2);

        [val_min_real, idx_min_real] = min(real_err);

        hold(ax3, 'on');
        plot(ax3, E(:,1), 'g:', 'DisplayName', 'Left Boundary Err');
        plot(ax3, E(:,2), 'b:', 'DisplayName', 'Right Boundary Err');
        plot(ax3, real_relative_err, '-', 'LineWidth', 1.5, 'DisplayName', 'Relative Error');
        plot(ax3, real_err, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Real Error');
        plot(ax3, idx_min_real, val_min_real, 'ro', 'DisplayName', 'Min Real Err');

        legend(ax3);
        xlabel(ax3, 'Итерации');
        ylabel(ax3, 'Ошибка');
        title(ax3, sprintf('Min Real Err: %.2e (iter %d)', val_min_real, idx_min_real));

        fprintf('Расчет завершен за %.4f сек. Итераций: %d\n', time_elapsed, size(E,1));
    end 
end
