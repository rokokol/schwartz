function schwartzapp2d
addpath('../utils');
addpath('../solvers');
addpath('../scripts');
    
fig = uifigure('Name', 'Метод Шварца 2D: Исследование', 'Position', [100 100 1300 700]); 

gl = uigridlayout(fig, [1, 2]);
gl.ColumnWidth = {350, '1x'};

panel = uipanel(gl, 'Title', 'Параметры задачи');

pLayout = uigridlayout(panel, [18, 1]);
pLayout.RowHeight = {22, 25, 22, 25, 22, 25, 22, 30, 22, 25, 22, 25, 22, 25, 22, 25, 22, 50};
pLayout.RowSpacing = 5;

uilabel(pLayout, 'Text', 'Узлы по X и Y (N, M):');
nmGrid = uigridlayout(pLayout, [1, 2]);
nmGrid.Padding = [0 0 0 0];
nmGrid.ColumnSpacing = 10;
nField = uieditfield(nmGrid, 'numeric', 'Value', 50, 'Limits', [5, Inf]);
mField = uieditfield(nmGrid, 'numeric', 'Value', 50, 'Limits', [5, Inf]);

uilabel(pLayout, 'Text', 'Точность (Accuracy):');
accField = uieditfield(pLayout, 'numeric', 'Value', 1e-3, 'ValueDisplayFormat', '%1.1e');

uilabel(pLayout, 'Text', 'Радиус пересечения r(n) (по X):');
radField = uieditfield(pLayout, 'text', 'Value', 'floor(n / 4)');

uilabel(pLayout, 'Text', 'Количество интервалов (по X):');
breakpointsField = uieditfield(pLayout, 'numeric', ...
  'Value', 2, ...             
  'Limits', [1, 30], ...      
  'RoundFractionalValues', 'on');

uilabel(pLayout, 'Text', 'Интервал X [a, b] и Y [c, d]:');
abcdGrid = uigridlayout(pLayout, [1, 4]);
abcdGrid.Padding = [0 0 0 0];
abcdGrid.ColumnSpacing = 5;
aField = uieditfield(abcdGrid, 'numeric', 'Value', 0);
aField.Placeholder = 'a';
bField = uieditfield(abcdGrid, 'numeric', 'Value', 1);
bField.Placeholder = 'b';
cField = uieditfield(abcdGrid, 'numeric', 'Value', 0);
cField.Placeholder = 'c';
dField = uieditfield(abcdGrid, 'numeric', 'Value', 1);
dField.Placeholder = 'd';

uilabel(pLayout, 'Text', 'Мультипликативный/Аддитивный метод');
method = uicheckbox(pLayout, 'Text', 'Использовать мультипликативный метод');

uilabel(pLayout, 'Text', 'Исходная функция f(x, y):');
funcField = uieditfield(pLayout, 'text', 'Value', 'sin(pi*x).*sin(pi*y)');

uilabel(pLayout, 'Text', 'Граничные(L[y], R[y], B[x], T[x]):');
bcGrid = uigridlayout(pLayout, [1, 4]);
bcGrid.Padding = [0 0 0 0];
bcGrid.ColumnSpacing = 5;
uLField = uieditfield(bcGrid, 'text', 'Value', '0.*y'); uLField.Placeholder = 'u_L(y)';
uRField = uieditfield(bcGrid, 'text', 'Value', '0.*y'); uRField.Placeholder = 'u_R(y)';
uBField = uieditfield(bcGrid, 'text', 'Value', '0.*x'); uBField.Placeholder = 'u_B(x)';
uTField = uieditfield(bcGrid, 'text', 'Value', '0.*x'); uTField.Placeholder = 'u_T(x)';

uilabel(pLayout, 'Text', 'Распределение узлов:');
distSlider = uislider(pLayout, 'Limits', [-2 2]);

timeLabel = uilabel(pLayout, 'Text', 'Время выполнения: --');

btn = uibutton(pLayout, 'Text', 'РАССЧИТАТЬ', ...
    'ButtonPushedFcn', @(btn, event) onRun());
btn.BackgroundColor = [0.30, 0.75, 0.93];
btn.FontWeight = 'bold';
btn.FontColor = 'white';

plotsGrid = uigridlayout(gl, [1, 3]);

ax1 = uiaxes(plotsGrid); title(ax1, 'Сравнение решений (3D)'); grid(ax1, 'on'); view(ax1, 3);
ax2 = uiaxes(plotsGrid); title(ax2, 'Ход итераций (3D)');      grid(ax2, 'on'); view(ax2, 3);
ax3 = uiaxes(plotsGrid); title(ax3, 'Ошибки');            grid(ax3, 'on');

function onRun()
    tStart = tic;
    n = round(nField.Value);
    m = round(mField.Value);

    if n <= 5 || m <= 5
        uialert(fig, 'Неправильное чилсло узлов n или m', 'Ошибка');
        return;
    end

    acc = accField.Value;

    a = aField.Value; b = bField.Value;
    c = cField.Value; d = dField.Value;

    if a >= b || c >= d
        uialert(fig, 'Начало отрезка должно быть меньше конца', 'Ошибка интервала');
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
    y = warped_linspace(c, d, m, distSlider.Value);
    bpoints = breakpointsField.Value;

    [Y_mesh, X_mesh] = meshgrid(y, x);

    try
        function_fcn = str2func(['@(x, y) ' funcField.Value]);
        F_vec = function_fcn(X_mesh, Y_mesh);
    catch
        uialert(fig, 'f(x, y) задана неверно', 'Неправильная функция');
        return;
    end

    try
        ul_fcn = str2func(['@(y) ' uLField.Value]); u_left = ul_fcn(y)';
        ur_fcn = str2func(['@(y) ' uRField.Value]); u_right = ur_fcn(y)';
        ub_fcn = str2func(['@(x) ' uBField.Value]); u_bottom = ub_fcn(x)';
        ut_fcn = str2func(['@(x) ' uTField.Value]); u_top = ut_fcn(x)';
        
        % Ensure column vectors match output of solve2d expectation
        if isrow(u_left), u_left = u_left'; end
        if isrow(u_right), u_right = u_right'; end
        if isrow(u_bottom), u_bottom = u_bottom'; end
        if isrow(u_top), u_top = u_top'; end
    catch
        uialert(fig, 'Граничные заданы неверно', 'Ошибка ГУ');
        return;
    end

    cla(ax1); cla(ax2); cla(ax3);

    % Exact solution via full 2D solver
    u_exact = solve2d(x, y, u_left, u_bottom, u_right, u_top, F_vec(2:end-1, 2:end-1));

    fStart = tic;
    [U_hist, e, iters] = schwartz_2d_symmetrical(x, y, bpoints, rad, F_vec, u_left, u_bottom, u_right, u_top, acc, method.Value, true, ax2);
    time_elapsed = toc(fStart);
    timeLabel.Text = sprintf('Время выполнения: %.4f сек.', time_elapsed);

    hold(ax1, 'off'); % Surf works better with hold off or clear
    
    if ~isempty(U_hist)
        U_final = U_hist{end};
    else
        U_final = u_exact; % fallback
    end
    
    % Plot exact and schwartz solutions as surfaces
    surf(ax1, X_mesh, Y_mesh, u_exact, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Точное');
    hold(ax1, 'on');
    surf(ax1, X_mesh, Y_mesh, U_final, 'FaceAlpha', 0.8, 'DisplayName', 'Шварц');
    
    legend(ax1);
    title(ax1, sprintf('N=%d, M=%d, Radius=%d', n, m, rad));

    real_err = zeros(size(e));
    for i=1:length(U_hist)
        real_err(i) = max(max(abs(U_hist{i}(2:end-1, 2:end-1) - u_exact(2:end-1, 2:end-1))));
    end

    [val_min_real, idx_min_real] = min(real_err);

    hold(ax3, 'on');
    plot(ax3, e, 'g-', 'DisplayName', 'Estimated Error');
    plot(ax3, real_err, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Real Error');
    plot(ax3, idx_min_real, val_min_real, 'ro', 'DisplayName', 'Min Real Err');

    legend(ax3);
    xlabel(ax3, 'Итерации');
    ylabel(ax3, 'Ошибка');
    title(ax3, sprintf('Min Real Err: %.2e (iter %d)', val_min_real, idx_min_real));

    fprintf('Расчет завершен за %.4f сек. Итераций: %d\n', time_elapsed, iters);
    toc(tStart)
end 
end
