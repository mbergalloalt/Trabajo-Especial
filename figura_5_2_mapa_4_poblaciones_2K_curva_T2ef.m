close all;
clear 

gamma = 42.58e6; % en H/2piT
G = 7; % en T/m
gamma = 2*pi*gamma;
tau1_min = 8e-6; % para calcular b
tau2_phys = 10e-3;
T2_phys = 600e-3;
D_vals   = [1e-9, 1e-10, 1e-11, 1e-12];
TE = 91e-6;
N = 50;
M = 8000; % numero de ecos

alpha = 1E0;
Nx = 100;
Ny = 100;
T1 = logspace(-13, -7, Ny);
T2 = logspace(-1, 1, Nx);

figure(1)
hold on
for k = 1:length(D_vals)
    D    = D_vals(k);

    % Calculo de tau1_max resolviendo la ecuacion cubica b*D = target_bD
    target_bD = 5;
    coefs = [2/3, tau2_phys, 0, -target_bD/(gamma^2 * G^2 * D)];
    r = roots(coefs);
    tau1_max = max(real(r(abs(imag(r)) < 1e-15)));
    fprintf('tau1_max = %.4e s\n', tau1_max);

    tau1_phys = logspace(log10(tau1_min), log10(tau1_max), N);

    % ---- DEFINIR VALORES DE b Y T ----
    b = gamma^2 * G^2 * tau1_phys.^2 .* (tau2_phys + (2/3)*tau1_phys); % Vector b-value
    t = linspace(TE, M*TE, M)'; % Vector temporal

    % ---- VERSION VECTORIZADA ----
    decay_b = exp(-b * D);
    decay_t = exp(-t/T2_phys - (1/12)*(gamma*G*TE)^2*D*t);
    Sig     = decay_b' * decay_t';

    % ---- INVERSION ----
    tau1 = b;
    tau2 = t;
    Z    = Sig;
    tau1 = tau1';

    K1 = exp(-tau1 * T1);       % Difusion
    K2 = exp(-tau2 * (1./T2));  % T2 relaxation data

    tic
    [S, resida] = flint(K1, K2, Z, alpha);
    inversion_time = toc;
    fprintf('Tiempo de inversion: %.4f s\n', inversion_time);
    contour(T2, T1, S, 20, 'HandleVisibility','off')
end

% ---- CURVA ANALÍTICA DE SUBESTIMACIÓN PARA T2=600ms ----
D_curve = logspace(-13, -7, 100);         
T2app   = 1 ./ (1/T2_phys + (gamma*G*TE)^2/12 * D_curve);

plot(T2app, D_curve, 'r-', 'LineWidth', 2, 'DisplayName', ...
    'T_2^{ef} = 1/(A + B)');
D_mark     = 1e-9;
T2ef_mark = 1 / (1/T2_phys + (gamma*G*TE)^2/12 * D_mark);

xline(T2_phys, '--k', 'LineWidth', 1.8, 'HandleVisibility', 'off');
text(T2_phys, 8e-14, '0,6s', ...
     'FontSize', 10, 'Color', 'k', 'Rotation', 0, ...
     'HorizontalAlignment', 'right', 'Clipping', 'off')
xline(T2ef_mark, '--r', 'LineWidth', 1.8, 'HandleVisibility', 'off');
text(T2ef_mark, 8e-14, strrep(sprintf('%.3fs', T2ef_mark), '.', ','), ...
     'FontSize', 10, 'Color', 'r', 'Rotation', 0, ...
     'HorizontalAlignment', 'right', 'Clipping', 'off')
plot(T2ef_mark, D_mark, 'ro', 'MarkerSize', 5, ...
     'MarkerFaceColor', 'r', 'HandleVisibility', 'off')
legend('Location', 'best', 'FontSize', 11);
hold off

set(gca, 'YScale', 'log', 'FontSize', 13)
set(gca, 'XScale', 'log', 'FontSize', 13)
xlabel('T_{2}(s)', 'FontSize', 14)
ylabel('D(m²/s)',  'FontSize', 14)
colorbar
xlim([1e-1, 1e1])
ylim([1e-13, 1e-7])
box on
axis square


