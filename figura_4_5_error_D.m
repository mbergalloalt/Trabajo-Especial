clear;
close all;

% Parametros fijos
gamma   = 2*pi*42.58e6;
G0      = 7;
tau2    = 10e-3;
nsteps  = 50;
tau1_min = 8e-6;

% Vectores de D y T2

D_vector  = logspace(-11, -8, 100);
T2_vector = logspace(-2, log10(2000e-3), 100);

% Matriz de error porcentual
error_pct = zeros(length(T2_vector), length(D_vector));

tau1_max_vector = zeros(1, length(D_vector));

for j = 1:length(D_vector)
    D = D_vector(j);
    % tau1_max para cada D
    coefs    = [2/3, tau2, 0, -1/(gamma^2 * G0^2 * D)];
    r        = roots(coefs);
    tau1_max = max(real(r(abs(imag(r)) < 1e-15)));
    tau1_max_vector(j) = tau1_max;  % guardar
    tau1     = logspace(log10(tau1_min), log10(tau1_max), nsteps)';
    b        = gamma^2 * G0^2 * tau1.^2 .* (tau2 + (2/3)*tau1);
    
    for k = 1:length(T2_vector)
        T2   = T2_vector(k);
        S_S0 = exp(-b*D - 2*tau1/T2);
        y    = log(S_S0);
        D_fit = -(b \ y);
        error_pct(k,j) = abs(D_fit - D) / D * 100;
    end
end

% Grilla para grafico
[X, Y] = meshgrid(D_vector, T2_vector*1e3);

figure(1);
pcolor(D_vector, T2_vector*1e3, error_pct);
shading flat;  % elimina la cuadrícula
xlabel('D (m^2/s)');
ylabel('T_2 (ms)');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
cb = colorbar;
cb.Label.String = '\DeltaD [%]';
cb.Label.FontSize = 12;
grid off;
axis square
colormap(turbo)


figure(2);
[C, h] = contour(D_vector, T2_vector*1e3, error_pct, [5, 5]);
% Extraer coordenadas del contorno
x_contorno = C(1, 2:end);   % valores de D
y_contorno = C(2, 2:end);   % valores de T2
semilogx(x_contorno, y_contorno, 'r-', 'LineWidth', 2);
xlabel('D [m^2/s]');
ylabel('T_2 [ms]');
title('Contorno donde error = 5%');
grid on;


figure(3);
pcolor(D_vector, T2_vector*1e3, error_pct);
shading flat;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
hold on;
plot(x_contorno, y_contorno, 'r-', 'LineWidth', 3);
hold off;
xlabel('D(m^2/s)', 'FontSize', 16);
ylabel('T_2(ms)', 'FontSize', 16);
cb = colorbar;
cb.Label.String = '\DeltaD [%]';
cb.Label.FontSize = 16;
% cb.Position(1) = cb.Position(1) + 0.115;
% cb.Position(3) = cb.Position(3) * 0.6;
grid off;
axis square
colormap(turbo)



