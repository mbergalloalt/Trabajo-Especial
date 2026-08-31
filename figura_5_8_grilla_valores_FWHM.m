close all;
clear

% ---- CARGAR ARCHIVOS ----

files = dir('C:\Users\marti\Desktop\Trabajo Final\Parte 3\Sim i9\Simulacion_20260601_120100\*.mat');
N_files = length(files);
D_list    = zeros(N_files, 1);
T2_list   = zeros(N_files, 1);
FWHM_list = zeros(N_files, 1);
for i = 1:N_files
    data = load(fullfile(files(i).folder, files(i).name), 'D', 'T2', 'FWHM_log', 'i2', 'Nx');
    D_list(i)    = data.D;
    T2_list(i)   = data.T2;
    FWHM_list(i) = data.FWHM_log;
end

% ---- RECONSTRUIR GRILLA ----

D_unique  = unique(D_list);
T2_unique = unique(T2_list);
Nd  = length(D_unique);
Nt2 = length(T2_unique);
fprintf('Valores unicos de D: %d\n', Nd)
fprintf('Valores unicos de T2: %d\n', Nt2)
FWHM_matrix = zeros(Nd, Nt2);
for i = 1:N_files
    k = find(D_unique  == D_list(i));
    j = find(T2_unique == T2_list(i));
    FWHM_matrix(k, j) = FWHM_list(i);
end

figure('Renderer', 'painters');
surf(T2_unique, D_unique, FWHM_matrix)
shading flat;
view(2)
set(gca, 'XScale', 'log', 'YScale', 'log')
set(gca, 'Box', 'on')
xlabel('T_2(s)', 'FontSize', 12)
ylabel('D(m^2/s)', 'FontSize', 12)
zlabel('FWHM', 'FontSize', 12)
xlim(sort([T2_unique(1) T2_unique(end)]))
ylim(sort([D_unique(1) D_unique(end)]))
clim([0 1])
colorbar
hold on
h1 = plot3(1.68579, 2.5e-9, 10, 's', 'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize',9);
h2 = plot3(2, 2.5e-9, 10, 's', 'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize',9);
h3 = plot3(2, 1.94236e-9, 10, 's', 'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize',9);
uistack([h1 h2 h3], 'top')
hold off
box on
axis square
colormap(turbo)

% ======================================================
% FIGURA APARTE: PROYECCION T2 + FWHM MEDIDO
% ======================================================
% Poné acá la ruta completa (o solo el nombre, si está en la misma
% carpeta que 'files') del archivo que quieras inspeccionar.

fileToPlot = 'C:\Users\marti\Desktop\Trabajo Final\Parte 3\Sim i9\Simulacion_20260601_120100\Data_D1p172e-9_T2605ms.mat';

dataProj = load(fileToPlot, ...
                'T2_axis', 'proj_T2', 'T2_left', 'T2_right', ...
                'D', 'T2', 'FWHM_log', 'S', 'D_axis');

figure;
semilogx(dataProj.T2_axis, dataProj.proj_T2, 'b', 'LineWidth', 2)
hold on
y_left  = interp1(dataProj.T2_axis, dataProj.proj_T2, dataProj.T2_left);
y_right = interp1(dataProj.T2_axis, dataProj.proj_T2, dataProj.T2_right);
plot(dataProj.T2_left,  y_left,  'ro', 'MarkerSize', 9, 'LineWidth', 2)
plot(dataProj.T2_right, y_right, 'ro', 'MarkerSize', 9, 'LineWidth', 2)
plot([dataProj.T2_left dataProj.T2_right], [0.5 0.5], 'r', 'LineWidth', 2)
yline(0.5, '--k', 'LineWidth', 1.2)
hold off
set(gca, 'XScale', 'log', 'FontSize', 12, 'Box', 'on')
xlim([min(dataProj.T2_axis) max(dataProj.T2_axis)])
ylim([0 1.05])
xlabel('T_2 [s]', 'FontSize', 14)
ylabel('Proyección normalizada', 'FontSize', 14)
title(sprintf('D = %.3g m^2/s, T_2 = %.3g s, FWHM_{log} = %.4g', ...
      dataProj.D, dataProj.T2, dataProj.FWHM_log), 'FontSize', 12)
legend({'Proyección','','','FWHM','Media altura'}, 'Location','best')
axis square

%======================================================
% FIGURA APARTE: MAPA D-T2 DE ESA MISMA POBLACION
%======================================================

figure;
contourf(dataProj.T2_axis, dataProj.D_axis, dataProj.S, 90, 'LineStyle', 'none')
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 14, 'Box', 'on')
xlabel('T_2 [s]', 'FontSize', 18)
ylabel('D [m^2/s]', 'FontSize', 18)
xlim([min(dataProj.T2_axis) max(dataProj.T2_axis)])
ylim([min(dataProj.D_axis) max(dataProj.D_axis)])
axis square

% Colormap con fondo blanco (misma logica que en el mapa D-T2 anterior)
nBase = 256;
base  = parula(nBase);
nRamp = 32;
ramp = [linspace(1,base(1,1),nRamp)', ...
        linspace(1,base(1,2),nRamp)', ...
        linspace(1,base(1,3),nRamp)'];
cmapBlanco = [ramp; base];
colormap(gca, cmapBlanco)

cb = colorbar;
cb.Label.String = 'Amplitud';
cb.FontSize = 12;
title(sprintf('D = %.3g m^2/s, T_2 = %.3g s', dataProj.D, dataProj.T2), ...
      'FontSize', 12)