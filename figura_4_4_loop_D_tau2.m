% Simulacion ecuacion 2.22 - PGSTE (Casanova)

clear;
close all;

%--------------------------------------------------------
% Constantes
%--------------------------------------------------------

gamma = 42.58e6;     % Hz/T
gamma = 2*pi*gamma;
G = 7;               % T/m
D = [1e-9 1e-10 1e-11];
nD = length(D);
T2 = 10e-3;
tau2_valores = [1e-3 10e-3 50e-3];
ntau2 = length(tau2_valores);
taumin = 8e-6;
nsteps = 32;

decaimiento_minimo = 0.1;

%--------------------------------------------------------
% Figura
%--------------------------------------------------------

figure('Position',[100 50 1200 950])
tl = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
tau2_labels = { ...
    '\tau_2 = 1ms',...
    '\tau_2 = 10ms',...
    '\tau_2 = 50ms'};

letters = 'abcdefghi';
contador = 0;

for j = 1:nD
    for k = 1:ntau2
        contador = contador + 1;
        tau2 = tau2_valores(k);

        %------------------------------------------------
        % Calculo de tau1max
        %------------------------------------------------

        b_deseado = -log(decaimiento_minimo)/D(j);
        fun = @(tau1) ...
            (gamma*G*tau1).^2 .* ...
            (tau2 + 2/3*tau1) ...
            - b_deseado;
        taumax = fzero(fun,5e-4);
        taumax = max(taumax,50e-6);

        %------------------------------------------------
        % Vector tau1
        %------------------------------------------------

        amin = log10(taumin);
        amax = log10(taumax);

        tau1_eje = logspace(amin,amax,nsteps);

        %------------------------------------------------
        % Calculo señal
        %------------------------------------------------

        b = (gamma*G.*tau1_eje).^2 .*(tau2 + 2/3.*tau1_eje);
        data_ideal = exp(-b*D(j));
        data_real = exp(-b*D(j)) .*exp(-2*tau1_eje/T2);

        %------------------------------------------------
        % Plot
        %------------------------------------------------

        ax = nexttile;
        plot(b*1e-9,...
             data_ideal,...
             '-',...
             'LineWidth',2,...
             'Color',[0.8500 0.3250 0.0980]);
        hold on
        plot(b*1e-9,...
             data_real,...
             '-',...
             'LineWidth',2,...
             'Color',[0 0.4470 0.7410]);
        hold off
        ylim([0 1.05])
        grid off
        box on
        set(gca,'FontSize',14,'LineWidth',1.1)

        %------------------------------------------------
        % Titulos de columnas
        %------------------------------------------------

        if j == 1
            title(tau2_labels{k},...
                'FontSize',14)
        end

        %------------------------------------------------
        % Etiquetas de filas
        %------------------------------------------------

        if k == 1
            exponent = round(log10(D(j)));
            ylabel({ ...
                sprintf('D = 10^{%d}m^2/s',exponent)
                ''
                'señal'},...
                'FontSize',14)
        end

        %------------------------------------------------
        % Ocultar etiquetas Y interiores
        %------------------------------------------------

        if k > 1
            ax.YTickLabel = [];
            ylabel('')
        end

        %------------------------------------------------
        % Xlabel solo en fila inferior
        %------------------------------------------------

        if j == 3
            xlabel('b (10^9 s/m^2)',...
                'FontSize',14)
        end

        %------------------------------------------------
        % Etiqueta tau1,max
        %------------------------------------------------

          text(0.98,0.92,...
    strrep(sprintf('\\tau_1^{max}=%.1fms',taumax*1000),'.', ','),...
'Units','normalized',...
'HorizontalAlignment','right',...
'VerticalAlignment','top',...
'BackgroundColor','white',...
'EdgeColor','black',...
'LineWidth',0.8,...
'Margin',4,...
'FontSize',12);
    end
end

%--------------------------------------------------------
% Leyenda única
%--------------------------------------------------------

% lgd = legend({'Difusión pura',...
%               'Difusión + T_2'},...
%               'Orientation','horizontal');
% 
% lgd.Layout.Tile = 'north';
% lgd.FontSize = 16;

%--------------------------------------------------------
% Etiquetas (a), (b), ... FUERA de los ejes
%--------------------------------------------------------

ax = findall(gcf,'Type','axes');

ax = ax(~arrayfun(@(a) isempty(a.Position),ax));

pos = cell2mat(get(ax,'Position'));

[~,idx] = sortrows([-pos(:,2) pos(:,1)]);

ax = ax(idx);

for i = 1:numel(ax)
    pos = ax(i).Position;
    annotation('textbox',...
        [pos(1)-0.035,...
         pos(2)+pos(4)+0.005,...
         0.03,...
         0.03],...
        'String',['(' letters(i) ')'],...
        'EdgeColor','none',...
        'FontWeight','bold',...
        'FontSize',13);
end