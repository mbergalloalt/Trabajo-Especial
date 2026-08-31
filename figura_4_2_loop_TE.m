% simulacion CPMG y ec. 2.19 del libro de Casanova

clear all
close all

gamma = 42.58e6;      % Hz/T
G = 7;                % T/m
gamma = 2*pi*gamma;

%-------------------------------------------------------
NE = 10000;
T2 = 200e-3;

TE = [100 120 140 160 180 200]*1e-6;

A = 1;
D = 2.3e-9;
%-------------------------------------------------------

figure(3)
clf
tiledlayout(3,2,'TileSpacing','compact','Padding','compact')
letters = 'abcdef';

for i = 1:length(TE)
    t(:,i) = linspace(TE(i),NE*TE(i),NE);

    data(:,i) = A*exp(-t(:,i)/T2);

    dataD(:,i) = exp(-(1/12)*(gamma*G*TE(i))^2*D*t(:,i));

    dataT(:,i) = data(:,i).*dataD(:,i);

    T2ef(i) = 1/(1/T2 + (1/12)*(gamma*G*TE(i))^2*D );

    ax = nexttile;

    plot(t(:,i),data(:,i),'LineWidth',2,'Color',[0 0.4470 0.7410])
    hold on
    plot(t(:,i),dataD(:,i),'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
    hold off
    xlim([0 1])
    ylim([0 1])
    xlabel('t(s)','FontSize',16)
    ylabel('señal','FontSize',16)

    box on

    %-------------------------------------
    % Etiqueta del panel
    %-------------------------------------

    text(-0.13,1.08,['(' letters(i) ')'],'Units','normalized','FontWeight','bold','FontSize',14)

    %-------------------------------------
    % T2 efectivo
    %-------------------------------------

    str = sprintf('T_2^{ef} = %.1fms',1000*T2ef(i));
    str = strrep(str,'.',',');

    text(0.98,0.90,str,...
        'Units','normalized',...
        'HorizontalAlignment','right',...
        'VerticalAlignment','top',...
        'BackgroundColor','white',...
        'EdgeColor',[0.4 0.4 0.4],...
        'Margin',4,...
        'FontSize',14);

end