clear 
close all

gamma = 42.58e6; % Hz/T
G = 7;           % T/m
gamma = 2*pi*gamma;

%-------------------------------------------------
NE = 10000;
TE = 100e-6;
T2 = 200e-3;
A  = 1;

D = [9 2.3 0.2 0.02]*1e-9;
t = linspace(TE,TE*NE,NE);

data = A*exp(-t/T2);

dataD = zeros(length(t),length(D));
T2ef  = zeros(size(D));

%-------------------------------------------------

figure(1)
clf

tiledlayout(2,2,...
    'TileSpacing','compact',...
    'Padding','compact')

letters = 'abcd';

for i = 1:4
    dataD(:,i) = exp(-(1/12)*(gamma*G*TE)^2*D(i)*t);
    T2ef(i) = 1/(1/T2 + (1/12)*(gamma*G*TE)^2*D(i));

    ax = nexttile;
    plot(t,data,'-','LineWidth',3.5,'Color',[0 0.4470 0.7410])
    hold on
    plot(t,dataD(:,i),'-','LineWidth',3.5,'Color',[0.8500 0.3250 0.0980])
    hold off
    xlabel('t(ms)','FontSize',16)
    ylabel('señal','FontSize',16)
    xlim([0 1])
    ylim([0 1])

    box on

    %-----------------------------------------
    % T2ef perfectamente alineado a la derecha
    %-----------------------------------------

    str = sprintf('T_2^{ef} = %.3fms',1000*T2ef(i));

    text(0.98,0.90,str,...
        'Units','normalized',...
        'HorizontalAlignment','right',...
        'VerticalAlignment','top',...
        'BackgroundColor','white',...
        'EdgeColor',[0.4 0.4 0.4],...
        'Margin',4,...
        'FontSize',12);

    %-----------------------------------------
    % Etiquetas (a), (b), (c), (d)
    %-----------------------------------------

    text(-0.13,1.08,...
        ['(' letters(i) ')'],...
        'Units','normalized',...
        'FontWeight','bold',...
        'FontSize',16);

end