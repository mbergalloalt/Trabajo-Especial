clear
close all;

gamma = 42.58e6; % en H/2piT
G = 7; % en Tm
gamma = 2*pi*gamma;
%---------------------------------------------
NE=10000;
TE=100e-6;  % en s
T2a=100e-3; % en s
T2b=50e-3;
A=1;
B=0;

%---------------------------------------------
% Vector temporal (columna)
t=linspace(TE,NE*TE,NE)';       % NE x 1

% Señal sin difusión
data=A*exp(-t/T2a) + B*exp(-t/T2b); 

% --- Lista de coeficientes de difusión que queremos barrer ---
D_list = logspace(-11, -8, 30);


% Inversion numerica
alpha = 1E-2;            
Nx=300;                  % número de puntos en el eje T2 (espacio de solución) 
T = logspace(-2, 0, Nx); % grid de tiempos T2 (0.01 s a 1 s)
 

% Construccion del kernel de Laplace
d=length(t);
for i = 1:d
      K1(i,:)=  exp(-t(i) * (1 ./ T));
end

%--- Prealocacion de resultados---
nD = numel(D_list);   %calcula cuantos valores de dif voy a barrer --> nD=4
S_all = zeros(Nx,nD); %prealoca matriz Nx x nD para guardar la distribucion
                      %de T2 obtenidas por inversion de Laplace, Cada columna 
                      %S_all(:,k) corresponderá a la distribución S(T) para
                      %el valor de difusión D_list(k).


% --- Bucle principal: para cada D, calcular atenuación por difusión, invertir ---
for k = 1:nD
    D = D_list(k);
    dataD = exp( - (1/12) * (gamma * G * TE)^2 * D .* t );
    dataT = data .* dataD;
    Z = dataT;

    [S, resida] = flint1D(K1, Z, alpha);
    S_all(:,k) = S(:);
end

figure;
[X, Y] = meshgrid(T, D_list);
waterfall(X,Y,S_all')
box off                % saca el marco del gráfico
ax = gca;
ax.Color = 'w';     % hace transparente el fondo
hx = xlabel('T_2(s)', 'FontSize', 16);
hy = ylabel('D(m^2/s)','FontSize', 16); 
% set(hx,'Rotation', 8)
% set(hy,'Rotation', -51)
set(gca, 'XScale', 'log');
set(gca,'YDir', 'reverse', 'YScale', 'log') % me falto ponerlo en escala log
xlim([0 0.2])
grid off;





