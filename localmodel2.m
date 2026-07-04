clear
clc

%% Parameters
D = 1; 
ni = 2; 
k2 = 1;
L = 15; % [-L,L] domain


%% Check if exactly one wave number falls into the TIR
figure()
hold on
box on
set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'GridColor', 'k');
set(gca, 'MinorGridColor', 'k');
plot([0 10],[0 0],'k:')

n = 1:20;  
[alphaT_local, k_loc] = critvaluelocal2(2*L, D, k2, ni, n);
alpha = alphaT_local+0.03;
DetMk=@(kk) D*kk.^4+(k2 - alpha*(ni-k2) + D*ni/k2).*kk.^2 - k2+ni;

kk=0:0.00001:5;
plot(kk,DetMk(kk));
axis([0 2 -0.45 0.8])
xlabel('k')
ylabel('Det$(M_k^l)$','interpreter','latex')
xline_h = yline(0, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Asse k (Det=0)');

% Position of some wave numbers
nn=0:20;
kk2_L=(2*pi*nn/(2*L)); 
plot(kk2_L,0*kk2_L,'k*', 'MarkerFaceColor', 'k')

%% Discretisation
dx = 0.05; % stepsize
S = 2*L/dx+1; % number of space steps
Tf = 4000; % final time
dt = 0.0005; % timestep
T = Tf/dt+1 ; % number of timesteps
xi=0:dx:2*L; 
Ti=0:dt:Tf;
save_interval = 500;
% Check CFL condition (CFL<1)
CFL1=2*D*dt/(dx^2);
fprintf('CFL1 %d \n', CFL1);


%% m and p initialization  

% Homogeneous steady state
m_eq = (-k2+ni)/ni;
p_eq = (-k2+ni)/k2;
par(1) =D;
par(2) = ni;
par(3) = k2;
par(4) = alpha;

fprintf('m_eq %d \n', m_eq);
fprintf('p_eq %d \n', p_eq);

% Current state vectors
m_curr = zeros(S,1);
p_curr = zeros(S,1);

% Initial conditions
for i = 1:S 
    m_curr(i) = m_eq + 0.05*rand; 
    p_curr(i) = p_eq ;
end

%% FTCS
[m_hist, p_hist, Ti_hist] = localFE2(m_curr,p_curr,par,dx,dt,T,save_interval);

%% Spatial-temporal plot
figure('Color','w');
subplot(1,2,1);
imagesc(xi, Ti_hist, m_hist');   
set(gca,'YDir','normal', 'XColor', 'k', 'YColor', 'k');
colormap(flipud(hot)); cb = colorbar; cb.Color = 'k';
xlabel('$x$ ','Interpreter','latex','FontSize',12);
ylabel('$t$ ','Interpreter','latex','FontSize',12);
title('$M(x,t)$','Interpreter','latex','FontSize',14, 'Color', 'k');
shading interp;

subplot(1,2,2);
imagesc(xi, Ti_hist, p_hist');   
set(gca,'YDir','normal', 'XColor', 'k', 'YColor', 'k');
colormap(flipud(hot)); cb = colorbar; cb.Color = 'k';
xlabel('$x$ ','Interpreter','latex','FontSize',12);
ylabel('$t$ ','Interpreter','latex','FontSize',12);
title('$P(x,t)$','Interpreter','latex','FontSize',14, 'Color', 'k');
shading interp;

