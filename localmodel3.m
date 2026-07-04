clear
clc

%% Parameters
D1 = 1;
D2 = 1;
k1 = 2;
k2 = 1;
k3 = 1;
k4 = 2; 
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

n=1:20;
[alphaT_local, k_loc] = critvaluelocal3(2*L, D1, D2, k1, k2, k3, k4, n);
alpha = alphaT_local + 0.1;
q6 = @(kk) (D1*D2).*kk.^6;
q4 = @(kk) (((k1*k4 + k2*k3)*D1*D2)/(k2*(k3+k4)) + k2*D2 + k4*D1-(alpha*D2*k4*(k1-k2))/(k3+k4)).*kk.^4;
q2 = @(kk) ( (k4*(k1-k2)*k2*D2)/(k2*(k3+k4))+ (k4*(k1*k4+k2*k3)*D1)/(k2*(k3+k4)) + k2*k4 - (alpha*k4^2*(k1-k2)*(k1*k4+k2*k3))/(k1*(k3+k4)^2)).*kk.^2;
q0 = (k4*(k1*k4 + k2*k3)*(k1-k2))/(k1*(k3+k4));
Qk = @(kk) q6(kk)+q4(kk)+q2(kk)+q0;

kk=0:0.00001:5; 
plot(kk,Qk(kk))
axis([0 10 -1.4 1])
xlabel('k')
ylabel('Q(k)','interpreter','latex')
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
Ti = 0:dt:Tf; 
save_interval = 500;
% check CFL condition (CFL<1)
CFL1=2*D1*dt/(dx^2);
fprintf('CFL1 %d \n', CFL1);


%% m, p, a initialisation 
% Homogeneous steady state
m_eq = k4 * (k1-k2)/(k1*k4 + k2*k3);
p_eq = (k4/k2)*(k1-k2)/(k4+k3);
a_eq = k3 * (k1-k2)/(k1*k4 + k2*k3);
par = zeros(7);
par(1) = D1;
par(2) = D2;
par(3) = k1;
par(4) = k2;
par(5) = k3;
par(6) = k4;
par(7) = alpha;

fprintf('m_eq %d \n', m_eq);
fprintf('p_eq %d \n', p_eq);
fprintf('a_eq %d \n', a_eq);

% Current state vectors
m_curr = zeros(S,1); 
p_curr = zeros(S,1); 
a_curr = zeros(S,1); 

% Initial conditions
for i = 1:S
    m_curr(i) = m_eq + 0.5*rand; 
    p_curr(i) = p_eq;
    a_curr(i) = a_eq;
end


%% FTCS
[m_hist, p_hist, a_hist, Ti_hist] = localFE3(m_curr,p_curr,a_curr,par,dx,dt,T,save_interval);


%% Spatial-temporal plot
figure('Color','w');
subplot(1,3,1);
imagesc(xi, Ti_hist, m_hist');   
set(gca,'YDir','normal', 'XColor', 'k', 'YColor', 'k');
colormap(flipud(hot)); cb = colorbar; cb.Color = 'k';
xlabel('$x$ ','Interpreter','latex','FontSize',12);
ylabel('$t$ ','Interpreter','latex','FontSize',12);
title('$M(x,t)$','Interpreter','latex','FontSize',14, 'Color', 'k');
shading interp;

subplot(1,3,2);
imagesc(xi, Ti_hist, p_hist');   
set(gca,'YDir','normal', 'XColor', 'k', 'YColor', 'k');
colormap(flipud(hot)); cb = colorbar; cb.Color = 'k';
xlabel('$x$ ','Interpreter','latex','FontSize',12);
ylabel('$t$ ','Interpreter','latex','FontSize',12);
title('$P(x,t)$','Interpreter','latex','FontSize',14, 'Color', 'k');
shading interp;

subplot(1,3,3);
imagesc(xi, Ti_hist, a_hist');   
set(gca,'YDir','normal', 'XColor', 'k', 'YColor', 'k');
colormap(flipud(hot)); cb = colorbar; cb.Color = 'k';
xlabel('$x$ ','Interpreter','latex','FontSize',12);
ylabel('$t$ ','Interpreter','latex','FontSize',12);
title('$A(x,t)$','Interpreter','latex','FontSize',14, 'Color', 'k');
shading interp;

