function [m_hist, p_hist, a_hist, Ti_hist] = localFE3(m_curr,p_curr,a_curr,par,dx,dt,T,save_interval)
D1 = par(1);
D2 = par(2);
k1 = par(3);
k2 = par(4);
k3 = par(5);
k4 = par(6);
alpha = par(7);
S = size(m_curr,1); 
m_next = zeros(S,1);
p_next = zeros(S,1);
a_next = zeros(S,1);
num_saved_steps = floor(T / save_interval) + 1;

% Historical matrix
m_hist = zeros(S, num_saved_steps);
p_hist = zeros(S, num_saved_steps);
a_hist =zeros(S, num_saved_steps);
Ti_hist = zeros(1, num_saved_steps); % Time saved vector

m_hist(:,1) = m_curr;
p_hist(:,1) = p_curr;
a_hist(:,1) = a_curr;
Ti_hist(1) = 0;
save_idx = 2;


tic;
for j = 1:T-1 % loop in time    
% left boundary space
% u(i-1) = u(S-1) when i = 1
    i = 1; 
    m_next(i) = m_curr(i) + dt*(m_curr(S-1)-2*m_curr(i)+m_curr(i+1))/(dx^2)... 
                    - alpha * dt * ((m_curr(i+1)+m_curr(i))*(p_curr(i+1)-p_curr(i))-(m_curr(i)+m_curr(S-1))*(p_curr(i)-p_curr(S-1)))/(2*(dx^2))...
                    + dt * p_curr(i)*(1-m_curr(i)) - dt * m_curr(i);
        p_next(i) = p_curr(i) + dt * D1 * (p_curr(S-1) - 2*p_curr(i) + p_curr(i+1))/(dx^2)... 
                    + dt * k1 * m_curr(i) / ( 1+ a_curr(i)) - dt * k2 * p_curr(i);
        a_next(i) = a_curr(i) + dt * D2 * (a_curr(S-1) -2*a_curr(i) +a_curr(i+1))/(dx^2)... 
                    + dt * k3 * m_curr(i) - dt* k4 * a_curr(i);

% central part space
    for i = 2:S-1
        m_next(i) = m_curr(i) + dt*(m_curr(i-1)-2*m_curr(i)+m_curr(i+1))/(dx^2)... 
                    - alpha * dt * ((m_curr(i+1)+m_curr(i))*(p_curr(i+1)-p_curr(i))-(m_curr(i)+m_curr(i-1))*(p_curr(i)-p_curr(i-1)))/(2*(dx^2))...
                    + dt * p_curr(i)*(1-m_curr(i)) - dt * m_curr(i);
        p_next(i) = p_curr(i) + dt * D1 * (p_curr(i-1) - 2*p_curr(i) + p_curr(i+1))/(dx^2)... 
                    + dt * k1 * m_curr(i) / ( 1+a_curr(i)) - dt * k2 * p_curr(i);
        a_next(i) = a_curr(i) + dt * D2 * (a_curr(i-1) -2*a_curr(i) +a_curr(i+1))/(dx^2)... 
                    + dt * k3 * m_curr(i) - dt* k4 * a_curr(i);
    end

% right boundary space
% u(i+1) = u(2) when i = S
    i = S; 
    m_next(i) = m_curr(i) + dt*(m_curr(i-1)-2*m_curr(i)+m_curr(2))/(dx^2)... 
                    - alpha * dt * ((m_curr(2)+m_curr(i))*(p_curr(2)-p_curr(i))-(m_curr(i)+m_curr(i-1))*(p_curr(i)-p_curr(i-1)))/(2*(dx^2))...
                    + dt * p_curr(i)*(1-m_curr(i)) - dt * m_curr(i);
        p_next(i) = p_curr(i) + dt * D1 * (p_curr(i-1) - 2*p_curr(i) + p_curr(2))/(dx^2)... 
                    + dt * k1 * m_curr(i) / ( 1+a_curr(i)) - dt * k2 * p_curr(i);
        a_next(i) = a_curr(i) + dt * D2 * (a_curr(i-1) -2*a_curr(i) +a_curr(2))/(dx^2)... 
                    + dt * k3 * m_curr(i) - dt* k4 * a_curr(i);


% State update
    m_curr = m_next;
    p_curr = p_next;
    a_curr = a_next;
   

% State saving every save_interval step
    if mod(j, save_interval) == 0 || j == T-1
        m_hist(:, save_idx) = m_curr;
        p_hist(:, save_idx) = p_curr;
        a_hist(:,save_idx) =a_curr;
        Ti_hist(save_idx) = j * dt;
        save_idx = save_idx + 1;
    end

end
toc;




