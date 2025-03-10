function numerical_ode_extended(h)
% Given ODE: dy/dt = -50*y + sin(t), y(0) = 1

t0 = 0;
tf = 2;
y0 = 1;

% Time vector
t = t0:h:tf;
n = length(t);

% Initialize solution vectors
y_FE = zeros(n,1);
y_ME = zeros(n,1);
y_BE = zeros(n,1);
y_RK4 = zeros(n,1);
y_AB2 = zeros(n,1);
y_AM2 = zeros(n,1);

y_FE(1) = y0;
y_ME(1) = y0;
y_BE(1) = y0;
y_RK4(1) = y0;
y_AB2(1) = y0;
y_AM2(1) = y0;

% Forward Euler
for i = 1:n-1
    y_FE(i+1) = y_FE(i) + h * (-50*y_FE(i) + sin(t(i)));
end

% Modified Euler
for i = 1:n-1
    f1 = -50*y_ME(i) + sin(t(i));
    y_pred = y_ME(i) + h * f1;
    f2 = -50*y_pred + sin(t(i+1));
    y_ME(i+1) = y_ME(i) + (h/2) * (f1 + f2);
end

% Backward Euler
for i = 1:n-1
    y_BE(i+1) = (y_BE(i) + h*sin(t(i+1))) / (1 + 50*h);
end

% Runge-Kutta 4th Order (RK4)
for i = 1:n-1
    k1 = -50*y_RK4(i) + sin(t(i));
    k2 = -50*(y_RK4(i) + 0.5*h*k1) + sin(t(i) + 0.5*h);
    k3 = -50*(y_RK4(i) + 0.5*h*k2) + sin(t(i) + 0.5*h);
    k4 = -50*(y_RK4(i) + h*k3) + sin(t(i) + h);
    y_RK4(i+1) = y_RK4(i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end

% Adams-Bashforth 2-step (AB2)
y_AB2(2) = y_RK4(2); % Use RK4 to get the second value
for i = 2:n-1
    f_n = -50*y_AB2(i) + sin(t(i));
    f_n1 = -50*y_AB2(i-1) + sin(t(i-1));
    y_AB2(i+1) = y_AB2(i) + (h/2) * (3*f_n - f_n1);
end

% Adams-Moulton 2-step (AM2)
y_AM2(2) = y_RK4(2); % Use RK4 to get the second value
for i = 2:n-1
    f_n = -50*y_AM2(i) + sin(t(i));
    f_n1 = -50*y_AM2(i-1) + sin(t(i-1));
    y_AM2(i+1) = y_AM2(i) + (h/12) * (5*f_n + 8*f_n1 - f_n1);
end

% Plot results
figure;
plot(t, y_FE, 'r', 'LineWidth', 1.5); hold on;
plot(t, y_ME, 'g', 'LineWidth', 1.5);
plot(t, y_BE, 'b', 'LineWidth', 1.5);
plot(t, y_RK4, 'm', 'LineWidth', 1.5);
plot(t, y_AB2, 'c', 'LineWidth', 1.5);
plot(t, y_AM2, 'k', 'LineWidth', 1.5);
legend('Forward Euler', 'Modified Euler', 'Backward Euler', 'RK4', 'Adams-Bashforth 2-step', 'Adams-Moulton 2-step');
xlabel('Time t'); ylabel('y(t)');
title('Numerical Solutions of ODE using Euler, Runge-Kutta, and Multi-step Methods');
grid on;
end
