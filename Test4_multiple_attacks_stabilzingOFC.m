clc; clear; close all;

% System Dynamics: x_dot = Ax + Bu, y = Cx
A = [-0.5 1; -1 -0.5];
B = [0; 1];
C = [1 0]; % Measuring only x1

% Solve LMI to find a stabilizing output-feedback gain L
setlmis([]); % Initialize LMI system
P = lmivar(1, [2 1]); % P: Positive definite matrix (2x2)
L = lmivar(2, [1 1]); % L: Output-feedback gain (1x1)

% Stability Constraint: A_cl = A - B*L*C should be Hurwitz
lmiterm([1 1 1 P], A, 1, 's');  % P*A + A'*P
lmiterm([1 1 1 L], -B, C, 's'); % -P*B*L*C - C'*L'*B'*P
lmiterm([1 1 1 0], -eye(2));    % -I (ensures P > 0)

% Solve the LMI problem
lmisys = getlmis;
[tmin, xfeas] = feasp(lmisys);
if tmin < 0
    P_sol = dec2mat(lmisys, xfeas, P);
    L_sol = dec2mat(lmisys, xfeas, L);
    L = L_sol; % Extract stabilizing gain
else
    error('LMI problem not feasible. Check system controllability/observability.');
end

% Observer Design (Luenberger Observer)
G = place(A', C', [-3 -4])'; % Observer poles (for stable estimation)

dt = 0.01; % Time step
T = 15; % Total simulation time
N = T/dt;
x = [1; 1]; % Initial true state
x_hat = [0; 0]; % Initial estimated state
u = 0; % Initialize control input

% Attack Time Windows
dos_start = 5; dos_end = 10; % DoS attack
dec_start = 8; dec_end = 12; % Deception attack

% Noise Parameters
noise_intensity_dos = 0.05;  % Increased state noise during DoS
noise_intensity_dec = 0.02;  % Increased sensor noise during deception

% Data Storage
x_hist = zeros(2, N);
x_hat_hist = zeros(2, N);
y_hist = zeros(1, N);
u_hist = zeros(1, N);
time = linspace(0, T, N);

% Simulation Loop
for k = 1:N
    t = time(k);
    
    % Compute true sensor output
    y = C*x; 
    
    % Apply deception attack to sensor output
    if t >= dec_start && t <= dec_end
        y = y + 0.5 + noise_intensity_dec * randn; % Bias + additional noise
    end
    
    % Observer update (Luenberger observer)
    x_hat = x_hat + dt * (A*x_hat + B*u + G*(y - C*x_hat));
    
    % Compute control input (output-feedback: u = -L * y)
    if t >= dos_start && t <= dos_end
        u = 0; % DoS attack: Control input is blocked
    else
        u = -L * y; % Uses measured output instead of full state
    end
    
    % Introduce disturbance (shaking effect during DoS)
    if t >= dos_start && t <= dos_end
        disturbance = noise_intensity_dos * randn(2,1); % More noise during DoS
    else
        disturbance = 0.01 * randn(2,1); % Small natural noise
    end
    
    % State update (Euler method)
    x = x + dt * (A*x + B*u) + disturbance;
    
    % Store data
    x_hist(:, k) = x;
    x_hat_hist(:, k) = x_hat;
    y_hist(k) = y;
    u_hist(k) = u;
end

% Plot State Trajectories with Attack Indications
figure;
hold on;

% Highlight DoS Attack Period (Red)
yLimits = ylim;
fill([dos_start dos_end dos_end dos_start], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none'); 

% Highlight Deception Attack Period (Green)
fill([dec_start dec_end dec_end dec_start], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
     'g', 'FaceAlpha', 0.15, 'EdgeColor', 'none');

% Plot True States with Noise
plot(time, x_hist(1,:), 'b', 'LineWidth', 1.5); % x1 (True)
plot(time, x_hist(2,:), 'r', 'LineWidth', 1.5); % x2

% Plot Estimated States
plot(time, x_hat_hist(1,:), 'b--', 'LineWidth', 1.2); % x1 (Estimated)
plot(time, x_hat_hist(2,:), 'r--', 'LineWidth', 1.2); % x2 (Estimated)

% Plot Sensor Output with Deception Attack
plot(time, y_hist, 'g-.', 'LineWidth', 1.5); % x1 (Manipulated sensor output)

% Labels and Legends
xlabel('Time (s)');
ylabel('State');
legend('DoS Attack', 'Deception Attack', ...
       'x_1 (True)', 'x_2 (True)', 'x_1 (Estimated)', 'x_2 (Estimated)', ...
       'x_1 (Sensor Output)');
title('LMI-Based Output-Feedback Stabilizing Controller');

% Add Attack Annotations
text((dos_start + dos_end)/2, yLimits(2)*0.8, 'DoS Attack', 'Color', 'r', 'FontSize', 12, ...
     'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text((dec_start + dec_end)/2, yLimits(2)*0.6, 'Deception Attack', 'Color', 'g', 'FontSize', 12, ...
     'FontWeight', 'bold', 'HorizontalAlignment', 'center');

hold off;