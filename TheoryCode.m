%% Code to model the theory behind the mechanism
clear; close all; clc;

%% Givens / Constants
r = 0.075;    % [m]
l = 0.26;     % [m]
d = 0.155;    % [m]

theta = linspace(0, 360*6, 1000);  % [deg]
w = [1,2,3,4,5,6];                 % [deg/s]

%% Modeling the Equations
v_mod = zeros(length(theta), length(w));

for i = 1:length(w)
    v_mod(:, i) = LCSMODEL(r, d, l, theta, w(i));
end

%% Plotting the Model Output
figure;
subplot(2,3,1)
hold on
plot(theta, v_mod(:,1));
xlabel('\theta (deg)');
ylabel('Collar Velocity (cm/s)');
title('V');
grid on;
hold off

subplot(2,3,2)
hold on
plot(theta, v_mod(:,2));
xlabel('\theta (deg)');
ylabel('Collar Velocity (cm/s)');
title('V');
grid on;
hold off

subplot(2,3,3)
hold on
plot(theta, v_mod(:,3));
xlabel('\theta (deg)');
ylabel('Collar Velocity (cm/s)');
title('V');
grid on;
hold off

subplot(2,3,4)
hold on
plot(theta, v_mod(:,4));
xlabel('\theta (deg)');
ylabel('Collar Velocity (cm/s)');
title('V');
grid on;
hold off

subplot(2,3,5)
hold on
plot(theta, v_mod(:,5));
xlabel('\theta (deg)');
ylabel('Collar Velocity (cm/s)');
title('V');
grid on;
hold off

subplot(2,3,6)
hold on
plot(theta, v_mod(:,6));
xlabel('\theta (deg)');
ylabel('Collar Velocity (cm/s)');
title('V');
grid on;
hold off
%% Function Definition
function v_B = LCSMODEL(r, d, l, theta, w)
    beta = asind((d - r * sind(theta)) / l);
    v_B = -w * ( r * (sind(theta) + cosd(theta) .* tand(beta)));
    
    % Convert the velocity from m/s to cm/s.
    v_B = v_B * 100;
end
