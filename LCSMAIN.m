%% Code to model the theory behind the mechanism
clear; close all; clc;

%% Givens / Constants
r = 0.075;    % [m]
l = 0.26;     % [m]
d = 0.155;    % [m]

theta = linspace(0, 360*7, 1000);  % [deg]

%% Modeling the Equations

n = 6; % Number of cycles

for i = 1:6
    % Data
    [theta_exp,w_exp,v_exp,time] = LCSDATA(i);

    % Model
    v_mod = LCSMODEL(r, d, l, theta_exp, w_exp);

    % Residual
    res_exp = transpose(v_exp) - v_mod;
    res_mod = v_mod - v_mod;
    std_res = std(res_exp);
    abs_res_exp = transpose(abs(v_exp)) - abs(v_mod);
    abs_res_mod = abs(v_mod) - abs(v_mod);
    abs_std_res = std(abs_res_exp);

    % Plotting
    PlotName = num2str(i+4) + ".5 V";

    figure(1)
    subplot(2,3,i)
    hold on
    plot(theta_exp,v_exp)
    plot(theta_exp, v_mod)
    xlim([0,360*n])
    ylim([-300 300])
    xlabel("\theta (deg)")
    ylabel("Collar Velocity (cm/s)")
    legend("Velocity From Data","Velocity From Model","Location","southoutside")
    title(PlotName)
    
    figure(2)
    subplot(2,3,i)
    hold on
    plot(time,res_exp)
    plot(time,res_mod)
    yline(std_res)
    xlim([min(time) max(time)])
    ylim([-40 40])
    xlabel("Time (s)")
    ylabel("Collar Velocity Residual(cm/s)")
    legend("Data Residual","Model (Baseline)","std. of Residual","Location","southoutside")
    title(PlotName)

    figure(3)
    subplot(2,3,i)
    hold on
    plot(time,abs_res_exp)
    plot(time,abs_res_mod)
    yline(abs_std_res)
    xlim([min(time) max(time)])
    ylim([-40 40])
    xlabel("Time (s)")
    ylabel("Abs. Collar Velocity Residual(cm/s)")
    legend("Data Abs. Residual","Model (Baseline)","std. of Abs. Residual","Location","southoutside")
    title(PlotName)

    figure(4)
    subplot(2,3,i)
    hold on
    plot(time,res_exp)
    plot(time,abs_res_exp)
    xlim([min(time) max(time)])
    ylim([-40 40])
    xlabel("Time (s)")
    ylabel("Collar Velocity Residual(cm/s)")
    legend("Data Residual","Data Abs. Residual","Location","southoutside")
    title(PlotName)

    figure(5)
    hold on
    plot(time,res_exp)
    if i == 1
    xlim([min(time) max(time)])
    end
    ylim([-40 40])
    xlabel("Time (s)")
    ylabel("Collar Velocity Residual(cm/s)")

    figure(6)
    hold on
    plot(theta_exp,res_exp)
    xlim([0,360*n])
    ylim([-40 40])
    xlabel("\theta (deg)")
end

figure(5)
legend("5.5v","6.5v","7.5v","8.5v","9.5v","10.5v","Location","eastoutside")
title("Time vs Residuals")

figure(6)
legend("5.5v","6.5v","7.5v","8.5v","9.5v","10.5v","Location","eastoutside")
title("Theta vs Residuals")


figure(1)
print('WheelPosCollarVel', '-dpng','-r300')

figure(2)
print('CollarVelResiduals', '-dpng','-r300')

figure(3)
print('AbsoluteCollarVelResiduals', '-dpng','-r300')

figure(4)
print('TimevBothSplot', '-dpng','-r300')

figure(5)
print('TimevBothTogether', '-dpng','-r300')

figure(6)
print('WheelPosvBoth', '-dpng','-r300')

%% Model Function
function v_B = LCSMODEL(r, d, l, theta, w)
    w = w * (pi/180); % converts w to rad/s

    % calculating beta & velocity of collar
    for j = 1:length(theta)
    beta(j) = asind((d - r * sind(theta(j))) / l);
    v_B(j) = ( -w(j) * ( r * (sind(theta(j)) + cosd(theta(j)) * tand(beta(j)))) ) * 100; % also converts m/s to cm/s
    end
end

%% Data Function
function [theta_exp,w_exp,v_exp,time] = LCSDATA(i)
    % Read in data
    TestName = "Test1_" + num2str(i+4) + "pt5V";
    Data = readmatrix(TestName);

    % Finds where slide position is 1st at a minimum
    TF = islocalmin(Data(:,3));    
    for j = 1:length(Data(:,1))
        if TF(j) == 1
            break
        end
    end

    % Re-index & zeroing data with post processing
    Data(:,2) = Data(:,2) - Data(j,2) + 152.5 - 360;
    TF2 = islocalmin(abs(Data(:,2)));    
    for k = 1:length(Data(:,1))
        if TF2(k) == 1
            break
        end
    end

    zeroedData = Data(k:length(Data(:,1)),:);   
    zeroedData(:,1) = zeroedData(:,1) - zeroedData(1,1);
    zeroedData(:,5) = zeroedData(:,5)/10; % mm to cm
    theta_exp = zeroedData(:,2);
    theta_exp = theta_exp(theta_exp <= 6*360); % deg
    zeroedData = zeroedData(1:length(theta_exp),:);
    v_exp = zeroedData(:,5); % cm/s
    time = zeroedData(:,1); % seconds
    w_exp = zeroedData(:,4); % deg/s
end