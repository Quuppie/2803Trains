%% Code to model the theory behind the mechanism
clear; close all; clc;

%% Givens / Constants
r = 0.075;    % [m]
l = 0.26;     % [m]
d = 0.155;    % [m]

theta = linspace(0, 360*7, 1000);  % [deg]
%w = [1,2,3,4,5,6]; % deg/s


%% Modeling the Equations

n = 6; % Number of cycles

for i = 1:6
    % Data
    [theta_exp,w_exp,v_exp,time] = LCSDATA(i);

    % Model
    v_mod = LCSMODEL(r, d, l, theta_exp, w_exp);

    % Residual
    res_exp = v_exp - v_exp;
    res_mod = v_mod - transpose(v_exp);
    std_res = std(res_mod);
    abs_res_exp = abs(v_exp) - abs(v_exp);
    abs_res_mod = abs(v_mod) - transpose(abs(v_exp));
    abs_std_res = std(abs_res_mod);

    % Plotting
    PlotName = num2str(i+4) + ".5 V";

    figure(1)
    subplot(2,3,i)
    hold on
    plot(theta_exp,v_exp)
    plot(theta_exp, v_mod)
    xlim([0,360*n])
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
    xlabel("Time (s)")
    ylabel("Collar Velocity Residual(cm/s)")
    legend("Data (Baseline)","Model Residual","std. of Residual","Location","southoutside")
    title(PlotName)

    figure(3)
    subplot(2,3,i)
    hold on
    plot(time,abs_res_exp)
    plot(time,abs_res_mod)
    yline(abs_std_res)
    xlim([min(time) max(time)])
    xlabel("Time (s)")
    ylabel("Abs. Collar Velocity Residual(cm/s)")
    legend("Data (Baseline)","Abs. Model Residual","std. of Abs. Residual","Location","southoutside")
    title(PlotName)

    figure(4)
    subplot(2,3,i)
    hold on
    plot(time,res_mod)
    plot(time,abs_res_mod)
    xlim([min(time) max(time)])
    xlabel("Time (s)")
    ylabel("Collar Velocity Residual(cm/s)")
    legend("Model Residual","Abs. Model Residual","Location","southoutside")
    title(PlotName)

    figure(5)
    hold on
    plot(time,res_mod)
    if i == 1
    xlim([min(time) max(time)])
    end
    xlabel("Time (s)")
    ylabel("Collar Velocity Residual(cm/s)")

    figure(6)
    hold on
    plot(theta_exp,res_mod)
    xlim([0,360*n])
    xlabel("\theta (deg)")
end

figure(5)
legend("5.5v","6.5v","7.5v","8.5v","9.5v","10.5v","Location","eastoutside")
title("Time vs Residuals")

figure(6)
legend("5.5v","6.5v","7.5v","8.5v","9.5v","10.5v","Location","eastoutside")
title("Theta vs Residuals")

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