%% Code to model the theory behind the mechanism
clear; close all; clc;

%% Givens / Constants
r = 0.075;    % [m]
l = 0.26;     % [m]
d = 0.155;    % [m]

theta = linspace(0, 360*6, 1000);  % [deg]
%w = [1,2,3,4,5,6]; % deg/s


%% Modeling the Equations

n = 6;

for i = 1:n
    % Data
    [theta_exp,w_exp,v_exp,time] = LCSDATA(i);
    % Model
    w = mean(w_exp)*(pi/180);
    v_mod(:, i) = LCSMODEL(r, d, l, theta, w);
    % Plotting
    figure(1)

    subplot(2,3,i)
    hold on
    plot(theta_exp,v_exp)
    plot(theta, v_mod(:,i))
    xlim([0,360*n])
    xlabel("\theta (deg)")
    ylabel("Collar Velocity (cm/s)")
    legend("Velocity From Data","Velocity From Model","Location","southoutside")
    PlotName = num2str(i+4) + ".5 V";
    title(PlotName)
    
end

%% Model Function
function v_B = LCSMODEL(r, d, l, theta, w)
    beta = asind((d - r * sind(theta)) / l);
    v_B = -w * ( r * (sind(theta) + cosd(theta) .* tand(beta)));
    
    % Convert the velocity from m/s to cm/s.
    v_B = v_B * 100;
end

%% Data Function
function [theta_exp,w_exp,v_exp,time] = LCSDATA(i)
    % Read in data
    TestName = "Test1_" + num2str(i+4) + "pt5V";
    StructName = ['T',num2str(i+4),'pt5'];
    Data.(StructName) = readmatrix(TestName);

    % Finds where slide position is 1st at a minimum
    TF = islocalmax(Data.(StructName)(:,3));    
    for j = 1:length(Data.(StructName)(:,1))
        if TF(j) == 1
            break
        end
    end

    % Re-index & zeroing data with post processing
    zeroedData.(StructName) = Data.(StructName)(j:length(Data.(StructName)(:,1)),:);   
    zeroedData.(StructName)(:,1) = zeroedData.(StructName)(:,1) - zeroedData.(StructName)(1,1);
    zeroedData.(StructName)(:,2) = zeroedData.(StructName)(:,2) - zeroedData.(StructName)(1,2);
    zeroedData.(StructName)(:,5) = zeroedData.(StructName)(:,5)/10; % mm to cm
    theta_exp = zeroedData.(StructName)(:,2);
    theta_exp = theta_exp(theta_exp <= 6*360); % deg
    zeroedData.(StructName) = zeroedData.(StructName)(1:length(theta_exp),:);
    v_exp = zeroedData.(StructName)(:,5); % cm/s
    time = zeroedData.(StructName)(:,1); % seconds
    w_exp = zeroedData.(StructName)(:,4); % deg/s
end