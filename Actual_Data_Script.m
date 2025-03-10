clear
clc
close all

%function [theta_exp,w_exp,v_exp,time] = LCSDATA(model)

for i = 1:6
    % Read in data
    TestName = "Test1_" + num2str(i+4) + "pt5V";
    StructName = ['T',num2str(i+4),'pt5'];
    Data.(StructName) = readmatrix(TestName);

    % Finds where slide position is 1st at a minimum
    TF = islocalmin(Data.(StructName)(:,3));    
    for j = 1:length(Data.(StructName)(:,1))
        if TF(j) == 1
            break
        end
    end

    n = 6; % number of cycles

    % Re-index & zeroing data with post processing
    zeroedData.(StructName) = Data.(StructName)(j:length(Data.(StructName)(:,1)),:);   
    zeroedData.(StructName)(:,1) = zeroedData.(StructName)(:,1) - zeroedData.(StructName)(1,1);
    zeroedData.(StructName)(:,2) = zeroedData.(StructName)(:,2) - zeroedData.(StructName)(1,2);
    zeroedData.(StructName)(:,5) = zeroedData.(StructName)(:,5)/1000;
    theta_exp = zeroedData.(StructName)(:,2);
    theta_exp = theta_exp(theta_exp <= n*360); % deg
    zeroedData.(StructName) = zeroedData.(StructName)(1:length(theta_exp),:);
    v_exp = zeroedData.(StructName)(:,5); % m/s
    time = zeroedData.(StructName)(:,1); % seconds
    w_exp = zeroedData.(StructName)(:,4); % deg/s
     
    % Graphing
    figure(1)
    hold on
    subplot(2,3,i)
    plot(theta_exp,v_exp)
    xlim([0,360*n])
    xlabel("\theta (deg)")
    ylabel("Collar Velocity (m/s)")
    legend("Velocity From Data","Velocity From Model","Location","southoutside")
    PlotName = num2str(i+4) + ".5 V";
    title(PlotName)
end

%end
