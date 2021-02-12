% least squares test

% x = [1;2;3;4];
% y = [5;6;7;8];
% 
% N = length(x);
% 
% 
% del = N*sum(x.^2)-(sum(x))^2;
% A = (sum(x.^2)+sum(y.^2)-sum(x)*sum(x.*y))/del;
% B = (N*sum(x.*y)-sum(x)*sum(y))/del;


% [p, S] = polyfit(input2, input1, 1);
% [y, delta] = polyval(p, input2, S);



%     [avgVelo1, errorEvlo1] = LSFunc(input2, input1);
    %[y, delta] = polyval(avgVelo, input2, errorEvlo);
%     AverageVelocity = [AverageVelocity; avgVelo1];
%     ErrorVelocity = [ErrorVelocity; errorEvlo1];



%% Pitot Static Probe Velocity calculations 
% for i = 2:1:60
%     input1 = PPPTVel(PPPTdata((500 * (i-1) + 1):500 * i, :));
%     input2 = PPPTdata((500 * (i-1) + 1):500 * i, end);
%     Velocity(i-1) = mean(input1);
%     Voltages1(i-1) = mean(input2);
% end
% 
% [bestFit,Error] = LSFunc(Voltages1',Velocity');
% 
% hold on
% plot(Voltages1, Velocity, 'o')
% plot(Voltages1,bestFit)
% errorbar(Voltages1,Velocity,Error);
% hold off


%% Ventury Tube Velocity calcuations 
% for i = 2:1:60
%     input1 = VTPTVel(VTPTdata((500 * (i-1) + 1):500 * i, :));
%     input2 = VTPTdata((500 * (i-1) + 1):500 * i, end);
%     Velocity_Venturi(i-1) = mean(input1);
%     Voltages_Venturi(i-1) = mean(input2);
% end


%% Old code
% OutVel1 = PPPTVel(PPPTdata(1:500,:));
% input0 = OutVel1;
% input3 = linspace(1,1.0001,500)';
% [avgVelo0, errorEvlo0] = LSFunc(input3, input0);
% AverageVelocity = [AverageVelocity; avgVelo0];
% ErrorVelocity = [ErrorVelocity; errorEvlo0];


%% Using Water Manometer 

water_manometer_data = readtable("ASEN 2002 Lab/Lab 2/Data collection/water_manometer_data.xlsx", 'VariableNamingRule', "preserve");
pitot_static_idx = ismember(water_manometer_data.("Was the manometer connected to the Venturi Tube or Pitot Static"), 'Pitot Static Probe'); % Logical
pitotStaticData = water_manometer_data(strcmp(water_manometer_data.("Was the manometer connected to the Venturi Tube or Pitot Static"), 'Pitot Static Probe'), :);

% Airspeed Calculation 

% function [airSpeed] = airSpeedCalc(airSpeedDiffPressure, R, Temp_Atm, Pressure_Atm, AreaRatio)
% 
% airSpeedDiffPressure = abs(airSpeedDiffPressure);
% airSpeed = sqrt((2*R*airSpeedDiffPressure.*Temp_Atm)./(PressureAtm*(1 - AreaRatio^2)));
% 
% end

%% Part 7 - Question 3

% temp = []; % some if statement to differentiate length between angles of attack
% 
% x  = ConfusedData;
% S  = numel(x);
% xx = reshape(x(1:S - mod(S, temp)), temp, []);
% y  = sum(xx, 1).' / temp;

