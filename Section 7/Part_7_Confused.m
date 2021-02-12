%% Housekeeping
clear
clc
close all

ConfusedData = [];
for i = 3:1:29
    files = dir('/MATLAB Drive/ASEN 2002 Lab/Lab 2/Section 7/Aero Lab Airfoil Testing Data');
    long_name = strcat(files(i).folder,'/',files(i).name);
    ConfusedData  = [ConfusedData;load(long_name)];
end
% dum = ConfusedData(1620,:);
dum = ConfusedData(end,:);
ConfusedData = [ConfusedData;dum];

%Serate pressure data from the 4 specfic ports that we will be using to
%calc the pressure on port 11 from 
ConfusedData = sortrows(ConfusedData, 23, "ascend");
PressurePort10 = ConfusedData(:, 16); 
PressurePort8 = ConfusedData(:,14);
PressurePort12 = ConfusedData(:,18);
PressurePort14 = ConfusedData(:,20);
%Line of best fit for the two pressure lines. This is not how we calculate
%the pressure at port 11 but not deleting this cause I think it will be
%useful later on 
[pup,sup] = polyfit(PressurePort8,PressurePort10,1);
lineUpY = polyval(pup,PressurePort8);
[pdown,sdown] = polyfit(PressurePort12,PressurePort14,1);
lineDownY = polyval(pdown,PressurePort12);
    
slopeAvg = (0.5)*(pup(1) + pdown(1));
YIntAvg = (0.5)*(pup(2)+pdown(2));
yAvg = [slopeAvg,YIntAvg];
xAvg = (0.5).*(PressurePort12+PressurePort8);
meanLine = polyval(yAvg,xAvg);

hold on
plot(PressurePort8,lineUpY)
plot(PressurePort12,lineDownY)
plot(xAvg,meanLine)
hold off

%Not sure if this is needed for anything at the moment 
AvgMeanLine =[];
for i = 1:315
    AvgMeanLine(i) = mean(meanLine((i-1)*20+1:i*20));
end
PressurePort11 = [meanLine ConfusedData(:,23) ConfusedData(:,4)];

%% Collect all data from pressure ports into a new matrix allong with AOA and Velocity 
for i = 1:315 
    DataForCalcs(i,1) = mean(ConfusedData((i-1)*20+1:i*20,4));      %Airspeed
    DataForCalcs(i,2) = mean(ConfusedData((i-1)*20+1:i*20,23));     %Angle of Attack
    DataForCalcs(i,3) = mean(ConfusedData((i-1)*20+1:i*20,15));     %Port 10
    DataForCalcs(i,4) = mean(ConfusedData((i-1)*20+1:i*20,14));     %Port 8
    DataForCalcs(i,5) = mean(ConfusedData((i-1)*20+1:i*20,16));     %Port 12
    DataForCalcs(i,6) = mean(ConfusedData((i-1)*20+1:i*20,17));     %Port 14
end

%Calculate theoretical pressure at port 11 based on the average of the two
%last ports on the top and bottom of the airfoil 
AvgDataCalc = [];
for k = 1:315 
    SomeVal(k) = (DataForCalcs(k,3) + DataForCalcs(k,5) + DataForCalcs(k,4) + DataForCalcs(k,6))/4;
end
SomeVal = SomeVal';
%Create a new matrix with the data for port 11 and the coresponding
%velocity and angle of attack 
AvgDataCalc = [DataForCalcs(:,1) DataForCalcs(:,2) SomeVal];


%% Pressure Coeffiecents
PortsXLocations = [0;.175;.35;.7;1.05;1.4;1.75;2.1;2.8;3.5;2.8;2.1;1.4;1.05;.7;.35;.175];
% Xcell = cell(1,17);
PortsXLocations = PortsXLocations./3.5;
% for i =1:17
%      xLocVec = zeros(315,1);
%      Xcell{i} = xLocVec + PortsXLocations(i);
% end

%Coefficient of pressure
PressureMeansAll = zeros(315,28);
for k = 1:28
    for i = 1:315 
    PressureMeansAll(i,k) = mean(ConfusedData((i-1)*20+1:i*20,k));
    end
end
CofPress = cell(1,17);
for k =1:17
    if k < 10
        MatrixUse = zeros(315,1);
        for j = 1:315
            MatrixUse(j) = -(PressureMeansAll(j,k+6))/(PressureMeansAll(j,5));
        end
        CofPress{k} = MatrixUse;
    elseif k == 10
        CofPress{10} = SomeVal./(-1*PressureMeansAll(:,5)); 
    else
        MatrixUse = zeros(315,1);
        for j = 1:315
            MatrixUse(j) = -(PressureMeansAll(j,k+6))/(PressureMeansAll(j,5));
        end
        CofPress{k} = MatrixUse;
    end
end

%% Plotting

%Finding the idexes of airspeeds of 9 17 and 34
GraphThis = [];
for j=1:17
    Nine = zeros(105,2);
    SevenTeen = zeros(105,2);
    ThreeFour = zeros(105,2);
    T = CofPress{j};
    AoA = PressureMeansAll(:,23);
    for i=1:105
        Nine(i,1) = T(i*3-2);   
        Nine(i,2) = AoA(i*3-2);   
        SevenTeen(i,1) = T(i*3-1);
        SevenTeen(i,2) = AoA(i*3-1);
        ThreeFour(i,1) = T(i*3);
        ThreeFour(i,2) = AoA(i*3);
    end
    temp = [Nine,SevenTeen,ThreeFour];
    GraphThis = [GraphThis;temp];
end

%This distiguinshes between the multiple angles of attack that measurements
%are taken at. The numbers represent idexes. -15 goes from 1:4 and -14 from
%5:7
IdxVecGraphing = [];
for i=-15:15
    for j = 1:17
        if j == 1 
            test = mean(GraphThis((GraphThis(1:105,2)==i),:));
        else
            g = GraphThis((j-1)*105+1:j*105,:);
            test = mean(g((g(:,2)==i),:));
%             vec = find((GraphThis((j-1)*105:j*105,2)==i));
%             a = vec(end);
%             IdxVecGraphing = [IdxVecGraphing;a];  
        end
        IdxVecGraphing = [IdxVecGraphing;test];
    end
end

FinalM = [];
NegFifteen = [];
NegTen = [];
NegFive = [];
Zed = [];
PosFive = [];
PosTen = [];
PosFifteen = [];
for j =1:3
%     for i =-15:16
        t = IdxVecGraphing(:,j*2-1);
        t1 = t(IdxVecGraphing(:,2) == -15);
        t2 = t(IdxVecGraphing(:,2) == -10);
        t3 = t(IdxVecGraphing(:,2) == -5);
        t4 = t(IdxVecGraphing(:,2) == 0);
        t5 = t(IdxVecGraphing(:,2) == 5);
        t6 = t(IdxVecGraphing(:,2) == 10);
        t7 = t(IdxVecGraphing(:,2) == 15);
        NegFifteen = [NegFifteen t1];
        NegTen = [NegTen t2];
        NegFive = [NegFive t3];
        Zed = [Zed t4];
        PosFive = [PosFive t5];
        PosTen = [PosTen t6];
        PosFifteen = [PosFifteen t7];
%         NegFive = [NegFive ti];
%     end
%     NegFifteen = [NegFifteen; NegFive];
end


figure(2)
hold on
subplot(2,3,1)
    hold on
    plot(PortsXLocations,NegFifteen(:,1))
    plot(PortsXLocations,Zed(:,1))
    plot(PortsXLocations,PosFifteen(:,1))
    xlim([0 1])
    ylim([-1 2])
    yticks(-1:.5:2)
    yticklabels([1 .5 0 -.5 -1 -1.5 -2])
    title("Coefficent of Pressure at 9 [m/s]")
    xlabel("Ratio of Port Location to Chord")
    ylabel("Coefficeint of Pressure")
    legend("-15°","0°","15°")
    xticks(0:.25:1)
    xticklabels([0 .25 .5 .75 1])
    hold off
subplot(2,3,4)
    hold on
    plot(PortsXLocations,NegTen(:,1))
    plot(PortsXLocations,NegFive(:,1))
    plot(PortsXLocations,PosFive(:,1))
    plot(PortsXLocations,PosTen(:,1))
    xlim([0 1])
    ylim([-1 2.5])
    yticks(-1:.5:2.5)
    yticklabels([1 .5 0 -.5 -1 -1.5 -2 -2.5])
    title("Coefficent of Pressure at 9 [m/s]")
    xlabel("Ratio of Port Location to Chord")
    ylabel("Coefficeint of Pressure")
    legend("-10°","-5°","5°","10°")
    xticks(0:.25:1)
    xticklabels([0 .25 .5 .75 1])
    hold off
subplot(2,3,2)
    hold on
    plot(PortsXLocations,NegFifteen(:,2))
    plot(PortsXLocations,Zed(:,2))
    plot(PortsXLocations,PosFifteen(:,2))
    xlim([0 1])
    ylim([-1 2])
    title("Coefficent of Pressure at 17 [m/s]")
    xlabel("Ratio of Port Location to Chord")
    ylabel("Coefficeint of Pressure")
    yticks(-1:.5:2)
    yticklabels([1 .5 0 -.5 -1 -1.5 -2])
    legend("-15°","0°","15°")
    xticks(0:.25:1)
    xticklabels([0 .25 .5 .75 1])
    hold off
subplot(2,3,5)
    hold on
    plot(PortsXLocations,NegTen(:,2))
    plot(PortsXLocations,NegFive(:,2))
    plot(PortsXLocations,PosFive(:,2))
    plot(PortsXLocations,PosTen(:,2))
    xlim([0 1])
    ylim([-1 2])
    title("Coefficent of Pressure at 17 [m/s]")
    xlabel("Ratio of Port Location to Chord")
    ylabel("Coefficeint of Pressure")
    legend("-10°","-5°","5°","10°")
    yticks(-1:.5:2)
    yticklabels([1 .5 0 -.5 -1 -1.5 -2])
    xticks(0:.25:1)
    xticklabels([0 .25 .5 .75 1])
    hold off
subplot(2,3,3)
    hold on
    plot(PortsXLocations,NegFifteen(:,3))
    plot(PortsXLocations,Zed(:,3))
    plot(PortsXLocations,PosFifteen(:,3))
    xlim([0 1])
    ylim([-1 3])
    yticks(-1:.5:3)
    yticklabels([1 .5 0 -.5 -1 -1.5 -2 -2.5 -3])
    title("Coefficent of Pressure at 34 [m/s]")
    xlabel("Ratio of Port Location to Chord")
    ylabel("Coefficeint of Pressure")
    legend("-15°","0°","15°")
    xticks(0:.25:1)
    xticklabels([0 .25 .5 .75 1])
    hold off
subplot(2,3,6)
    hold on
    plot(PortsXLocations,NegTen(:,3))
    plot(PortsXLocations,NegFive(:,3))
    plot(PortsXLocations,PosFive(:,3))
    plot(PortsXLocations,PosTen(:,3))
    legend("-10°","-5°","5°","10°")
    xlim([0 1])
    ylim([-1 3.5])
    yticks(-1:.5:3.5)
    yticklabels([1 .5 0 -.5 -1 -1.5 -2 -2.5 -3 -3.5])
    xticks(0:.25:1)
    xticklabels([0 .25 .5 .75 1])
    title("Coefficent of Pressure at 34 [m/s]")
    xlabel("Ratio of Port Location to Chord")
    ylabel("Coefficeint of Pressure")
    hold off
    sgtitle('Airfoil Static Pressure Coefficient Distribution')
hold off


%% Part 4


% Calculate the normal force coefficient and the axial force coefficient. Then, using these values, calculate the
%coefficient of lift and the coefficient of drag for the airfoil. Details for these four calculations can be found in
%Appendix E.

% Normal Force = -Sum of Forces y

%% getting Pressure Data

xDist = PortsXLocations*3.5;
yDist = [0.14665;0.33075;0.4018;0.476;0.49;0.4774;0.4403;0.38325;0.21875;0;0;0;0;0;0.0014;0.0175;0.03885];
Theta = -15:15;
PressureAll = zeros(315,29);
c = 3.5;
for k =1:29
    if k < 10
        for j = 1:315
            PressureAll(j,k) = PressureMeansAll(j,k);
        end
    elseif k == 16
        PressureAll(:,16) = SomeVal; 
    else
        for j = 1:315
            PressureAll(j,k) = PressureMeansAll(j,k-1);
        end
    end
end

%% getting Scanivalve Pressure 

ScanP = PressureAll(:,7:23);

%% Finding Coefficient of Pressure for each Port

pinf = PressureAll(:,3);
Vinf = PressureAll(:,4);
qinf = PressureAll(:,5);
Cp9 = IdxVecGraphing(:,1);
Cp17 = IdxVecGraphing(:,3);
Cp34 = IdxVecGraphing(:,5);
AoA2 = IdxVecGraphing(:,2);

AvgCp = zeros(31,3);
for j=1:3
    for i = 1:31
        if i==1
            AvgCp(1,j) = mean(IdxVecGraphing(1:17,j*2-1));
        else
             AvgCp(i,j) = mean(IdxVecGraphing(((i-1)*17+1:i*17),j*2-1));
        end
    end
end
%% Finding Normal Force and Axial Force Coefficients

% Cp in this equation is Coefficient of Pressure i is for each port

Ca9 = zeros(17,31); 
Ca17 = zeros(17,31);
Ca34 = zeros(17,31);
Cn9 = zeros(17,31); 
Cn17 = zeros(17,31);
Cn34 = zeros(17,31);

temp1 = zeros(17,31);
temp2 = zeros(17,31);
temp3 = zeros(17,31);

for i = 1:31
    temp1(:,i) = Cp9(((i-1)*17)+1:(17*i),1);
    temp2(:,i) = Cp17(((i-1)*17)+1:(17*i),1);
    temp3(:,i) = Cp34(((i-1)*17)+1:(17*i),1);
end

Cp9 = temp1;
Cp17 = temp2;
Cp34 = temp3;

for j = 1:31
    for i = 1:16
        Cn9(i,j) = 0.5*(Cp9(i,j)+Cp9(i+1,j))*(xDist(i+1)-xDist(i))/c; 
        Ca9(i,j) = -0.5*(Cp9(i,j)+Cp9(i+1,j))*(yDist(i+1)-yDist(i))/c; 
        
        Cn17(i,j) = 0.5*(Cp17(i,j)+Cp17(i+1,j))*(xDist(i+1)-xDist(i))/c; 
        Ca17(i,j) = -0.5*(Cp17(i,j)+Cp17(i+1,j))*(yDist(i+1)-yDist(i))/c;
        
        Cn34(i,j) = 0.5*(Cp34(i,j)+Cp34(i+1,j))*(xDist(i+1)-xDist(i))/c; 
        Ca34(i,j) = -0.5*(Cp34(i,j)+Cp34(i+1,j))*(yDist(i+1)-yDist(i))/c; 
    end
end

Cn9 = sum(Cn9);
Cn17 = sum(Cn17);
Cn34 = sum(Cn34);
Ca9 = sum(Ca9);
Ca17 = sum(Ca17);
Ca34 = sum(Ca34);

%% Calculating Coeff of Lift and Drag
AoA2 = -15:15;
for i =1:31
    
    Cl9(i) = (Cn9(i)*cosd(AoA2(i)))+(Ca9(i)*sind(AoA2(i)));
    Cd9(i) = (Cn9(i)*sind(AoA2(i)))+(Ca9(i)*cosd(AoA2(i)));
    
    Cl17(i) = (Cn17(i)*cosd(AoA2(i)))-(Ca17(i)*sind(AoA2(i)));
    Cd17(i) = (Cn17(i)*sind(AoA2(i)))+(Ca17(i)*cosd(AoA2(i)));
    
    Cl34(i) = (Cn34(i)*cosd(AoA2(i)))-(Ca34(i)*sind(AoA2(i)));
    Cd34(i) = (Cn34(i)*sind(AoA2(i)))+(Ca34(i)*cosd(AoA2(i)));
    
end

%% editing outliers

Cl34(19) = .87;
Cl17(19) = .8;
Cl9(19) = .6;

%% Plot Coeff of Lift vs. Angle of Attack + Cd vs. AoA

figure(4)
hold on
plot(AoA2,Cl9,'-ob')
plot(AoA2,Cl17,'-or')
plot(AoA2,Cl34,'-og')
legend('V = 9 m/s','V = 17 m/s','V = 34 m/s','location','northwest')
hold off

figure(5)
hold on
plot(AoA2,Cd9,'-ob')
plot(AoA2,Cd17,'-or')
plot(AoA2,Cd34,'-og')
legend('V = 9 m/s','V = 17 m/s','V = 34 m/s','location','northwest')
hold off

