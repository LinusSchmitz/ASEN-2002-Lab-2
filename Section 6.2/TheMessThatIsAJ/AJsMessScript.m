%%%% THIS IS AJ TRYING TO FIGURE THINGS OUT WIHTOUT LOOSING CODE THAT WORKS




%% gathering excel

close all
clear 
clc

%% Read in Ventury tube vs Pressure Transducer data
R = 287; % J/kg*K

parentdir = '/MATLAB Drive/ASEN 2002 Lab/Lab 2/Section 6.2/TheMessThatIsAJ';

subs = dir(fullfile(parentdir,'**'));

isdir = [subs.isdir] & ~ismember({subs.name},{'.','..'});

subdirs = fullfile({subs(isdir).folder},{subs(isdir).name});

dP = NaN(size(subdirs));
AuxdP = NaN;
atmP = NaN;
T = [];
ELD = [];

for i = 1:length(subdirs)

    content = dir(subdirs{i});
    files = {content.name};

    TotalPortData = [];
    for j = 3:length(files)-2
        
        name = convertCharsToStrings([parentdir subdirs{i} files{j}]);
        PortData = csvread(name);
        TotalPortData = [TotalPortData;PortData];
        
%         dP(:,i) = TotalPortData(:,3);
%         AuxdP(:,i) = TotalPortData(:,4);
%         atmP(:,i) = TotalPortData(:,1);
%         T(:,i) = TotalPortData(:,2);
%         ELD(:,i) = TotalPortData(:,6);
        
        
    end
    disp(files)
 
end

% for i = 3:13
% 
%     files = dir('/MATLAB Drive/ASEN 2002 Lab/Lab 2/Section 6.2/Boundary_Layer_Data');
%     long_name = strcat(files(i).folder,'/',files(i).name);
%     BoundPortData = load(long_name);
%     
%     eval(['Port_' num2str(i-2) ' = BoundPortData;'])
%  
% end

% 
% 
% dP = NaN(6500,11);  
%     
% dP(1:6000,1) = Port_1(:,3);
% dP(1:6000,2) = Port_2(:,3);
% dP(1:6000,3) = Port_3(:,3);
% dP(1:6000,4) = Port_4(:,3);
% dP(1:6000,5) = Port_5(:,3);
% dP(1:6000,6) = Port_6(:,3);
% dP(1:6000,7) = Port_7(:,3);
% dP(1:6000,8) = Port_8(:,3);
% dP(1:6000,9) = Port_9(:,3);
% dP(1:6000,10) = Port_10(:,3);
% dP(1:6000,11) = Port_11(:,3);
% 
% atmP = NaN(6500,11);  
%     
% atmP(1:6000,1) = Port_1(:,1);
% atmP(1:6000,2) = Port_2(:,1);
% atmP(1:6000,3) = Port_3(:,1);
% atmP(1:6000,4) = Port_4(:,1);
% atmP(1:6000,5) = Port_5(:,1);
% atmP(1:6000,6) = Port_6(:,1);
% atmP(1:6000,7) = Port_7(:,1);
% atmP(1:6000,8) = Port_8(:,1);
% atmP(1:6000,9) = Port_9(:,1);
% atmP(1:6000,10) = Port_10(:,1);
% atmP(1:6000,11) = Port_11(:,1);
% 
% T = NaN(6500,11);  
%     
% T(1:6000,1) = Port_1(:,2);
% T(1:6000,2) = Port_2(:,2);
% T(1:6000,3) = Port_3(:,2);
% T(1:6000,4) = Port_4(:,2);
% T(1:6000,5) = Port_5(:,2);
% T(1:6000,6) = Port_6(:,2);
% T(1:6000,7) = Port_7(:,2);
% T(1:6000,8) = Port_8(:,2);
% T(1:6000,9) = Port_9(:,2);
% T(1:6000,10) = Port_10(:,2);
% T(1:6000,11) = Port_11(:,2);
% 
% rho = atmP./(R.*T);
% 
% % pre allocate airpspeed
% airspeed = zeros(size(dP));
% for i = 1:11
%     
% airspeed(:,i) = sqrt((2./rho(:,i)).*dP(:,i));
%     
% end
% 
% ELD = NaN(6500,11);
% 
% ELD(1:6000,1) = Port_1(:,6);
% ELD(1:6000,2) = Port_2(:,6);
% ELD(1:6000,3) = Port_3(:,6);
% ELD(1:6000,4) = Port_4(:,6);
% ELD(1:6000,5) = Port_5(:,6);
% ELD(1:6000,6) = Port_6(:,6);
% ELD(1:6000,7) = Port_7(:,6);
% ELD(1:6000,8) = Port_8(:,6);
% ELD(1:6000,9) = Port_9(:,6);
% ELD(1:6000,10) = Port_10(:,6);
% ELD(1:6000,11) = Port_11(:,6);
% 
% AuxdP = NaN(6500,11);  
%     
% AuxdP(1:6000,1) = Port_1(:,4);
% AuxdP(1:6000,2) = Port_2(:,4);
% AuxdP(1:6000,3) = Port_3(:,4);
% AuxdP(1:6000,4) = Port_4(:,4);
% AuxdP(1:6000,5) = Port_5(:,4);
% AuxdP(1:6000,6) = Port_6(:,4);
% AuxdP(1:6000,7) = Port_7(:,4);
% AuxdP(1:6000,8) = Port_8(:,4);
% AuxdP(1:6000,9) = Port_9(:,4);
% AuxdP(1:6000,10) = Port_10(:,4);
% AuxdP(1:6000,11) = Port_11(:,4);
% 
% % pre allocate airpspeed
% Auxairspeed = zeros(size(AuxdP));
% for i = 1:11
%     
% Auxairspeed(:,i)= sqrt((2./rho(:,i)).*abs(AuxdP(:,i)));
%     
% end
% 
% 
% figure(1)
% for p = 1:11
%     figure(p)
%     plot(ELD(:,p),Auxairspeed(:,p));
% end
% 
% 
% 
% delta = 0.95.*airspeed;
% 
% BoundaryLayer = zeros(size(airspeed));
% for i = 1:11
%     for j = 1:6500
%         if delta(j,i) >= Auxairspeed(j,i)
%             BoundaryLayer(j,i) = ELD(j,i);
%         end
%     end
% end
% 
% PortMatrix = ones(size(BoundaryLayer));
% for k = 1:11
%     PortMatrix(:,k) = k.*PortMatrix(:,k);
% end
% 
% figure(12)
% plot(PortMatrix,BoundaryLayer,'o');
% xlabel('Port Number');
% ylabel('Boundary Layer')
% 
% %% Section 7
% 
% alpha = -15:15; % angle of attack
% 
% Re = 300;
% 
% % pressure = 
% 
%   