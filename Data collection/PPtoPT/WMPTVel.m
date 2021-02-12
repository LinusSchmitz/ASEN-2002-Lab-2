function [output1] = WMPTVel(input1, input2)
    deltaP = input1(:,4);
    R = 287;
    temp = input2(:,2);
    Pressure_atm = input2(:,1);
    A2 = 1;
    A1 = 9.5;
    output1 = sqrt(2.*deltaP .* (R .* temp)./(Pressure_atm .* (1- (A2/A1)^2)));
end