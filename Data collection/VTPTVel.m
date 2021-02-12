function [output1] = VTPTVel(input1)
    deltaP = input1(:,3);
    R = 287;
    temp = input1(:,2);
    Pressure_atm = input1(:,1);
    A2 = 1;
    A1 = 9.5;
    output1 = sqrt(2.*deltaP .* (R .* temp)./(Pressure_atm .* (1- (A2/A1)^2)));
end