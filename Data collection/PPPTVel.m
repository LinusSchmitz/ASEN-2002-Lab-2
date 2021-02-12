function [output1] = PPPTVel(input1)
    deltaP = input1(:,3);
    R = 287;
    temp = input1(:,2);
    Pressure_atm = input1(:,1);
    output1 = sqrt(2.*deltaP .* (R .* temp)./Pressure_atm);
end