%Error in Pitot Static Measurement of Airspeed
%Del_P is measured from data
%T is atmospheric temperature
%P_atm is atmospheric pressure
%R is the universal gas constant
%Anything starting with "un" is uncertainty of the given variable
%dv terms are partial derivates of airspeed of the pitot static with
%   respect to the given variable. 


dvdDel_P = (R*T*unDelP)/(P_atm*sqrt(2*DelP*(R*T)/(P_atm)));
dvdP_atm = (-DelP*R*T*unP_atm)/((P_atm)^(2)*sqrt(2*DelP*(R*T)/(P_atm)));
dvdT = (DelP*R*unT)/(P_atm*sqrt(2*DelP*(R*T)/(P_atm)));

unV = sqrt(dvdDel_P^2+dvdP_atm^2+dvdT^2);