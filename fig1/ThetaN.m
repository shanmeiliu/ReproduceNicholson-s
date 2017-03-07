function dydt = ThetaN( tau,y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gamma=1;
dydt=gamma*y/(y^2+gamma^2*tau^2+gamma*tau);

end

