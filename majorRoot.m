function [ root ] = majorRoot( Rohm, Rct,Rd, Capacitance, Tau )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

significantRoots(1) = -1/sqrt(3)*sqrt((Capacitance * Rohm *...
    (-Capacitance * Rd * (Rct + Rd) * Rohm + 2 * Rct * (Rct + Rohm) * Tau) +...
    sqrt( Capacitance^2 * Rohm^2 * (Capacitance^2 * Rd^2 * (Rct + Rd)^2 * Rohm^2 + ...
    Capacitance * Rct * Rd * Rohm * (2* Rct * (Rct + Rd) - ( 7 * Rct + 4 * Rd ) * Rohm )...
    * Tau + Rct^2 * (Rct + Rohm)^2 * Tau^2)))/(Capacitance^2 * Rct^2 * Rohm^2));

significantRoots(2) = 1/sqrt(3)*sqrt((Capacitance * Rohm *...
    (-Capacitance * Rd * (Rct + Rd) * Rohm + 2 * Rct * (Rct + Rohm) * Tau) +...
    sqrt( Capacitance^2 * Rohm^2 * (Capacitance^2 * Rd^2 * (Rct + Rd)^2 * Rohm^2 + ...
    Capacitance * Rct * Rd * Rohm * (2* Rct * (Rct + Rd) - ( 7 * Rct + 4 * Rd ) * Rohm )...
    * Tau + Rct^2 * (Rct + Rohm)^2 * Tau^2)))/(Capacitance^2 * Rct^2 * Rohm^2));

significantRoots(3) = - 1/sqrt(3)*sqrt((Capacitance * Rohm *...
    (-Capacitance * Rd * (Rct + Rd) * Rohm + 2 * Rct * (Rct + Rohm) * Tau) -...
    sqrt( Capacitance^2 * Rohm^2 * (Capacitance^2 * Rd^2 * (Rct + Rd)^2 * Rohm^2 + ...
    Capacitance * Rct * Rd * Rohm * (2* Rct * (Rct + Rd) - ( 7 * Rct + 4 * Rd ) * Rohm )...
    * Tau + Rct^2 * (Rct + Rohm)^2 * Tau^2)))/(Capacitance^2 * Rct^2 * Rohm^2));

significantRoots(4) = 1/sqrt(3)*sqrt((Capacitance * Rohm *...
    (-Capacitance * Rd * (Rct + Rd) * Rohm + 2 * Rct * (Rct + Rohm) * Tau) -...
    sqrt( Capacitance^2 * Rohm^2 * (Capacitance^2 * Rd^2 * (Rct + Rd)^2 * Rohm^2 + ...
    Capacitance * Rct * Rd * Rohm * (2* Rct * (Rct + Rd) - ( 7 * Rct + 4 * Rd ) * Rohm )...
    * Tau + Rct^2 * (Rct + Rohm)^2 * Tau^2)))/(Capacitance^2 * Rct^2 * Rohm^2));

significantRoots = significantRoots(significantRoots==real(significantRoots));
significantRoots = significantRoots(significantRoots>0);


[~,indx] = max(2.*Rd.*Tau.^2./...
        (Capacitance.^2.*Rohm.^2.*significantRoots.^4.*(Rd.*(Rct+Rd)+Rct.^2.*significantRoots.^2)-...
        Capacitance.*Rohm.*significantRoots.^2.*(2.*Rd.*(Rct+Rd)-Rd.*Rohm + 2.*Rct.*(Rct+Rohm).*significantRoots.^2).*Tau +...
        (Rd.*(Rct+Rohm+Rd)+(Rct+Rohm).^2.*significantRoots.^2).*Tau.^2));   

root = significantRoots(indx);


end

