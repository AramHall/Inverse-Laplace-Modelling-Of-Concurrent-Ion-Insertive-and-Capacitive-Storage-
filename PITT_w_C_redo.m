function [  sumCurrent ] = PITT_w_C_redo( Rohm, Rct,Rd, Cpt, Tau, Estep, time  )
% Concurrent Ion-Insertive and Capacitive Storage PITT model to PITT data

% Estep is the stepped potential, as a single value. This should be the
% actual potential step, not the absolute value, as this is what determines
% if the curret is +ve or -ve.

% Rohm, Rct,Rd, Cpt and Tau are parameters, and should be single values. 

% time is a time series vector, for which the current will be modelled with
% the given parameters Rohm, Rct,Rd, Cpt, Tau and Estep

% substituted to make the equation more readable    
t = time;

% using the parameters Rohm/Rct/Rd/Cpt/Tau, finds the requested roots
bRoots  = PITT_root_finder_cmplx_redo_asymptotesGuess( Rohm, Rct,Rd, Cpt, Tau );
   

% initialise the matrix for the current transients. 
Currents = zeros(length(time),length(bRoots));
   

% evaluates the current from every sum term (the length of which is set by the number of roots)
for k = 1:length(bRoots)
    
    % substituted to make the equation more readable    
    x = bRoots(k);  
    
    tempCurrent =  exp(-x.^2.*t./Tau).*2.*Rd.*Estep.*Tau.^2./...
        (Cpt.^2.*Rohm.^2.*x.^4.*(Rd.*(Rct+Rd)+Rct.^2.*x.^2)-...
        Cpt.*Rohm.*x.^2.*(2.*Rd.*(Rct+Rd)-Rd.*Rohm + 2.*Rct.*(Rct+Rohm).*x.^2).*Tau +...
        (Rd.*(Rct+Rohm+Rd)+(Rct+Rohm).^2.*x.^2).*Tau.^2) ;
    
    tempCurrent(isnan(tempCurrent)) = 0;
    
    Currents(:,k) =tempCurrent;
    


end

% figure
% mesh(bRoots,time,Currents)

% sums the currents at every time step to form the modelled current
% transient
sumCurrent = sum(Currents,2);


end

