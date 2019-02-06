function [  sumCurrent ] = PITT_w_C_redo( Rohm, Rct,Rd, Cpt, Tau, bLength, Estep, time  )
% Concurrent Ion-Insertive and Capacitive Storage PITT model to PITT data

% bLength is the number roots used in the solution. The greater the number,
% the more accurate, however the more computationally intensive. At bLength
% values > ~5, improvement is really only seen in the short time domain. It
% is recommended to use values between 5 (faster fitting) and 15(greater accuracy).

% Estep is the stepped potential, as a single value. This should be the
% actual potential step, not the absolute value, as this is what determines
% if the curret is +ve or -ve.

% initialise the matrix for the current transients. 
Currents = zeros(length(time),bLength);

% substituted to make the equation more readable    
t = time;

% using the parameters Rohm/Rct/Rd/Cpt/Tau, finds the requested roots
bRoots  = PITT_root_finder_cmplx_redo( Rohm, Rct,Rd, Cpt, Tau, bLength );

%%%%    in some cases the quick root rinder above might miss a root. In this
%%%%    case the more rigerous root finder can be used instead:
%%%%    bRoots  = PITT_root_finder_cmplx_redo_rigourous( Rohm, Rct,Rd, Cpt, Tau, bLength );

% evaluates the current from every sum term (the length of which is set by the number of roots)
for k = 1:bLength
    
    % substituted to make the equation more readable    
    x = bRoots(k);  
    
    Currents(:,k) = exp(-x.^2.*t./Tau).*2.*Rd.*Estep.*Tau.^2./...
        (Cpt.^2.*Rohm.^2.*x.^4.*(Rd.*(Rct+Rd)+Rct.^2.*x.^2)-...
        Cpt.*Rohm.*x.^2.*(2.*Rd.*(Rct+Rd)-Rd.*Rohm + 2.*Rct.*(Rct+Rohm).*x.^2).*Tau +...
        (Rd.*(Rct+Rohm+Rd)+(Rct+Rohm).^2.*x.^2).*Tau.^2) ;

end

% sums the currents at every time step to form the modelled current
% transient
sumCurrent = sum(Currents,2);


end

