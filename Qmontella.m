function [ qDiffussion ] = Qmontella( Rohm, Rct,Rd, Tau, bLength, Estep  )
%UNTITLED3 Summary of this function goes here
%   Figures out the capacity from insertive storafe based on montella model


Lambda = Rd/(Rohm+ Rct);

bRoots  = PITT_root_finder_redo(Lambda, bLength );

qDiffussion = 2 * Estep / (Rct + Rohm) * sum( Lambda .* Tau ./ ((Lambda^2 + Lambda + bRoots.^2 ) .* bRoots.^2) );

end

