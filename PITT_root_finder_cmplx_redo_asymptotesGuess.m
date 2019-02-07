function [ bRoots ] = PITT_root_finder_cmplx_redo_asymptotesGuess( Rohm, Rct,Rd, Capacitance, Tau, bLength )
% Solves transendental equation for the roots of the Concurrent
% Ion-Insertive and Capacitive Storage PITT model. This version is the
% more rigerous approach, and is more likely to find a root if it lies on a
% veritcal asymptote. The function is evaluated for a range of values,
% and roots are found via indexing a change of sign.

% transcendental equation
transcendentalEq=@(b,Rohm, Rct,Rd, Capacitance, Tau) cot(b) - b.*(Rct + Rohm - Capacitance .* Rct .* Rohm .* ( b.^2./Tau))./(Rd.*(1-Capacitance .* Rohm .*( b.^2 ./ Tau)));

% Sensativity of the root finding. As a rule there tends to be (at least) a
% root every period of pi, hence the length is set as pi * the desired
% number of roots.
step = [1e-10:0.001:(bLength+1)*pi()];

painfulAsymptote = sqrt(Tau/(Rohm*Capacitance));
babySteps = 1/Rd/100;

step = [step, [ painfulAsymptote*(1-babySteps*100) : babySteps : painfulAsymptote*(1+babySteps*100) ] ];
step = sort(step);


% Evaluates the transendental equation over the desired range
bValTemp = transcendentalEq(step,Rohm, Rct,Rd, Capacitance, Tau);

% finds changes of sign of the evaluated transendental equation
bRoots = step(find(bValTemp(1:end-1)>0 & bValTemp(2:end) <0));

% check whether the number of roots founds is equal or less than the
% requested number
if bRoots < bLength+1
    
    figure
    plot(step(1:1000:end),bValTemp(1:1000:end),'o')
    ylim([-1e2,1e2])
    hold on
    plot(bRoots, zeros(1,length(bRoots)),'o','MarkerSize',12)
    
    f = warndlg('Not enough roots found!','Warning')
    
    
end

% % %     given the interval has to be fairly fine to not miss any roots,
% % %     plotting is restricted to every 1000th step.

    figure
    plot(step(1:1:end),bValTemp(1:1:end),'o')
    ylim([-1e2,1e2])
    hold on
    plot(bRoots, zeros(1,length(bRoots)),'o','MarkerSize',12)

bRoots = bRoots(1:bLength);
