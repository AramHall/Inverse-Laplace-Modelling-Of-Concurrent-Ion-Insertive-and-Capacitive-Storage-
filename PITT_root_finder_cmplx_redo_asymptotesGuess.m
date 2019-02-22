function [ bRoots ] = PITT_root_finder_cmplx_redo_asymptotesGuess( Rohm, Rct,Rd, Capacitance, Tau )
% Solves transendental equation for the roots of the Concurrent
% Ion-Insertive and Capacitive Storage PITT model. This version is the
% more rigerous approach, and is more likely to find a root if it lies on a
% veritcal asymptote. The function is evaluated for a range of values,
% and roots are found via indexing a change of sign.

%%%%% NB: I break if Rd > ~ 1e15

% transcendental equation
transcendentalEq=@(b,Rohm, Rct,Rd, Capacitance, Tau) cot(b) - b.*(Rct + Rohm - Capacitance .* Rct .* Rohm .* ( b.^2./Tau))./(Rd.*(1-Capacitance .* Rohm .*( b.^2 ./ Tau)));

% Sensativity of the root finding. As a rule there tends to be (at least) a
% root every period of pi, hence the length is set as pi * the desired
% number of roots.

[ mostSignifRoot ] = majorRoot( Rohm, Rct,Rd, Capacitance, Tau );

if isempty(mostSignifRoot)
    mostSignifRoot = 0;
end
% evaluates to at least 10 roots, even if the major one is one of the first
% few

painfulAsymptote = sqrt(Tau/(Rohm*Capacitance));

step = [ [1e-10:0.001:(20*pi)], ...
    [painfulAsymptote-4*pi : 0.001 : painfulAsymptote + 4*pi],...
    [mostSignifRoot-4*pi : 0.001 : mostSignifRoot + 4*pi] ];


babySteps = min(1/Rd/100 , sqrt(Tau/(Rohm*Capacitance))/100000);
babySteps = max(1e-6, babySteps);

step = [ step, [ (painfulAsymptote-babySteps*10000) : babySteps : (painfulAsymptote+babySteps*1000) ] , [ (mostSignifRoot-babySteps*10000) : babySteps : (mostSignifRoot+babySteps*1000) ] ];

if Rd > 1e8 
    step = [step,  [ (painfulAsymptote-(1/Rd)*10) : 1/Rd/10000 : (painfulAsymptote+(1/Rd)*10) ] ];
elseif Capacitance < 1e-4
    step = [step,  [ (painfulAsymptote-Capacitance) : Capacitance/100000 : (painfulAsymptote+Capacitance) ] ];
end

step = sort(step);
step = step(step>0);


% Evaluates the transendental equation over the desired range
bValTemp = transcendentalEq(step,Rohm, Rct,Rd, Capacitance, Tau);


% finds changes of sign of the evaluated transendental equation
bRootsFirst = step(find(bValTemp(1:end-1)>0 & bValTemp(2:end) <0));

opts = optimset('Display','off');

for jjj = 1:length(bRootsFirst)
    
    
    bRootsSecond(jjj) = fsolve(@(b)transcendentalEq(b,Rohm, Rct,Rd, Capacitance, Tau), bRootsFirst(jjj),opts);
    
    
end

if Rd > 1e8
    asymptoteCheck = abs(bRootsFirst./ bRootsSecond-1);
    asymptoteRoot = find(asymptoteCheck > 1e-2);
    bRootsSecond(asymptoteRoot) = bRootsFirst(asymptoteRoot);
elseif Capacitance < 1e-4
    asymptoteCheck = abs(bRootsFirst./ bRootsSecond-1);
    asymptoteRoot = find(asymptoteCheck > 1e-2);
    bRootsSecond(asymptoteRoot) = bRootsFirst(asymptoteRoot);
end
    


bRoots = bRootsSecond;

% check whether the number of roots founds is equal or less than the
% requested number

% % %     given the interval has to be fairly fine to not miss any roots,
% % %     plotting is restricted to every 1000th step.

%     figure
%     plot(step(1:1:end),bValTemp(1:1:end),'o')
%     ylim([-1e2,1e2])
%     hold on
%     plot(bRoots, zeros(1,length(bRoots)),'o','MarkerSize',12)


end

