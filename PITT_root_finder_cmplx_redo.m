function [ bRoots ] = PITT_root_finder_cmplx_redo( Rohm, Rct,Rd, Capacitance, Tau, bLength )
% Solves transendental equation for the roots of the Concurrent
% Ion-Insertive and Capacitive Storage PITT model. This version is the
% default, quick method, and in some cases will miss a root if a veritcal
% asymptote. If the root lies on a veritcal asymptote, solver functions
% seem ot struggle. Hence the function is evaluated for a range of values,
% and roots are found via indexing a change of sign.

% bLength is the number of roots to find
 
% transcendental equation
transcendentalEq=@(b,Rohm, Rct,Rd, Capacitance, Tau) cot(b) - b.*(Rct + Rohm - Capacitance .* Rct .* Rohm .* ( b.^2./Tau))./(Rd.*(1-Capacitance .* Rohm .*( b.^2 ./ Tau)));

% Sensativity of the prelim root finding. As a rule there tends to be a
% root every period of pi, hence the length is set as pi * the desired
% number of roots.
step = [1e-10:0.001:bLength*pi()];

% This controls the sensativity of the finer root finding
interval = 0.001;

% evaluates the transendental equation over the rough step size
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
   
    
    
% does a finer pass at the previously identified roots     
else
    
    %runs through all the previously identified roots
    for ii = 1:bLength
        
        clear step bValTemp
        
        % creates a finer resultion evaluation vector
        step = [bRoots(ii)*(1-interval):0.0000001:bRoots(ii)*(1+interval)];
        
        % reevaluates the transendental equation around a root
        bValTemp = transcendentalEq(step,Rohm, Rct,Rd, Capacitance, Tau);
        
        % finds changes of sign of the evaluated transendental equation
        bRootsFinePassTemp = step(find(bValTemp(1:end-1)>0 & bValTemp(2:end) <0));
        
        % checks whether the interval was big enough to find the root, and
        % if not expands the interval
        if isempty(bRootsFinePassTemp)
            
            % increases the interval by an order of magnitude
            step = [bRoots(ii)*(1-interval*10):0.0000001:bRoots(ii)*(1+interval*10)];
            
            bValTemp = transcendentalEq(step,Rohm, Rct,Rd, Capacitance, Tau);
            
            bRootsFinePassTemp = step(find(bValTemp(1:end-1)>0 & bValTemp(2:end) <0));
            
            % if the root still isn't on the interval, sets the root such
            % so that it'll throw an error when checked later. 
            if isempty(bRootsFinePassTemp)
                bRootsFinePass(ii) = 99999;
            else
                bRootsFinePass(ii) = bRootsFinePassTemp;
            end
            
        else
             bRootsFinePass(ii) = bRootsFinePassTemp;
        end
        
        
    end
end

bRoots = bRootsFinePass;
