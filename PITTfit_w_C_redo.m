function [ fittedCurrent , fittedParams,ci] = PITTfit_w_C_redo( current, time, Estep,initParams )
% Function for fitting Concurrent Ion-Insertive and Capacitive Storage PITT
% model to PITT data

% bLength is the number roots used in the solution. The greater the number,
% the more accurate, however the more computationally intensive. At bLength
% values > ~5, improvement is really only seen in the short time domain. It
% is recommended to use values between 5 (faster fitting) and 15(greater accuracy). 

% Estep is the stepped potential, as a single value. This should be the
% actual potential step, not the absolute value, as this is what determines
% if the curret is +ve or -ve.

% sets the initial parameter guesses 
% if none are provided, the default is a fast faradaic process, with a
% double layer capacity ~1/20th of the faradaic process.
if isempty(initParams)
   pInx = [5;25;60;1e-5;30];
else
    pInx = initParams;
end

% Fitted parameters have structure:
% fittedParams(1) = Rohm
% fittedParams(2) = Rct
% fittedParams(3) = Rd
% fittedParams(4) = Cpt
% fittedParams(5) = Tau

% sets solver options to more strict conditions
opts = optimoptions('lsqnonlin','FunctionTolerance',1e-20,'OptimalityTolerance',1e-20,...
    'MaxIterations', 500,'MaxFunctionEvaluations',500,'StepTolerance',1e-20);


% sets solver model, based on the PITT_w_C_redo function
iModel01 =@(fittedParams) PITT_w_C_redo( fittedParams(1), fittedParams(2),...
    fittedParams(3), fittedParams(4), fittedParams(5), Estep, time ) - current;

% nonlinear least squares to fit the model to PITT data
%[fittedParams,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(iModel01,pInx,[1e-5;1e-5;1e-5;1e-8;1e-1],[100;1000;Inf;1;1000],opts);
[fittedParams,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(iModel01,pInx,[1e-5,1e-5,1e-5,1e-10,0.01],[1000;1e5;1e6;100;1e4],opts);
%[fittedParams,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(iModel01,pInx,[0,0,0,0,0],[Inf;Inf;Inf;Inf;Inf]);

% Calculates a 95% CI for each of the fitted parameters
ci = nlparci(fittedParams,resid,'jacobian',J);

% uses the fitted parameters to generate the current transient
fittedCurrent = PITT_w_C_redo( fittedParams(1), fittedParams(2),...
fittedParams(3), fittedParams(4), fittedParams(5), Estep, time  );



end










