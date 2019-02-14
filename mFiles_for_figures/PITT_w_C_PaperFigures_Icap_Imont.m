%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I(t)/Imont
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% I(t)/Imont     Varied Lambda          Qc = Qd/1000

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc = Qd x  1000
% Imont taken from the values in I(t) but no cap

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = [0.01,0.1,1,10,100];
%Cap set @ Qd*1000
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 20;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));  

for j = 1:length(RdVals)
    
    insCurrent = PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  );
    Qd(j) = trapz(Time,insCurrent);
    Cap(j) =  Qd(j)/1000/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, Estep, Time);
    
    j 
    plotCurrent = sumCurrent./insCurrent;
    plotCurrent(isnan(plotCurrent)) = [];
    plotCurrent(plotCurrent==1) = [];
    plot( Time(1:length(plotCurrent))./Tau, plotCurrent, 'color', col(j,:))
    legEnt{j} = ['\Lambda = ' num2str(RdVals(j))];
    
end

ylabel('I(t) / I_I_n_s')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
ylim([1 1.01])




%% I(t)/Imont     Varied R/R          Qc = Qd/1000

% Lambda constant = 1
% Rohm:Rct = 0.01 0.1 1 10 100
% This for Qc<<Qd, say 1/1000
% Imont taken from the values in I(t) but no cap

clear 


Rohm = 1;
RctVals = [0.01,0.1,1,10,100];
RdVals = Rohm + RctVals;
%Cap set @ Qd/1000
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 20;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    insCurrent = PITT_montella_redo( Rohm, RctVals(j),RdVals(j), Tau, bLength, Estep, Time  );
    Qd(j) = trapz(Time,insCurrent);
    Cap(j) =  Qd(j)/1000/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, RctVals(j) ,RdVals(j), Cap(j) , Tau, Estep, Time);
    j 
    plotCurrent = sumCurrent./insCurrent;
    plotCurrent(isnan(plotCurrent)) = [];
    plotCurrent(plotCurrent==1) = [];
    plot( Time(1:length(plotCurrent))./Tau, plotCurrent, 'color', col(j,:))
    legEnt{j} = ['R_c_t/R_\Omega = ' num2str(RctVals(j))];
    
end

ylabel('I(t) / I_I_n_s')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([1e-5 1e5])
ylim([1e-4 1])



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I(t)/Icap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% I(t)/I(0)     Varied ?          Qc = Qd x 1000

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc = Qd x  1000
% Icap taken from the values in I(t) Rohm and Cap

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = [0.01,0.1,1,10,100];
%Cap set @ Qd*1000
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 40;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Cap(j) =  Qd(j)*1000/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, Estep, Time);
    capCurrent = Estep/Rohm .* exp(-Time./(Cap(j).*Rohm));
    j 
    plot( Time./Tau, sumCurrent./capCurrent', 'color', col(j,:))
    legEnt{j} = ['\Lambda = ' num2str(RdVals(j))];
    
end

ylabel('I(t) / I_E_D_L_C')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
xlim([1e-5 1e5])
ylim([1 1.01])











%% I(t)/Icap     Varied R/R          Qc = Qd*1000

% Lambda constant = 1
% Rohm:Rct = 0.01 0.1 1 10 100
% This for Qc>>Qd, say 1000/1
% Icap taken from the values in I(t) Rohm and Cap

clear 


Rohm = 1;
RctVals = [0.01,0.1,1,10,100];
RdVals = Rohm + RctVals;
%Cap set @ Qd/1000
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 40;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, RctVals(j),RdVals(j), Tau, bLength, Estep, Time  ));
    Cap(j) =  Qd(j)*1000/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, RctVals(j) ,RdVals(j), Cap(j) , Tau, Estep, Time);
    capCurrent = Estep/Rohm .* exp(-Time./(Cap(j).*Rohm));
    j 
    plot( Time./Tau, sumCurrent./capCurrent', 'color', col(j,:))
    legEnt{j} = ['R_c_t/R_\Omega = ' num2str(RctVals(j))];
    
end

ylabel('I(t) / I_E_D_L_C')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([1e-5 1e5])
ylim([1 1.01])



