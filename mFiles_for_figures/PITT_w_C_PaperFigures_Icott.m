%% I(t)/I_cottrell     Varied ?          Edge Case 1: c=0

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for C = 0

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = [0.01,0.1,1,10,100];
Cap = 0;
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap , Tau, bLength, Estep, Time);
    Icott = trapz(Time,sumCurrent)./sqrt(pi().* Time .* Tau);
    j 
    plot( Time./Tau, sumCurrent./Icott', 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I_c_o_t_t_r_e_l_l')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
xlim([1E-5 100])
ylim([0 1.1])



%% I(t)/I_c_o_t_t_r_e_l_l     Varied ?          Qc = Qd/1000

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc<<Qd, say 1/1000

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = [0.01,0.1,1,10,100];
%Cap set @ Qd/1000
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)/(1000*Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    Icott = trapz(Time,sumCurrent)./sqrt(pi().* Time .* Tau);
    j 
    plot( Time./Tau, sumCurrent./Icott', 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I_c_o_t_t_r_e_l_l')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
xlim([1E-5 100])
ylim([0 1.1])



%% I(t)/I_c_o_t_t_r_e_l_l     Varied ?          Qc = Qd

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc = Qd

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = [0.01,0.1,1,10,100];
%Cap set @ Qd
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    Icott = trapz(Time,sumCurrent)./sqrt(pi().* Time .* Tau);
    j 
    plot( Time./Tau, sumCurrent./Icott', 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I_c_o_t_t_r_e_l_l')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
xlim([1E-5 100])
ylim([0 1.1])


















%% I(t)/I_c_o_t_t_r_e_l_l     Varied ?          Qc = Qd x 1000

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc = Qd x  1000

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = [0.01,0.1,1,10,100];
%Cap set @ Qd*1000
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)*1000/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    Icott = trapz(Time,sumCurrent)./sqrt(pi().* Time .* Tau);
    j 
    plot( Time./Tau, sumCurrent./Icott', 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I_c_o_t_t_r_e_l_l')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
xlim([1E-5 100])
ylim([0 1.1])





%% I(t)/I_c_o_t_t_r_e_l_l     Varied ?          Qc only

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc only
% Uses the cap value from the previous graph

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = [0.01,0.1,1,10,100];   
%Cap set @ Qd*1000
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)*1000/(Estep);
    sumCurrent = Estep./Rohm.*exp(-Time./(Cap(j).*Rohm));
    Icott = trapz(Time,sumCurrent)./sqrt(pi().* Time .* Tau);
    j 
    plot( Time./Tau, sumCurrent./Icott, 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I_c_o_t_t_r_e_l_l')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
xlim([1E-5 100])
ylim([0 1.1])






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   New bit for normalising Icott with only Qd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% I(t)/I_c_o_t_t_r_e_l_l     Varied ?          Qc = Qd/1000 Qmontella=Qcott

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc<<Qd, say 1/1000

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = [0.01,0.1,1,10,100];
%Cap set @ Qd/1000
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)/(1000*Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    Icott = Qd(j)./sqrt(pi().* Time .* Tau);
    j 
    plot( Time./Tau, sumCurrent./Icott', 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I_c_o_t_t_r_e_l_l')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
xlim([1E-5 100])
ylim([0 1.1])



%% I(t)/I_c_o_t_t_r_e_l_l     Varied ?          Qc = Qd Qmontella=Qcott

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc = Qd

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = [0.01,0.1,1,10,100];
%Cap set @ Qd
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    Icott = Qd(j)./sqrt(pi().* Time .* Tau);
    j 
    plot( Time./Tau, sumCurrent./Icott', 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I_c_o_t_t_r_e_l_l')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
xlim([1E-5 100])
ylim([0 1.1])


















%% I(t)/I_c_o_t_t_r_e_l_l     Varied ?          Qc = Qd x 1000 Qmontella=Qcott

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc = Qd x  1000

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = [0.01,0.1,1,10,100];
%Cap set @ Qd*1000
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)*1000/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    Icott = Qd(j)./sqrt(pi().* Time .* Tau);
    j 
    plot( Time./Tau, sumCurrent./Icott', 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I_c_o_t_t_r_e_l_l')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
xlim([1E-5 100])
ylim([0 1.1])




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Just a wee test to make sure the normalisation is working
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc = Qd x  1000

clear 


Rohm = 5;
Rct = 5;
RdVals = [0.01,0.1,1,10,100].*10;
%Cap set @ Qd*1000
Tau = 10;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)*1000/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    Icott = Qd(j)./sqrt(pi().* Time .* Tau);
    j 
    plot( Time./Tau, sumCurrent./Icott', 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I_c_o_t_t_r_e_l_l')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
xlim([1E-5 100])
ylim([0 1.1])

