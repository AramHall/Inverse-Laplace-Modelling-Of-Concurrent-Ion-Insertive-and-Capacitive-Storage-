%% I(t)/I(0)     Varied ?          Edge Case 1: c=0

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for C = 0
% Note I(0) in this case is Estep/(Rohm+ Rct)

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
    j 
    plot( Time./Tau, sumCurrent./(Estep./(Rohm+Rct)), 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')




%% I(t)/I(0)     Varied ?          Qc = Qd/1000

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
bLength = 40;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Cap(j) =  Qd(j)/(1000*Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    j 
    plot( Time./Tau, (sumCurrent./(Estep/Rohm)),'color', col(j,:))
    legEnt{j} = ['\Lambda = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([1e-4 1])



%% I(t)/I(0)     Varied ?          Qc = Qd

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
bLength = 40;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    j 
    plot( Time./Tau, sumCurrent./(Estep/Rohm), 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ylim([1e-4 1])



















%% I(t)/I(0)     Varied ?          Qc = Qd x 1000

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
    j 
    plot( Time./Tau, sumCurrent./(Estep/Rohm), 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')












%% I(t)/I(0)     Varied ?          Qc only

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
    j 
    plot( Time./Tau, sumCurrent./(Estep/Rohm), 'color', col(j,:))
    legEnt{j} = ['? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Varied Rct/Rohm
%%%%%%%%%%%%%%%%%%%%%%%%%%



%% I(t)/I(0)     Varied R/R          Qc = Qd/1000

% Lambda constant = 1
% Rohm:Rct = 0.01 0.1 1 10 100
% This for Qc<<Qd, say 1/1000

clear 


Rohm = 1;
RctVals = [0.01,0.1,1,10,100];
RdVals = Rohm + RctVals;
%Cap set @ Qd/1000
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, RctVals(j),RdVals(j), Tau, bLength, Estep, Time  ));
    Cap(j) =  Qd(j)/(1000*Estep);
    sumCurrent = PITT_w_C_redo( Rohm, RctVals(j) ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    j 
    plot( Time./Tau, sumCurrent./(Estep/Rohm), 'color', col(j,:))
    legEnt{j} = ['R_c_t/R_\Omega = ' num2str(RctVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([1e-5 1e5])
ylim([1e-4 1])



%% I(t)/I(0)     Varied R/R          Qc = Qd

% Lambda constant = 1
% Rohm:Rct = 0.01 0.1 1 10 100
% This for Qc=Qd

clear 


Rohm = 1;
RctVals = [0.01,0.1,1,10,100];
RdVals = Rohm + RctVals;
%Cap set @ Qd/1000
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, RctVals(j),RdVals(j), Tau, bLength, Estep, Time  ));
    Cap(j) =  Qd(j)/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, RctVals(j) ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    j 
    plot( Time./Tau, sumCurrent./(Estep/Rohm), 'color', col(j,:))
    legEnt{j} = ['R_c_t/R_\Omega = ' num2str(RctVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([1e-5 1e5])
ylim([1e-4 1])





%% I(t)/I(0)     Varied R/R          Qc = Qd*1000

% Lambda constant = 1
% Rohm:Rct = 0.01 0.1 1 10 100
% This for Qc>>Qd, say 1000/1

clear 


Rohm = 1;
RctVals = [0.01,0.1,1,10,100];
RdVals = Rohm + RctVals;
%Cap set @ Qd/1000
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, RctVals(j),RdVals(j), Tau, bLength, Estep, Time  ));
    Cap(j) =  Qd(j)*1000/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, RctVals(j) ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    j 
    plot( Time./Tau, sumCurrent./(Estep/Rohm), 'color', col(j,:))
    legEnt{j} = ['R_c_t/R_\Omega = ' num2str(RctVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([1e-5 1e5])
ylim([1e-4 1])



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brickabrack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% I(t)/I(0)     Varied ?          Qc = Qd

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
Time = logspace(-5,10,1000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(2*length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    guessCurrent = PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  )' + Estep./Rohm.*exp(-Time./(Rohm.*Cap(j)));
    j 
%     plot( Time./Tau, sumCurrent, 'color', col(2*j-1,:))
%     plot( Time./Tau, guessCurrent, 'color', col(2*j,:))
    plot( Time./Tau, abs(sumCurrent'-guessCurrent)./sumCurrent', 'color', col(2*j,:))
    legEnt{2*j-1} = ['Real ? = ' num2str(RdVals(j))];
    legEnt{2*j} = ['Est ? = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')










%% I(t)/I(0)    Fixed Lambda         Varied Qc = Qd

%  Lambda =1
% This for Qc = Qd

clear 


Rohm = [0.5];
Rct = Rohm*0.1;
RdVals = (Rohm + Rct)*100;
mult = logspace(-3,3,7);
%Cap set @ Qd
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(mult));

for j = 1:length(mult)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals, Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals, Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)*mult(j)/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm, Rct ,RdVals, Cap(j) , Tau, bLength, Estep, Time);
    j 
    plot( Time./Tau, sumCurrent./(Estep/Rohm), 'color', col(j,:))
    legEnt{j} = ['Q_E_D_L_C / Q_L_i = ' num2str(mult(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')












