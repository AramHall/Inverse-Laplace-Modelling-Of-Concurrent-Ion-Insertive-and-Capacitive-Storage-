



%% I(t)/I(0)     Varied ?          Qc = Qd

% Varied Lambda =1
%Fixed Q/Q
%Fixed Lambda
%Varied Rd
%Fixed Rct/Rohm
% This for Qc = Qd

clear 


RdVals = [0.01,0.1,1,10,100];
Rohm = RdVals./2;
Rct = RdVals./2;
%Cap set @ Qd
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm(j), Rct(j),RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm(j), Rct(j),RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm(j), Rct(j) ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    j 
    plot( Time./Tau, sumCurrent./(Estep/Rohm(j)), 'color', col(j,:))
    legEnt{j} = ['Rd = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')














%% I(t)/I(0)     Varied ?          Qc = Qd x 1000

%  Lambda =1
%Fixed Q/Q
%Fixed Lambda
%Varied Romh
%Fixed Rct
% This for Qc = Qd

clear 


Rohm = [0.01,0.1,1,10,100];
Rct = Rohm*0.1;
RdVals = Rohm + Rct;
%Cap set @ Qd
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm(j), Rct(j),RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm(j), Rct(j),RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)*1000/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm(j), Rct(j) ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    j 
    plot( Time./Tau, sumCurrent./(Estep/Rohm(j)), 'color', col(j,:))
    legEnt{j} = ['Rd = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')









%% I(t)/I(0)     Varied ?          Qc = Qd

%  Lambda =1
%Fixed Q/Q
%Fixed Lambda
%Varied Romh
%Fixed Rct
% This for Qc = Qd

clear 


Rohm = [0.01,0.1,1,10,100];
Rct = Rohm;
RdVals = Rohm + Rct;
%Cap set @ Qd
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm(j), Rct(j),RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm(j), Rct(j),RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm(j), Rct(j) ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    j 
    plot( Time./Tau, sumCurrent./(Estep/Rohm(j)), 'color', col(j,:))
    legEnt{j} = ['Rd = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')



















%% I(t)/I(0)     Varied ?          Qc = Qd x 0.001

%  Lambda =1
%Fixed Q/Q
%Fixed Lambda
%Varied Romh
%Fixed Rct
% This for Qc = Qd

clear 


Rohm = [0.01,0.1,1,10,100];
Rct = Rohm;
RdVals = Rohm + Rct;
%Cap set @ Qd
Tau = 1;
Time = logspace(-5,10,10000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(length(RdVals));

for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm(j), Rct(j),RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm(j), Rct(j),RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)*0.001/(Estep);
    sumCurrent = PITT_w_C_redo( Rohm(j), Rct(j) ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    j 
    plot( Time./Tau, sumCurrent./(Estep/Rohm(j)), 'color', col(j,:))
    legEnt{j} = ['Rd = ' num2str(RdVals(j))];
    
end

ylabel('I_C_+_D(t) / I(0)')
xlabel('t/\tau')
set(gca,'fontsize', 12);
legend(legEnt)
set(gca, 'XScale', 'log')



















