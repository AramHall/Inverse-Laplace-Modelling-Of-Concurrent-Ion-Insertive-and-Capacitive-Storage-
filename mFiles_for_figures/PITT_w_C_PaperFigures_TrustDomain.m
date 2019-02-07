%% Trust Domains    Varied ?, t/tau          Sanity check C = 0

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc = Qd

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = logspace(-2,3,11);
%Cap set @ Qd
Tau = 1;
Time = logspace(-5,10,1000);
bLength = 10;
Estep = 0.025;

figure 
hold on
col = jet(2*length(RdVals));

for j = 1:length(RdVals)
    
    Cap(j) = 0;
    
    sumCurrent(:,j) = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    
    

    guessMontella(:,j) = PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  )';
    guessCott(:,j) = trapz(Time,sumCurrent(:,j))./sqrt(pi().* Time .* Tau);
       
    
    errorMontella(:,j) = abs(sumCurrent(:,j)-guessMontella(:,j))./sumCurrent(:,j);
    errorCott(:,j) = abs(sumCurrent(:,j)-guessCott(:,j))./sumCurrent(:,j);
    j 

    
end




%% Trust Domains    Varied ?, t/tau          Qc = Qd * 0.001

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc = Qd*0.001

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = logspace(-2,3,11);
%Cap set @ Qd*0.001
Tau = 1;
Time = logspace(-5,10,1000);
bLength = 10;
Estep = 0.025;



for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)/(1000*Estep);
    
    sumCurrent(:,j) = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    
    
    guessCapMontSum(:,j) = PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  )' + Estep./Rohm.*exp(-Time./(Rohm.*Cap(j)));
    
    guessMontella(:,j) = PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  )';
    
    guessRC(:,j) = Estep./Rohm.*exp(-Time./(Rohm.*Cap(j)));
    
    guessCott(:,j) = trapz(Time,sumCurrent(:,j))./sqrt(pi().* Time .* Tau);
       
    
    errorCapMontSum(:,j) = abs(sumCurrent(:,j)-guessCapMontSum(:,j))./sumCurrent(:,j);
    errorMontella(:,j) = abs(sumCurrent(:,j)-guessMontella(:,j))./sumCurrent(:,j);
    errorRC(:,j) = abs(sumCurrent(:,j)-guessRC(:,j))./sumCurrent(:,j);
    errorCott(:,j) = abs(sumCurrent(:,j)-guessCott(:,j))./sumCurrent(:,j);
    j 

    
end




%% Trust Domains    Varied ?, t/tau          Qc = Qd 

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc = Qd

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = logspace(-2,3,11);
%Cap set @ Qd
Tau = 1;
Time = logspace(-5,10,1000);
bLength = 10;
Estep = 0.025;


for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)/(Estep);
    
    sumCurrent(:,j) = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    
    
    guessCapMontSum(:,j) = PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  )' + Estep./Rohm.*exp(-Time./(Rohm.*Cap(j)));
    
    guessMontella(:,j) = PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  )';
    
    guessRC(:,j) = Estep./Rohm.*exp(-Time./(Rohm.*Cap(j)));
    
    guessCott(:,j) = trapz(Time,sumCurrent(:,j))./sqrt(pi().* Time .* Tau);
       
    
    errorCapMontSum(:,j) = abs(sumCurrent(:,j)-guessCapMontSum(:,j))./sumCurrent(:,j);
    errorMontella(:,j) = abs(sumCurrent(:,j)-guessMontella(:,j))./sumCurrent(:,j);
    errorRC(:,j) = abs(sumCurrent(:,j)-guessRC(:,j))./sumCurrent(:,j);
    errorCott(:,j) = abs(sumCurrent(:,j)-guessCott(:,j))./sumCurrent(:,j);
    j 

    
end





%% Trust Domains    Varied ?, t/tau          Qc = Qd *1000

% Varied Lambda from 0.01,0.1,1,10,100
% Rohm = Rct = 0.5
% Rd for various Lambda = the Lambda Values (bc Lambda = Rd/(Rohm+Rct)
% This for Qc = Qd*1000

clear 


Rohm = 0.5;
Rct = 0.5;
RdVals = logspace(-2,3,11);
%Cap set @ Qd*1000
Tau = 1;
Time = logspace(-5,10,1000);
bLength = 10;
Estep = 0.025;


for j = 1:length(RdVals)
    
    
    Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
    Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
    Cap(j) =  Qd(j)*1000/(Estep);
    
    sumCurrent(:,j) = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
    
    
    guessCapMontSum(:,j) = PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  )' + Estep./Rohm.*exp(-Time./(Rohm.*Cap(j)));
    
    guessMontella(:,j) = PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  )';
    
    guessRC(:,j) = Estep./Rohm.*exp(-Time./(Rohm.*Cap(j)));
    
    guessCott(:,j) = trapz(Time,sumCurrent(:,j))./sqrt(pi().* Time .* Tau);
       
    
    errorCapMontSum(:,j) = abs(sumCurrent(:,j)-guessCapMontSum(:,j))./sumCurrent(:,j);
    errorMontella(:,j) = abs(sumCurrent(:,j)-guessMontella(:,j))./sumCurrent(:,j);
    errorRC(:,j) = abs(sumCurrent(:,j)-guessRC(:,j))./sumCurrent(:,j);
    errorCott(:,j) = abs(sumCurrent(:,j)-guessCott(:,j))./sumCurrent(:,j);
    j 

    
end






%%

clear error1
error1 = errorCapMontSum;
error1(error1 >0.01) = 10;
error1(error1 <0.01) = 1;
error1(error1 >1) = 0;
figure
mesh(Time,RdVals, error1')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('t/\tau')
ylabel('\Lambda')
zlabel('Validity domain')

clear error1
error1 = errorMontella;
error1(error1 >0.01) = 10;
error1(error1 <0.01) = 1;
error1(error1 >1) = 0;
figure
mesh(Time,RdVals, error1')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('t/\tau')
ylabel('\Lambda')
zlabel('Validity domain')

clear error1
error1 = errorRC;
error1(error1 >0.01) = 10;
error1(error1 <0.01) = 1;
error1(error1 >1) = 0;
figure
mesh(Time,RdVals, error1')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('t/\tau')
ylabel('\Lambda')
zlabel('Validity domain')

clear error1
error1 = errorCott;
error1(error1 >0.01) = 10;
error1(error1 <0.01) = 1;
error1(error1 >1) = 0;
figure
mesh(Time,RdVals, error1')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('t/\tau')
ylabel('\Lambda')
zlabel('Validity domain')

