


%% Trust Domains    

% Multilevel graph for trust domains
% varied Lambda and Rct/Rohm
% at each, plot trust domains for Q/Q vs t/Tau
% at the moment only checks Iedlc, Iins, Iedlc + Iins

clear 

% constants
Rohm = 1;
Tau = 1;
Time = logspace(-5,5,1000);
RctByRohmRatio = logspace(-2,2,5);
Lambda = logspace(-2,2,5);
QedlcByQinsRatio = logspace(-3,3,31);
bLength = 10;
Estep = 0.025;

% preallocate trustDomain
trustDomain = cell(length(RctByRohmRatio),length(Lambda)) ;

for j = 1:length(RctByRohmRatio)
    
    
    for k = 1:length(Lambda)
        
        
        for l = 1:length(QedlcByQinsRatio)
            
            Rct = QedlcByQinsRatio(l) * Rohm;
            Rd = Lambda(k) * (Rct + Rohm);
            Qd = trapz(Time,PITT_montella_redo( Rohm, Rct,Rd, Tau, bLength, Estep, Time  ));
            Cap = Qd / Estep * QedlcByQinsRatio(l);
            
            sumCurrent = PITT_w_C_redo( Rohm, Rct ,Rd, Cap , Tau, bLength, Estep, Time);
            
            guessMontella = PITT_montella_redo( Rohm, Rct,Rd, Tau, bLength, Estep, Time  )';
            
            guessRC = Estep./Rohm.*exp(-Time./(Rohm.*Cap));
            
            guessCapMontSum = guessMontella + guessRC;
            
            
            errorCapMontSum = abs(sumCurrent-guessCapMontSum)./sumCurrent;
            errorCapMontSum(errorCapMontSum >0.01) = 10;
            errorCapMontSum(errorCapMontSum <0.01) = 1;
            errorCapMontSum(errorCapMontSum >1) = 0;
            
            errorMontella = abs(sumCurrent-guessMontella)./sumCurrent;
            errorMontella(errorMontella >0.01) = 10;
            errorMontella(errorMontella <0.01) = 1;
            errorMontella(errorMontella >1) = 0;
            
            errorRC = abs(sumCurrent-guessRC)./sumCurrent;
            errorRC(errorRC >0.01) = 10;
            errorRC(errorRC <0.01) = 1;
            errorRC(errorRC >1) = 0;
           
            
            trustDomain{j,k}{1}(:,l) = errorCapMontSum;
            trustDomain{j,k}{2}(:,l) = errorMontella;
            trustDomain{j,k}{3}(:,l) = errorRC;

            
            l
            
        end
        
        stem3(QedlcByQinsRatio,Time,trustDomain{j,k}{2},'color','b')
        zlim([0.5,2])
        hold on
        stem3(QedlcByQinsRatio,Time,trustDomain{j,k}{3},'color','b')
        k
    end

    j
end





% Rohm = 0.5;
% Rct = 0.5;
% RdVals = logspace(-2,3,11);
% %Cap set @ Qd*0.001
% Tau = 1;
% Time = logspace(-5,10,1000);
% bLength = 10;
% Estep = 0.025;
% 
% 
% 
% for j = 1:length(RdVals)
%     
%     
%     Qd(j) = trapz(Time,PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  ));
%     Qd2(j) = sum(PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  (1:end-1)  ).*diff(Time'));
%     Cap(j) =  Qd(j)/(1000*Estep);
%     
%     sumCurrent(:,j) = PITT_w_C_redo( Rohm, Rct ,RdVals(j), Cap(j) , Tau, bLength, Estep, Time);
%     
%     
%     guessCapMontSum(:,j) = PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  )' + Estep./Rohm.*exp(-Time./(Rohm.*Cap(j)));
%     
%     guessMontella(:,j) = PITT_montella_redo( Rohm, Rct,RdVals(j), Tau, bLength, Estep, Time  )';
%     
%     guessRC(:,j) = Estep./Rohm.*exp(-Time./(Rohm.*Cap(j)));
%     
%     guessCott(:,j) = trapz(Time,sumCurrent(:,j))./sqrt(pi().* Time .* Tau);
%        
%     
%     errorCapMontSum(:,j) = abs(sumCurrent(:,j)-guessCapMontSum(:,j))./sumCurrent(:,j);
%     errorMontella(:,j) = abs(sumCurrent(:,j)-guessMontella(:,j))./sumCurrent(:,j);
%     errorRC(:,j) = abs(sumCurrent(:,j)-guessRC(:,j))./sumCurrent(:,j);
%     errorCott(:,j) = abs(sumCurrent(:,j)-guessCott(:,j))./sumCurrent(:,j);
%     j 
% 
%     
% end
% 
% 
% 
% 



% 
% 
% %%
% 
% clear error1
% error1 = errorCapMontSum;
% error1(error1 >0.01) = 10;
% error1(error1 <0.01) = 1;
% error1(error1 >1) = 0;
% figure
% mesh(Time,RdVals, error1')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlabel('t/\tau')
% ylabel('\Lambda')
% zlabel('Validity domain')
% 
% clear error1
% error1 = errorMontella;
% error1(error1 >0.01) = 10;
% error1(error1 <0.01) = 1;
% error1(error1 >1) = 0;
% figure
% mesh(Time,RdVals, error1')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlabel('t/\tau')
% ylabel('\Lambda')
% zlabel('Validity domain')
% 
% clear error1
% error1 = errorRC;
% error1(error1 >0.01) = 10;
% error1(error1 <0.01) = 1;
% error1(error1 >1) = 0;
% figure
% mesh(Time,RdVals, error1')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlabel('t/\tau')
% ylabel('\Lambda')
% zlabel('Validity domain')
% 
% clear error1
% error1 = errorCott;
% error1(error1 >0.01) = 10;
% error1(error1 <0.01) = 1;
% error1(error1 >1) = 0;
% figure
% mesh(Time,RdVals, error1')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlabel('t/\tau')
% ylabel('\Lambda')
% zlabel('Validity domain')
% 
