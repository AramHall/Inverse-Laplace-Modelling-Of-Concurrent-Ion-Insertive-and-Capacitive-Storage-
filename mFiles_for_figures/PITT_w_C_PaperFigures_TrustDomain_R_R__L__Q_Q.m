


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
bLength = 5;
Estep = 0.025;

% preallocate trustDomain
trustDomain = cell(length(RctByRohmRatio),length(Lambda)) ;
trustDomainErr = cell(length(RctByRohmRatio),length(Lambda)) ;

for j = 1:length(RctByRohmRatio)
    
    
    for k = 1:length(Lambda)
        tic
        
        for l = 1:length(QedlcByQinsRatio)
            
            Rct = RctByRohmRatio(j) * Rohm;
            Rd = Lambda(k) * (Rct + Rohm);
            Qd = trapz(Time,PITT_montella_redo( Rohm, Rct,Rd, Tau, bLength, Estep, Time  ));
            Cap = Qd / Estep * QedlcByQinsRatio(l);
            
            sumCurrent = PITT_w_C_redo( Rohm, Rct ,Rd, Cap , Tau, bLength, Estep, Time);
            
            guessMontella = PITT_montella_redo( Rohm, Rct,Rd, Tau, bLength, Estep, Time  )';
            
            guessRC = Estep./Rohm.*exp(-Time./(Rohm.*Cap));
            
            guessCapMontSum = guessMontella + guessRC;
            
            
            
            
            errorCapMontSum = abs(sumCurrent-guessCapMontSum')./sumCurrent;
            trustDomainErr{j,k}{1}(:,l) = errorCapMontSum;
            errorCapMontSum(errorCapMontSum >0.01) = 9999;
            errorCapMontSum(errorCapMontSum <0.01) = 1;
            errorCapMontSum(errorCapMontSum >1) = 0;
            
            errorMontella = abs(sumCurrent-guessMontella')./sumCurrent;
            trustDomainErr{j,k}{2}(:,l) = errorMontella;
            errorMontella(errorMontella >0.01) = 9999;
            errorMontella(errorMontella <0.01) = 10;
            errorMontella(errorMontella >10) = 0;
            
            errorRC = abs(sumCurrent-guessRC')./sumCurrent;
            trustDomainErr{j,k}{3}(:,l) = errorRC;
            errorRC(errorRC >0.01) = 9999;
            errorRC(errorRC <0.01) = 100;
            errorRC(errorRC >100) = 0;
           
            
            trustDomain{j,k}{1}(:,l) = errorCapMontSum;
            trustDomain{j,k}{2}(:,l) = errorMontella;
            trustDomain{j,k}{3}(:,l) = errorRC;

        end
        k
        toc
    end

    j
end


%%
% for j = 1:length(RctByRohmRatio)
%     
%     
%     for k = 1:length(Lambda)
%             
%         subplot(length(RctByRohmRatio),length(Lambda),(j-1)*length(Lambda)+k)
%         contour(Time,QedlcByQinsRatio,[trustDomain{j,k}{2}]',[10,10])
%         hold on
%         contour(Time,QedlcByQinsRatio,[trustDomain{j,k}{3}]',[100,100])
%         set(gca, 'YScale', 'log')
%         set(gca, 'XScale', 'log')
%         view(0,90)
% 
%     end
% end

for j = 1:length(RctByRohmRatio)
    
    
    for k = 1:length(Lambda)
            
            errorCapMontSum = trustDomainErr{j,k}{1};
            errorCapMontSum(errorCapMontSum >0.01) = 9999;
            errorCapMontSum(errorCapMontSum <0.01) = 1;
            errorCapMontSum(errorCapMontSum >1) = 0;
            
            
            errorMontella = trustDomainErr{j,k}{2};
            errorMontella(errorMontella >0.01) = 9999;
            errorMontella(errorMontella <0.01) = 1.0001;
            errorMontella(errorMontella >1.0001) = 0;
            
            
            errorRC = trustDomainErr{j,k}{3};
            errorRC(errorRC >0.01) = 9999;
            errorRC(errorRC <0.01) = 1.0002;
            errorRC(errorRC >1.0002) = 0;
            
            trustDomain{j,k}{1} = errorCapMontSum;
            trustDomain{j,k}{2} = errorMontella;
            trustDomain{j,k}{3} = errorRC;
    end
end



% figure
% 
% for j = 1:3
%     
%     
%     for k = 1:3
%             
%         subplot(3,3,(j-1)*(3)+k)
%         stem3(Time,QedlcByQinsRatio,[trustDomain{j,k}{2}]','color','b')
%         hold on
%         stem3(Time,QedlcByQinsRatio,[trustDomain{j,k}{3}]','color','r','marker','.')
%         set(gca, 'YScale', 'log')
%         set(gca, 'XScale', 'log')
%         zlim([0.5,150])
%         view(0,90)
%         if k> 1
%            set(gca,'YTick',[]);
%         end
%         if j< 3
%            set(gca,'XTick',[]);
%         end
% 
%     end
% end

%%
figure

for j = 1:3
    
    
    for k = 1:3
            
        subplot(3,3,(j-1)*(3)+k)
       
        surf(Time,QedlcByQinsRatio,[trustDomain{(j*2-1),(k*2-1)}{2}]','FaceAlpha',0.1,'EdgeAlpha', 0.5,'EdgeColor', 'blue','FaceColor','blue')
        hold on
        surf(Time,QedlcByQinsRatio,[trustDomain{(j*2-1),(k*2-1)}{3}]','FaceAlpha',0.1,'EdgeAlpha', 0.5,'EdgeColor', 'red','FaceColor','red')
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
        zlim([1 1.0003])
        xlim([1e-5 1e5])
        ylim([1e-3 1e3])
        view(0,90)
        grid off
        if k> 1
           set(gca,'YTick',[]);
        else
            ylabel('Q_E_D_L_C / Q_I_n_s')
        end
        if j< 3
           set(gca,'XTick',[]);
        else
            xlabel('t/\tau')
        end
        if j==1
            if k==3
                legend('I_I_n_s Valid','I_E_D_L_C Valid')
            end
        end
        

    end
end



% figure
% 
% for j = 1:2:length(RctByRohmRatio)
%     
%     
%     for k = 1:2:length(Lambda)
%             
%         subplot(length(RctByRohmRatio)-2,length(Lambda)-2,(j-1)*(length(Lambda)-2)+k)
%         stem3(Time,QedlcByQinsRatio,[trustDomain{j,k}{2}]','color','b')
%         hold on
%         stem3(Time,QedlcByQinsRatio,[trustDomain{j,k}{3}]','color','r','marker','.')
%         set(gca, 'YScale', 'log')
%         set(gca, 'XScale', 'log')
%         zlim([0.5,150])
%         view(0,90)
%         if k> 1
%            set(gca,'YTick',[]);
%         end
%         if j< 5
%            set(gca,'XTick',[]);
%         end
% 
%     end
% end


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
