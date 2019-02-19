function [ bRoots ] = PITT_root_finder_redo( Lambda, bLength )
%UNTITLED Solves b.Tan(b)=Lambda for PITT fitting

transendentalFunction=@(b,Lambda) b.*tan(b) - Lambda ;

%This controls the sensativity of the rootfinding
step = [[1e-10:0.0001:1],[(1+1e-10):0.001:(bLength+1)*pi()]];

bValTemp = transendentalFunction(step,Lambda) ;

bRootsFirst = step(find(bValTemp(1:end-1)<0 & bValTemp(2:end) >0));

opts = optimset('Display','off');

if length(bRootsFirst)<bLength
        
        diffVal = diff(bValTemp);
        bRootsOnePointFive = step(find(diffVal<0));
        
        for jjj = 1:length(bRootsOnePointFive)
            bRootsSecond(jjj) = fsolve(@(b)transendentalFunction(b,Lambda), bRootsOnePointFive(jjj),opts);
        end
        
else

        for jjj = 1:length(bRootsFirst)
            bRootsSecond(jjj) = fsolve(@(b)transendentalFunction(b,Lambda), bRootsFirst(jjj),opts);
            
        end
        
end



if length(bRootsSecond) < bLength
    
    figure
    plot(step(1:1000:end),bValTemp(1:1000:end),'o')
    ylim([-1e2,1e2])
    hold on
    plot(bRoots, zeros(1,length(bRoots)),'o','MarkerSize',12)
    
    f = warndlg('Not enough roots found!','Warning')
    
    

end

%     figure
%     plot(step(1:1000:end),bValTemp(1:1000:end),'o')
%     ylim([-1e2,1e2])
%     hold on
%     plot(bRoots, zeros(1,length(bRoots)),'o','MarkerSize',12)

bRoots = bRootsSecond(1:bLength);
end

