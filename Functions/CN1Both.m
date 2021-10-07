function [ VectorConn, VectorConnstr ] = CN1Both( FilteredFFT, Powers, FilteredFFT1, Powers1, EpochTime, Thereshold_Corr, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures)
%function [ VectorConn, VectorConnstr ] = CN( FilteredFFT, Powers, EpochTime, Thereshold_Corr, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures)
%   Plot the Correlation Network & Find the Average Connection and the
%   Strength of the Connections of the network



    P = (Powers) / (max(Powers));

    S(:,:) = FilteredFFT;
    CoefCorr = zeros(length(Electrodes));
    for j = 1: length(Electrodes)
    for k =j+1:length(Electrodes)
        R = corrcoef(S(j,:),S(k,:));
        CoefCorr(j,k) = sqrt(((R(1,2)+1)/2).^3+mean([P(j) P(k)]).^2);
        CoefCorr(k,j) = sqrt(((R(1,2)+1)/2).^3+mean([P(j) P(k)]).^2);
    end
    end

    Cmin = min(CoefCorr(:));
    Cmax = max(CoefCorr(:));
    CoefCorr = (CoefCorr) / (Cmax);
     
    
    Cat = [];
    ConnStr = zeros(1,length(Regions));
     EstimatedSource = zeros(length(Regions),2);
     for r = 1:length(Regions)
         max1 = 0;
         xy1 = [ 0, 0];
         max2 = 0;
         xy2 = [ 0, 0];
        R = cell2mat(Regions(r));
        for k=1:length(R)
            for j=k+1:length(R)
                if (CoefCorr(R(k),R(j)) > Thereshold_Corr) 
                    Cat = [Cat r];
%                     ConnStr(r) = ConnStr(r) + (1 - (1 - CoefCorr(R(k),R(j)))/(1 - Thereshold_Corr));
                    tmpedge = (1 - (1 - CoefCorr(R(k),R(j)))/(1 - Thereshold_Corr));
                    if tmpedge > max1
                        max2 = max1;
                        xy2 = xy1;
                        max1 = tmpedge;
                        xy1 = mean([Centers(R(j),:);Centers(R(k),:)]);
                    elseif tmpedge > max2 
                        max2 = tmpedge;
                        xy2 = mean([Centers(R(j),:);Centers(R(k),:)]);
                    end
                end
            end
        end
         ConnStr(r) = mean([max1,max2]);
         EstimatedSource(r,:) = mean([xy1;xy2]);
    end
    
    C = categorical(Cat,[1:length(Regions)],RegionLabels');
    [Ns,Categories] = histcounts(C);
    for r = 1:length(ConnStr)
        if Ns(r) == 0
            VectorConnstr(r) = 0;
        else
%          VectorConnstrDist(r) = ConnStr(r) ./ Ns(r);
            VectorConnstr(r) = ConnStr(r);
        end
    end
    
    for r = 1:length(Regions)
        Ns(r) = Ns(r)/nchoosek(length(cell2mat(Regions(r))),2);
    end
    VectorConn = Ns;

    
    P1 = (Powers1) / (max(Powers1));

    S1(:,:) = FilteredFFT1;
    CoefCorr1 = zeros(length(Electrodes));
    for j = 1: length(Electrodes)
    for k =j+1:length(Electrodes)
        R = corrcoef(S1(j,:),S1(k,:));
        CoefCorr1(j,k) = sqrt(((R(1,2)+1)/2).^3+mean([P1(j) P1(k)]).^2);
        CoefCorr1(k,j) = sqrt(((R(1,2)+1)/2).^3+mean([P1(j) P1(k)]).^2);
    end
    end

    Cmin1 = min(CoefCorr1(:));
    Cmax1 = max(CoefCorr1(:));
    CoefCorr1 = (CoefCorr1) / (Cmax1);
     
    
    Cat1 = [];
    ConnStr1 = zeros(1,length(Regions));
     EstimatedSource1 = zeros(length(Regions),2);
     for r = 1:length(Regions)
         max11 = 0;
         xy11 = [ 0, 0];
         max21 = 0;
         xy21 = [ 0, 0];
        R = cell2mat(Regions(r));
        for k=1:length(R)
            for j=k+1:length(R)
                if (CoefCorr1(R(k),R(j)) > Thereshold_Corr) 
                    Cat1 = [Cat1 r];
%                     ConnStr(r) = ConnStr(r) + (1 - (1 - CoefCorr(R(k),R(j)))/(1 - Thereshold_Corr));
                    tmpedge1 = (1 - (1 - CoefCorr(R(k),R(j)))/(1 - Thereshold_Corr));
                    if tmpedge1 > max11
                        max21 = max11;
                        xy21 = xy11;
                        max11 = tmpedge1;
                        xy11 = mean([Centers(R(j),:);Centers(R(k),:)]);
                    elseif tmpedge1 > max21 
                        max21 = tmpedge1;
                        xy21 = mean([Centers(R(j),:);Centers(R(k),:)]);
                    end
                end
            end
        end
         ConnStr1(r) = mean([max11,max21]);
         EstimatedSource1(r,:) = mean([xy11;xy21]);
    end
    
    C1 = categorical(Cat1,[1:length(Regions)],RegionLabels');
    [Ns1,Categories1] = histcounts(C1);
    for r = 1:length(ConnStr1)
        if Ns1(r) == 0
            VectorConnstr1(r) = 0;
        else
%          VectorConnstrDist(r) = ConnStr(r) ./ Ns(r);
            VectorConnstr1(r) = ConnStr1(r);
        end
    end
    
    for r = 1:length(Regions)
        Ns1(r) = Ns1(r)/nchoosek(length(cell2mat(Regions(r))),2);
    end
    VectorConn1 = Ns1;    
    
    
    if strcmp(Figures,'on')
    
        figure('Visible',Visibility)
        subplot(2,2,1)
        
%         BrainMap( Centers, Type );
%         BrainMapScalpPowers( Centers, FilteredFFT );
        BrainMapScalpPowers( Centers, CoefCorr );
        
        hold on
        for k=1:length(Labels)-1
           for j=k+1:length(Labels)
              if (CoefCorr(k,j) > Thereshold_Corr) 
                 Connstrength = 1 - (1 - CoefCorr(k,j))/(1 - Thereshold_Corr);
                 line([Centers(k,1) Centers(j,1)],[Centers(k,2) Centers(j,2)],'color',1-[Connstrength Connstrength Connstrength],'LineWidth',(2*Connstrength+0.5));
              end
           end
        end
        
        tlable=[ num2str(EpochTime(1) )  ' -- ' num2str(EpochTime(2))];

        axis equal
        set(gcf, 'Position', get(0, 'Screensize'))
        for t=1:length(Labels)
            text(Centers(t,1),Centers(t,2),Labels{t},'HorizontalAlignment','center','FontSize',10);
        end
        
     [~,I] = max(ConnStr(3:end-1));
     circles(3,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','k','linewidth',2)
%      circles(6,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','r','linewidth',2)
     circles(9,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','k','linewidth',2)
%      circles(12,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','r','linewidth',2)
     circles(15,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','k','linewidth',2)
     circles(20,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','k','linewidth',2)
     
        title({Name,'Correlation Network - Filtered',['\color{blue}' tlable]});
        for s = 1:Seizures(1)
            if (Seizures(2*s)<= EpochTime(1) ) && (Seizures((2*s)+1) >= EpochTime(2) )
                title({Name,'Average Connections',['\color{red}' tlable]});
            end
        end
		
		subplot(2,2,3)
        
%         BrainMap( Centers, Type );
%         BrainMapScalpPowers( Centers, FilteredFFT );
        BrainMapScalpPowers( Centers, CoefCorr1 );
        
        hold on
        for k=1:length(Labels)-1
           for j=k+1:length(Labels)
              if (CoefCorr1(k,j) > Thereshold_Corr) 
                 Connstrength1 = 1 - (1 - CoefCorr1(k,j))/(1 - Thereshold_Corr);
                 line([Centers(k,1) Centers(j,1)],[Centers(k,2) Centers(j,2)],'color',1-[Connstrength1 Connstrength1 Connstrength1],'LineWidth',(2*Connstrength1+0.5));
              end
           end
        end
        
        tlable=[ num2str(EpochTime(1) )  ' -- ' num2str(EpochTime(2))];

        axis equal
        set(gcf, 'Position', get(0, 'Screensize'))
        for t=1:length(Labels)
            text(Centers(t,1),Centers(t,2),Labels{t},'HorizontalAlignment','center','FontSize',10);
        end
        
     [~,I] = max(ConnStr1(3:end-1));
     circles(3,[EstimatedSource1(I+2,1) EstimatedSource1(I+2,2)],'color','k','linewidth',2)
%      circles(6,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','r','linewidth',2)
     circles(9,[EstimatedSource1(I+2,1) EstimatedSource1(I+2,2)],'color','k','linewidth',2)
%      circles(12,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','r','linewidth',2)
     circles(15,[EstimatedSource1(I+2,1) EstimatedSource1(I+2,2)],'color','k','linewidth',2)
     circles(20,[EstimatedSource1(I+2,1) EstimatedSource1(I+2,2)],'color','k','linewidth',2)
     
        title({Name,'Correlation Network - Raw',['\color{blue}' tlable]});
        for s = 1:Seizures(1)
            if (Seizures(2*s)<= EpochTime(1) ) && (Seizures((2*s)+1) >= EpochTime(2) )
                title({Name,'Average Connections',['\color{red}' tlable]});
            end
        end
		
		
        subplot(2,2,2)

        bar([Ns;Ns1]')
        set(gca,'xticklabel',RegionLabels)
        ylim([0 1]);
        title({Name,'Average Connections',['\color{blue}' tlable]});
        for s = 1:Seizures(1)
            if (Seizures(2*s)<= EpochTime(1) ) && (Seizures((2*s)+1) >= EpochTime(2) )
                title({Name,'Average Connections',['\color{red}' tlable]});
            end
        end
        legend('Filtered','Raw')
        
        subplot(2,2,4)
        bar([VectorConnstr;VectorConnstr1]')
        set(gca,'xticklabel',RegionLabels)
        ylim([0 1]);
        title({Name,'Connections Strength',['\color{blue}' tlable]});
        for s = 1:Seizures(1)
            if (Seizures(2*s)<= EpochTime(1) ) && (Seizures((2*s)+1) >= EpochTime(2) )
                title({Name,'Average Connections',['\color{red}' tlable]});
            end
        end 
        legend('Filtered','Raw')
        set(gcf, 'Position', get(0, 'Screensize'))

    end
    
end