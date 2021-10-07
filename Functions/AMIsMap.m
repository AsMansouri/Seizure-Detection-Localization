function [] = AMIsMap( AMIsPairs, EpochTime, Thereshold_Corr, Centers, ChLabels, Electrodes, Seizures, Name, Visibility)
%function [] = AMIsMap( AMIsPairs, EpochTime, Centers, Labels, Electrodes, Seizures, Name, Visibility)
%   Plot the Pairwise AMIs Networks between channels
%   The Background indicates the Entropies of channels


    Entropies = diag(AMIsPairs);
    
    Cmax = max(AMIsPairs(:));
    AMIsPairs = (AMIsPairs) / (Cmax);
    
    
    figure('Visible',Visibility)

    testX = Centers(:,1);
    testY = Centers(:,2);
    testZ = Entropies;

    % Create polar data
    [t,r] = meshgrid((0:0.5:360)*pi/180,[0:max(abs(Centers(:)))/10:(max(abs(Centers(:)))+max(abs(Centers(:)))/10)]);
    % Convert to Cartesian
    [X,Y] = pol2cart(t,r);

    F = scatteredInterpolant(testX,testY,testZ);
    Z = F(X,Y);
    smoothZ = smoothdata(Z);
    contourf(X,Y,smoothZ,10,'LineColor', 'none');
%     hold on
    daspect([1 1 1])
    colorbar('Location','westoutside')
    colormap(jet)
%     axis image;
    set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
        
    hold on
    for k=1:length(Electrodes)
       for j=k+1:length(Electrodes)
            if (CoefCorr(k,j) > Thereshold_Corr) 
                Connstrength = 1 - (1 - AMIsPairs(k,j))/(1 - Thereshold_Corr);
                line([Centers(k,1) Centers(j,1)],[Centers(k,2) Centers(j,2)],'color',[Connstrength Connstrength Connstrength],'LineWidth',(2*Connstrength+0.5));
            end
       end
    end
    
    tlable=[ num2str(EpochTime(1) )  ' -- ' num2str(EpochTime(2))];

    axis equal
    set(gcf, 'Position', get(0, 'Screensize'))
    
    for t=1:length(ChLabels)
        text(Centers(t,1),Centers(t,2),ChLabels{t},'HorizontalAlignment','center','FontSize',10);
    end
        
    title({Name,'AMI Network',['\color{blue}' tlable],'Non-ictal'});
    for s = 1:Seizures(1)
        if (Seizures(2*s)<= EpochTime(1) ) && (Seizures((2*s)+1) >= EpochTime(2) )
            title({Name,'AMI Connections',['\color{red}' tlable],'Ictal period'});
        end
    end
    
%      [~,I] = max(ConnStr(3:end-1));
%      circles(3,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','k','linewidth',2)
% %      circles(6,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','r','linewidth',2)
%      circles(9,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','k','linewidth',2)
% %      circles(12,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','r','linewidth',2)
%      circles(15,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','k','linewidth',2)
%      circles(20,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','k','linewidth',2)
     

		
    
%     title({Name,'Connections Strength',['\color{blue}' tlable]});
%         
%         
%         subplot(2,2,2)
% 
%         bar(Ns)
%         set(gca,'xticklabel',RegionLabels)
%         ylim([0 1]);
%         title({Name,'Average Connections',['\color{blue}' tlable]});
%         for s = 1:Seizures(1)
%             if (Seizures(2*s)<= EpochTime(1) ) && (Seizures((2*s)+1) >= EpochTime(2) )
%                 title({Name,'Average Connections',['\color{red}' tlable]});
%             end
%         end
%         subplot(2,2,4)
%         bar(VectorConnstr)
%         set(gca,'xticklabel',RegionLabels)
%         ylim([0 1]);
%         for s = 1:Seizures(1)
%             if (Seizures(2*s)<= EpochTime(1) ) && (Seizures((2*s)+1) >= EpochTime(2) )
%                 title({Name,'Average Connections',['\color{red}' tlable]});
%             end
%         end            
%         set(gcf, 'Position', get(0, 'Screensize'))


    
end