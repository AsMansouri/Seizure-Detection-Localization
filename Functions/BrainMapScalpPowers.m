function [  ] = BrainMapScalpPowers( Centers, FilteredFFT )
%function [  ] = BrainMap( Centers, Type )
%   Draw the Background Map of a brain and Electrode Locations


%      Rt = max(max(Centers(:,1))-min(Centers(:,1)),max(Centers(:,2))-min(Centers(:,2)))/2;
%      if Type == 10
%          set(gca,'ydir','reverse');
%          circles(Rt+25,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
%          circles(Rt/3,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
%          circles(25,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
%          circles(5,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
%          circles(5,Centers,'color','b','linewidth',0.5)
%      else
%          if Type ~=1 || Type ~=11
%              circles(Rt+25,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
%              circles(Rt/3,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
%              circles(25,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
%              circles(5,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
%              circles(5,Centers,'color','b','linewidth',0.5)
%          else
%              circles(Rt+30,[0 0],'color','y','linewidth',2)
%              circles(Rt/2,[0 0],'color','y','linewidth',2)
%              circles(50,[0 0],'color','y','linewidth',2)
%              circles(10,[0 0],'color','y','linewidth',2)
%              circles(20,Centers,'color','b','linewidth',0.5)
%          end
%      end
     

    S(:,:) = FilteredFFT;
    ElecPSD = sum(S,2);
    
    testX = Centers(:,1);
    testY = Centers(:,2);
    testZ = (ElecPSD-min(ElecPSD(:)))/(max(ElecPSD(:))-min(ElecPSD(:)));

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
     
end

