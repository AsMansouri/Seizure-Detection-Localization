function [  ] = BrainMap( Centers, Type )
%function [  ] = BrainMap( Centers, Type )
%   Draw the Background Map of a brain and Electrode Locations

     Rt = max(max(Centers(:,1))-min(Centers(:,1)),max(Centers(:,2))-min(Centers(:,2)))/2;
     if Type == 10
         set(gca,'ydir','reverse');
         circles(Rt+25,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
         circles(Rt/3,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
         circles(25,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
         circles(5,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
         circles(5,Centers,'color','b','linewidth',0.5)
     else
         if Type ~=1
             circles(Rt+25,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
             circles(Rt/3,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
             circles(25,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
             circles(5,[mean(Centers(:,1)) mean(Centers(:,2))],'color','y','linewidth',2)
             circles(5,Centers,'color','b','linewidth',0.5)
         else
             circles(Rt+30,[0 0],'color','y','linewidth',2)
             circles(Rt/2,[0 0],'color','y','linewidth',2)
             circles(50,[0 0],'color','y','linewidth',2)
             circles(10,[0 0],'color','y','linewidth',2)
             circles(20,Centers,'color','b','linewidth',0.5)
         end
     end
end

