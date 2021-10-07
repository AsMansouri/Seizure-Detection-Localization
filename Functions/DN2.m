function [ VectorConnDist, VectorConnstrDist ] = DN2( FilteredFFT, EpochTime, Thereshold_Dist, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures, FlagSource)
%function [ VectorConnDist, VectorConnstrDist ] = DN2( FilteredFFT, EpochTime, Thereshold_Dist, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures, FlagSource)
%   Detailed explanation goes here
                

 S(:,:) = FilteredFFT;
 DistMat = zeros(length(Electrodes));
 for k = 1:length(Electrodes)
     for j = k+1:length(Electrodes)
        DistMat(k,j) = sum((S(k,:)-S(j,:)).^2);
        DistMat(j,k) = sum((S(k,:)-S(j,:)).^2);
     end
 end 


 Dmin = min(DistMat(:));
 Dmax = max(DistMat(:));
 tmp = (DistMat - Dmin) / (Dmax - Dmin);
 DistMat = tmp;
 

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
%                         ConnStr(r) = ConnStr(r) + CoefCorr(R(k),R(j));
            if (DistMat(R(k),R(j)) < Thereshold_Dist) 
%                        Cat(r) = Cat(r) + (1 - Distn(R(k),R(j)));
               Cat = [Cat r];
%                ConnStr(r) = ConnStr(r) + ((Thereshold_Dist-DistMat(R(k),R(j)))/Thereshold_Dist);
                tmpedge = (Thereshold_Dist-DistMat(R(k),R(j)))/Thereshold_Dist;
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
     if (xy1(1) ~= 0 || xy1(2) ~= 0) && (xy2(1) ~= 0 || xy2(2) ~= 0)
         EstimatedSource(r,:) = mean([xy1;xy2]);
     elseif (xy1(1) ~= 0 || xy1(2) ~= 0)
         EstimatedSource(r,:) = xy1;
     elseif (xy2(1) ~= 0 || xy2(2) ~= 0)
         EstimatedSource(r,:) = xy2;
     else
         EstimatedSource(r,:) = [0,0];
     end
%      ConnStr(r) = ConnStr(r)/nchoosek(length(cell2mat(Regions(r))),2);
 end

 C = categorical(Cat,[1:length(Regions)],RegionLabels');
 [Ns,Categories] = histcounts(C);

 for r = 1:length(ConnStr)
     if Ns(r) == 0
         VectorConnstrDist(r) = 0;
     else
%          VectorConnstrDist(r) = ConnStr(r) ./ Ns(r);
          VectorConnstrDist(r) = ConnStr(r);
     end
 end


%          histogram(C,'BarWidth',0.5)
 for r = 1:length(Regions)
%              Ns(r) = Cat(r)/sum(1:length(cell2mat(Regions(r))));
     Ns(r) = Ns(r)/nchoosek(length(cell2mat(Regions(r))),2);
 end
 VectorConnDist = Ns;

 if strcmp(Figures,'on')

     figure('Visible',Visibility)
     subplot(2,2,[1 3])
     
%      BrainMap( Centers, Type );
%      BrainMapScalpPowers( Centers, FilteredFFT );
     BrainMapScalpPowers( Centers, DistMat );     
     
     hold on
     for k=1:length(Labels)-1
         for j=k+1:length(Labels)
            if (DistMat(k,j) < Thereshold_Dist) 
               Connstrength = 1 - ((Thereshold_Dist-DistMat(k,j))/Thereshold_Dist);
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
     if FlagSource == 1
         if EstimatedSource(I+2,1) ~= 0 || EstimatedSource(I+2,2) ~= 0
             circles(3,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','k','linewidth',2)
             circles(9,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','k','linewidth',2)
             circles(15,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','k','linewidth',2)
             circles(20,[EstimatedSource(I+2,1) EstimatedSource(I+2,2)],'color','k','linewidth',2)
         end   
     end
%      axis image;

     title({Name,'Distances Network',['\color{blue}' tlable]});
     for s = 1:Seizures(1)
         if (Seizures(2*s)<= EpochTime(1) ) && (Seizures((2*s)+1) >= EpochTime(2) )
             title({Name,'Distances Network',['\color{red}' tlable]});
         end
     end

     subplot(2,2,2)

     bar(Ns)
     set(gca,'xticklabel',RegionLabels)
     ylim([0 1]);
     title({Name,'Average Connections',['\color{blue}' tlable]});
     for s = 1:Seizures(1)
         if (Seizures(2*s)<= EpochTime(1) ) && (Seizures((2*s)+1) >= EpochTime(2) )
             title({Name,'Average Connections',['\color{red}' tlable]});
         end
     end
     subplot(2,2,4)
     bar(VectorConnstrDist)
     set(gca,'xticklabel',RegionLabels)
     ylim([0 1]);
     title({Name,'Connection Strength',['\color{blue}' tlable]});
     for s = 1:Seizures(1)
         if (Seizures(2*s)<= EpochTime(1) ) && (Seizures((2*s)+1) >= EpochTime(2) )
             title({Name,'Connection Strength',['\color{red}' tlable]});
         end
     end

     set(gcf, 'Position', get(0, 'Screensize'))
 end

end

