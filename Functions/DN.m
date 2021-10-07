function [ VectorConnDist, VectorConnstrDist ] = DN( FilteredFFT, EpochTime, Thereshold_Dist, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures)
%function [ VectorConnDist, VectorConnstrDist ] = DN( FilteredFFT, EpochTime, Thereshold_Dist, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures)
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
 for r = 1:length(Regions)
     R = cell2mat(Regions(r));
     for k=1:length(R)
         for j=k+1:length(R)
%                         ConnStr(r) = ConnStr(r) + CoefCorr(R(k),R(j));
            if (DistMat(R(k),R(j)) < Thereshold_Dist) 
%                        Cat(r) = Cat(r) + (1 - Distn(R(k),R(j)));
               Cat = [Cat r];
               ConnStr(r) = ConnStr(r) + ((Thereshold_Dist-DistMat(R(k),R(j)))/Thereshold_Dist);
            end
         end
     end
%      ConnStr(r) = ConnStr(r)/nchoosek(length(cell2mat(Regions(r))),2);
 end

 C = categorical(Cat,[1:length(Regions)],RegionLabels');
 [Ns,Categories] = histcounts(C);

 for r = 1:length(ConnStr)
     if Ns(r) == 0
         VectorConnstrDist(r) = 0;
     else
         VectorConnstrDist(r) = ConnStr(r) ./ Ns(r);
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

