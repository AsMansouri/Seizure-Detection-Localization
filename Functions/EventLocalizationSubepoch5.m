function [DetectedType,DetectedRegion] = EventLocalizationSubepoch5(FreqBand, mainepoch, EventFreqBand, RegionLabels, ff)
%function [DetectedType, DetectedRegion] = EventLocalization(FreqBand, FeatureNum, RegionLabels , NNextEpochs)
%   Detailed explanation goes here

    DetectedType = [];
    DetectedRegion = [];
    
    
    %% Generalized or Focal
    f = 6;
% f = ff;
    DetectedType = 'Focal Seizure';
    if FreqBand(f).VectorConn(mainepoch,end) > 0.6
        DetectedType = 'Generalized';
    end
    
   
    %% Sub Epoch 
    PickedEvent = 0;
    if PickedEvent == 0
        for e = 1:size(EventFreqBand(f).VectorConn,1)
            if sum(EventFreqBand(f).VectorConn(e,3:end-1)>0.4)
                PickedEvent = e;
                break;
            end
        end
    end
    
    if PickedEvent == 0
        potentialSubEpochs = find(diff(EventFreqBand(f).VectorConn(:,end))>0.1)+1;
        if ~isempty(potentialSubEpochs)
            for i = 1:length(potentialSubEpochs)
                PickedEvent = potentialSubEpochs(i);
                [~,I] = sort(EventFreqBand(f).VectorConn(PickedEvent,3:end-1),'descend');
                for Ii = 1:length(I)
                    if EventFreqBand(f).VectorConn(PickedEvent,I(Ii)+2) > 0.4
                       if EventFreqBand(f).VectorConnstr(PickedEvent,I(Ii)+2) > EventFreqBand(f).VectorConn(PickedEvent,I(Ii)+2)
                           break;
                       end
                    end
                end
            end
        end
    end
    
    if PickedEvent == 0   
        [~,PickedEvent] = max(max(EventFreqBand(f).VectorConn(:,3:end-1),[],2));
    end
    
    if PickedEvent == 0
        epoch = 1;
    else
        DetectedRegion = [DetectedRegion num2str(PickedEvent) ' - '];
        epoch = PickedEvent;
    end
    
    
    %% Right or Left or Bilateral
    f = 6;
% f = ff;
    RoL1 = (EventFreqBand(f).VectorConn(epoch,1)*EventFreqBand(f).VectorConnstr(epoch,1)) / ...
           (EventFreqBand(f).VectorConn(epoch,2)*EventFreqBand(f).VectorConnstr(epoch,2));
       
%     f = 8;
%     RoL2 = (FreqBand(f).VectorConnDist(mainepoch,1) / FreqBand(f).VectorConnDist(mainepoch,2));
%     RoL3 = (FreqBand(f).VectorConnstrDist(mainepoch,1) / FreqBand(f).VectorConnstrDist(mainepoch,2));
      
%     RoLmeans = mean([RoL1 RoL2 RoL3]);
    
    f = 6;
% f = ff;

    if ((EventFreqBand(f).VectorConn(epoch,1)*EventFreqBand(f).VectorConnstr(epoch,1)) / ...
        (EventFreqBand(f).VectorConn(epoch,2)*EventFreqBand(f).VectorConnstr(epoch,2)) > 0.9)  && ...
       ((EventFreqBand(f).VectorConn(epoch,1)*EventFreqBand(f).VectorConnstr(epoch,1)) / ...
        (EventFreqBand(f).VectorConn(epoch,2)*EventFreqBand(f).VectorConnstr(epoch,2)) < 1.1)

        DetectedRegion = [DetectedRegion 'Bilateral - '];
    else
        if EventFreqBand(f).VectorConn(epoch,1)*EventFreqBand(f).VectorConnstr(epoch,1) > EventFreqBand(f).VectorConn(epoch,2)*EventFreqBand(f).VectorConnstr(epoch,2)  %-Right
            DetectedRegion = [DetectedRegion cell2mat(RegionLabels(1)) ' - '];
        else
            DetectedRegion = [DetectedRegion cell2mat(RegionLabels(2)) ' - '];
        end
    end
    
   
    %% Source Region
    
    f = 6;
% f = ff;
    if not(strcmp(DetectedRegion,'Generalized'))
        [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)).*EventFreqBand(f).VectorConn(epoch,(3:end-1)));
        r = I+2;
        if EventFreqBand(f).VectorConnstr(epoch,r) > 0.4
            if EventFreqBand(f).VectorConnstr(epoch,r) > EventFreqBand(f).VectorConn(epoch,r)
                DetectedRegion = [DetectedRegion ' -1- ' cell2mat(RegionLabels(r))];
            else
                [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)));
                r = I+2;
                DetectedRegion = [DetectedRegion ' -2- ' cell2mat(RegionLabels(r))];
            end
        else
            [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)));
            r = I+2;
            DetectedRegion = [DetectedRegion ' -3- ' cell2mat(RegionLabels(r))];
        end
    end

    
end

