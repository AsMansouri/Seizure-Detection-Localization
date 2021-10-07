function [DetectedType1,DetectedType2,DetectedType3,DetectedType4,DetectedRegion1,DetectedRegion2,DetectedRegion3,DetectedRegion4,DetectedRegion5,DetectedRegion6,DetectedRegion7,DetectedRegion8,DetectedRegion9] = EventLocalizationSubepoch4(FreqBand, mainepoch, EventFreqBand, RegionLabels)
%function [DetectedType, DetectedRegion] = EventLocalization(FreqBand, FeatureNum, RegionLabels , NNextEpochs)
%   Detailed explanation goes here

    DetectedRegion1 = [];
    DetectedRegion2 = [];
    DetectedRegion3 = [];
    DetectedRegion4 = [];
    DetectedRegion5 = [];
    DetectedRegion6 = [];
    DetectedRegion7 = [];
    DetectedRegion8 = [];
    DetectedRegion9 = [];
        
    f = 6;
    
    
    DetectedType1 = 'Focal Seizure';
    DetectedType2 = 'Focal Seizure';
    DetectedType3 = 'Focal Seizure';
    DetectedType4 = 'Focal Seizure';
    
    if FreqBand(f).VectorConn(mainepoch,end) > 0.6
        DetectedType1 = 'Generalized';
    end
    if FreqBand(f).VectorConn(mainepoch,end) > 0.5
        if FreqBand(f).VectorConnstr(mainepoch,1) > 0 && FreqBand(f).VectorConnstr(mainepoch,2) > 0
            if FreqBand(f).VectorConnstr(mainepoch,1)/FreqBand(f).VectorConnstr(mainepoch,2) < 2 
                DetectedType2 = 'Generalized';
            end
        end
    end
    
    f=8;
    if FreqBand(f).VectorConnDist(mainepoch,end) > 0.6
        DetectedType3 = 'Generalized';
    end
    if FreqBand(f).VectorConnstrDist(mainepoch,end) > 0.6
        DetectedType4 = 'Generalized';
    end
    
    
    
    
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
        DetectedRegion1 = [DetectedRegion1 '0->1 - '];
        epoch = 1;
    else
        DetectedRegion1 = [DetectedRegion1 num2str(PickedEvent-1) '->' num2str(PickedEvent) ' - '];
        epoch = PickedEvent;
    end
    f = 6;
    
    if not(strcmp(DetectedRegion1,'Generalized'))
        [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)).*EventFreqBand(f).VectorConn(epoch,(3:end-1)));
        r = I+2;
        if EventFreqBand(f).VectorConnstr(epoch,r) > 0.4
            if EventFreqBand(f).VectorConnstr(epoch,r) > EventFreqBand(f).VectorConn(epoch,r)
                DetectedRegion1 = [DetectedRegion1 ' -1- ' cell2mat(RegionLabels(r))];
            else
                [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)));
                r = I+2;
                DetectedRegion1 = [DetectedRegion1 ' -2- ' cell2mat(RegionLabels(r))];
            end
        else
            [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)));
            r = I+2;
            DetectedRegion1 = [DetectedRegion1 ' -3- ' cell2mat(RegionLabels(r))];
        end
    end

    PickedEvent = 0;
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
        DetectedRegion2 = [DetectedRegion2 '0->1 - '];
        epoch = 1;
    else
        DetectedRegion2 = [DetectedRegion2 num2str(PickedEvent) ' - '];
        epoch = PickedEvent;
    end
    f = 6;

    if not(strcmp(DetectedRegion2,'Generalized'))
        [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)).*EventFreqBand(f).VectorConn(epoch,(3:end-1)));
        r = I+2;
        if EventFreqBand(f).VectorConnstr(epoch,r) > 0.4
            if EventFreqBand(f).VectorConnstr(epoch,r) > EventFreqBand(f).VectorConn(epoch,r)
                DetectedRegion2 = [DetectedRegion2 ' -1- ' cell2mat(RegionLabels(r))];
            else
                [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)));
                r = I+2;
                DetectedRegion2 = [DetectedRegion2 ' -2- ' cell2mat(RegionLabels(r))];
            end
        else
            [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)));
            r = I+2;
            DetectedRegion2 = [DetectedRegion2 ' -3- ' cell2mat(RegionLabels(r))];
        end
    end
    
    
    
    
    
    
    PickedEvent = 0;
    if PickedEvent == 0   
        [~,PickedEvent] = max(max(EventFreqBand(f).VectorConn(:,3:end-1),[],2));
    end
    if PickedEvent == 0
        DetectedRegion3 = [DetectedRegion3 '0->1 - '];
        epoch = 1;
    else
        DetectedRegion3 = [DetectedRegion3 num2str(PickedEvent) ' - '];
        epoch = PickedEvent;
    end
    f = 6;

    if not(strcmp(DetectedRegion3,'Generalized'))
        [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)).*EventFreqBand(f).VectorConn(epoch,(3:end-1)));
        r = I+2;
        if EventFreqBand(f).VectorConnstr(epoch,r) > 0.4
            if EventFreqBand(f).VectorConnstr(epoch,r) > EventFreqBand(f).VectorConn(epoch,r)
                DetectedRegion3 = [DetectedRegion3 ' -1- ' cell2mat(RegionLabels(r))];
            else
                [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)));
                r = I+2;
                DetectedRegion3 = [DetectedRegion3 ' -2- ' cell2mat(RegionLabels(r))];
            end
        else
            [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)));
            r = I+2;
            DetectedRegion3 = [DetectedRegion3 ' -3- ' cell2mat(RegionLabels(r))];
        end
    end
    
    
    
    f = 6;

    if ((EventFreqBand(f).VectorConn(epoch,1)*EventFreqBand(f).VectorConnstr(epoch,1)) / ...
        (EventFreqBand(f).VectorConn(epoch,2)*EventFreqBand(f).VectorConnstr(epoch,2)) > 0.9)  && ...
       ((EventFreqBand(f).VectorConn(epoch,1)*EventFreqBand(f).VectorConnstr(epoch,1)) / ...
        (EventFreqBand(f).VectorConn(epoch,2)*EventFreqBand(f).VectorConnstr(epoch,2)) < 1.1)

        DetectedRegion4 = [DetectedRegion4 'Bilateral - '];
    else
        if EventFreqBand(f).VectorConn(epoch,1)*EventFreqBand(f).VectorConnstr(epoch,1) > EventFreqBand(f).VectorConn(epoch,2)*EventFreqBand(f).VectorConnstr(epoch,2)  %-Right
            DetectedRegion4 = [DetectedRegion4 cell2mat(RegionLabels(1)) ' - '];
        else
            DetectedRegion4 = [DetectedRegion4 cell2mat(RegionLabels(2)) ' - '];
        end
    end


    
        f = 6;

    if ((EventFreqBand(f).VectorConn(epoch,1)*EventFreqBand(f).VectorConnstr(epoch,1)) / ...
        (EventFreqBand(f).VectorConn(epoch,2)*EventFreqBand(f).VectorConnstr(epoch,2)) > 0.8)  && ...
       ((EventFreqBand(f).VectorConn(epoch,1)*EventFreqBand(f).VectorConnstr(epoch,1)) / ...
        (EventFreqBand(f).VectorConn(epoch,2)*EventFreqBand(f).VectorConnstr(epoch,2)) < 1.25)

        DetectedRegion4 = [DetectedRegion4 'Bilateral - '];
    else
        if EventFreqBand(f).VectorConn(epoch,1)*EventFreqBand(f).VectorConnstr(epoch,1) > EventFreqBand(f).VectorConn(epoch,2)*EventFreqBand(f).VectorConnstr(epoch,2)  %-Right
            DetectedRegion4 = [DetectedRegion4 cell2mat(RegionLabels(1)) ' - '];
        else
            DetectedRegion4 = [DetectedRegion4 cell2mat(RegionLabels(2)) ' - '];
        end
    end
    
    
    
    f = 6;
    if (FreqBand(f).VectorConn(mainepoch,1) / FreqBand(f).VectorConn(mainepoch,2)) > 0.83  && ...
       (FreqBand(f).VectorConn(mainepoch,1) / FreqBand(f).VectorConn(mainepoch,2)) < 1.2

        DetectedRegion5 = [DetectedRegion5 'Bilateral - '];
    else
        if FreqBand(f).VectorConn(mainepoch,1) > FreqBand(f).VectorConn(mainepoch,2)  %-Right
            DetectedRegion5 = [DetectedRegion5 cell2mat(RegionLabels(1)) ' - '];
        else
            DetectedRegion5 = [DetectedRegion5 cell2mat(RegionLabels(2)) ' - '];
        end
    end

        f = 6;
    if (FreqBand(f).VectorConnstr(mainepoch,1) / FreqBand(f).VectorConnstr(mainepoch,2)) > 0.83  && ...
       (FreqBand(f).VectorConnstr(mainepoch,1) / FreqBand(f).VectorConnstr(mainepoch,2)) < 1.2

        DetectedRegion6 = [DetectedRegion6 'Bilateral - '];
    else
        if FreqBand(f).VectorConnstr(mainepoch,1) > FreqBand(f).VectorConnstr(mainepoch,2)  %-Right
            DetectedRegion6 = [DetectedRegion6 cell2mat(RegionLabels(1)) ' - '];
        else
            DetectedRegion6 = [DetectedRegion6 cell2mat(RegionLabels(2)) ' - '];
        end
    end
    
    f = 8;
    if (FreqBand(f).VectorConnDist(mainepoch,1) / FreqBand(f).VectorConnDist(mainepoch,2)) > 0.83  && ...
       (FreqBand(f).VectorConnDist(mainepoch,1) / FreqBand(f).VectorConnDist(mainepoch,2)) < 1.2

        DetectedRegion7 = [DetectedRegion7 'Bilateral - '];
    else
        if FreqBand(f).VectorConnDist(mainepoch,1) > FreqBand(f).VectorConnDist(mainepoch,2)  %-Right
            DetectedRegion7 = [DetectedRegion7 cell2mat(RegionLabels(1)) ' - '];
        else
            DetectedRegion7 = [DetectedRegion7 cell2mat(RegionLabels(2)) ' - '];
        end
    end
    
    f = 8;
    if (FreqBand(f).VectorConnstrDist(mainepoch,1) / FreqBand(f).VectorConnstrDist(mainepoch,2)) > 0.83  && ...
       (FreqBand(f).VectorConnstrDist(mainepoch,1) / FreqBand(f).VectorConnstrDist(mainepoch,2)) < 1.2

        DetectedRegion8 = [DetectedRegion8 'Bilateral - '];
    else
        if FreqBand(f).VectorConnstrDist(mainepoch,1) > FreqBand(f).VectorConnstrDist(mainepoch,2)  %-Right
            DetectedRegion8 = [DetectedRegion8 cell2mat(RegionLabels(1)) ' - '];
        else
            DetectedRegion8 = [DetectedRegion8 cell2mat(RegionLabels(2)) ' - '];
        end
    end
    
    f = 8;
    if ((FreqBand(f).VectorConnDist(mainepoch,1)*FreqBand(f).VectorConnstrDist(mainepoch,1)) / ...
        (FreqBand(f).VectorConnDist(mainepoch,2)*FreqBand(f).VectorConnstrDist(mainepoch,2)) > 0.8)  && ...
       ((FreqBand(f).VectorConnDist(mainepoch,1)*FreqBand(f).VectorConnstrDist(mainepoch,1)) / ...
        (FreqBand(f).VectorConnDist(mainepoch,2)*FreqBand(f).VectorConnstrDist(mainepoch,2)) < 1.25)

        DetectedRegion9 = [DetectedRegion9 'Bilateral - '];
    else
        if FreqBand(f).VectorConnDist(mainepoch,1)*FreqBand(f).VectorConnstrDist(mainepoch,1) > FreqBand(f).VectorConnDist(mainepoch,2)*FreqBand(f).VectorConnstrDist(mainepoch,2)  %-Right
            DetectedRegion9 = [DetectedRegion9 cell2mat(RegionLabels(1)) ' - '];
        else
            DetectedRegion9 = [DetectedRegion9 cell2mat(RegionLabels(2)) ' - '];
        end
    end


end

