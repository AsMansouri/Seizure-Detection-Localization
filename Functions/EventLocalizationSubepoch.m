function [DetectedType, DetectedRegion] = EventLocalizationSubepoch(EventFreqBand, RegionLabels)
%function [DetectedType, DetectedRegion] = EventLocalization(FreqBand, FeatureNum, RegionLabels , NNextEpochs)
%   Detailed explanation goes here

    DetectedType = [];
    DetectedRegion = [];
    f = 6;
    
    
    PickedEvent = 0;
    for e = 1:size(EventFreqBand(f).VectorConn,1)
        if sum(EventFreqBand(f).VectorConn(e,3:end-1)>0.4)
            PickedEvent = e;
            break;
        end
    end
    
    if PickedEvent == 0
        [~,PickedEvent] = max(max(EventFreqBand(f).VectorConn(:,3:end-1),[],2));
    end

    epoch = PickedEvent;
    

    if EventFreqBand(f).VectorConn(epoch,end) > 0.5
        DetectedType = 'Generalized';
        DetectedRegion = 'Generalized';
    else
        DetectedType = 'Focal Seizure';
        if EventFreqBand(f).VectorConn(epoch,1) / EventFreqBand(f).VectorConn(epoch,2) > 1.5   %-Right
            DetectedRegion = [DetectedRegion cell2mat(RegionLabels(1)) ' - '];
        else
            if EventFreqBand(f).VectorConn(epoch,2) / EventFreqBand(f).VectorConn(epoch,1) > 1.5   %-Left
                DetectedRegion = [DetectedRegion cell2mat(RegionLabels(2)) ' - '];
            else
                if EventFreqBand(f).VectorConn(epoch,1) == EventFreqBand(f).VectorConn(epoch,2)    %-Bilateral
                    DetectedRegion = [DetectedRegion 'Bilateral - '];
    %                                 DetectedType = ['Bilateral'];
                else
                    if EventFreqBand(f).VectorConn(epoch,1) > EventFreqBand(f).VectorConn(epoch,2)
                        DetectedRegion = [DetectedRegion cell2mat(RegionLabels(1)) ',' cell2mat(RegionLabels(2)) ' - '];                           
                    else
                        DetectedRegion = [DetectedRegion cell2mat(RegionLabels(2)) ',' cell2mat(RegionLabels(1)) ' - '];        
                    end               
                end                
            end                      
        end
    end

    if not(strcmp(DetectedRegion,'Generalized'))
%             [~,I] = max(FreqBand(f).VectorConnstr(epoch,(3:end-1)));
        [~,I] = max(EventFreqBand(f).VectorConnstr(epoch,(3:end-1)).*EventFreqBand(f).VectorConn(epoch,(3:end-1)));
        r = I+2;
        if EventFreqBand(f).VectorConnstr(epoch,r) > 0.4
            if EventFreqBand(f).VectorConnstr(epoch,r) > EventFreqBand(f).VectorConn(epoch,r)
                DetectedRegion = [DetectedRegion ' ' cell2mat(RegionLabels(r))];
            end
        end
    end

end

