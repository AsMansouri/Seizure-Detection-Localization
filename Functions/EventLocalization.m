function [DetectedType, DetectedRegion] = EventLocalization(FreqBand, FeatureNum, RegionLabels , NNextEpochs)
%function [DetectedType, DetectedRegion] = EventLocalization(FreqBand, FeatureNum, RegionLabels , NNextEpochs)
%   Detailed explanation goes here

    DetectedType = [];
    DetectedRegion = [];
    f = 8;
    epoch = FreqBand(f).FeatureVectors(FeatureNum).time;
    
    if FeatureNum>1 && ((FreqBand(f).FeatureVectors(FeatureNum).time - FreqBand(f).FeatureVectors(FeatureNum-1).time)> NNextEpochs+1)

        f = 6;
        %                 ConnVector = FreqBand(f).VectorConn(epoch-NPrevEpochs:epoch+NNextEpochs,:);

        if FreqBand(f).VectorConn(epoch,end) > 0.5
            DetectedType = 'Generalized';
            DetectedRegion = 'Generalized';
        else
            DetectedType = 'Focal Seizure';
            if FreqBand(f).VectorConn(epoch,1) / FreqBand(f).VectorConn(epoch,2) > 1.5   %-Right
                DetectedRegion = [DetectedRegion cell2mat(RegionLabels(1)) ' - '];
            else
                if FreqBand(f).VectorConn(epoch,2) / FreqBand(f).VectorConn(epoch,1) > 1.5   %-Left
                    DetectedRegion = [DetectedRegion cell2mat(RegionLabels(2)) ' - '];
                else
                    if FreqBand(f).VectorConn(epoch,1) == FreqBand(f).VectorConn(epoch,2)    %-Bilateral
                        DetectedRegion = [DetectedRegion 'Bilateral - '];
        %                                 DetectedType = ['Bilateral'];
                    else
                        if FreqBand(f).VectorConn(epoch,1) > FreqBand(f).VectorConn(epoch,2)
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
            [~,I] = max(FreqBand(f).VectorConnstr(epoch,(3:end-1)).*FreqBand(f).VectorConn(epoch,(3:end-1)));
            r = I+2;
            if FreqBand(f).VectorConnstr(epoch,r) > 0.4
                if FreqBand(f).VectorConnstr(epoch,r) > FreqBand(f).VectorConn(epoch,r)
                    DetectedRegion = [DetectedRegion ' ' cell2mat(RegionLabels(r))];
                end
            end
        end
%         if not(strcmp(DetectedRegion,'Generalized'))
%             for r = 3:size(FreqBand(f).VectorConn,2)-1
%                 if FreqBand(f).VectorConnstr(epoch,r) > 0.4
%                     if FreqBand(f).VectorConnstr(epoch,r) > FreqBand(f).VectorConn(epoch,r)
%                         DetectedRegion = [DetectedRegion ' ' cell2mat(RegionLabels(r))];
%                     end
%                 end
%             end
%         end

    elseif FeatureNum == 1
        
        f = 6;

        if FreqBand(f).VectorConn(epoch,end) > 0.5
            DetectedType = 'Generalized';
            DetectedRegion = 'Generalized';
        else
            DetectedType = 'Focal Seizure';
            if FreqBand(f).VectorConn(epoch,1) / FreqBand(f).VectorConn(epoch,2) > 1.5   %-Right
                DetectedRegion = [DetectedRegion cell2mat(RegionLabels(1)) ' - '];
            else
                if FreqBand(f).VectorConn(epoch,2) / FreqBand(f).VectorConn(epoch,1) > 1.5   %-Left
                    DetectedRegion = [DetectedRegion cell2mat(RegionLabels(2)) ' - '];
                else
                    if FreqBand(f).VectorConn(epoch,1) == FreqBand(f).VectorConn(epoch,2)
                        DetectedRegion = [DetectedRegion 'Bilateral - '];
                    else
                        if FreqBand(f).VectorConn(epoch,1) > FreqBand(f).VectorConn(epoch,2)
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
            [~,I] = max(FreqBand(f).VectorConnstr(epoch,(3:end-1)).*FreqBand(f).VectorConn(epoch,(3:end-1)));
            r = I+2;
            if FreqBand(f).VectorConnstr(epoch,r) > 0.4
                if FreqBand(f).VectorConnstr(epoch,r) > FreqBand(f).VectorConn(epoch,r)
                    DetectedRegion = [DetectedRegion ' ' cell2mat(RegionLabels(r))];
                end
            end
        end
%         if not(strcmp(DetectedRegion,'Generalized'))
%             for r = 3:size(FreqBand(f).VectorConn,2)-1
%                 if FreqBand(f).VectorConnstr(epoch,r) > 0.4
%                     if FreqBand(f).VectorConnstr(epoch,r) > FreqBand(f).VectorConn(epoch,r)
%                         DetectedRegion = [DetectedRegion ' ' cell2mat(RegionLabels(r))];
%                     end
%                 end
%             end
%         end

    end
 
end

