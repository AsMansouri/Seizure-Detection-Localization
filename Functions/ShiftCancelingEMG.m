function [EEGDataShifted,EEGDataNoEMG] = ShiftCancelingEMG(EEGData, Type, BlockSize, WindowMove, Regions, RegionRef)
%Canceling the EMG by shifting and removing the average of signals
%   Detailed explanation goes here


it = floor((length(EEGData)-BlockSize)/WindowMove);

EEGDataShifted = zeros(size(EEGData));
EEGDataNoEMG = zeros(size(EEGData));
for i = 0:it
    clear tmp
    if i == 0 || i == it
        EEGDataShifted(:,(i*WindowMove)+1:((i*WindowMove)+BlockSize)) = EEGData(:,(i*WindowMove)+1:((i*WindowMove)+BlockSize));    
        EEGDataNoEMG(:,(i*WindowMove)+1:((i*WindowMove)+BlockSize)) = EEGData(:,(i*WindowMove)+1:((i*WindowMove)+BlockSize));    
    else
    if Type == 1 || 11
        tmpchann = [];
        tmp = EEGData(:,(i*WindowMove)+1:((i*WindowMove)+BlockSize));
        for r = 3:length(Regions)-1
            clear tmpRegionEEG
            tmpchann = [tmpchann Regions{r}(RegionRef(r))];
            EEGDataShifted(Regions{r}(RegionRef(r)),(i*WindowMove)+1:((i*WindowMove)+BlockSize)) = EEGData(Regions{r}(RegionRef(r)),(i*WindowMove)+1:((i*WindowMove)+BlockSize));
            for ch = 1:length(Regions{r})
                if ~ismember(Regions{r}(ch),tmpchann)
                    D = finddelay(tmp(Regions{r}(RegionRef(r)),:),tmp(ch,:),WindowMove);
                    EEGDataShifted(Regions{r}(ch),(i*WindowMove)+1:((i*WindowMove)+BlockSize)) = EEGData(Regions{r}(ch),(i*WindowMove)+1+D:((i*WindowMove)+BlockSize)+D);
                end
            end
            tmpch = [Regions{r}(RegionRef(r)) setdiff(Regions{r},tmpchann)];
            if length(tmpch) > 1
                tmpRegionEEG = EEGDataShifted(tmpch,(i*WindowMove)+1:((i*WindowMove)+BlockSize));
                tmpAvrg = mean(tmpRegionEEG);
                EEGDataNoEMG(tmpch,(i*WindowMove)+1:((i*WindowMove)+BlockSize)) = tmpRegionEEG - tmpAvrg;
            else
                EEGDataNoEMG(tmpch,(i*WindowMove)+1:((i*WindowMove)+BlockSize)) = EEGDataShifted(tmpch,(i*WindowMove)+1:((i*WindowMove)+BlockSize));
            end
            tmpchann = [tmpchann tmpch];
        end
    end
    end
end

end

