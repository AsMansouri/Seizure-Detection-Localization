function [ ] = PlotAllChann( EEGData,Fs,Name,Labels,Seizures,Visibility,varargin)
%function [ ] = PlotAllChann( EEGData,Fs,Name,Labels,Seizures,Visibility)
%function [ ] = PlotAllChann( EEGData,Fs,Name,Labels,Seizures,Visibility,header.starttime)
%function [ ] = PlotAllChann( EEGData,Fs,Name,Labels,Seizures,Visibility,header.starttime, Number of Time points)
%   Plot all the channels together and indicate the seizure periods with
%   dashed red lines. 
%   Given starttime, on x axis the time is shown.
%   The number of time points can be indicated after the starttime
   
     if ~isempty(varargin)
        if isa(varargin{1},'char')
            if length(varargin) == 2 && isa(varargin{2},'double')
                NumTimePonits = cell2mat(varargin(2))-1;
            else
                NumTimePonits = 10;   % default is 11 time points
            end
            ST1 = datetime(varargin(1),'Format','HH.mm.ss');
            Timeseries = seconds(0:(round(length(EEGData(1,:)))/Fs)/NumTimePonits:round(length(EEGData(1,:)))/Fs);
            TimeStringC = cell(NumTimePonits+1,1);
            for i=1:NumTimePonits+1
                TimeStringC{i} = datestr(ST1+Timeseries(i),'HH:MM:SS');
            end 
        elseif isa(varargin{1},'double')
            NumTimePonits = cell2mat(varargin(1))-1;
            Timeseries = seconds(0:(round(length(EEGData(1,:)))/Fs)/NumTimePonits:round(length(EEGData(1,:)))/Fs);
            TimeStringC = cell(NumTimePonits+1,1);
            for i=1:NumTimePonits+1
               TimeStringC{i} = string(Timeseries(i));
            end
        else
            NumTimePonits = 10;
            Timeseries = seconds(0:(round(length(EEGData(1,:)))/Fs)/NumTimePonits:round(length(EEGData(1,:)))/Fs);
            TimeStringC = cell(NumTimePonits+1,1);
            for i=1:NumTimePonits+1
               TimeStringC{i} = string(Timeseries(i));
            end
        end
     else
         NumTimePonits = 10;
         Timeseries = seconds(0:(round(length(EEGData(1,:)))/Fs)/NumTimePonits:round(length(EEGData(1,:)))/Fs);
         TimeStringC = cell(NumTimePonits+1,1);
         for i=1:NumTimePonits+1
            TimeStringC{i} = string(Timeseries(i));
         end
     end

     sig = flipud(EEGData);
     T = size(sig,2);
     t = linspace(0,round(T/Fs),T);

     % calculate shift
     mi = min(sig,[],2);
     ma = max(sig,[],2);
     shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]);
     shift = repmat(shift,1,T);  

     figure('Visible',Visibility)
     plot(t,sig+shift)
     title({Name})
        % Edit Axes
     set(gca,'ytick',mean(sig+shift,2),'yticklabel',fliplr(Labels))
     grid on
     ylim([mi(1) max(max(shift+sig))])      
        % Identify Seizures
     for S = 1: Seizures(1)    
        yL = get(gca,'YLim');
        line([Seizures(2*S) Seizures(2*S)],yL,'Color','r','LineStyle','--');
        line([Seizures(2*S+1) Seizures(2*S+1)],yL,'Color','r','LineStyle','--');
     end
     set(gca, 'Xtick',linspace(t(1),t(end),length(TimeStringC)),'XTickLabel',TimeStringC);
     xlabel('Time')
     %ylabel('Channels')
     xlim([t(1) t(end)])
     set(gcf, 'Position', get(0, 'Screensize'))
end

