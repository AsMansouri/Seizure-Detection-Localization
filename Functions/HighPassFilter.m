function [ FilteredData ] = HighPassFilter( PureData,Fs,varargin )
%function [ FilteredData ] = HighPassFilter( PureData,Fs )
%   10th Order Cheby IIR filter  HighPass filter with cutoff freq (0.1)Hz
%   Also in this function the 60Hz (59.5 to 60.5) noise is filtered
%   by an Elliptic filter
%
%function [ FilteredData ] = HighPassFilter( PureData,Fs, CutOff Freq )
%   10th Order Cheby IIR filter  HighPass filter with cutoff freq
%   of "CutOff Freq"
%   Also in this function the 60Hz (59.05 to 60.05) noise is filtered 
%   by an Elliptic filter


 %%---Filtering Data
order    = 6;
fpass  = (0.1); % default
if ~isempty(varargin) && isa(varargin{1},'double')
    fpass  = varargin{1};
end

% [b,a] = butter(order,fpass, 'high');			   

hpFilt = designfilt('highpassiir', 'FilterOrder', order, 'PassbandFrequency', fpass, ...
                    'StopbandAttenuation', 40, 'PassbandRipple', 0.05, ...
                    'SampleRate', Fs, 'DesignMethod', 'ellip');
			   

NotchFilt = designfilt('bandstopiir', 'FilterOrder', order, 'StopbandAttenuation', 40, ...
               'StopbandFrequency1', 59, 'StopbandFrequency2', 61,  ...
               'SampleRate', Fs, 'DesignMethod', 'cheby2');
			   
s=size(PureData);
for i=1:s(1)
    tmp = PureData(i,:);
    tmpHp = filtfilt(hpFilt,tmp);
    FilteredData(i,:) = filtfilt(NotchFilt,tmpHp);
end
end

