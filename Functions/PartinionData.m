function [ Epoch ,FreqBand  ] = PartinionData( EEGData, Fs, BlockSize, WindowMove, FreqBand, Format,varargin )
%function [ Epoch ,FreqBand  ] = PartinionData( EEGData, Fs, BlockSize, WindowMove, FreqBand, Format )
%   Partinion Epochs, filtered data, and energies of bands. 
%   Formats are 1 = FFT, 2 = DCT, 3 = FFT & DCT.
%   
 if ~isempty(varargin)
     FlagWindow = 1;
 else
     FlagWindow = 0;
 end

 df=Fs/BlockSize;
 it = floor((length(EEGData)-BlockSize)/WindowMove);
 
 switch Format
     case 1   
       for i = 0:it
         clear tmp tmp1 tmp2 Wtmp
         tmp = EEGData(:,(i*WindowMove)+1:((i*WindowMove)+BlockSize));
         Epoch.RawData(i+1,:,:) = tmp;
%          tmpn = normr(tmp);
         if FlagWindow 
%              Wtmp = tmp .* repmat(hamming(size(tmp,2))',size(tmp,1),1);
             Wtmp = tmp .* repmat(tukeywin(size(tmp,2))',size(tmp,1),1);
             tmp1 = abs(fft(Wtmp,size(Wtmp,2),2));  
         else
             tmp1 = abs(fft(tmp,size(tmp,2),2));
         end
         
         Epoch.FFTData(i+1,:,:) = tmp1;
                  
         for f = 1:size(FreqBand,2)            
             FreqRange=FreqBand(f).Range;           
             tmp2 = tmp1(:,round((FreqRange(1)/df)):round((FreqRange(2)/df)));
%              tmp2n = normr(tmp2);
             FreqBand(f).FilteredDataFFT(i+1,:,:) = tmp2;          
             FreqBand(f).TimeEnergies(i+1,:) = sum(tmp.^2,2);         
             FreqBand(f).FFTEnergies(i+1,:) = sum(tmp2.^2,2);         
             Epoch.Period(i+1,:) = [((i)*WindowMove/Fs) ((i*WindowMove+BlockSize)/Fs)];
             
         end
       end
       
     case 2
       for i = 0:it
         clear tmp tmp1dct tmp2dct
         tmp = EEGData(:,(i*WindowMove)+1:((i*WindowMove)+BlockSize));
         Epoch.RawData(i+1,:,:) = tmp;
%          tmpn = normr(tmp);
         for o = 1:size(EEGData,1)
            tmp1dct(o,:) = abs(dct(tmp(o,:)));
         end
         Epoch.DCTData(i+1,:,:) = tmp1dct;
                  
         for f = 1:size(FreqBand,2)
             FreqRange=FreqBand(f).Range;
             tmp2dct = tmp1dct(:,2*round((FreqRange(1)/df)):2*round((FreqRange(2)/df)));
%              tmp2n = normr(tmp2);
             FreqBand(f).FilteredDataDCT(i+1,:,:) = tmp2dct;
             FreqBand(f).TimeEnergies(i+1,:) = sum(tmp.^2,2);
             FreqBand(f).DCTEnergies(i+1,:) = sum(tmp2dct.^2,2);
             Epoch.Period(i+1,:) = [((i)*WindowMove/Fs) ((i*WindowMove+BlockSize)/Fs)];
         end
       end
     case 3
       for i = 0:it
         clear tmp tmp1 tmp2 tmp1dct tmp2dct
         tmp = EEGData(:,(i*WindowMove)+1:((i*WindowMove)+BlockSize));
         Epoch.RawData(i+1,:,:) = tmp;
%          tmpn = normr(tmp);

         tmp1 = abs(fft(tmp,size(tmp,2),2));  
         Epoch.FFTData(i+1,:,:) = tmp1;
         
         for o = 1:size(EEGData,1)
            tmp1dct(o,:) = abs(dct(tmp(o,:)));
         end
         Epoch.DCTData(i+1,:,:) = tmp1dct;
                  
         for f = 1:size(FreqBand,2)
             
             FreqRange=FreqBand(f).Range;
             tmp2 = tmp1(:,round((FreqRange(1)/df)):round((FreqRange(2)/df)));
             tmp2dct = tmp1dct(:,2*round((FreqRange(1)/df)):2*round((FreqRange(2)/df)));
%              tmp2n = normr(tmp2);

             FreqBand(f).FilteredDATAFFT(i+1,:,:) = tmp2;          
             FreqBand(f).FilteredDataDCT(i+1,:,:) = tmp2dct;
             FreqBand(f).TimeEnergies(i+1,:) = sum(tmp.^2,2);
             FreqBand(f).FFTEnergies(i+1,:) = sum(tmp2.^2,2);
             FreqBand(f).DCTEnergies(i+1,:) = sum(tmp2dct.^2,2);

             Epoch.Period(i+1,:) = [((i)*WindowMove/Fs) ((i*WindowMove+BlockSize)/Fs)];
              
         end  
       end
 end
 
end

