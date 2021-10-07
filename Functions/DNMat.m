function [ DistMatrix ] = DNMat( FilteredFFT, Electrodes)
%function [ DistMat ] = DNMat( FilteredFFT, Electrodes)
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
     DistMatrix = tmp;

end

