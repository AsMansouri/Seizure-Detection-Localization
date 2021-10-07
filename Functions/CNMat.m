function [ MatrixConn] = CNMat( FilteredFFT, Powers, Electrodes)
%function [ MatrixConn] = CNMat( FilteredFFT, Powers, Electrodes)
%   Find the correlation connections



    P = (Powers) / (max(Powers));

    S(:,:) = FilteredFFT;
    CoefCorr = zeros(length(Electrodes));
    for j = 1: length(Electrodes)
    for k =j+1:length(Electrodes)
        R = corrcoef(S(j,:),S(k,:));
        CoefCorr(j,k) = sqrt(((R(1,2)+1)/2).^3+mean([P(j) P(k)]).^2);
        CoefCorr(k,j) = sqrt(((R(1,2)+1)/2).^3+mean([P(j) P(k)]).^2);
    end
    end

    Cmin = min(CoefCorr(:));
    Cmax = max(CoefCorr(:));
    CoefCorr = (CoefCorr) / (Cmax);
     
    MatrixConn = CoefCorr;
             
  
end