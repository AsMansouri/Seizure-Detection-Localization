function [AMIs] = AMI(Data,NLevel,MAXAMIK)
%[AMIs] = AMI(Data,NLevel,MAXAMIK)
%   Calculate the AMI by uniformly quantization. The Data should be a row
%   vector or a row matrix.
    
   Mins = min(Data,[],2);
   Maxs = max(Data,[],2);
   
   NormData = ((2.*Data)-(Maxs+Mins))./(Maxs-Mins);
   Data = NormData;
   
   NCol = size(Data,1);
   TotalSamp = size(Data,2);
   AMIs = zeros(NCol, MAXAMIK);


for i = 1: NCol
   Level = (Maxs(i)-Mins(i))/NLevel;
%    Levels = [-inf Mins(i)+Level:Level:Maxs(i)-Level inf];
   Levels = [Mins(i):Level:Maxs(i)];
   JointProb = zeros(NLevel);
   Prob = zeros(1,NLevel);
   Loc = cell(1,NLevel);
   for m = 1:NLevel
       tmp = find(Data(i,:) >= Levels(m) & Data(i,:) < Levels(m+1));
       Loc{m} = tmp;
%            Counters(i,m) = length(tmp);
       Prob(m) = length(tmp)/TotalSamp;
   end
   for k = 1:MAXAMIK
       for p = 1:NLevel    
           for q = 1:NLevel
               [Lia,~] = ismember(Loc{p},Loc{q}-k);
               JointProb(p,q) = sum(Lia)/(TotalSamp-k);
               if Prob(p)~= 0 && Prob(q)~=0 && JointProb(p,q)~=0
                   AMIs(i,k) = AMIs(i,k) + JointProb(p,q)*log2(JointProb(p,q)/(Prob(p)*Prob(q)));
               end
           end
       end
   end
   
end

end

