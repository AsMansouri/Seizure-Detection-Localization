%% EEG Analysis

clear;
clc;

addpath(genpath('..\Functions'));


%% -------------Initiation 
 
 Figures = 'off';
 Visibility = 'off';
 OrigionFlag = false;
 
 Thereshold_Corr = 0.85;
 Thereshold_Dist = 0.15;
 Thereshold_PSD  = 5;
 
 EpochSize = 10;
 WindowsOverlap = 0.5 ;
 N = 17;                   % Number of epochs for normalization the Normalized Energies
 Format = 1;               % Formats are 1 = FFT, 2 = DCT, 3 = FFT & DCT.
 
 NPrevEpochs = 2;
 NNextEpochs = 2;
 
 MarginalDetectionTime = 2;
 
 XLLine = 1;
 XLName = ['Result\Results-',datestr(now,'HH-MM-SS'),'.xlsx'];
 
 BannedSequences = [437:478 518:535 594:618 654:710];      %--Their Age ranges 3<x<19
 BannedSequences2 = [437:478 518:535 654:677];             %--Their Age ranges 3 =< x <= 19
%% -------------------------

 Directory = 'Result\';
 if ~exist(Directory,'dir')
     mkdir('Result');
 end
 
 ConfigTxt([XLName(1:end-5) '-Config.txt'],Thereshold_Corr,Thereshold_Dist,Thereshold_PSD,EpochSize,WindowsOverlap,N,Format,NPrevEpochs,NNextEpochs,MarginalDetectionTime);
 
 
initialVars = who;
initialVars(end+1:end+3) = {'initialVars','datanum','f'};

clearvars('-except',initialVars{:})

%---Automaticly
for datanum = 5:145
    tic;
%     if ismember(datanum,[22:24 30:36 48:54 58:70 71:78])
%         continue;
%     end
    if ismember(datanum,BannedSequences2)
        continue;
    end
    if ismember(datanum,[22:24 30:36 58:70])
        continue;
    end
    if ismember(datanum,[22:24 30:36 40:44 48:54 58:70 71:78])
        continue;
    end
    
    %% --- Info of the Data ( Name, Type, Seizures , ...) 
    
    DataInfo = ReadDataInfo( datanum , '../Datasets/MITinfo.xlsx');
    Name = cell2mat(DataInfo.Name);
    Type = cell2mat(DataInfo.EEGMAPType); % DataSet & Channels Type
    Seizures = cell2mat(DataInfo.Seizures);
    PatientID = extractBefore(Name,'_');
    
    if Type == 4 || Type == 5
        continue;
    end
    
    %%--- Types Info: Electrodes, Their Locations & Regions 
    [Electrodes, Centers, Regions, RegionLabels]  = TypeInfo(Type);

    %%--- Reading EEG data
        % -- Including Preprocessing Filtering Data
    if Type ~= 10  %edf format
        [Data, Fs, Labels, recorddata, header] = ReadEDFformatEEG(DataInfo.Address, Electrodes);
    else
        FileAddress = ['..\dataset\Karunya\' Name '\' Name '.xlsx'];
        [Data, Fs, Labels] = ReadExcelformatEEG(FileAddress);
    end


     %%---Filtering Non EEG Channels
     EEGData = Data(Electrodes(:),:);
     
     BlockSize=round(EpochSize*Fs);
     WindowMove=round(EpochSize*(1-WindowsOverlap)*Fs);
 
    
     for i=2:length(Seizures)
        SeizureBlock(i-1) = Seizures(i)*Fs;
        Seizure(i-1) = round(SeizureBlock(i-1)/WindowMove);
     end
     

% -------Freuency Bands
FreqBand = FreqBand_init();


ST1 = datetime(header.starttime,'Format','HH.mm.ss');
NewST = datestr(ST1,'HH.MM.SS');

  %---Make folders
Directory = ['Result\' Name '\'];
if ~exist(Directory,'dir')
    mkdir(Directory);
end

%% Saving Raw EEG signals

%  % ---- Plot All channel EEG data
%  PlotAllChann( EEGData,Fs,Name,Labels,Seizures,Visibility,NewST );
%  saveas(gcf,[Directory,Name,'- EEG Signals.tif']);
%  close all force
% 
%  % ---- Raw Signal of Brain Regions
%  if ~exist([Directory '\Regions Raw Data'],'dir')
%      mkdir([Directory '\Regions Raw Data']);
%      for r = 1:length(Regions)
%          PlotAllChann( EEGData(cell2mat(Regions(r)),:),Fs,[Name ' - ' cell2mat(RegionLabels(r))],Labels(cell2mat(Regions(r))),Seizures,Visibility,NewST)
%          saveas(gcf,[Directory,'\Regions Raw Data\',Name,' - ',cell2mat(RegionLabels(r)),' Signals.tif']);
%          close all force
%      end
%  end
%% Partitioning the epochs

tmpPartitionName = [Name '-' num2str(EpochSize) '-' num2str(WindowsOverlap) '.mat'];
tmpPartitionFile = ['J:\EEG Project 2\Datasets\Partitions\' tmpPartitionName];

if exist(tmpPartitionFile,'file')
    disp(['------------Loading ' tmpPartitionFile '  Partitions--------------------']);
    load(tmpPartitionFile)
else
    disp('--------------Partitioning epochs ---------------------');
    [ Epoch, FreqBand ] = PartinionData( EEGData, Fs, BlockSize, WindowMove, FreqBand, Format);
end
clear tmpPartitionName tmpPartitionFile

%% Analysis  
FreqBand = PickVariables(FreqBand, Format);
FreqBand = Connection_init(FreqBand, RegionLabels);
FreqBand = TotalAvgSqrEngery1(FreqBand);

 f = 8; 
 FreqBand(f).TotEnergiesNorm = FreqBand(f).TotEnergies(1:min(2*N,length(FreqBand(f).TotEnergies)))/mean(FreqBand(f).TotEnergies(1:min(N,length(FreqBand(f).TotEnergies))));
 FreqBand(f).TotEnergiesNorm11 = FreqBand(f).TotEnergiesNorm;
 f = 6;
 FreqBand(f).TotEnergiesNorm = FreqBand(f).TotEnergies(1:min(2*N,length(FreqBand(f).TotEnergies)))/mean(FreqBand(f).TotEnergies(1:min(N,length(FreqBand(f).TotEnergies))));
 FreqBand(f).TotEnergiesNorm11 = FreqBand(f).TotEnergiesNorm;
 
 f =8;
 Directory = ['Result\' Name '\' 'Bands' '\' FreqBand(f).Name '\'];
 if ~exist(Directory,'dir')
     mkdir(Directory);
 end

 
 FeatureVectors = [];
 DetectedSeizure = [];
 features = 0;

 
 Thereshold_PSD_Adap = zeros(size(FreqBand(f).FilteredDATA,1),1);
 Thresh_Adapt_epochs = [1:2*N];
 Thereshold_PSD_Adap(1:2*N)= Thereshold_PSD;
 
%  for i = 2*N+1:length(Thereshold_PSD_Adap)
%     Thereshold_PSD_Adap(i) = mean(FreqBand(8).TotEnergiesNorm(1:i-1))*5;
%  end
%  figure
% plot(Thereshold_PSD_Adap)
% hold on
% plot(FreqBand(f).TotEnergiesNorm)
% hold off



 e = 2*N+1;
 while e <= length(FreqBand(f).TotEnergies)

     f = 8;
     [FreqBand] = NormalizedEnergies(FreqBand, f, e, N);
     
%      if FreqBand(f).TotEnergiesNorm(e) * FreqBand(f).TotEnergiesNorm2(e) > Thereshold_PSD
     if FreqBand(f).TotEnergiesNorm(e) > Thereshold_PSD_Adap(e-1)
    
         features = features + 1;
         FeatureVectors(features).time = e;
         FeatureVectors(features).type = 'FP';
         FeatureVectors(features).SeizureNum = [];
         FeatureVectors(features).PBI = FreqBand(f).TotEnergiesNorm(e);
         FeatureVectors(features).PBI1 = FreqBand(f).TotEnergiesNorm11(e);
         
%          Thereshold_PSD_Adap(e) = FreqBand(f).TotEnergiesNorm(e)*0.1 + Thereshold_PSD_Adap(e-1)*0.9;
%          Thereshold_PSD_Adap(e) = (mean(FreqBand(8).TotEnergiesNorm(1:e-(2*N)))*0.5 + mean(FreqBand(8).TotEnergiesNorm(e-(2*N)+1:e-(N+1)))*0.25 + mean(FreqBand(8).TotEnergiesNorm(e-(N):e))*0.25)*5;
         [Thereshold_PSD_Adap] = ThreshAdaptation(FreqBand, Thereshold_PSD_Adap, 8, e, N, 10);
		 
         for s = 1:Seizures(1)
%                  tmp = ismember([e-NPrevEpochs:e+NNextEpochs] , [Seizure((2*s)-1):Seizure(2*s)] );
             %---TP : before the seizure ending and in not too far from the beginning
             if (e >= Seizure((2*s)-1) - MarginalDetectionTime) && (e <= Seizure(2*s))   
                 FeatureVectors(features).type = 'TP';
                 FeatureVectors(features).SeizureNum = s;
                 DetectedSeizure = [DetectedSeizure s];
             end
         end
    
         giffrsnum = 0; 
         for i = max(e-NPrevEpochs,1):min(length(FreqBand(f).TotEnergies),e+NNextEpochs)
             
             
             
             giffrsnum = giffrsnum + 1;    
             if i > e
                 [FreqBand] = NormalizedEnergies(FreqBand, f, i, N);
                 [Thereshold_PSD_Adap] = ThreshAdaptation(FreqBand, Thereshold_PSD_Adap, 8, i, N, 10);
%                  Thereshold_PSD_Adap(i) = Thereshold_PSD_Adap(i-1);
             end
            %%------------DIST Network
            %--HighGammas
            f = 6;
            [FreqBand(f).VectorConnDist(i,:), FreqBand(f).VectorConnstrDist(i,:)] = DN1( FreqBand(f).FilteredDATA(i,:,:), Epoch.Period(i,:), Thereshold_Dist, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures);

            FeatureVectors(features).HGDN((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConnDist(i,:);
            FeatureVectors(features).HGDNstr((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConnstrDist(i,:);  
            
            if strcmp(Figures,'on')
                tlable=[ num2str(Epoch.Period(i,1) )  ' -- ' num2str(Epoch.Period(i,2))];
                saveas(gcf,[Directory,Name,'- ',FreqBand(f).Name,'- Correlation Network-',tlable,'.tif']);
                close all force
            end
            %%------------Correlation Network
            %--HighGammas
            f = 6;
            [FreqBand(f).VectorConn(i,:), FreqBand(f).VectorConnstr(i,:)] = CN1( FreqBand(f).FilteredDATA(i,:,:), FreqBand(f).Energies(i,:), Epoch.Period(i,:), Thereshold_Corr, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures);

            FeatureVectors(features).HGCN((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConn(i,:);
            FeatureVectors(features).HGCNstr((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConnstr(i,:);
            
            if strcmp(Figures,'on')
                tlable=[ num2str(Epoch.Period(i,1) )  ' -- ' num2str(Epoch.Period(i,2))];
                saveas(gcf,[Directory,Name,'- ',FreqBand(f).Name,'- Correlation Network-',tlable,'.tif']);
                close all force
            end
            %%------------DIST Network
          %--Theta-Alpha 
            f = 8;
            [FreqBand(f).VectorConnDist(i,:), FreqBand(f).VectorConnstrDist(i,:)] = DN1( FreqBand(f).FilteredDATA(i,:,:), Epoch.Period(i,:), Thereshold_Dist, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures);

            FeatureVectors(features).DN((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConnDist(i,:);
            FeatureVectors(features).DNstr((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConnstrDist(i,:);

%             GIFmaker( gcf, giffrsnum, Directory, [ 'DN - ' num2str(features),'- ',FreqBand(f).Name], 0.5);
            if strcmp(Figures,'on')
                tlable=[ num2str(Epoch.Period(i,1) )  ' -- ' num2str(Epoch.Period(i,2))];
                saveas(gcf,[Directory,Name,'- ',FreqBand(f).Name,'- Distances Network-',tlable,'.tif']);
                close all force
            end

            %%------------Correlation Network
          %--Theta-Alpha
            f = 8;
            [FreqBand(f).VectorConn(i,:), FreqBand(f).VectorConnstr(i,:)] = CN1( FreqBand(f).FilteredDATA(i,:,:), FreqBand(f).Energies(i,:), Epoch.Period(i,:), Thereshold_Corr, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures);

            FeatureVectors(features).CN((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConn(i,:);
            FeatureVectors(features).CNstr((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConnstr(i,:);

%             GIFmaker( gcf, giffrsnum, Directory, [ 'CN - ' num2str(features),'- ',FreqBand(f).Name], 0.5);

            if strcmp(Figures,'on')
                tlable=[ num2str(Epoch.Period(i,1) )  ' -- ' num2str(Epoch.Period(i,2))];
                saveas(gcf,[Directory,Name,'- ',FreqBand(f).Name,'- Correlation Network-',tlable,'.tif']);
                close all force
            end
         end

         %--DN verification of TP/FP
         f = 8;
         figure('Visible',Visibility)
         subplot(211)
         plot(FreqBand(f).VectorConnstrDist(max(e-NPrevEpochs,1):min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end));
         title('DN Theta Alpha- STR - General Connections');
         subplot(212)
         plot(FreqBand(f).VectorConnDist(max(e-NPrevEpochs,1):min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end));
         title('DN Theta Alpha- #Connections- General Connections');
         set(gcf, 'Position', get(0, 'Screensize'))
         saveas(gcf,[Directory,Name,'- ',FreqBand(f).Name,'- DN verification- ',num2str(features),' .tif']);
         close all force  
         
         f = 8;
%          tmp = FreqBand(f).VectorConnDist(max(e-NPrevEpochs,1):min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);
         tmp = FreqBand(f).VectorConnDist(e:min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);
         if sum(diff(tmp)>0) > (2/3*(length(tmp)-1))
             FeatureVectors(features).DNVerif1 = 'Yes';
         else
             FeatureVectors(features).DNVerif1 = 'No';
         end
         
         if sum(tmp>0.2) > (2/3*(length(tmp)))
             FeatureVectors(features).DNVerif2 = 'Yes';
         else
             FeatureVectors(features).DNVerif2 = 'No';
         end
         
         if sum(diff(tmp)>0) > (1/2*(length(tmp)-1))
             FeatureVectors(features).DNVerif3 = 'Yes';
         else
             FeatureVectors(features).DNVerif3 = 'No';
         end
         if sum(tmp>0.2) > (1/2*(length(tmp)))
             FeatureVectors(features).DNVerif4 = 'Yes';
         else
             FeatureVectors(features).DNVerif4 = 'No';
         end
         clear tmp

         f = 8;
         tmp = FreqBand(f).VectorConnDist(e-1:min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);         
         if sum(tmp>0.2) > (2/3*(length(tmp)))
             FeatureVectors(features).DNVerif222 = 'Yes';
         else
             FeatureVectors(features).DNVerif222 = 'No';
         end
         clear tmp;

         f = 8;
         tmp = FreqBand(f).VectorConn(e:min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);
         if sum(diff(tmp)>0) > (2/3*(length(tmp)-1))
             FeatureVectors(features).DNVerif5 = 'Yes';
         else
             FeatureVectors(features).DNVerif5 = 'No';
         end
         if sum(tmp>0.2) > (2/3*(length(tmp)))
             FeatureVectors(features).DNVerif6 = 'Yes';
         else
             FeatureVectors(features).DNVerif6 = 'No';
         end
         if sum(diff(tmp)>0) > (1/2*(length(tmp)-1))
             FeatureVectors(features).DNVerif7 = 'Yes';
         else
             FeatureVectors(features).DNVerif7 = 'No';
         end
         if sum(tmp>0.2) > (1/2*(length(tmp)))
             FeatureVectors(features).DNVerif8 = 'Yes';
         else
             FeatureVectors(features).DNVerif8 = 'No';
         end
         clear tmp

         f = 6;
         tmp = FreqBand(f).VectorConnDist(e:min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);
         if sum(diff(tmp)>0) > (2/3*(length(tmp)-1))
             FeatureVectors(features).DNVerif21 = 'Yes';
         else
             FeatureVectors(features).DNVerif21 = 'No';
         end
         
         if sum(tmp>0.2) > (2/3*(length(tmp)))
             FeatureVectors(features).DNVerif22 = 'Yes';
         else
             FeatureVectors(features).DNVerif22 = 'No';
         end
         
         if sum(diff(tmp)>0) > (1/2*(length(tmp)-1))
             FeatureVectors(features).DNVerif23 = 'Yes';
         else
             FeatureVectors(features).DNVerif23 = 'No';
         end
         if sum(tmp>0.2) > (1/2*(length(tmp)))
             FeatureVectors(features).DNVerif24 = 'Yes';
         else
             FeatureVectors(features).DNVerif24 = 'No';
         end
         clear tmp


         tmp = FreqBand(f).VectorConn(e:min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);
         if sum(diff(tmp)>0) > (2/3*(length(tmp)-1))
             FeatureVectors(features).DNVerif25 = 'Yes';
         else
             FeatureVectors(features).DNVerif25 = 'No';
         end
         if sum(tmp>0.2) > (2/3*(length(tmp)))
             FeatureVectors(features).DNVerif26 = 'Yes';
         else
             FeatureVectors(features).DNVerif26 = 'No';
         end
         if sum(diff(tmp)>0) > (1/2*(length(tmp)-1))
             FeatureVectors(features).DNVerif27 = 'Yes';
         else
             FeatureVectors(features).DNVerif27 = 'No';
         end
         if sum(tmp>0.2) > (1/2*(length(tmp)))
             FeatureVectors(features).DNVerif28 = 'Yes';
         else
             FeatureVectors(features).DNVerif28 = 'No';
         end
         clear tmp
         
         f = 6;
         tmp = FreqBand(f).VectorConnDist(e-1:min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);         
         if sum(tmp>0.2) > (2/3*(length(tmp)))
             FeatureVectors(features).DNVerif2222 = 'Yes';
         else
             FeatureVectors(features).DNVerif2222 = 'No';
         end
         clear tmp;
         
         
%          Test = [];
%          for ff = 1:12
%             Test(ff,:,:) = FreqBand(ff).FFTEnergies(max(1,e-NPrevEpochs):min(e+NNextEpochs,length(FreqBand(ff).FFTEnergies)),:);
%          end
%          X = tensor(Test);
%          for cp = 1:10
%             M1 = cp_als(X,cp);
%             lamdas = M1.lambda./M1.lambda(1);
%             if lamdas(end) <= 0.2
%                 break;
%             end
%          end
%          vizopts = {'PlotCommands',{'bar','line','bar'},...
%                     'ModeTitles',{'Bands','Epochs','Channels'},...
%                     'BottomSpace',0.10,'HorzSpace',0.04,'Normalize',0};
%          info1 = viz(M1,'Figure',1,vizopts{:});
%          set(gcf, 'Position', get(0, 'Screensize'), 'Visible',Visibility)
%          saveas(gcf,[Directory Name '- Tesnors - ' num2str(features) ' .tif']);
%          close all force  
% 
%          EpochFeat = M1.u{2};
%          if (EpochFeat(end-NNextEpochs+2,1)- EpochFeat(end-NNextEpochs+1,1)) > 0
%              FeatureVectors(features).Tensor1 = 'Yes';
%          else
%              FeatureVectors(features).Tensor1 = 'No';
%          end
%          if any(EpochFeat(:,end)<0)
%              FeatureVectors(features).TensorLast = 'No';
%          else
%              FeatureVectors(features).TensorLast = 'Yes';
%          end
% 
%          if (EpochFeat(end-NNextEpochs+2,1)- EpochFeat(end-NNextEpochs+1,1)) > 0
%              [B,I] = sort(abs(lamdas-.2),'ascend');
%              if any(EpochFeat(:,I(1))<0)
%                  FeatureVectors(features).TensorComb = 'No';
%              else
%                  FeatureVectors(features).TensorComb = 'Yes';
%              end
%          else
%              FeatureVectors(features).TensorComb = 'No';
%          end
% 
%          
%          tmpr = M1.u{3};
%          figure('Visible',Visibility)
%          for subs = 1:size(tmpr,2)
%              subplot(size(tmpr,2),2,(2*subs)-1)
%              bar(tmpr(:,subs));
%              subplot(size(tmpr,2),2,2*subs)
%              tmprr = [];
%              for r = 1:length(Regions)
%                  tmprr(r,subs) = sum(tmpr(Regions{r},subs))/length(Regions{r});
%              end
%              bar(tmprr(:,subs));
%              set(gca,'xticklabel',RegionLabels)
%              ylim([0 1]);
%          end
%          set(gcf, 'Position', get(0, 'Screensize'), 'Visible',Visibility)
%          saveas(gcf,[Directory Name '- Tesnors - ' num2str(features) ' - Regions.tif']);
%          close all force 

%              
%              tmpEpoch = [];
%              for o = -NNextEpochs:NPrevEpochs
%                  tmpEpoch = [tmpEpoch Epoch.RawData(e+o,:,:)];
%              end
%              tmpEpoch1(:,:) = tmpEpoch;
%              tmpAMIs = AMI(tmpEpoch1,16,100);
%              clear tmpEpoch tmpEpoch1
%              figure('Visible',Visibility)
%              plot(tmpAMIs')
%              title({Name,['AMIs' num2str(features)]});
%              set(gcf, 'Position', get(0, 'Screensize'))
%              saveas(gcf,[Directory Name '- AMIs - ' num2str(features) '.tif']);
%              close all force 

         e = e + NNextEpochs;
     else
%          Thereshold_PSD_Adap(e) = (mean(FreqBand(8).TotEnergiesNorm(1:e-(N+1)))*0.75 + mean(FreqBand(8).TotEnergiesNorm(e-(N+1):e-1))*0.25)*5;
%          Thresh_Adapt_epochs(end+1) = e;
         
         [Thereshold_PSD_Adap] = ThreshAdaptation(FreqBand, Thereshold_PSD_Adap, 8, e, N, 10);
%          Thereshold_PSD_Adap(e) = (mean(FreqBand(8).TotEnergiesNorm(1:e-(2*N)))*0.5 + mean(FreqBand(8).TotEnergiesNorm(e-(2*N)+1:e-(N+1)))*0.25 + mean(FreqBand(8).TotEnergiesNorm(e-(N):e))*0.25)*5;
         Thresh_Adapt_epochs(end+1) = e;
         
%          Thereshold_PSD_Adap(e) = mean(FreqBand(8).TotEnergiesNorm(Thresh_Adapt_epochs))*5;
%          Thresh_Adapt_epochs(end+1) = e;
     end
     e = e+1;
 end

 %---- FN

 tmpFN = ismember([1:Seizures(1)] , DetectedSeizure );
 for t = 1:length(tmpFN)
     if tmpFN(t) == 0
         e = round((Seizure((2*t)-1)+Seizure(2*t))/2);

         f = 8;
         features = features + 1;
         FeatureVectors(features).type = 'FN'; 
         FeatureVectors(features).time = e;
         FeatureVectors(features).SeizureNum = t;
         if e <= length(FreqBand(f).TotEnergiesNorm)
             FeatureVectors(features).PBI = FreqBand(f).TotEnergiesNorm(e);
             FeatureVectors(features).PBI1 = FreqBand(f).TotEnergiesNorm11(e);
         else
             FeatureVectors(features).PBI = -100;
             FeatureVectors(features).PBI1 = -100;
         end
             
         for i = max(e-NPrevEpochs,1):min(length(FreqBand(f).TotEnergies),e+NNextEpochs)

            %%------------DIST Network
            f = 6;
            [FreqBand(f).VectorConnDist(i,:), FreqBand(6).VectorConnstrDist(i,:)] = DN( FreqBand(f).FilteredDATA(i,:,:), Epoch.Period(i,:), Thereshold_Dist, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures);
            FeatureVectors(features).HGDN((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConnDist(i,:);
            FeatureVectors(features).HGDNstr((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConnstrDist(i,:);
            if strcmp(Figures,'on')
                tlable=[ num2str(Epoch.Period(i,1) )  ' -- ' num2str(Epoch.Period(i,2))];
                saveas(gcf,[Directory,Name,'- ',FreqBand(f).Name,'- Correlation Network-',tlable,'.tif']);
                close all force
            end
            %%------------Correlation Network
            f = 6;
            [FreqBand(f).VectorConn(i,:), FreqBand(f).VectorConnstr(i,:)] = CN( FreqBand(f).FilteredDATA(i,:,:), FreqBand(f).Energies(i,:), Epoch.Period(i,:), Thereshold_Corr, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures);
            FeatureVectors(features).HGCN((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConn(i,:);
            FeatureVectors(features).HGCNstr((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConnstr(i,:);
            if strcmp(Figures,'on')
                tlable=[ num2str(Epoch.Period(i,1) )  ' -- ' num2str(Epoch.Period(i,2))];
                saveas(gcf,[Directory,Name,'- ',FreqBand(f).Name,'- Correlation Network-',tlable,'.tif']);
                close all force
            end  

            %%------------DIST Network
            f = 8;
            [FreqBand(f).VectorConnDist(i,:), FreqBand(f).VectorConnstrDist(i,:)] = DN( FreqBand(f).FilteredDATA(i,:,:), Epoch.Period(i,:), Thereshold_Dist, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures);
            FeatureVectors(features).DN((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConnDist(i,:);
            FeatureVectors(features).DNstr((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConnstrDist(i,:);
%                     FeatureVectors(f,features,(i- max(e-5,1)+1),:) = FreqBand(f).VectorConnDist(i,:);
            if strcmp(Figures,'on')
                tlable=[ num2str(Epoch.Period(i,1) )  ' -- ' num2str(Epoch.Period(i,2))];
                saveas(gcf,[Directory,Name,'- ',FreqBand(f).Name,'- Distances Network-',tlable,'.tif']);
                close all force
            end                    
            %%------------Correlation Network
            f = 8;
            [FreqBand(f).VectorConn(i,:), FreqBand(f).VectorConnstr(i,:)] = CN( FreqBand(f).FilteredDATA(i,:,:), FreqBand(f).Energies(i,:), Epoch.Period(i,:), Thereshold_Corr, Regions, RegionLabels, Centers, Labels, Type, Electrodes, Seizures, Name, Visibility, Figures);
            FeatureVectors(features).CN((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConn(i,:);
            FeatureVectors(features).CNstr((i- max(e-NPrevEpochs,1)+1),:) = FreqBand(f).VectorConnstr(i,:);
%                   
            if strcmp(Figures,'on')
                tlable=[ num2str(Epoch.Period(i,1) )  ' -- ' num2str(Epoch.Period(i,2))];
                saveas(gcf,[Directory,Name,'- ',FreqBand(f).Name,'- Correlation Network-',tlable,'.tif']);
                close all force
            end

         end

         f = 8;
         figure('Visible',Visibility)
         subplot(211)
         plot(FreqBand(f).VectorConnstrDist(max(e-NPrevEpochs,1):min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end));
         title('DN Theta Alpha- STR - General Connections');
         subplot(212)
         plot(FreqBand(f).VectorConnDist(max(e-NPrevEpochs,1):min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end));
         title('DN Theta Alpha- #Connections- General Connections');
         set(gcf, 'Position', get(0, 'Screensize'))
         saveas(gcf,[Directory,Name,'- ',FreqBand(f).Name,'- DN verification- FN- ',num2str(features),' .tif']);
         close all force  

         
         f = 8;
%          tmp = FreqBand(f).VectorConnDist(max(e-NPrevEpochs,1):min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);
         tmp = FreqBand(f).VectorConnDist(e:min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);
         if sum(diff(tmp)>0) > (2/3*(length(tmp)-1))
             FeatureVectors(features).DNVerif1 = 'Yes';
         else
             FeatureVectors(features).DNVerif1 = 'No';
         end
         
         if sum(tmp>0.2) > (2/3*(length(tmp)))
             FeatureVectors(features).DNVerif2 = 'Yes';
         else
             FeatureVectors(features).DNVerif2 = 'No';
         end
         
         if sum(diff(tmp)>0) > (1/2*(length(tmp)-1))
             FeatureVectors(features).DNVerif3 = 'Yes';
         else
             FeatureVectors(features).DNVerif3 = 'No';
         end
         if sum(tmp>0.2) > (1/2*(length(tmp)))
             FeatureVectors(features).DNVerif4 = 'Yes';
         else
             FeatureVectors(features).DNVerif4 = 'No';
         end
         clear tmp

         f = 8;
         tmp = FreqBand(f).VectorConnDist(e-1:min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);         
         if sum(tmp>0.2) > (2/3*(length(tmp)))
             FeatureVectors(features).DNVerif222 = 'Yes';
         else
             FeatureVectors(features).DNVerif222 = 'No';
         end
         clear tmp;

         f = 8;
         tmp = FreqBand(f).VectorConn(e:min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);
         if sum(diff(tmp)>0) > (2/3*(length(tmp)-1))
             FeatureVectors(features).DNVerif5 = 'Yes';
         else
             FeatureVectors(features).DNVerif5 = 'No';
         end
         if sum(tmp>0.2) > (2/3*(length(tmp)))
             FeatureVectors(features).DNVerif6 = 'Yes';
         else
             FeatureVectors(features).DNVerif6 = 'No';
         end
         if sum(diff(tmp)>0) > (1/2*(length(tmp)-1))
             FeatureVectors(features).DNVerif7 = 'Yes';
         else
             FeatureVectors(features).DNVerif7 = 'No';
         end
         if sum(tmp>0.2) > (1/2*(length(tmp)))
             FeatureVectors(features).DNVerif8 = 'Yes';
         else
             FeatureVectors(features).DNVerif8 = 'No';
         end
         clear tmp

         f = 6;
         tmp = FreqBand(f).VectorConnDist(e:min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);
         if sum(diff(tmp)>0) > (2/3*(length(tmp)-1))
             FeatureVectors(features).DNVerif21 = 'Yes';
         else
             FeatureVectors(features).DNVerif21 = 'No';
         end
         
         if sum(tmp>0.2) > (2/3*(length(tmp)))
             FeatureVectors(features).DNVerif22 = 'Yes';
         else
             FeatureVectors(features).DNVerif22 = 'No';
         end
         
         if sum(diff(tmp)>0) > (1/2*(length(tmp)-1))
             FeatureVectors(features).DNVerif23 = 'Yes';
         else
             FeatureVectors(features).DNVerif23 = 'No';
         end
         if sum(tmp>0.2) > (1/2*(length(tmp)))
             FeatureVectors(features).DNVerif24 = 'Yes';
         else
             FeatureVectors(features).DNVerif24 = 'No';
         end
         clear tmp


         tmp = FreqBand(f).VectorConn(e:min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);
         if sum(diff(tmp)>0) > (2/3*(length(tmp)-1))
             FeatureVectors(features).DNVerif25 = 'Yes';
         else
             FeatureVectors(features).DNVerif25 = 'No';
         end
         if sum(tmp>0.2) > (2/3*(length(tmp)))
             FeatureVectors(features).DNVerif26 = 'Yes';
         else
             FeatureVectors(features).DNVerif26 = 'No';
         end
         if sum(diff(tmp)>0) > (1/2*(length(tmp)-1))
             FeatureVectors(features).DNVerif27 = 'Yes';
         else
             FeatureVectors(features).DNVerif27 = 'No';
         end
         if sum(tmp>0.2) > (1/2*(length(tmp)))
             FeatureVectors(features).DNVerif28 = 'Yes';
         else
             FeatureVectors(features).DNVerif28 = 'No';
         end
         clear tmp
         
         f = 6;
         tmp = FreqBand(f).VectorConnDist(e-1:min(length(FreqBand(f).TotEnergies),e+NNextEpochs),end);         
         if sum(tmp>0.2) > (2/3*(length(tmp)))
             FeatureVectors(features).DNVerif2222 = 'Yes';
         else
             FeatureVectors(features).DNVerif2222 = 'No';
         end
         clear tmp;

%          Test = [];
%          for ff = 1:12
%             Test(ff,:,:) = FreqBand(ff).FFTEnergies(e-NPrevEpochs:e+NNextEpochs,:);
%          end
%          X = tensor(Test);
%          for cp = 1:10
%             M1 = cp_als(X,cp);
%             lamdas = M1.lambda./M1.lambda(1);
%             if lamdas(end) <= 0.2
%                 break;
%             end
%          end
%          vizopts = {'PlotCommands',{'bar','line','bar'},...
%                     'ModeTitles',{'Bands','Epochs','Channels'},...
%                     'BottomSpace',0.10,'HorzSpace',0.04,'Normalize',0};
%          info1 = viz(M1,'Figure',1,vizopts{:});
%          set(gcf, 'Position', get(0, 'Screensize'), 'Visible',Visibility)
%          saveas(gcf,[Directory Name '- Tesnors - ' num2str(features) ' .tif']);
%          close all force  
% 
%          EpochFeat = M1.u{2};
%          if (EpochFeat(end-NNextEpochs+2,1)- EpochFeat(end-NNextEpochs+1,1)) > 0
%              FeatureVectors(features).Tensor1 = 'Yes';
%          else
%              FeatureVectors(features).Tensor1 = 'No';
%          end
%          if any(EpochFeat(:,end)<0)
%              FeatureVectors(features).TensorLast = 'No';
%          else
%              FeatureVectors(features).TensorLast = 'Yes';
%          end
% 
%          if (EpochFeat(end-NNextEpochs+2,1)- EpochFeat(end-NNextEpochs+1,1)) > 0
%              [B,I] = sort(abs(lamdas-.2),'ascend');
%              if any(EpochFeat(:,I(1))<0)
%                  FeatureVectors(features).TensorComb = 'No';
%              else
%                  FeatureVectors(features).TensorComb = 'Yes';
%              end
%          else
%              FeatureVectors(features).TensorComb = 'No';
%          end
% 
%          tmpr = M1.u{3};
%          figure('Visible',Visibility)
%          for subs = 1:size(tmpr,2)
%              subplot(size(tmpr,2),2,(2*subs)-1)
%              bar(tmpr(:,subs));
%              subplot(size(tmpr,2),2,2*subs)
%              tmprr = [];
%              for r = 1:length(Regions)
%                  tmprr(r,subs) = sum(tmpr(Regions{r},subs))/length(Regions{r});
%              end
%              bar(tmprr(:,subs));
%              set(gca,'xticklabel',RegionLabels)
%              ylim([0 1]);
%          end
%          set(gcf, 'Position', get(0, 'Screensize'), 'Visible',Visibility)
%          saveas(gcf,[Directory Name '- Tesnors - ' num2str(features) ' - Regions.tif']);
%          close all force 


%              tmpEpoch = [];
%              for o = -NNextEpochs:NPrevEpochs
%                  tmptmpEpoch(:,:) = Epoch.RawData(e+o,:,:);
%                  tmpEpoch = [tmpEpoch tmptmpEpoch];
%                  clear tmptmpEpoch
%              end
%              tmpAMIs = AMI(tmpEpoch,16,100);
%              
%              figure('Visible',Visibility)
%              plot(tmpAMIs')
%              title({Name,['AMIs' num2str(features)]});
%              set(gcf, 'Position', get(0, 'Screensize'))
%              saveas(gcf,[Directory Name '- AMIs - ' num2str(features) '.tif']);
%              close all force 



     end
 end

 
 f = 8;
 FreqBand(f).FeatureVectors = FeatureVectors;

 if ~isempty(FreqBand(f).FeatureVectors)
 for i = 1:length(FreqBand(f).FeatureVectors)
     detect = 'No'; % other criteria for enhancing the accuracy
%              if (abs(XYsvds(i,12))+1.2) > (abs(XYsvds(i,11)))
%                 detect = 'Yes';
%              else
%                  detect = 'No';
%              end
     ST1 = datetime(NewST,'Format','HH.mm.ss');
     ST = seconds(mean(Epoch.Period(FreqBand(f).FeatureVectors(i).time,1:2))) + ST1 ;
     Time = datestr(ST,'HH.MM.SS');
     BeginLatency = [];
     if ~isempty(FeatureVectors(i).SeizureNum)
         BeginLatency =  abs(Seizures(2*FeatureVectors(i).SeizureNum) - mean(Epoch.Period(FreqBand(f).FeatureVectors(i).time,1)));
     end
     
     DetectedType = [];
     DetectedRegion = [];
     
     if i>1 && ((FreqBand(f).FeatureVectors(i).time - FreqBand(f).FeatureVectors(i-1).time)~= NNextEpochs+1)
        epoch = FreqBand(f).FeatureVectors(i).time;
        f = 6;
        ConnVector = FreqBand(f).VectorConn(max(1,epoch-NPrevEpochs):min(epoch+NNextEpochs,length(FreqBand(f).VectorConn)),:);

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
            for r = 3:size(FreqBand(f).VectorConn,2)-1
                if FreqBand(f).VectorConnstr(epoch,r) > 0.4
                    if FreqBand(f).VectorConnstr(epoch,r) > FreqBand(f).VectorConn(epoch,r)
                        DetectedRegion = [DetectedRegion ' ' cell2mat(RegionLabels(r))];
                    end
                end
            end
        end

     elseif i == 1
        epoch = FreqBand(f).FeatureVectors(i).time;
        f = 6;
        ConnVector = FreqBand(f).VectorConn(max(1,epoch-NPrevEpochs):min(epoch+NNextEpochs,length(FreqBand(f).VectorConn)),:);

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
            for r = 3:size(FreqBand(f).VectorConn,2)-1
                if FreqBand(f).VectorConnstr(epoch,r) > 0.4
                    if FreqBand(f).VectorConnstr(epoch,r) > FreqBand(f).VectorConn(epoch,r)
                        DetectedRegion = [DetectedRegion ' ' cell2mat(RegionLabels(r))];
                    end
                end
            end
        end

     end

     f = 8;
     if XLLine == 1           
%                  xlswrite(XLName,{'ID','#','Time','Estimated Start (sec)','Estimated Start (#Epoch)','State', },FreqBand(f).Name,'A1');
        xlswrite(XLName,{'ID','Patient','#','Time','PBI','PBI11','Latency','Type','Source location','State','Type', 'Verifications:' },FreqBand(f).Name,'A1');                 
        XLLine = XLLine + 1;
     end
    
     xlswrite(XLName,{Name,PatientID,num2str(i),Time,num2str(FeatureVectors(i).PBI),num2str(FeatureVectors(i).PBI1),...
         int2str(BeginLatency),DetectedType,DetectedRegion,detect,FreqBand(f).FeatureVectors(i).type,...
         FeatureVectors(i).DNVerif1,FeatureVectors(i).DNVerif2,FeatureVectors(i).DNVerif3,FeatureVectors(i).DNVerif4,...
         FeatureVectors(i).DNVerif5,FeatureVectors(i).DNVerif7,FeatureVectors(i).DNVerif8,...
         FeatureVectors(i).DNVerif21,FeatureVectors(i).DNVerif22,FeatureVectors(i).DNVerif23,FeatureVectors(i).DNVerif24,...
         FeatureVectors(i).DNVerif27,FeatureVectors(i).DNVerif222,FeatureVectors(i).DNVerif2222},...
         FreqBand(f).Name,['A' int2str(XLLine)]);
     xlswrite(XLName,{Name,num2str(i),num2str((size(EEGData,2)/Fs)/60),num2str(Seizures(1))},'Info',['A' int2str(XLLine)]);
     XLLine = XLLine + 1;
 end
 else
     if XLLine == 1           
%                  xlswrite(XLName,{'ID','#','Time','Estimated Start (sec)','Estimated Start (#Epoch)','State', },FreqBand(f).Name,'A1');
        xlswrite(XLName,{'ID','Patient','#','Time','PBI','PBI11','Latency','Type','Source location','State','Type', 'Verifications:' },FreqBand(f).Name,'A1');                 
        XLLine = XLLine + 1;
        xlswrite(XLName,{Name,PatientID,num2str(0)},FreqBand(f).Name,['A' int2str(XLLine)]);  
        xlswrite(XLName,{Name,num2str(0),num2str((size(EEGData,2)/Fs)/60),num2str(Seizures(1))},'Info',['A' int2str(XLLine)]);
        XLLine = XLLine + 1;
     else
        xlswrite(XLName,{Name,PatientID,num2str(0)},FreqBand(f).Name,['A' int2str(XLLine)]);    
        xlswrite(XLName,{Name,num2str(0),num2str((size(EEGData,2)/Fs)/60),num2str(Seizures(1))},'Info',['A' int2str(XLLine)]);
        XLLine = XLLine + 1;
     end
 end
 XLLine = XLLine + 1;

 %% visualization
 f = 8;
 figure('Visible',Visibility)
 plot(FreqBand(f).TotEnergiesNorm)
 title([Name ' - PSD']);
 yL = get(gca,'YLim');
 for s = 1:Seizures(1)
   line([Seizure((2*s)-1) Seizure((2*s)-1)],yL,'Color','r','LineStyle',':');
   line([Seizure(2*s) Seizure(2*s)],yL,'Color','r','LineStyle',':');
 end
 xlim([1 length(FreqBand(f).TotEnergiesNorm)])
 line([1 length(FreqBand(f).TotEnergiesNorm)],[Thereshold_PSD Thereshold_PSD] ,'Color','k','LineStyle','--');
 ylabel([FreqBand(f).Name '-PSD']);
 xlabel('Epochs');
 set(gcf, 'Position', get(0, 'Screensize'))
 saveas(gcf,[Directory Name '- Powers.tif']);
 close all force  

  f = 8;
 figure('Visible',Visibility)
 plot(FreqBand(f).TotEnergiesNorm)
 hold on
 plot(Thereshold_PSD_Adap)
 title([Name ' - PSD & Adaptive Thresh']);
 yL = get(gca,'YLim');
 for s = 1:Seizures(1)
   line([Seizure((2*s)-1) Seizure((2*s)-1)],yL,'Color','r','LineStyle',':');
   line([Seizure(2*s) Seizure(2*s)],yL,'Color','r','LineStyle',':');
 end
 xlim([1 length(FreqBand(f).TotEnergiesNorm)])
 line([1 length(FreqBand(f).TotEnergiesNorm)],[Thereshold_PSD Thereshold_PSD] ,'Color','k','LineStyle','--');
 ylabel([FreqBand(f).Name '-PSD']);
 xlabel('Epochs');
 set(gcf, 'Position', get(0, 'Screensize'))
 saveas(gcf,[Directory Name '- Powers - AdaptiveThresh.tif']);
 close all force 
 
   f = 8;
 figure('Visible',Visibility)
 plot(FreqBand(f).TotEnergiesNorm11)
 hold on
 plot(Thereshold_PSD_Adap)
 title([Name ' - PSD & Adaptive Thresh1']);
 yL = get(gca,'YLim');
 for s = 1:Seizures(1)
   line([Seizure((2*s)-1) Seizure((2*s)-1)],yL,'Color','r','LineStyle',':');
   line([Seizure(2*s) Seizure(2*s)],yL,'Color','r','LineStyle',':');
 end
 xlim([1 length(FreqBand(f).TotEnergiesNorm11)])
 line([1 length(FreqBand(f).TotEnergiesNorm11)],[Thereshold_PSD Thereshold_PSD] ,'Color','k','LineStyle','--');
 ylabel([FreqBand(f).Name '-PSD']);
 xlabel('Epochs');
 set(gcf, 'Position', get(0, 'Screensize'))
 saveas(gcf,[Directory Name '- Powers - AdaptiveThresh1.tif']);
 close all force 
 
  f = 8;
 figure('Visible',Visibility)
 plot(FreqBand(f).TotEnergiesNorm11)
 title([Name ' - PSD 1']);
 yL = get(gca,'YLim');
 for s = 1:Seizures(1)
   line([Seizure((2*s)-1) Seizure((2*s)-1)],yL,'Color','r','LineStyle',':');
   line([Seizure(2*s) Seizure(2*s)],yL,'Color','r','LineStyle',':');
 end
 xlim([1 length(FreqBand(f).TotEnergiesNorm11)])
 line([1 length(FreqBand(f).TotEnergiesNorm11)],[Thereshold_PSD Thereshold_PSD] ,'Color','k','LineStyle','--');
 ylabel([FreqBand(f).Name '-PSD 1']);
 xlabel('Epochs');
 set(gcf, 'Position', get(0, 'Screensize'))
 saveas(gcf,[Directory Name '- Powers1.tif']);
 close all force  
 
 
%  Test = [];
%  for i = 1:12
%     Test(i,:,:) = FreqBand(i).FFTEnergies;
%  end
%  X = tensor(Test);
%  M1 = cp_als(X,4);
%  vizopts = {'PlotCommands',{'bar','line','bar'},...
%             'ModeTitles',{'Bands','Epochs','Channels'},...
%             'BottomSpace',0.10,'HorzSpace',0.04,'Normalize',0};
%  info1 = viz(M1,'Figure',1,vizopts{:});
%  set(gcf, 'Position', get(0, 'Screensize'), 'Visible',Visibility)
%  saveas(gcf,[Directory Name '- Tesnors - All.tif']);
%  close all force  

 figure('Visible',Visibility)
 plot(FreqBand(f).TotEnergiesNorm2)
 title([Name ' - PSD - Normalized']);
 yL = get(gca,'YLim');
 for s = 1:Seizures(1)
   line([Seizure((2*s)-1) Seizure((2*s)-1)],yL,'Color','r','LineStyle',':');
   line([Seizure(2*s) Seizure(2*s)],yL,'Color','r','LineStyle',':');
 end
 xlim([1 length(FreqBand(f).TotEnergiesNorm2)])
 line([1 length(FreqBand(f).TotEnergiesNorm2)],[.5 .5] ,'Color','k','LineStyle','--');
 ylabel([FreqBand(f).Name '-PSD - Norm']);
 xlabel('Epochs');
 set(gcf, 'Position', get(0, 'Screensize'))
 saveas(gcf,[Directory Name '- Powers - Norm.tif']);
 close all force  


 figure('Visible',Visibility)
 plot(FreqBand(f).TotEnergiesNorm2.*FreqBand(f).TotEnergiesNorm')
 title([Name ' - PSD - Norm * Norm2']);
 yL = get(gca,'YLim');
 for s = 1:Seizures(1)
   line([Seizure((2*s)-1) Seizure((2*s)-1)],yL,'Color','r','LineStyle',':');
   line([Seizure(2*s) Seizure(2*s)],yL,'Color','r','LineStyle',':');
 end
 xlim([1 length(FreqBand(f).TotEnergiesNorm2)])
 line([1 length(FreqBand(f).TotEnergiesNorm2)],[Thereshold_PSD Thereshold_PSD] ,'Color','k','LineStyle','--');
 ylabel([FreqBand(f).Name '-PSD - Norm*Norm2']);
 xlabel('Epochs');
 set(gcf, 'Position', get(0, 'Screensize'))
 saveas(gcf,[Directory Name '- Powers - Norm Norm2.tif']);
 close all force  

toc
end

disp("---------DONE---------------");

% end


%% Functions
% Functions inside this code
% * 1) FreqBand initials
function FreqBand = FreqBand_init()

    FreqBand(1).Range = [0.5 4];
    FreqBand(1).Name = 'Delta';
    FreqBand(2).Range = [4 8];
    FreqBand(2).Name = 'Theta';
    FreqBand(3).Range = [8 14];
    FreqBand(3).Name = 'Alpha';
    FreqBand(4).Range = [14 30];
    FreqBand(4).Name = 'Beta';
    FreqBand(5).Range = [30 80];
    FreqBand(5).Name = 'Gamma';
    FreqBand(6).Range = [80 128];
    FreqBand(6).Name = 'High-Gamma';
    FreqBand(7).Range = [0.5 8];
    FreqBand(7).Name = 'Delta-Theta'; 
    FreqBand(8).Range = [4 14];
    FreqBand(8).Name = 'Theta-Alpha'; 
    FreqBand(9).Range = [8 30];
    FreqBand(9).Name = 'Alpha-Beta';
    FreqBand(10).Range = [14 80];
    FreqBand(10).Name = 'Beta-Gamma';
    FreqBand(11).Range = [30 128];
    FreqBand(11).Name = 'Gammas';
    FreqBand(12).Range = [0.5 128];
    FreqBand(12).Name = 'All';
    
end

% * 2) Reading EEG file in ".edf" format
function [Data, Fs, Labels, recorddata, header] = ReadEDFformatEEG(FileAddress, Electrodes)
    [header, recorddata] = edfread(FileAddress);
    Fs=(header.samples(1))/(header.duration);
    %%---Labels
    for i=1:length(Electrodes)
        Labels{i} = header.label{Electrodes(i)};
    end
    Data = HighPassFilter( recorddata,Fs);
end

% * 3) Reading EEG file in ".xlsx" format
function [Data, Fs, Labels] = ReadExcelformatEEG(FileAddress)
    Fs = 256;  
    [num,txt,~] = xlsread(FileAddress);
    %%---Labels
    Labels = txt;
    Data = HighPassFilter( num',Fs );
end

% * 4) Initiate the connection network matrices (Zero matrices)
function [FreqBand] = Connection_init(FreqBand, RegionLabels)
    for f = 1:length(FreqBand)
        FreqBand(f).VectorConn = zeros(length(FreqBand(f).Energies),length(RegionLabels));
        FreqBand(f).VectorConnstr = zeros(length(FreqBand(f).Energies),length(RegionLabels));
        FreqBand(f).VectorConnDist = zeros(length(FreqBand(f).Energies),length(RegionLabels));
        FreqBand(f).VectorConnstrDist = zeros(length(FreqBand(f).Energies),length(RegionLabels));
    end
end

% * 5) Calculating the Normalized epoch Total energies of an epoch
%      base on the [2N to N] previous epochs 
function [FreqBand] = NormalizedEnergies(FreqBand, f, i, N)

FreqBand(f).TotEnergiesNorm11(i) =FreqBand(f).TotEnergies(i)/mean(FreqBand(f).TotEnergies(max(1,i-2*N):i-N));
% FreqBand(f).TotEnergiesNorm(i) =((FreqBand(f).TotEnergies(i)/mean(FreqBand(f).TotEnergies(i-2*N:i-N)))/mean(FreqBand(f).TotEnergiesNorm(1:2*N)));
% FreqBand(f).TotEnergiesNorm(i) =((FreqBand(f).TotEnergies(i)/mean(FreqBand(f).TotEnergies(i-2*N:i-N)))- min(FreqBand(f).TotEnergiesNorm(max(1,i-N):i-1)))/(max(FreqBand(f).TotEnergiesNorm(max(1,i-N):i-1))-min(FreqBand(f).TotEnergiesNorm(max(1,i-N):i-1)));
% FreqBand(f).TotEnergiesNorm(i) =((FreqBand(f).TotEnergies(i)/mean(FreqBand(f).TotEnergies(i-2*N:i-N)))- min(FreqBand(f).TotEnergiesNorm(max(1,i-2*N):i-1)))/(max(FreqBand(f).TotEnergiesNorm(max(1,i-2*N):i-1))-min(FreqBand(f).TotEnergiesNorm(max(1,i-2*N):i-1)));
% FreqBand(f).TotEnergiesNorm(i) =((FreqBand(f).TotEnergies(i)/mean(FreqBand(f).TotEnergies(i-2*N:i-N)))- min(FreqBand(f).TotEnergies(max(1,i-2*N):i-1)))/(max(FreqBand(f).TotEnergies(max(1,i-2*N):i-1))-min(FreqBand(f).TotEnergies(max(1,i-2*N):i-1)));
FreqBand(f).TotEnergiesNorm(i) =(FreqBand(f).TotEnergies(i)- min(FreqBand(f).TotEnergies(max(1,i-2*N):i-N)))/(max(FreqBand(f).TotEnergies(max(1,i-2*N):i-N))-min(FreqBand(f).TotEnergies(max(1,i-2*N):i-N)));

if FreqBand(f).TotEnergiesNorm(i) < 0
    FreqBand(f).TotEnergiesNorm(i) = 0;
end

% FreqBand(f).TotEnergiesNorm(i) = mean(FreqBand(f).TotEnergiesNorm1(i)+FreqBand(f).TotEnergiesNorm11(i));

tmpN = 50;
if i == 1 || i == 2
    FreqBand(f).TotEnergiesNorm2(i) = 0.1;
else
    tmpX = FreqBand(f).TotEnergies(max(1,i-tmpN):i);
    minim = min(tmpX);
    maxim = max(tmpX);
    FreqBand(f).TotEnergiesNorm2(i) = (tmpX(end) - minim) / (maxim - minim);
end

end

% * 6) Calculating the Average Total Energies of epochs 
%      base on the half of electrodes with less energies in each Freq band
% function [FreqBand] = TotalAvgSqrEngery(FreqBand)
% %          FreqBand(f).TotEnergies = mean(FreqBand(f).Energies.^2,2);
%     for f = 1:length(FreqBand)
%         tmp = sort(FreqBand(f).Energies.^2,2);
%         FreqBand(f).TotEnergies = mean(tmp(:,1:ceil(end/2)),2);  
%     end
% end

% * 6) Calculating the Average Total Energies of epochs 
%      base on the half of electrodes with less energies in each Freq band
function [FreqBand] = TotalAvgSqrEngery1(FreqBand)
%          FreqBand(f).TotEnergies = mean(FreqBand(f).Energies.^2,2);
    for f = 1:length(FreqBand)
        tmp = sort(FreqBand(f).Energies,2);
        FreqBand(f).TotEnergies = mean(tmp(:,1:ceil(end/(4/3))),2);  
    end
end


% * 7) Picking the FFT or DCT energies and the Filtered bands 
%      base on the defined "Format"
function [FreqBand] = PickVariables(FreqBand, Format)
    if Format == 1 || Format == 3
        for f = 1:length(FreqBand)
            FreqBand(f).Energies = FreqBand(f).FFTEnergies;
            FreqBand(f).FilteredDATA = FreqBand(f).FilteredDataFFT;
        end
    elseif Fromat == 2
        for f = 1:length(FreqBand)
            FreqBand(f).Energies = FreqBand(f).DCTEnergies;
            FreqBand(f).FilteredDATA = FreqBand(f).FilteredDataDCT;
        end
    end
end

% * 8) Write the config of the analysis to a txt file 
%      Shows the default values of the parameters and thresholds
function [] = ConfigTxt(Name,Thereshold_Corr,Thereshold_Dist,Thereshold_PSD,EpochSize,WindowsOverlap,N,Format,NPrevEpochs,NNextEpochs,MarginalDetectionTime)

fileID = fopen(Name,'w');

fprintf(fileID,'Thereshold_Corr = %.2f \r\n',Thereshold_Corr);
fprintf(fileID,'Thereshold_Dist = %.2f \r\n',Thereshold_Dist);
fprintf(fileID,'Thereshold_PSD = %.2f \r\n',Thereshold_PSD);
 
fprintf(fileID,'EpochSize = %.2f \r\n',EpochSize);
fprintf(fileID,'WindowsOverlap = %.2f \r\n',WindowsOverlap);
fprintf(fileID,'N = %.2f\r\n',N);
fprintf(fileID,'Format = %.2f\r\n',Format);

fprintf(fileID,'NPrevEpochs = %.2f\r\n',NPrevEpochs);
fprintf(fileID,'NNextEpochs = %.2f\r\n',NNextEpochs);

fprintf(fileID,'MarginalDetectionTime = %.2f\r\n',MarginalDetectionTime);

fclose(fileID);

end

% * 9) Draw the 2D representation of the AMIs of EEG signals
%      With SVD dimensionality reduction
function [] = AMISVDs (Features, Epoch,  NNextEpochs, NPrevEpochs, Path, Visibility)
    
    c=zeros(3,3);
    c(1, :) = [0,0,0]; % FN
    c(2, :) = [0,1,0]; % TP
    c(3, :) = [1,0,0]; % FP

    
figure('Visible',Visibility,'Position', get(0, 'Screensize'))
for f = 1:length(Features)
    e = Features(f).time;
    
    switch Features(f).type
        case 'FN'
            tmpc = [0,0,0];
        case 'TP'
            tmpc = [0,1,0];        
        case 'FP'
            tmpc = [1,0,0];
    end
    tmpEpoch = [];
    for o = -NNextEpochs:NPrevEpochs
        tmptmpEpoch(:,:) = Epoch.RawData(e+o,:,:);
        tmpEpoch = [tmpEpoch tmptmpEpoch];
        clear tmptmpEpoch
    end
    tmpAMIs = AMI(tmpEpoch,32,100);
    Features(f).AMIs = tmpAMIs;
    
    [U,S,V] = svd(tmpAMIs);
    
    %--2D
    SVDrec = U(:,1:2)*S(1:2,1:2)*V(1:2,1:2)';
    % scatter(SVDrec(:,1),SVDrec(:,2),25)
    scatter(SVDrec(:,1),SVDrec(:,2),25,tmpc,'filled')
    hold on   
end
    title('AMIs SVD - ' )
    xlabel('Lambda 1 ')
    ylabel('Lambda 1 ')
    saveas(gcf,[Path ' AMIs SVD.jpg']);
    close all force;
    
    for f = 1:length(Features)
        figure('Visible',Visibility,'Position', get(0, 'Screensize'))
        plot(Features(f).AMIs')
        title(['AMIs' num2str(f)]);
        saveas(gcf,[Path ' - AMIs - ' num2str(f) '.jpg']);
        close all force 
    end
end


function [Thereshold_PSD_Adap] = ThreshAdaptation(FreqBand, Thereshold_PSD_Adap, f, i, N, K)

%     Thereshold_PSD_Adap(i) = Thereshold_PSD_Adap(i-1);
    Thereshold_PSD_Adap(i) = ((mean(Thereshold_PSD_Adap(1:i-(2*N)))/K)*0.5 + ...
                              mean(FreqBand(f).TotEnergiesNorm(i-(2*N)+1:i-(N+1)))*0.25 +...
                              mean(FreqBand(f).TotEnergiesNorm(i-(N):i))*0.25)*K;

end
