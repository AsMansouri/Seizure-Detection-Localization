function [ Electrodes ] = EEGChannels( Type )
%function [ Electrodes ] = EEGChannels( Type )
%   Give the Avialable Electrodes

switch Type
     case 1  % Our Data   % ------ 24 & 25 are LOC and EKG
         Electrodes=[2:23 26:32];
     case 2  % ---- PhysioNet Data
         Electrodes=[1:18,20:22];
     case 3
         Electrodes = [[1:4] [6:9] [11:12] [14:17] [19:22] [24:28]];
     case 4
         Electrodes = [[1:3] [5:9] [11:13] [15:19] [21:23] [25:29]];
     case 5
         Electrodes = [[1:3] [5:9] [11:13] [15:19] [21:23] [25:29]];
     case 6
         Electrodes = [[1:4] [6:9] [11:12] [14:17] [19:22]];
     case 7
         Electrodes = [[1:4] [6:9] [11:13] [15:18] [20:23] [25:29] [31:38]];
     case 10
         Electrodes=[1:16];   % ---- Karunya Data
end

end

