function [ DataInfo ] = ReadDataInfo( datanum , varargin )
%function [ DataInfo ] = ReadDataInfo( datanum )
%   Get the info of a specific data from info.xlsx (one row)

if ~isempty(varargin)
    InfoAddress = varargin{1};
    if ~exist(InfoAddress, 'file')
        InfoAddress = '..\info.xlsx';
        warningMessage = sprintf('Warning: file does not exist:\n%s\n Using the default ''Info.xlsx'' file ...', InfoAddress);
        uiwait(msgbox(warningMessage));
    end
else
    InfoAddress = '..\info.xlsx';
end

[num,txt,raw] = xlsread(InfoAddress);
i =  datanum + 1;
DataInfo.Name = raw(i,1);
DataInfo.EEGMAPType = raw(i,7);
if cell2mat(DataInfo.EEGMAPType) ~= 10
%     DataInfo.Address = ['..\dataset\' cell2mat(raw(i,3)) '\' cell2mat(DataInfo.Name) '.edf'];
    DataInfo.Address = cell2mat(raw(i,3));
    DataInfo.Seizures(1) = raw(i,4);
    DataInfo.Seizures(2) = raw(i,5);
    DataInfo.Seizures(3) = raw(i,6);
    if cell2mat(DataInfo.Seizures(1)) > 1
        for j = 2:cell2mat(DataInfo.Seizures(1))
            DataInfo.Seizures(2*j) = raw(i,6+((j-1)*2));
            DataInfo.Seizures((2*j)+1) = raw(i,6+((j-1)*2)+1);
        end
    end
else
    [num1,txt1,info] = xlsread('D:\EEG Project\dataset\Karunya\Regions Info.xlsx');
    DataInfo.Seizures(1) = raw(i,4);
    DataInfo.Seizures(2) = raw(i,5);
    DataInfo.Seizures(3) = raw(i,6);
    DataInfo.Name = info(i,1);
    DataInfo.Age = info(i,2);
    DataInfo.Sex = info(i,3);
    DataInfo.ProvisionalDiagnosis = info(i,4);
    DataInfo.Disorder = info(i,5);
    DataInfo.SeizureType = info(i,6);
    DataInfo.Region = info(i,7);

end

end

