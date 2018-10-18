% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,logFile] = UiO_pca(data_struct,subj_name,EEG,logFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the last processed data (if availeble)
% subj_name: subject name according to csvfile
% logFile: logFile of previous function. If empty [] this function will
%       load the last processed logFile (if availeble)
%
% This function will compress the data to 99% of explained variance. This
% step is important for later ICA-analysis
% 
% by questions:
% b.e.juel@medisin.uio.no
% sevenius.nilsen@gmail.com
% benjamin.thuerer@kit.edu
% 
function [EEG,logFile] = UiO_pca(data_struct,subj_name,EEG,logFile)

if nargin < 2
    error('provide at least data_struct and subject name. See help UiO_pca')
end


% check if EEG structure is provided. If not, load previous data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,'epoched');   
    else
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end


% change to double precision and check if continous or epoched data
ObsData = double(EEG.data);
if ndims(EEG.data) == 3
    ObsData = reshape(ObsData,size(ObsData,1),size(ObsData,2)*size(ObsData,3));
end

% remove the mean (zero mean), compute single value decomposition
% with the single values
ZeroData = bsxfun(@minus,ObsData,mean(ObsData,2));
[VecData,ValData] = svd(ZeroData); %single vectors and single values
ValData = diag(ValData*ValData'); %reak single values

% decompress to 99.9% of the variance
perVar = ones(1,size(ValData,1));
k = 100;
i = size(ValData,1);
while k > 99.99
    k = (sum(ValData(1:i))/sum(ValData(:)))*100;
    i = i-1;
end

EEG.lastPC = i+1;

% keep the 99.9% components and multiply them with the original data
% then decompress the data
PostPCAData = VecData(:,1:EEG.lastPC)' * ZeroData;
DecompData = (PostPCAData' * VecData(:,1:EEG.lastPC)')';

if ndims(EEG.data) == 3
    postCompData = reshape(DecompData,size(EEG.data,1),size(EEG.data,2),size(EEG.data,3));
else
    postCompData = DecompData;
end

EEG.data = postCompData;

disp(['Data compressed to ' num2str(EEG.lastPC) ' principal components']);

% loc file entry
logFile{end+1} = {'after_pca',['data is compressed to ' num2str(EEG.lastPC) ' principal components which explain ' ...
    '99.9% of the variance']};

if str2double(data_struct.plot_always)==1
    UiO_plots(data_struct,subj_name,EEG,logFile);
end

disp('data PCA is done')

end