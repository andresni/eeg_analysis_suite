% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,locFile] = UiO_pca(data_struct,subj_name,EEG,locFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the last processed data (if availeble)
% subj_name: subject name according to csvfile
% locFile: locFile of previous function. If empty [] this function will
%       load the last processed locFile (if availeble)
%
% This function will compress the data to 99% of explained variance. This
% step is important for later ICA-analysis
% 
% by questions: benjamin.thuerer@kit.edu
% 
function [EEG,locFile] = UiO_pca(data_struct,subj_name,EEG,locFile)

if nargin < 2
    error('provide at least data_struct and subject name. See help UiO_pca')
end


% check if EEG structure is provided. If not, load previous data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,locFile] = UiO_load_data(data_struct,subj_name,'epoched');   
    else
        [EEG,locFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end


% change to double precision and check if continous or epoched data
ObsData = double(EEG.data);
if ndims(EEG.data) == 3
    ObsData = reshape(ObsData,size(ObsData,1),size(ObsData,2)*size(ObsData,3));
end

% remove the mean (zero mean), compute the eigenvectors of the covariance
% matrix and sort eigenvectors and values in descending order
ZeroData = bsxfun(@minus,ObsData,mean(ObsData,2));
CovData = cov(ZeroData');
[VecData,ValData] = eig(CovData);
VecData = VecData(:,end:-1:1); %eigenvectors
ValData = diag(ValData);
ValData = ValData(end:-1:1); %eigenvalues

% decompress to 99.9% of the variance
perVar = ones(1,size(ValData,1));
for k = 1:size(ValData,1)
    perVar(k) = (sum(ValData(1:k))/sum(ValData(:)))*100;
end
lastPC = find(diff(perVar > 99.9))+1;
EEG.lastPC = lastPC;

% keep the 99% components and multiply them with the original data
% then decompress the data
PostPCAData = VecData(:,1:lastPC)' * ObsData;
DecompData = (PostPCAData' * VecData(:,1:lastPC)')';

if ndims(EEG.data) == 3
    postCompData = reshape(DecompData,size(EEG.data,1),size(EEG.data,2),size(EEG.data,3));
end

EEG.data = postCompData;

disp('Data compressed to ' num2str(lastPC) ' principal components');

% loc file entry
locFile{end+1} = {'after_pca',['data is compressed to ' num2str(lastPC) ' principal components which explain ' ...
    '99.9% of the variance']};

if str2double(data_struct.plot_always)==1
    UiO_plots(data_struct,subj_name,EEG,locFile);
end

end