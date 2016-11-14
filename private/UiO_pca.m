% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% Processing steps (according to Makotos preprocessing pipeline):
% 1: load subject data
% 3: high-pass filter (1hz)
% 4: downsample (1000hz)
% 5: import channel info
% 6: remove bad channels + subspace reconstruction (clean_rawdata)
% 7: interpolate bad channels
% 8: rereference to the average (add the reference channel?)
% 9: remove line noise (cleanline)
% 10: epoch data (-1.5 to 1.5)
% 11: cut TMS artifacts and interpolate by mean of baseline
% 12: rejecting bad epoches (excluding EOG channels)
% 13: centering + compression (subtracting the mean and use PCA)
% 14: Data with 99% of Variance goes into ICA with reduced rank
% 15: save Data after ICA
% 
%%
function EEG = UiO_pca(data_struct,sbj_name,locFile)

EEG = pop_epoch( EEG, {  'R128'  }, epochRange, 'newname', ' resampled epochs', 'epochinfo', 'yes');

%remove bad epochs
EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan-2],eegThreshEpoch(1),eegThreshEpoch(2),EEG.times(1)/1000,EEG.times(end)/1000,0,1);
EEG.accBadEpochs = [];
EEG.accBadEpochs = find(EEG.reject.rejthresh); %store rejected trials

% performing PCA
ObsData = double(EEG.data);
ObsData = reshape(ObsData,size(ObsData,1),size(ObsData,2)*size(ObsData,3));
ZeroData = bsxfun(@minus,ObsData,mean(ObsData,2));
CovData = cov(ZeroData');
[VecData,ValData] = eig(CovData);
VecData = VecData(:,end:-1:1); %eigenvectors
ValData = diag(ValData);
ValData = ValData(end:-1:1); %eigenvalues

% decompress to 99% of the variance
perVar = ones(1,size(ValData,1));
for k = 1:size(ValData,1)
    perVar(k) = (sum(ValData(1:k))/sum(ValData(:)))*100;
end
lastPC = find(diff(perVar > 99))+1;
EEG.lastPC = lastPC;

PostPCAData = VecData(:,1:lastPC)' * ObsData;
DecompData = (PostPCAData' * VecData(:,1:lastPC)')';

postCompData = reshape(DecompData,size(EEG.data,1),size(EEG.data,2),size(EEG.data,3));
EEG.data = postCompData;

save_file = [data_path data_name '_after_ICA'];
save(save_file,'EEG');


end