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
function EEG = UiO_ica(data_struct,sbj_name,locFile)

% performing ICA
EEG = pop_runica(EEG,'pca',lastPC,'extended',1,'interupt','on');
EEG.data = single(EEG.data);
% EEG.icaact = single(EEG.icaact);