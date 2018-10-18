% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,logFile] = UiO_norm_power(data_struct,subj_name,EEG,logFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the 'after_wavelet' data (if available)
% subj_name: subject name according to csvfile
% logFile: logFile of previous function. If empty [] this function will
%       load the 'after_wavelet' logFile (if availeble)
%
% This function normalizes the EEG data (after wavelet decomposition) according
% to a baseline period (set in csv-file). Several options of normalization
% are possible: (1) ERD/ERS, (2) db, (3)
% 
% by questions:
% b.e.juel@medisin.uio.no
% sevenius.nilsen@gmail.com
% benjamin.thuerer@kit.edu
%

function [EEG,logFile] = UiO_norm_power(data_struct,subj_name,EEG,logFile)


if nargin < 2
    error('provide at least data_struct and subject name. See help UiO_norm_power')
end

% check if EEG structure is provided. If not, load previous data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,'after_wavelet');   
    else
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end

% get baseline range from csv-file. Note: length and start must be in ms
bl_range = [str2double(data_struct.baseline_start) (str2double(data_struct.baseline_start) + str2double(data_struct.baseline_length))];

[~,baselineidx(1)] = min(abs(EEG.times-bl_range(1)));
[~,baselineidx(2)] = min(abs(EEG.times-bl_range(2)));

%% compute baseline normalization according to csv options using median

baseline_power = median(EEG.data(:,:,baselineidx(1):baselineidx(2)),3);

%% new baseline normalization for DAVOS study
% baseline_power = EEG.bl;

%%

if str2double(data_struct.baseline_method) == 1
    
    % first option ERD/ERS using median    
    normChange = 100 * bsxfun(@minus,EEG.data,repmat(baseline_power,1,1,size(EEG.data,3)))./ repmat(baseline_power,1,1,size(EEG.data,3));
    norm_method = 'ERD';
elseif str2double(data_struct.baseline_method) == 2
    
    % second option decibel
    normChange = 10*log10(EEG.data./repmat(baseline_power,1,1,size(EEG.data,3)));
    norm_method = 'db';
end

EEG.data = abs(normChange);

logFile{end+1} = {'after_norm_power',['Power normalization is done using ' norm_method]};

if str2double(data_struct.plot_always)==1
    UiO_plots(data_struct,subj_name,EEG,logFile);
end

end