% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,locFile] = UiO_ica(data_struct,subj_name,EEG,locFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the last processed data (if availeble)
% subj_name: subject name according to csvfile
% locFile: locFile of previous function. If empty [] this function will
%       load the last processed locFile (if availeble)
%
% This function will perform ICA either on continous or epoched data.
% Please note that this may take a while to run. ICA tries to work with
% double-data. Please make sure that eeglab is running in double-precision
% mode: eeglab-->memory options
% 
% by questions: benjamin.thuerer@kit.edu
%
function [EEG,locFile] = UiO_ica(data_struct,subj_name,EEG,locFile)

if nargin < 2
    error('provide at least data_struct and subject name. See help UiO_ica')
end


% check if EEG structure is provided. If not, load previous data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,locFile] = UiO_load_data(data_struct,subj_name,'after_pca');   
    else
        [EEG,locFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end

% change to double precision
EEG.data = double(EEG.data);

% perform ICA on reduced rank after pca
if isfield(EEG,'lastPC')
    lastPC = EEG.lastPC;
    EEG = pop_runica(EEG,'pca',lastPC,'extended',1,'interupt','on');
else
    EEG = pop_runica(EEG,'extended',1,'interupt','on');
end


% loc file entry
locFile{end+1} = {'after_ica',['Independent components are computed and stored in the EEG struct']};

if str2double(data_struct.plot_always)==1
    UiO_plots(data_struct,subj_name,EEG,locFile);
end

end
