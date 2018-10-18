% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,logFile] = UiO_load_data(data_struct,subj_name,file_step)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% subj_name: subject name according to csvfile
% file_step: last step in EEG-analysis to load that file
%
% This function will load the eeg file and logFile of the last processing
% step in the EEG-analysis.
% 
% by questions:
% b.e.juel@medisin.uio.no
% sevenius.nilsen@gmail.com
% benjamin.thuerer@kit.edu
%
function [EEGF,LOCF] = UiO_load_data(data_struct,subj_name,file_step,specific_name)

% create file and loc name
if nargin < 4
    load_name = [subj_name{1} '_' data_struct.session '_' file_step];
    loc_name = [subj_name{1} '_' data_struct.session '_logFile'];
    loc_error = 0;
else
    load_name = data_struct.load_data;
    if strfind(load_name,'.mat')
        load_name = load_name(1:end-4);
    end
    % can not read loc file because loc file might not available or
    % different name then eeg-file...
    loc_error = 1;
    loc_name = [subj_name{1} '_' data_struct.session '_logFile'];
end

% check if file_save path is provided and check for / or \
if str2double(data_struct.save_folder) == 0
    if isempty(strfind(data_struct.vhdrsource,'\'))
        char_idx = strfind(data_struct.vhdrsource,'/'); 
    else
        char_idx = strfind(data_struct.vhdrsource,'\');
    end
    data_path = data_struct.vhdrsource(1:char_idx(end));
    load_file = [data_path load_name];
    load_loc = [data_path loc_name];
else
    if isempty(strfind(data_struct.save_folder,'\'))
        char_idx = strfind(data_struct.save_folder,'/'); 
    else
        char_idx = strfind(data_struct.save_folder,'\');
    end
    data_path = [data_struct.save_folder(1:char_idx(end))];
    load_file = [data_path subj_name{1} '\' data_struct.session '\'  load_name];
    load_loc = [data_path subj_name{1} '\' data_struct.session '\' loc_name];
end



disp(['load data: ' load_name]);
load([load_file '.mat']);

exist EEG;
if ans == 0
    loc_error = loc_error+1;
end

if loc_error == 0
    load([load_loc '.mat']);
    EEGF = EEG;
    LOCF = logFile;
elseif loc_error == 1
    EEGF = EEG;
    LOCF = [];
elseif loc_error == 2 
    % if file comes from SSP, variable is not called EEG
    sprt_name = strfind(load_name,'_');
    vrbl_name = load_name(1:sprt_name(1));
    data_ssp = eval('vrbl_name');
    
    data_ssp = rearrange(data_ssp,[3,2,1]);
    EEGF = struct;
    EEGF.data = data_ssp;
    EEGF.srate = 1000;
    EEGF.times = -size(data_ssp,2)/2:size(data_ssp,2)/2-1;
    LOCF = [];
end
end