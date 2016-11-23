% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,locFile] = UiO_save(data_struct,subj_name,EEG,locFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% subj_name: subject name to store the file
% EEG: EEG structure of previous function.
% locFile: locFile of previous function.
% 
% This function will save the eeg file (.mat) and the loc file (.txt) in the raw
% data folder (if save_folder: 0 ) or in the provided folder path
% (save_folder)
% 
%
% by questions: benjamin.thuerer@kit.edu
% 
function [EEG,locFile] = UiO_save(data_struct,subj_name,EEG,locFile)

if nargin < 4
    error('this function needs data structure, subject name, processed EEG data, and locFile')
elseif isempty(EEG)
    error('this function needs processed EEG data. Run a processing step before this')
end

ending_name = locFile{end};
endName = ending_name{1};

% create file and loc name
save_name = [subj_name{1} '_' data_struct.session '_' endName];
loc_name = [subj_name{1} '_' data_struct.session '_locFile'];

% save file and loc to the same folder as the raw data if no save_folder
% path provided in csvfile
if str2double(data_struct.save_folder) == 0
    if isempty(strfind(data_struct.vhdrsource,'\'))
        char_idx = strfind(data_struct.vhdrsource,'/'); 
    else
        char_idx = strfind(data_struct.vhdrsource,'\');
    end
    data_path = data_struct.vhdrsource(1:char_idx(end));
    save_file = [data_path  save_name];
    save_loc = [data_path loc_name];
else
    % if save_folder path provided, check if / or \ is used and save file
    % and loc-file. If subfolders (name, session) are nonexistent, create
    % them
    subfolder1 = subj_name{1};
    subfolder2 = data_struct.session;
    
    if isempty(strfind(data_struct.save_folder,'\'))
        save_folder1 = [data_struct.save_folder '/' subfolder1];
        save_folder2 = [data_struct.save_folder '/' subfolder1 '/' subfolder2];
        if ~exist(save_folder1, 'dir')
          mkdir(save_folder1);
        end
        if ~exist(save_folder2, 'dir')
          mkdir(save_folder2);
        end
        save_file = [save_folder2 '/' save_name];
        save_loc = [save_folder2 '/' loc_name]; 
    else
        save_folder1 = [data_struct.save_folder '\' subfolder1];
        save_folder2 = [data_struct.save_folder '\' subfolder1 '\' subfolder2];
        if ~exist(save_folder1, 'dir')
          mkdir(save_folder1);
        end
        if ~exist(save_folder2, 'dir')
          mkdir(save_folder2);
        end
        save_file = [save_folder2 '\' save_name];
        save_loc = [save_folder2 '\' loc_name];
    end    
end

disp(['save file: ' save_name]);
% EEG.data = single(EEG.data);
save([save_file '.mat'],'EEG','-v7.3');
save([save_loc '.mat'],'locFile');

pause(0.5)

end