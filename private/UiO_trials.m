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
function [EEG,locFile] = UiO_trials(data_struct,subj_name,EEG,locFile)

if nargin < 2
    error('provide at least data_struct and subject name. See help UiO_trials')
end

% check if EEG structure is provided. If not, load preprocessed data
if isempty(EEG)
    [EEG,locFile] = UiO_load_data(data_struct,subj_name);
end
  
% epoch data to trials (EEG*times*trials)
EEG = pop_epoch( EEG, {  'R128'  }, [str2double(data_struct.trial_start)/EEG.srate str2double(data_struct.trial_end)/EEG.srate] ...
    , 'newname', ' resampled epochs', 'epochinfo', 'yes');

%remove bad epochs either automized (1) or by inspection (2)
if str2double(data_struct.trial_rejection) == 1
    EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan-2],-str2double(data_struct.trial_threshold),str2double(data_struct.trial_threshold),EEG.times(1)/1000,EEG.times(end)/1000,0,1);
elseif str2double(data_struct.trial_rejection) == 2
    MElec = squeeze(mean(EEG.data,1));
    [~,idx(1)] = min(abs(EEG.times-(-500)));
    [~,idx(2)] = min(abs(EEG.times-(500)));
    
    goodTrial = [];
    badTrial = [];

    % plot the average over channels for each trial and distinguish between good
    % and bad trials according to the mouse / space click
    for Ei = 1:size(MElec,2)
        h = figure;
        plot(EEG.times(idx(1):idx(2)),MElec(idx(1):idx(2),Ei));
        title(['Trial ' int2str(Ei) ': press space to remove trial']);
        grid
        axis([-500 500 -15 15]);
        ylabel('Amplitude (\muV)'); xlabel('Time (ms)');
        button = waitforbuttonpress;
        if button == 0
            goodTrial(end+1) = Ei;
        else
            badTrial(end+1) = Ei;
        end
        close(h)
    end

    % plot the average for each trial again and show the actual
    % marking of the trial. Press space if you changed your mind
    for Ei = 1:size(MElec,2)
        if find(goodTrial==Ei)
            Stigma = 'good';
        else
            Stigma = 'bad';
        end
        h = figure;
        plot(EEG.times(idx(1):idx(2)),MElec(idx(1):idx(2),Ei));
        title({['Trial ' int2str(Ei) ' marked as ' Stigma] ; ...
            ['if you want to change the mark (regardless of direction g-->b; b-->g) press space']});
        grid
        axis([-500 500 -15 15]);
        ylabel('Amplitude (\muV)'); xlabel('Time (ms)');
        button = waitforbuttonpress;
        if button ~= 0
            if find(badTrial==Ei)
                idx_Ei = badTrial == Ei;
                badTrial(idx_Ei) = [];
            else
                badTrial(end+1) = Ei;
            end
        end
        close(h)
    end

    % sort bad trials and remove bad trials from the dataset
    badTrial = sort(badTrial);
    EEG.nbchan = EEG.nbchan-length(badChan);
    EEG.data(badChan,:) = [];
else
    warning('trial rejection is not accuratly provided in csvfile. Rejection is done automatically')
end
    
    
EEG.accBadEpochs = [];
EEG.accBadEpochs = find(EEG.reject.rejthresh); %store rejected trials

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function will load, if EEG not provided, the relevant file and loc
% paths + names to load these data for subsequent processing.
function [EEG,LOCF] = UiO_load_data(data_struct,subj_name)

% create file and loc name
load_name = [subj_name{1} '_' data_struct.session '_preprocessed'];
loc_name = [subj_name{1} '_' data_struct.session '_locFile'];

% check if file_save path is provided and check for / or \
if str2double(data_struct.save_folder) == 0
    if isempty(strfind(data_struct.vhdrsource,'\'))
        char_idx = strfind(data_struct.vhdrsource,'/'); 
    else
        char_idx = strfind(data_struct.vhdrsource,'\');
    end
    data_path = data_struct.vhdrsource(1:char_idx(end));
    load_file = [data_path  load_name];
    load_loc = [data_path loc_name];
else
    % if save_folder path provided, check if / or \ is used
    subfolder1 = subj_name{1};
    subfolder2 = data_struct.session;
    
    if isempty(strfind(data_struct.save_folder,'\'))
        load_folder2 = [data_struct.save_folder '/' subfolder1 '/' subfolder2];
        load_file = [load_folder2 '/' load_name];
        load_loc = [load_folder2 '/' loc_name]; 
    else
        load_folder2 = [data_struct.save_folder '\' subfolder1 '\' subfolder2];
        load_file = [load_folder2 '\' load_name];
        load_loc = [load_folder2 '\' loc_name];
    end    
end

load([load_file '.mat']);
load([load_loc '.mat']);

EEG.data = double(EEG.data);
LOCF = locFile;

end
