    % EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,locFile] = UiO_TMS_correction(data_struct,~,~,locFile)
% 
% data_struct: structure of the csv-file specified for subject and
%           experiment
% ~: sbj_name and EEG not necessary for this function
% locFile: not necessary for this function as it is run first
% 
% this function interpolates stimulus artifacts induced by TMS (or other
% stimuli). According to the csv-file, the interpolation is done by
% 1: automatic for each trial according to the value (x > mean + std)
% 2: manually by inspection of the gran average plot over all trials and
% channels
% 3: manually by the values provided in the csv-file
%
% by questions: benjamin.thuerer@kit.edu
% 
function [EEG,locFile] = UiO_TMS_correction(data_struct,subj_name,EEG,locFile)

if nargin < 1
    error('provide at least data_struct. See help UiO_TMS_correction')
end


% check if EEG structure is provided. If not, load raw data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        
        % check if the filepath is seperated by / or \ and and seperate file-path
        % from file-name
        if isempty(strfind(data_struct.vhdrsource,'\'))
            char_idx = strfind(data_struct.vhdrsource,'/'); 
        else
            char_idx = strfind(data_struct.vhdrsource,'\');
        end

        data_name = data_struct.vhdrsource(char_idx(end)+1:end);
        data_path = data_struct.vhdrsource(1:char_idx(end));

        % check if header ending is provided, otherwise add .vhdr to the file-name
        if strcmp(data_name(end-4:end),'.vhdr')
            EEG = pop_loadbv(data_path, data_name, [], []);
        else
            EEG = pop_loadbv(data_path, [data_name '.vhdr'], [], []);
        end
    
    else
        [EEG,locFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end

    
EEG.data = bsxfun(@minus,EEG.data,mean(EEG.data,2));

% check which kind of interpolation is required
if str2double(data_struct.tms_interpolation) == 1 % 1: automatic
    % epoch data
    EEGN = pop_epoch( EEG, {  'R128'  }, [-1 0.1], 'newname', ' resampled epochs', 'epochinfo', 'yes');

    % find indices in the time-domain
    [~,idx(1)] = min(abs(EEGN.times-(-10)));
    [~,idx(2)] = min(abs(EEGN.times-(-10)));
    [~,idx(3)] = min(abs(EEGN.times-(12)));
    [~,idx(4)] = min(abs(EEGN.times-(-100)));
    
    % average over electrodes and trials and take the absolute
    MElec = squeeze(mean(EEGN.data,1));
    MTrial = squeeze(mean(MElec(idx(4):idx(1),:),1));
    STrial = squeeze(std(MElec(idx(4):idx(1),:),1));
    MElec = abs(MElec);
    % define threshold (mean + std) and find indices of exeeding threshold
    % for every trial
    CIdx = zeros(2,size(MTrial,2));
    for i = 1:size(MTrial,2)
        TMSCut(i,:) = MElec(idx(2):idx(3),i) > MTrial(i)+(0.2*STrial(i));
        CIdx(1,i) = min(find(TMSCut(i,:))-1)+idx(2)-1;
        CIdx(2,i) = max(find(TMSCut(i,:))+1)+idx(2)-1;
    end
    
    % find event/trial indices of the TMS-stimuli
    getStimuli = strcmp({EEG.event.type},'boundary');
    getLatency = {EEG.event.latency};
    getLatency(getStimuli) = [];
    Tsample = 1/EEG.srate;
    
    % interpolate for every event/trial the date inbetween the threshold
    % indices using a baseline which is set -0.05sec to -0.01sec befor stimuli
    for ISi = 1:size(getLatency,2)
        accS = [getLatency{ISi}+(CIdx(1,ISi)-5000):getLatency{ISi}+(CIdx(2,ISi)-5000)];
        EEG.data(:,getLatency{ISi}+(CIdx(1,ISi)-5000):getLatency{ISi}+(CIdx(2,ISi)-5000)) = repmat(mean(EEG.data(:,getLatency{ISi}-(0.05/Tsample): ...
            getLatency{ISi}-(0.01/Tsample)),2),[1 length(accS)]);
    end
    
    % store mean tms cut for locFile
    tms_cut_start = mean(CIdx(1,:),2);
    tms_cut_end = mean(CIdx(2,:),2);
    disp(['TMS-artifact corrected between: ' num2str(tms_cut_start) ' and ' num2str(tms_cut_end) ' ms (mean over trials)']);
elseif str2double(data_struct.tms_interpolation) == 2 % 2: manually by visual inspection
    % epoch data
    EEGN = pop_epoch( EEG, {  'R128'  }, [-1 0.1], 'newname', ' resampled epochs', 'epochinfo', 'yes');
    
    % find indices in the time domain close to the stimulus
    [~,idx(1)] = min(abs(EEGN.times-(-5)));
    [~,idx(2)] = min(abs(EEGN.times-(15)));
    
    % compute the grand average mean over trials and electrodes
    GAMean = squeeze(mean(mean(EEGN.data,1),3));
    GAMean = GAMean-mean(GAMean);
    
    % open figure which stores the x-coordinates of two mouse clicks
    h = figure;
    plot(EEGN.times(idx(1):idx(2)),GAMean(idx(1):idx(2)));
    axis([-5 15 -50 50])
    title('Grand Average: please click where the TMS-artifact starts (1. click) and ends (2. click)')
    [x_I,~] = ginput(2);
    close(h)
    
    % use the x-coordinates of the mous clicks and find their indices
    [~,x_ipt(1)] = min(abs(EEGN.times-(x_I(1))));
    [~,x_ipt(2)] = min(abs(EEGN.times-(x_I(2))));
    
    % find event/trial indices of the TMS-stimuli
    getStimuli = strcmp({EEG.event.type},'boundary');
    getLatency = {EEG.event.latency};
    getLatency(getStimuli) = [];
    Tsample = 1/EEG.srate;
    
    % interpolate between indices (mous clickes) from -0.05sec to -0.01sec
    % before stimuli
    for ISi = 1:size(getLatency,2)
        accS = [getLatency{ISi}+(x_ipt(1)-5000):getLatency{ISi}+(x_ipt(2)-5000)];
        EEG.data(:,getLatency{ISi}+(x_ipt(1)-5000):getLatency{ISi}+(x_ipt(2)-5000)) = repmat(mean(EEG.data(:,getLatency{ISi}-(0.05/Tsample): ...
            getLatency{ISi}-(0.01/Tsample)),2),[1 length(accS)]);
    end
    
    % store mean tms cut for locFile
    tms_cut_start = x_I(1);
    tms_cut_end = x_I(2);
    disp(['TMS-artifact corrected between: ' num2str(tms_cut_start) ' and ' num2str(tms_cut_end) ' ms']);
elseif str2double(data_struct.tms_interpolation) == 3 % 3: manually by csv-file
    % transform csv time indices
    idx(1) = str2double(data_struct.interpolate_lower_tms)/EEG.srate;
    idx(2) = str2double(data_struct.interpolate_higher_tms)/EEG.srate;
    
    % find event/trial indices of the TMS-stimuli
    getStimuli = strcmp({EEG.event.type},'boundary');
    getLatency = {EEG.event.latency};
    getLatency(getStimuli) = [];
    Tsample = 1/EEG.srate;
    
    % interpolate between indices (csv-file) from -0.05sec to -0.01sec
    % before stimuli
    for ISi = 1:size(getLatency,2)
        accS = [getLatency{ISi}+(idx(1)/Tsample):getLatency{ISi}+(idx(2)/Tsample)];
        EEG.data(:,getLatency{ISi}+(idx(1)/Tsample):getLatency{ISi}+(idx(2)/Tsample)) = repmat(mean(EEG.data(:,getLatency{ISi}-(0.05/Tsample): ...
            getLatency{ISi}-(0.01/Tsample)),2),[1 length(accS)]);
    end
    
    % store mean tms cut for locFile
    tms_cut_start = data_struct.interpolate_lower_tms;
    tms_cut_end = data_struct.interpolate_higher_tms;
    disp(['TMS-artifact corrected between: ' num2str(tms_cut_start) ' and ' num2str(tms_cut_end) ' ms']);
else
    error('TMS interpolation correction not defined')
end
   
locFile{end+1} = {'after_tms',['EEG data is now TMS-artifact correctet using ' ...
    data_struct.tms_interpolation ' (1=automatic, 2=by visual inspection, 3=manually by csvfile).' ...
    ' Artifact is cutted between ' num2str(tms_cut_start) ' and ' num2str(tms_cut_start) ' ms.']};

if str2double(data_struct.plot_always)==1
    UiO_plots(data_struct,subj_name,EEG,locFile);
end
end