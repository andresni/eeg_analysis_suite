% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,locFile] = UiO_preprocessing(data_struct,~,EEG,locFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the raw data (without TMS-artifact interpolation)
% ~: subj_name not necessary for this function
% locFile: locFile of previous function. If empty [] this function will be
%       the first entry to the locFile
% 
% This function does the necessary preprocessing steps for EEG analysis.
% According to the csv-file, preprocessing is done using the specified parameters.
% this function does:
% 1: load subject data
% 3: high-pass filter
% 4: downsample
% 5: import channel info
% 6: remove bad channels + subspace reconstruction (clean_rawdata)
% 7: interpolate bad channels
% 8: rereference to the average
% 9: remove line noise (cleanline)
% 
%
% by questions: benjamin.thuerer@kit.edu
% 
function [EEG,locFile] = UiO_preprocessing(data_struct,~,EEG,locFile)

if nargin < 1
    error('provide at least data_struct. See help UiO_preprocessing')
end

% check if EEG structure is provided. If not, load raw data
if isempty(EEG)
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
end
  
% get all relevant parameters for preprocessing (see function at the buttom)
[Nsrate,HpassF,LpassF,BurstC,LNFreq] = UiO_getParameters(data_struct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change channel names for brainamp (this is specific to the recording
% system and the cap. If you use a different cap or recording system make
% shure you change this code!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RField1 = strcmp({EEG.chanlocs.labels},'AF8'); %change Af8 to AFz
RField2 = strcmp({EEG.chanlocs.labels},'AF7'); %change AF7 to FCz
EEG.chanlocs(RField1).labels = 'AFz';
EEG.chanlocs(RField2).labels = 'FCz';

% resample with mild anti-aliasing filter for possible connectivity analysis
if ~isempty(Nsrate)
    EEG = pop_resample(EEG, Nsrate, 0.8, 0.4);
end

% filter the data
if ~isempty(HpassF)
    EEG = pop_eegfiltnew(EEG, [], HpassF, [], true, [], 0);
end
if ~isempty(LpassF)
    EEG = pop_eegfiltnew(EEG, LpassF, [], [], true, [], 0);
end

% read channel locations according to Label in EEG.chanlocs
EEG = pop_chanedit(EEG, 'lookup','.\elec\standard_1005.elc');
chLocs = EEG.chanlocs;

% check if cleaning and channel rejection should be done automatically
% if yes, clean data using automatic subspace reconstruction (ASR)
if str2double(data_struct.cleaning_artifacts) == 1
    if str2double(data_struct.channel_rejection) == 2
        % clean data without removing any channel or time-bins
        EEG = clean_rawdata(EEG, 'off', 'off', 'off', 'off', BurstC, 'off');
        
        % clean channels manually by visual inspection
        % epoch the data and compute the average over trials
        EEGN = pop_epoch( EEG, {  'R128'  }, [-1 1], 'newname', ' resampled epochs', 'epochinfo', 'yes');
        [~,idx(1)] = min(abs(EEGN.times-(-1000)));
        [~,idx(2)] = min(abs(EEGN.times-(1000)));
        MElec = squeeze(mean(EEGN.data,3));
        goodChan = [];
        badChan = [];
        
        % plot the average for each channel and distinguish between good
        % and bad channels according to the mouse / space click
        for Ei = 1:size(MElec,1)
            h = figure;
            plot(EEGN.times(idx(1):idx(2)),MElec(Ei,idx(1):idx(2)));
            title(['Channel ' int2str(Ei) ': press space to remove channel']);
            grid
            axis([-1000 1000 -15 15]);
            ylabel('Amplitude (\muV)'); xlabel('Time (ms)');
            button = waitforbuttonpress;
            if button == 0
                goodChan(end+1) = Ei;
            else
                badChan(end+1) = Ei;
            end
            close(h)
        end
        
        % plot the average for each channel again and show the actual
        % marking of the channel. Press space if you changed your mind
        for Ei = 1:size(MElec,1)
            if find(goodChan==Ei)
                Stigma = 'good';
            else
                Stigma = 'bad';
            end
            h = figure;
            plot(EEGN.times(idx(1):idx(2)),MElec(Ei,idx(1):idx(2)));
            title({['Channel ' int2str(Ei) ' marked as ' Stigma] ; ...
                ['if you want to change the mark (regardless of direction g-->b; b-->g) press space']});
            grid
            axis([-1000 1000 -15 15]);
            ylabel('Amplitude (\muV)'); xlabel('Time (ms)');
            button = waitforbuttonpress;
            if button ~= 0
                if find(badChan==Ei)
                    idx_Ei = badChan == Ei;
                    badChan(idx_Ei) = [];
                else
                    badChan(end+1) = Ei;
                end
            end
            close(h)
        end
        
        % sort bad channels and remove bad channels from the dataset
        badChan = sort(badChan);
        EEG.nbchan = EEG.nbchan-length(badChan);
        EEG.data(badChan,:) = [];
        EEG.chanlocs(:,badChan) = [];
    elseif str2double(data_struct.channel_rejection) == 1
        % clean the data but do not filter or remove time-bins
        EEG = clean_rawdata(EEG, [], 'off', [], [], BurstC, 'off');       
    else
        warning('no correct channel rejection type provided. default cleaning with channel rejection on')
    end
elseif str2double(data_struct.cleaning_artifacts) == 0
    if str2double(data_struct.channel_rejection) == 1
        % clean channel manually by visual inspection and do not clean data
        % epoch the data and compute the average over trials
        EEGN = pop_epoch( EEG, {  'R128'  }, [-1 1], 'newname', ' resampled epochs', 'epochinfo', 'yes');
        [~,idx(1)] = min(abs(EEGN.times-(-1000)));
        [~,idx(2)] = min(abs(EEGN.times-(1000)));
        MElec = squeeze(mean(EEGN.data,3));
        goodChan = [];
        badChan = [];
        
        % plot the average for each channel and distinguish between good
        % and bad channels according to the mouse / space click
        for Ei = 1:size(MElec,1)
            h = figure;
            plot(EEGN.times(idx(1):idx(2)),MElec(Ei,idx(1):idx(2)));
            title(['Channel ' int2str(Ei) ': click space to remove channel']);
            grid
            axis([-1000 1000 -15 15]);
            ylabel('Amplitude (\muV)'); xlabel('Time (ms)');
            button = waitforbuttonpress;
            if button == 0
                goodChan(end+1) = Ei;
            else
                badChan(end+1) = Ei;
            end
            close(h)
        end
        
        % plot the average for each channel again and show the actual
        % marking of the channel. Press space if you changed your mind
        for Ei = 1:size(MElec,1)
            if find(goodChan==Ei)
                Stigma = 'good';
            else
                Stigma = 'bad';
            end
            h = figure;
            plot(EEGN.times(idx(1):idx(2)),MElec(Ei,idx(1):idx(2)));
            title({['Channel ' int2str(Ei) ' marked as ' Stigma] ; ...
                ['if you want to change the mark (regardless of direction g-->b; b-->g) click space']});
            grid
            axis([-1000 1000 -15 15]);
            ylabel('Amplitude (\muV)'); xlabel('Time (ms)');
            button = waitforbuttonpress;
            if button ~= 0
                if find(badChan==Ei)
                    idx_Ei = badChan == Ei;
                    badChan(idx_Ei) = [];
                else
                    badChan(end+1) = Ei;
                end
            end
            close(h)
        end
        
        % sort bad channels and remove bad channels from the dataset
        badChan = sort(badChan);
        EEG.nbchan = EEG.nbchan-length(badChan);
        EEG.data(badChan,:) = [];
        EEG.chanlocs(:,badChan) = [];
    else
        warning('no cleaning and nor channel rejection done. This might mess up a later re-referencing')
    end
elseif str2double(data_struct.cleaning_artifacts) > 1
    warning('no correct cleaning artifacts type provided. No cleaning is done!')
end

%store removed channels
EEG.CHremoved = [];
EEG.CHremoved = setdiff({chLocs.labels},{EEG.chanlocs.labels}); 

% interpolate bad channels
EEG = pop_interp(EEG, chLocs,'spherical'); %interpolate removed channels

% re-reference the data either to a channel or average
if isempty(str2num(data_struct.re_referencing))
    ChNum = strcmp({EEG.chanlocs.labels},data_struct.re_referencing);
    EEG = pop_reref( EEG, ChNum);
    EEG.ref = data_struct.re_referencing;
elseif str2double(data_struct.re_referencing) == 1
    EEG.data = bsxfun(@minus, EEG.data, sum(EEG.data,1)/EEG.nbchan);
    EEG.ref = 'average';
end

% correct data for line noise without notch-filter
EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:EEG.trials] ,'computepower',0,'linefreqs',[LNFreq LNFreq*2],'normSpectrum',0,'p',0.01,'pad',2,'plotfigures' ...
    ,0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);
   
locFile{end+1} = {'preprocessed',['preprocessed with (sampling rate; ' ...
    'highpassfilter; lowpassfilter; burst criterion; line noise correction: ' ...
    num2str(Nsrate) '; ' num2str(HpassF) '; ' num2str(LpassF) '; ' ...
    num2str(BurstC) '; ' num2str(LNFreq)]};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes the provided parameters (of the csv file) out of the
% struct (data_struct) and makes them availeble for the preprocessing

function [Nsrate,HpassF,LpassF,BurstC,LNFreq] = UiO_getParameters(data_struct)
% this function defines the relevant parameteres according to the csvfile

% sampling rate
if str2double(data_struct.downsample_rate) > 1
    Nsrate = str2double(data_struct.downsample_rate); %which sampling rate for resample?
elseif str2double(data_struct.downsample_rate) == 0
    Nsrate = [];
else
    warning('sampling rate to low. New sampling rate: 1000 hz')
    Nsrate = 1000;
end

% highpassfilter bandpass edge
if str2double(data_struct.highpassfilter) == 0
    warning('are you sure you do not want to high-pass filter the data?')
    HpassF = [];
else
    HpassF = str2double(data_struct.highpassfilter); %threshold for high-pass filter (this means cut-off at HpassF - 0.5 hz!)
end

% lowpassfilter bandpass edge
if str2double(data_struct.lowpassfilter) == 0  
    LpassF = [];
else
    LpassF = str2double(data_struct.lowpassfilter);
    if LpassF < 20
        warning('are you sure you want to low pass below 20 hz?')
    end
end

% Burst-criterion: defines how strict cleaning will be done
if str2double(data_struct.correction_sensitivity) == 0  
    BurstC = [];
    warning('No Burst criterion (cleaning sensitivity) provided. Data will not be cleaned!')
else
    BurstC = str2double(data_struct.correction_sensitivity);
end

% line noise frequency
if str2double(data_struct.linenoise_frequency) == 0
    warning('line noise frequency not defined. changing to 50 Hz')
    LNFreq = 50;
else
    LNFreq = str2double(data_struct.linenoise_frequency);
end

end