% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,logFile] = UiO_preprocessing(data_struct,~,EEG,logFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the raw data (without TMS-artifact interpolation)
% ~: subj_name not necessary for this function
% logFile: logFile of previous function. If empty [] this function will be
%       the first entry to the logFile
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
% by questions:
% benjamin.thuerer@kit.edu
% 
function [EEG,logFile] = KIT_preprocessing(data_struct,subj_name,EEG,logFile)

if nargin < 1
    error('provide at least data_struct. See help UiO_preprocessing')
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
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end
  
% get all relevant parameters for preprocessing (see function at the buttom)
[Nsrate,HpassF,LpassF,BurstC,LNFreq] = UiO_getParameters(data_struct);



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
EEG = pop_chanedit(EEG, 'lookup','C:\\Users\\Promotion\\Matlab_Scripts\\eeglab13_5_4b\\plugins\\dipfit2.3\\standard_BEM\\elec\\standard_1005.elc', ...
    'load',{'C:\Users\Promotion\03_Studie-III_DAVOS\data\Raw_Data\06_Task-EEG-Daten\Chan31_RefCz.locs' 'filetype' 'autodetect'},'eval', ...
    'chans = pop_chancenter( chans, [],[]);','append',31,'changefield',{32 'labels' 'Cz'},'changefield',{32 'theta' '90'},'changefield',{32 'radius' '0'});

chLocs = EEG.chanlocs;

% remove mean over channel
EEG.data = bsxfun(@minus,EEG.data,mean(EEG.data,2));

EEG.mainEventTrigger = 'S  7';
multPl = 50; % y-axis distance for plot
CHlabels = {EEG.chanlocs.labels};

% check if cleaning and channel rejection should be done automatically
% if yes, clean data using automatic subspace reconstruction (ASR)
if str2double(data_struct.cleaning_artifacts) == 1
    if str2double(data_struct.channel_rejection) == 2
        % clean data without removing any channel or time-bins
        EEGT = EEG;
        EEG = clean_rawdata(EEG, 'off', [0.25 0.75], 'off', 'off', BurstC, 'off');
        
        % run clean_rawdata only to find bad channels
        disp('cleaning channels again to compare bad channels');
        EEGT = clean_rawdata(EEGT, [], [0.25 0.75], 0.85, [], 'off', 'off');
        bad_chans_clean = setdiff({chLocs.labels},{EEGT.chanlocs.labels});

        % clean channels manually by visual inspection
        % epoch the data and compute the average over trials
        EEGN = pop_epoch( EEG, {  EEG.mainEventTrigger  }, [-1 1], 'newname', ' resampled epochs', 'epochinfo', 'yes');
        [~,idx(1)] = min(abs(EEGN.times-(-500)));
        [~,idx(2)] = min(abs(EEGN.times-(500)));
        goodChan = [];
        badChan = [];
        
        % plot for each channel all trials and distinguish between good
        % and bad channels according to the mouse / space click
        for Ei = 1:size(EEGN.data,1)
            % check if bad channel in clean_rawdata
            if find(strcmp(bad_chans_clean,CHlabels{Ei}))
                bad_chan_mark = 'bad channel';
            else
                bad_chan_mark = ' ';
            end
            
            CutTrial = squeeze(EEGN.data(Ei,:,:));
            CutTrial = bsxfun(@minus,CutTrial,mean(CutTrial,1));
            
            max_value = 1.1*(mean(max(CutTrial)) + mean(abs(min(CutTrial))));
            
            n = 1:max_value:size(CutTrial,2)*max_value;
            CutTrial = bsxfun(@plus,CutTrial,n);
            trial_cut = size(CutTrial,2);

            % load figure with bigger size. Then do either four or three subplots each
            % with 1/4 (or 3) of the trials according to the total amount of trials.
            if size(EEGN.data,3) > 96
                devTrials = 4;
                h = figure('units','normalized','position',[.05 .05 .9 .9]);
                subplot(1,4,1); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),1:floor(trial_cut/devTrials)),'color',[0,0,0]);
                grid
                title(['trials ' int2str(1) ' to ' int2str(floor(trial_cut/devTrials))]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) -multPl n(floor(trial_cut/devTrials))+multPl]);
                xlabel('Time (ms)')

                subplot(1,4,2); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),floor(trial_cut/devTrials)+1:(floor(trial_cut/devTrials))*2),'color',[0,0,0]);
                grid
                title(['trials ' int2str(floor(trial_cut/devTrials)+1) ' to ' int2str((floor(trial_cut/devTrials))*2)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n(floor(trial_cut/devTrials)) n((floor(trial_cut/devTrials))*2)+multPl]);
                xlabel('Time (ms)')

                subplot(1,4,3); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),(floor(trial_cut/devTrials))*2+1:(floor(trial_cut/devTrials))*3),'color',[0,0,0]);
                grid
                title(['trials ' int2str((floor(trial_cut/devTrials))*2+1) ' to ' int2str((floor(trial_cut/devTrials))*3)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n((floor(trial_cut/devTrials))*2) n((floor(trial_cut/devTrials))*3)+multPl]);
                xlabel('Time (ms)')
                
                subplot(1,4,4); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),(floor(trial_cut/devTrials))*3+1:end),'color',[0,0,0]);
                grid
                title(['trials ' int2str((floor(trial_cut/devTrials))*3+1) ' to ' int2str(trial_cut)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n((floor(trial_cut/devTrials))*3) n(end)+multPl]);
                xlabel('Time (ms)')

                suptitle({['channel ' CHlabels{Ei} ', \color{red}' bad_chan_mark] ; ...
                    ['\color{black} if you want to reject: press space; if you want to keep: mouse click']})
            else
                devTrials = 3;
                h = figure('units','normalized','position',[.05 .05 .9 .9]);
                subplot(1,3,1); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),1:floor(trial_cut/devTrials)),'color',[0,0,0]);
                grid
                title(['trials ' int2str(1) ' to ' int2str(floor(trial_cut/devTrials))]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) -multPl n(floor(trial_cut/devTrials))+multPl]);
                xlabel('Time (ms)')

                subplot(1,3,2); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),floor(trial_cut/devTrials)+1:(floor(trial_cut/devTrials))*2),'color',[0,0,0]);
                grid
                title(['trials ' int2str(floor(trial_cut/devTrials)+1) ' to ' int2str((floor(trial_cut/devTrials))*2)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n(floor(trial_cut/devTrials)) n((floor(trial_cut/devTrials))*2)+multPl]);
                xlabel('Time (ms)')

                subplot(1,3,3); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),(floor(trial_cut/devTrials))*2+1:end),'color',[0,0,0]);
                grid
                title(['trials ' int2str((floor(trial_cut/devTrials))*2+1) ' to ' int2str(trial_cut)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n((floor(trial_cut/devTrials))*2) n(end)+multPl]);
                xlabel('Time (ms)')

                suptitle({['channel ' CHlabels{Ei} ', \color{red}' bad_chan_mark] ; ...
                    ['\color{black} if you want to reject: press space; if you want to keep: mouse click']})
            end
            
            button = waitforbuttonpress;
            if button == 0
                goodChan(end+1) = Ei;
            else
                badChan(end+1) = Ei;
            end
            close(h)
        end
        
        % plot trials for each channel again and show the actual
        % marking of the channel. Press space if you changed your mind
        val = [];
        for Ei = 1:size(EEGN.data,1)
            if find(goodChan==Ei)
                Stigma = 'good';
            else
                Stigma = 'bad';
            end
            CutTrial = squeeze(EEGN.data(Ei,:,:));
            CutTrial = bsxfun(@minus,CutTrial,mean(CutTrial,1));
        
            max_value = 1.1*(mean(max(CutTrial)) + mean(abs(min(CutTrial))));
            
            n = 1:max_value:size(CutTrial,2)*max_value;
            CutTrial = bsxfun(@plus,CutTrial,n);
            trial_cut = size(CutTrial,2);

            % load figure with bigger size. Then do four or three subplots each
            % with 1/3 of the trials.
            if size(EEGN.data,3) > 96
                devTrials = 4;
                h = figure('units','normalized','position',[.05 .05 .9 .9]);
                subplot(1,4,1); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),1:floor(trial_cut/devTrials)),'color',[0,0,0]);
                grid
                title(['trials ' int2str(1) ' to ' int2str(floor(trial_cut/devTrials))]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) -multPl n(floor(trial_cut/devTrials))+multPl]);
                xlabel('Time (ms)')

                subplot(1,4,2); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),floor(trial_cut/devTrials)+1:(floor(trial_cut/devTrials))*2),'color',[0,0,0]);
                grid
                title(['trials ' int2str(floor(trial_cut/devTrials)+1) ' to ' int2str((floor(trial_cut/devTrials))*2)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n(floor(trial_cut/devTrials)) n((floor(trial_cut/devTrials))*2)+multPl]);
                xlabel('Time (ms)')

                subplot(1,4,3); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),(floor(trial_cut/devTrials))*2+1:(floor(trial_cut/devTrials))*3),'color',[0,0,0]);
                grid
                title(['trials ' int2str((floor(trial_cut/devTrials))*2+1) ' to ' int2str((floor(trial_cut/devTrials))*3)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n((floor(trial_cut/devTrials))*2) n((floor(trial_cut/devTrials))*3)+multPl]);
                xlabel('Time (ms)')
                
                subplot(1,4,4); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),(floor(trial_cut/devTrials))*3+1:end),'color',[0,0,0]);
                grid
                title(['trials ' int2str((floor(trial_cut/devTrials))*3+1) ' to ' int2str(trial_cut)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n((floor(trial_cut/devTrials))*3) n(end)+multPl]);
                xlabel('Time (ms)')

                suptitle({['channel ' CHlabels{Ei} ' is marked as ' Stigma] ; ...
                    ['if you want to change the mark (regardless of direction g-->b; b-->g) press space'] ; ...
                    ['if you want to skip this additional loop press: e']})
            else
                devTrials = 3;
                h = figure('units','normalized','position',[.05 .05 .9 .9]);
                subplot(1,3,1); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),1:floor(trial_cut/devTrials)),'color',[0,0,0]);
                grid
                title(['trials ' int2str(1) ' to ' int2str(floor(trial_cut/devTrials))]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) -multPl n(floor(trial_cut/devTrials))+multPl]);
                xlabel('Time (ms)')

                subplot(1,3,2); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),floor(trial_cut/devTrials)+1:(floor(trial_cut/devTrials))*2),'color',[0,0,0]);
                grid
                title(['trials ' int2str(floor(trial_cut/devTrials)+1) ' to ' int2str((floor(trial_cut/devTrials))*2)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n(floor(trial_cut/devTrials)) n((floor(trial_cut/devTrials))*2)+multPl]);
                xlabel('Time (ms)')

                subplot(1,3,3); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),(floor(trial_cut/devTrials))*2+1:end),'color',[0,0,0]);
                grid
                title(['trials ' int2str((floor(trial_cut/devTrials))*2+1) ' to ' int2str(trial_cut)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n((floor(trial_cut/devTrials))*2) n(end)+multPl]);
                xlabel('Time (ms)')

                suptitle({['channel ' CHlabels{Ei} ' is marked as ' Stigma] ; ...
                    ['if you want to change the mark (regardless of direction g-->b; b-->g) press space'] ; ...
                    ['if you want to skip this additional loop press: e']})
            end      
            
            button = waitforbuttonpress;
            % breaks loop if "e" is pressed
            val=double(get(h,'CurrentCharacter'));
            if val == 101
                close(h)
                break 
            end
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
        EEG = clean_rawdata(EEG, [], [0.25 0.75], 0.85, [], BurstC, 'off');       
    else
        warning('no correct channel rejection type provided. default cleaning with channel rejection on')
    end
    
% if no cleaning should be done but channel rejection manually:
elseif str2double(data_struct.cleaning_artifacts) == 0
    if str2double(data_struct.channel_rejection) == 1
        % clean channel manually by visual inspection and do not clean data
        % epoch the data and compute the average over trials
        EEGN = pop_epoch( EEG, {  EEG.mainEventTrigger  }, [-1 1], 'newname', ' resampled epochs', 'epochinfo', 'yes');
        [~,idx(1)] = min(abs(EEGN.times-(-500)));
        [~,idx(2)] = min(abs(EEGN.times-(500)));
        goodChan = [];
        badChan = [];
        
        % run clean_rawdata only to find bad channels
        disp('cleaning channels again to compare bad channels');
        EEGT = EEG;
        EEGT = clean_rawdata(EEGT, [], [0.25 0.75], 0.85, [], 'off', 'off');
        bad_chans_clean = setdiff({chLocs.labels},{EEGT.chanlocs.labels});
        
        % plot for each channel plot all trials and distinguish between good
        % and bad channels according to the mouse / space click
        for Ei = 1:size(EEGN.data,1)
            if find(strcmp(bad_chans_clean,CHlabels{Ei}))
                bad_chan_mark = 'bad channel';
            else
                bad_chan_mark = ' ';
            end
            
            CutTrial = squeeze(EEGN.data(Ei,:,:));
            CutTrial = bsxfun(@minus,CutTrial,mean(CutTrial,1));
        
            max_value = 1.1*(mean(max(CutTrial)) + mean(abs(min(CutTrial))));
            
            n = 1:max_value:size(CutTrial,2)*max_value;
            CutTrial = bsxfun(@plus,CutTrial,n);
            trial_cut = size(CutTrial,2);

             % load figure with bigger size. Then do either four or three subplots each
            % with 1/4 (or 3) of the trials according to the total amount of trials.
            if size(EEGN.data,3) > 96
                devTrials = 4;
                h = figure('units','normalized','position',[.1 .1 .9 .9]);
                subplot(1,4,1); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),1:floor(trial_cut/devTrials)),'color',[0,0,0]);
                grid
                title(['trials ' int2str(1) ' to ' int2str(floor(trial_cut/devTrials))]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) -multPl n(floor(trial_cut/devTrials))+multPl]);
                xlabel('Time (ms)')

                subplot(1,4,2); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),floor(trial_cut/devTrials)+1:(floor(trial_cut/devTrials))*2),'color',[0,0,0]);
                grid
                title(['trials ' int2str(floor(trial_cut/devTrials)+1) ' to ' int2str((floor(trial_cut/devTrials))*2)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n(floor(trial_cut/devTrials)) n((floor(trial_cut/devTrials))*2)+multPl]);
                xlabel('Time (ms)')

                subplot(1,4,3); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),(floor(trial_cut/devTrials))*2+1:(floor(trial_cut/devTrials))*3),'color',[0,0,0]);
                grid
                title(['trials ' int2str((floor(trial_cut/devTrials))*2+1) ' to ' int2str((floor(trial_cut/devTrials))*3)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n((floor(trial_cut/devTrials))*2) n((floor(trial_cut/devTrials))*3)+multPl]);
                xlabel('Time (ms)')
                
                subplot(1,4,4); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),(floor(trial_cut/devTrials))*3+1:end),'color',[0,0,0]);
                grid
                title(['trials ' int2str((floor(trial_cut/devTrials))*3+1) ' to ' int2str(trial_cut)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n((floor(trial_cut/devTrials))*3) n(end)+multPl]);
                xlabel('Time (ms)')

                suptitle({['channel ' CHlabels{Ei} ', \color{red}' bad_chan_mark] ; ...
                    ['\color{black} if you want to reject: press space; if you want to keep: mouse click']})
            else
                devTrials = 3;
                h = figure('units','normalized','position',[.1 .1 .9 .9]);
                subplot(1,3,1); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),1:floor(trial_cut/devTrials)),'color',[0,0,0]);
                grid
                title(['trials ' int2str(1) ' to ' int2str(floor(trial_cut/devTrials))]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) -multPl n(floor(trial_cut/devTrials))+multPl]);
                xlabel('Time (ms)')

                subplot(1,3,2); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),floor(trial_cut/devTrials)+1:(floor(trial_cut/devTrials))*2),'color',[0,0,0]);
                grid
                title(['trials ' int2str(floor(trial_cut/devTrials)+1) ' to ' int2str((floor(trial_cut/devTrials))*2)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n(floor(trial_cut/devTrials)) n((floor(trial_cut/devTrials))*2)+multPl]);
                xlabel('Time (ms)')

                subplot(1,3,3); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),(floor(trial_cut/devTrials))*2+1:end),'color',[0,0,0]);
                grid
                title(['trials ' int2str((floor(trial_cut/devTrials))*2+1) ' to ' int2str(trial_cut)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n((floor(trial_cut/devTrials))*2) n(end)+multPl]);
                xlabel('Time (ms)')

                suptitle({['channel ' CHlabels{Ei} ', \color{red}' bad_chan_mark] ; ...
                    ['\color{black} if you want to reject: press space; if you want to keep: mouse click']})
            end
            
            button = waitforbuttonpress;
            if button == 0
                goodChan(end+1) = Ei;
            else
                badChan(end+1) = Ei;
            end
            close(h)
        end
        
        % plot trials for each channel again and show the actual
        % marking of the channel. Press space if you changed your mind
        val = [];
        for Ei = 1:size(EEGN.data,1)
            if find(goodChan==Ei)
                Stigma = 'good';
            else
                Stigma = 'bad';
            end
            CutTrial = squeeze(EEGN.data(Ei,:,:));
            CutTrial = bsxfun(@minus,CutTrial,mean(CutTrial,1));
        
            max_value = 1.1*(mean(max(CutTrial)) + mean(abs(min(CutTrial))));
            
            n = 1:max_value:size(CutTrial,2)*max_value;
            CutTrial = bsxfun(@plus,CutTrial,n);
            trial_cut = size(CutTrial,2);

            % load figure with bigger size. Then do four or three subplots each
            % with 1/3 of the trials.
            if size(EEGN.data,3) > 96
                devTrials = 4;
                h = figure('units','normalized','position',[.1 .1 .9 .9]);
                subplot(1,4,1); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),1:floor(trial_cut/devTrials)),'color',[0,0,0]);
                grid
                title(['trials ' int2str(1) ' to ' int2str(floor(trial_cut/devTrials))]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) -multPl n(floor(trial_cut/devTrials))+multPl]);
                xlabel('Time (ms)')

                subplot(1,4,2); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),floor(trial_cut/devTrials)+1:(floor(trial_cut/devTrials))*2),'color',[0,0,0]);
                grid
                title(['trials ' int2str(floor(trial_cut/devTrials)+1) ' to ' int2str((floor(trial_cut/devTrials))*2)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n(floor(trial_cut/devTrials)) n((floor(trial_cut/devTrials))*2)+multPl]);
                xlabel('Time (ms)')

                subplot(1,4,3); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),(floor(trial_cut/devTrials))*2+1:(floor(trial_cut/devTrials))*3),'color',[0,0,0]);
                grid
                title(['trials ' int2str((floor(trial_cut/devTrials))*2+1) ' to ' int2str((floor(trial_cut/devTrials))*3)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n((floor(trial_cut/devTrials))*2) n((floor(trial_cut/devTrials))*3)+multPl]);
                xlabel('Time (ms)')
                
                subplot(1,4,4); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),(floor(trial_cut/devTrials))*3+1:end),'color',[0,0,0]);
                grid
                title(['trials ' int2str((floor(trial_cut/devTrials))*3+1) ' to ' int2str(trial_cut)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n((floor(trial_cut/devTrials))*3) n(end)+multPl]);
                xlabel('Time (ms)')

                suptitle({['channel ' CHlabels{Ei} ' is marked as ' Stigma] ; ...
                    ['if you want to change the mark (regardless of direction g-->b; b-->g) press space'] ; ...
                    ['if you want to skip this additional loop press: e']})
            else
                devTrials = 3;
                h = figure('units','normalized','position',[.1 .1 .9 .9]);
                subplot(1,3,1); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),1:floor(trial_cut/devTrials)),'color',[0,0,0]);
                grid
                title(['trials ' int2str(1) ' to ' int2str(floor(trial_cut/devTrials))]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) -multPl n(floor(trial_cut/devTrials))+multPl]);
                xlabel('Time (ms)')

                subplot(1,3,2); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),floor(trial_cut/devTrials)+1:(floor(trial_cut/devTrials))*2),'color',[0,0,0]);
                grid
                title(['trials ' int2str(floor(trial_cut/devTrials)+1) ' to ' int2str((floor(trial_cut/devTrials))*2)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n(floor(trial_cut/devTrials)) n((floor(trial_cut/devTrials))*2)+multPl]);
                xlabel('Time (ms)')

                subplot(1,3,3); plot(EEGN.times(idx(1):idx(2)),CutTrial(idx(1):idx(2),(floor(trial_cut/devTrials))*2+1:end),'color',[0,0,0]);
                grid
                title(['trials ' int2str((floor(trial_cut/devTrials))*2+1) ' to ' int2str(trial_cut)]);
                set(gca,'YTick',[]);
                axis([min(EEGN.times(idx(1):idx(2))) max(EEGN.times(idx(1):idx(2))) n((floor(trial_cut/devTrials))*2) n(end)+multPl]);
                xlabel('Time (ms)')

                suptitle({['channel ' CHlabels{Ei} ' is marked as ' Stigma] ; ...
                    ['if you want to change the mark (regardless of direction g-->b; b-->g) press space'] ; ...
                    ['if you want to skip this additional loop press: e']})
            end 
            
            % breaks loop if "e" is pressed
            val=double(get(h,'CurrentCharacter'));
            if val == 101
                break
            end
            
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

% save original channel locations for interpolation (see UiO_ica_cleaning)
EEG.chLocs = chLocs;

% re-reference the data either to a channel or average
if isempty(str2num(data_struct.re_referencing))
    ChNum = strcmp({EEG.chanlocs.labels},data_struct.re_referencing);
    EEG = pop_reref( EEG, ChNum);
    EEG.ref = data_struct.re_referencing;
elseif str2double(data_struct.re_referencing) == 1
    EEG = pop_reref( EEG, [],'refloc',struct('theta',{90},'radius',{0},'labels',{'Cz'},'sph_theta',{[-90]}, ...
        'sph_phi',{[90]},'X',{[3.75e-33]},'Y',{[-6.12e-17]},'Z',{[1]},'sph_radius',{[1]},'type',{''},'ref',{''},'urchan',{[]},'datachan',{0}));
    EEG.ref = 'average';
end

% remove mean over channel
EEG.data = bsxfun(@minus,EEG.data,mean(EEG.data,2));

% reduce line noise either by notch-filter or cleanline
if str2double(data_struct.notch_filter) == 0
    EEG = pop_cleanline(EEG,'ChanCompIndices',[1:EEG.nbchan],'SignalType','Channels','computepower',0,'LineFrequencies',[LNFreq LNFreq*2 LNFreq*4],'normSpectrum',0,'p',0.01,'pad',2,'plotfigures' ...
    ,0,'scanforlines',1,'tau',100,'verb',1,'winsize',4,'winstep',4);
elseif str2double(data_struct.notch_filter) == 1
    filter_deg = 3;
    LNFnotch = [(LNFreq-5)*2/EEG.srate, (LNFreq+5)*2/EEG.srate];
    [b,a] = butter(filter_deg,LNFnotch,'stop');
    EEG.data = filter(b,a,EEG.data);
    disp(['line noise reduced by notch filter']);
else
    error('no notch_filter stated in the csv file')
end
 
% loc file entry
logFile{end+1} = {'preprocessed',['preprocessed with (sampling rate; ' ...
    'highpassfilter; lowpassfilter; burst criterion; line noise correction: ' ...
    num2str(Nsrate) '; ' num2str(HpassF) '; ' num2str(LpassF) '; ' ...
    num2str(BurstC) '; ' num2str(LNFreq)]};

if str2double(data_struct.plot_always)==1
    UiO_plots(data_struct,subj_name,EEG,logFile);
end

disp('data preprocessing is done')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will crop the data according to the first and last trial or
% to the provided time in the csvfile. The idea is that at some recordings
% at the beginning and end some things are recorded (talking) which should
% not be considered for preprocessing.
% 
% 
function EEG = UiO_crop_data(data_struct,EEG)

if isempty(strfind(data_struct.crop_data,','))
    % find first and last stimuli and define crop borders as stimuli 1 -
    % 10seconds and stimuli end + 10 seconds
    getStimuli = strcmp({EEG.event.type},'boundary');
    getLatency = {EEG.event.latency};
    getLatency(getStimuli) = [];
    x = getLatency{1}-10*EEG.srate;
    y = getLatency{end}+10*EEG.srate;
elseif str2souble(data_struct.crop_data) == 0
    x = 1;
    y = size(EEG.data,2);
else
    %define stimuli 1 and stimuli end according to csvfile
    n_comma = strfind(data_struct.crop_data,',');
    x_str = data_struct.crop_data(1:n_comma-1);
    y_str = data_struct.crop_data(n_comma+1:end);
    [~,x] = min(abs(EEG.times-str2double(data_struct.crop_data(1)))); 
    [~,y] = min(abs(EEG.times-str2double(data_struct.crop_data(2)))); 
end

EEG.data = EEG.data(:,x:y);
EEG.times = EEG.times(x:y);

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