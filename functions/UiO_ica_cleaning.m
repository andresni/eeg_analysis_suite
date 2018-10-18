% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,logFile] = UiO_ica_cleaning(data_struct,subj_name,EEG,logFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the 'after_ica' data (if availeble)
% subj_name: subject name according to csvfile
% logFile: logFile of previous function. If empty [] this function will
%       load the 'after_ica' logFile (if availeble)
%
% This function will clean independent components automized (1) or manually
% by visual inspection (2). Manually by inspection is highly recommendet
% for ICs as the automized method works not very good.
% ICA cleaning works on both, continous and epoched data.
% 
% by questions:
% b.e.juel@medisin.uio.no
% sevenius.nilsen@gmail.com
% benjamin.thuerer@kit.edu
%
function [EEG,logFile] = UiO_ica_cleaning(data_struct,subj_name,EEG,logFile)

if nargin < 2
    error('provide at least data_struct and subject name. See help UiO_ica_cleaning')
end


% check if EEG structure is provided. If not, load previous data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,'after_ica');   
    else
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end

% if eeglab_options not saving properly, icaact might be deleted. This must
% be reconstructed now
if isempty(EEG.icaact)
    if ndims(EEG.data)==3
        for i=1:EEG.trials
            EEG.icaact(:,:,i)=(EEG.icaweights*EEG.icasphere)*squeeze(EEG.data(EEG.icachansind,:,i));
        end
    else
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    end
end

% time period (ms) for epochs in ICA (just for plotting) and adjust output
% file
timePeriod = [-500 500];
out = 'adjust_output_can_be_deleted';


%% ADJUST needs epoched data. Fake epochs if necessary
if size(EEG.data,3) == 1
    lag=5; %epoch duration (in seconds)
    disp(['Continuous dataset epoched in ' num2str(lag) ' sec long epochs to compute feature statistics over time']);
    % check whether '5sec' events are present
    si=0; 
    for i=1:length(EEG.event) 
        if(strcmp(EEG.event(1,i).type, '5sec')==1) 
            si=1; 
        end
    end
    if si==0 %add events
        ntrials=floor((EEG.xmax-EEG.xmin)/lag);
        nevents=length(EEG.event);
        for index=1:ntrials
            EEG.event(index+nevents).type='5sec';
            EEG.event(index+nevents).latency=1+(index-1)*lag*EEG.srate; %EEG.srate is the sampling frequency
        end

        EEG=eeg_checkset(EEG,'eventconsistency');
    end

    EEGep = pop_epoch( EEG, {  '5sec'  }, [0 lag], 'newname', [EEG.setname '_ep5'] , 'epochinfo', 'yes');
    %         % removing baseline
    %         EEGep = pop_rmbase( EEGep, [0  0]);
    EEGep = eeg_checkset(EEGep);

    % collects ICA data from EEG
    if isempty(EEGep.icaact)
    warning('EEG.icaact missing! Recomputed from EEG.icaweights, EEG.icasphere and EEG.data');
    % Next instruction: see eeg_checkset
        for i=1:EEGep.trials
            EEGep.icaact(:,:,i) = (EEGep.icaweights*EEGep.icasphere)*squeeze(EEGep.data(EEGep.icachansind,:,i));
        end
    end
else
    EEGep = EEG;
end


%%  clean ICs using ADJUST

    
% check if IC-rejection is automatic (1) or manually (2)
if str2double(data_struct.ica_rejection) == 1
    % run ADJUST und mark bad ICs according to the outcome of ADJUST
    [art] = ADJUST(EEGep,out);
    badICs = art;
    goodICs = setdiff(1:size(EEGep.icaact,1),badICs);
elseif str2double(data_struct.ica_rejection) == 2   
    % run ADJUST only to mark ICs in the plotting
    [~, horiz, vert, blink, disc] = ADJUST(EEGep,out);

    % set good and bad ic variables and find time indices in the data
    [~,idx_t(1)] = min(abs(EEGep.times-timePeriod(1)));
    [~,idx_t(2)] = min(abs(EEGep.times-timePeriod(2)));
    goodICs = [];
    badICs = [];
    
    % for every component: plot the component with ERP, power-spectra,
    % power over trials, and topographic plot    
    for Ci = 1:size(EEGep.icaact,1)
        artHE = [];
        artVE = [];
        artB = [];
        artD = [];
        
        % test if ADJUST would exclude this IC and why
        if find(horiz==Ci)
            artHE = 'HEye'; % horizontal eye movement
        end
        if find(vert==Ci)
            artVE = 'VEye'; % vertical eye movement
        end
        if find(blink==Ci)
            artB = 'Blink'; % blink artifact
        end
        if find(disc==Ci)
            artD = 'Discont'; % data drift
        end

        % check if ICA is done on continous or epoched data and plot
        % accordingly
        if ndims(EEGep.icaact) == 3
            % open figure with 4 subplots (1 to 4)
            h = figure('units','normalized','position',[.1 .1 .8 .8]);

            % 1: mean time (ERP) plot from -500 to 500 ms
            subplot(2,2,1), plot(EEGep.times(idx_t(1):idx_t(2)),mean(EEGep.icaact(Ci,idx_t(1):idx_t(2),:),3));
            grid
            title('ERP response');
            ylabel('Arbitrary Units'); xlabel('ms');

            % 2: power spectra from 0 to 60 hz
            icasmthg = (EEGep.icaweights(Ci,:)*EEGep.icasphere)*reshape(EEGep.data(EEG.icachansind,:,:), length(EEGep.icachansind), EEGep.trials*EEGep.pnts); 
            [spectra,freq] = spectopo( icasmthg, EEGep.pnts, EEGep.srate, 'mapnorm', EEGep.icawinv(:,Ci), 'plot','off');
            freq = freq';
            [~,idx_f] = min(abs(freq-60));
%             spectra = 10*log10(spectra./freq);
            subplot(2,2,2),plot(freq(1:idx_f)',spectra(1:idx_f),'LineWidth',2) %,axis([0 60 min(abs(real(spectra(1:idx_f)))) max(abs(real(spectra(1:idx_f))))]);
            ylabel('dezibel'); xlabel('frequency');
            title('power spectra');

            % 3: single trial values from -500 to 500ms
            dataT = squeeze(EEGep.icaact(Ci,idx_t(1):idx_t(2),:));
            dataT = bsxfun(@minus,dataT,mean(dataT(:)));
            dataT = dataT.^2;
            clims = [0 mean(dataT(:))+4*std(dataT(:))];
            xlim = timePeriod;
            ylim = [1 size(EEGep.icaact,3)];
            subplot(2,2,3),imagesc(xlim,ylim,flipud(dataT'),clims);
            colorbar; set(gca,'YDir','normal');
            title('single trial power values');

            % 4: Topographie plot 2d
            subplot(2,2,4),topoplot(EEGep.icawinv(:,Ci),EEGep.chanlocs,'electrodes','on'); colorbar
            title('2d headplot');
            suptitle({['IC ' int2str(Ci) ' : \color{red}' artHE ' ' artVE ' ' artB ' ' artD] ; ...
                ['\color{black} press space to reject and left mouse click to keep the component']})

            % seperate good and bad ic's by mouse click / space press
            button = waitforbuttonpress; 
            if button==0  
                goodICs = [goodICs, Ci];
            else
                badICs = [badICs, Ci];
            end
            close(h)
        else
            % open figure with 3 subplots (1 to 3)
            h = figure('units','normalized','position',[.1 .1 .8 .8]);

            % 1: power spectra from 0 to 100 hz
            icasmthg = (EEGep.icaweights(Ci,:)*EEGep.icasphere)*EEGep.data(EEGep.icachansind,:);
            [spectra,freq] = spectopo( icasmthg, EEGep.pnts, EEGep.srate, 'mapnorm', EEGep.icawinv(:,Ci), 'plot','off');
            [~,idx_f] = min(abs(freq-60));
            subplot(2,2,2),plot(freq(1:idx_f),spectra(1:idx_f),'LineWidth',2) %,axis([0 60 min(spectra(1:idx_f))-1 max(spectra(1:idx_f))+1]);
            title('power spectra');
            
            % 4: Topographie plot 2d
            subplot(2,2,2),topoplot(EEGep.icawinv(:,Ci),EEGep.chanlocs,'electrodes','on'); colorbar
            title('2d headplot');
            suptitle(['CI ' int2str(Ci) ' : ' artHE ' ' artVE ' ' artB ' ' artD])
            
            % 1: ongoin time plot
            subplot(2,2,[3,4]), plot(EEGep.times,EEGep.icaact(Ci,:));
            grid
            title('continous time');
            ylabel('Arbitrary Units'); xlabel('ms');

            % seperate good and bad ic's by mouse click / space press
            button = waitforbuttonpress; 
            if button==0  
                goodICs = [goodICs, Ci];
            else
                badICs = [badICs, Ci];
            end
            close(h)
        end
    end
else
    error('ica_rejection is not 1 or 2 in the csvfile')
end

% reject bad components on single trial or continous data
decompProj = EEG.icawinv(:, goodICs)*eeg_getdatact(EEG, 'component', goodICs, 'reshape', '2d');

if ndims(decompProj) == 3
    decompProj = reshape(decompProj,size(decompProj,1),EEG.pnts,EEG.trials);
    EEG.data(EEG.icachansind,:,:) = decompProj;
else
    EEG.data(EEG.icachansind,:) = decompProj;
end

% convert data to single and delete unnecissary matrices (save disc space)
EEG.icaact  = [];
EEG.icawinv     = EEG.icawinv(:,goodICs);
EEG.icaweights  = EEG.icaweights(goodICs,:);
EEG.specicaact  = [];
EEG.specdata    = [];
EEG.reject      = [];

% remove EOG channels if exist

if find(cell2mat(strfind({EEG.chanlocs.labels},'EOG')))
    label_cell = strfind({EEG.chanlocs.labels},'EOG');
    for i = 1:length(label_cell)
        if ~isempty(label_cell{i} == 1)
            label_idx(i) = 1;
        else
            label_idx(i) = 0;
        end
    end
    
    label_idx = label_idx == 1;
    if ndims(EEG.data)==3
        EEG.data(label_idx,:,:) = [];
    elseif ndims(EEG.data)==2
        EEG.data(label_idx,:) = [];
    else
        error('not 2 and not 3 dimensions. What else could it be?')
    end
    EEG.chanlocs(label_idx) = [];
end

disp([num2str(length(badICs)) ' ICs rejected']);

% clean up a bit
EEG.nbchan = size(EEG.data,1);


% loc file entry
logFile{end+1} = {'ica_cleaned',['data is cleaned by ICA. In total ' num2str(length(badICs)) ...
    ' were removed.']};

if str2double(data_struct.plot_always)==1
    UiO_plots(data_struct,subj_name,EEG,logFile);
end

disp('data ICA cleaning is done')

end