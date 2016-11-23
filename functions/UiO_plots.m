% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,locFile] = UiO_plots(data_struct,subj_name,EEG,locFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the last processed data (if availeble)
% subj_name: subject name according to csvfile
% locFile: locFile of previous function. If empty [] this function will
%       load the last processed locFile (if availeble)
%
% This function will plot the EEG at the step provided in the csvfile or at
% the step which is processed before this plot function
% 
% by questions: benjamin.thuerer@kit.edu
% 
function [EEG,locFile] = UiO_plots(data_struct,subj_name,EEG,locFile)

if nargin < 2
    error('provide at least data_struct and subject name. See help UiO_trials')
end

exist data_struct.load_data;

if ans == 0
    data_struct.load_data = '0';
end

% check if EEG structure is provided. If not, load previous data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,locFile] = UiO_load_data(data_struct,subj_name,'epoched');   
    else
        [EEG,locFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end

% plot according to the data and what make sense
if strcmp(locFile{end}{1},'after_tms')
    eegplot(EEG.data);
    
    % compute power spectral density
    concat_data = reshape(EEG.data,1,numel(EEG.data));
    transform = fft(concat_data,EEG.srate/2);
    PSD = transform.*conj(transform)/(EEG.srate/2);
    f = 1000/(EEG.srate/2)*(0:127);
    figure;
    plot(f,PSD(1:128))
    title('Power spectral density over all channels')
    xlabel('Frequency (Hz)')
    
elseif strcmp(locFile{end}{1},'preprocessed')
    eegplot(EEG.data);
    
    % compute power spectral density
    concat_data = reshape(EEG.data,1,numel(EEG.data));
    transform = fft(concat_data,EEG.srate/2);
    PSD = transform.*conj(transform)/(EEG.srate/2);
    f = 1000/(EEG.srate/2)*(0:127);
    figure;
    plot(f,PSD(1:128))
    title('Power spectral density over all channels')
    xlabel('Frequency (Hz)')
    
elseif strcmp(locFile{end}{1},'after_pca')
    eegplot(EEG.data);
    
    % compute power spectral density
    concat_data = reshape(EEG.data,1,numel(EEG.data));
    transform = fft(concat_data,EEG.srate/2);
    PSD = transform.*conj(transform)/(EEG.srate/2);
    f = 1000/(EEG.srate/2)*(0:127);
    figure;
    plot(f,PSD(1:128))
    title('Power spectral density over all channels')
    xlabel('Frequency (Hz)')
    
elseif strcmp(locFile{end}{1},'after_ica')
    eegplot(EEG.data);
    
    % compute power spectral density
    concat_data = reshape(EEG.data,1,numel(EEG.data));
    transform = fft(concat_data,EEG.srate/2);
    PSD = transform.*conj(transform)/(EEG.srate/2);
    f = 1000/(EEG.srate/2)*(0:127);
    figure;
    plot(f,PSD(1:128))
    title('Power spectral density over all channels')
    xlabel('Frequency (Hz)')
    
    % if icaact empty
    if isempty(EEG.icaact)
        if ndims(EEG.data)==3
            for i=1:EEG.trials
                EEG.icaact(:,:,i)=(EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:,i);
            end
        else
            EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data(EEG.icachansind,:);
        end
    end
    
    % Check if ICA done over trials and then plot ICs with different
    % informations
    for Ci = 1:size(EEG.icaact,1)
        if ndims(EEG.icaact) == 3
            % open figure with 4 subplots (1 to 4)
            h = figure;

            % 1: mean time (ERP) plot from -500 to 500 ms
            subplot(2,2,1), plot(EEG.times(idx_t(1):idx_t(2)),mean(EEG.icaact(Ci,idx_t(1):idx_t(2),:),3));
            grid
            title('ERP response');
            ylabel('Arbitrary Units'); xlabel('ms');

            % 2: power spectra from 0 to 100 hz
            icasmthg = (EEG.icaweights(Ci,:)*EEG.icasphere)*reshape(EEG.data(EEG.icachansind,:,:), length(EEG.icachansind), EEG.trials*EEG.pnts); 
            [spectra,freq] = spectopo( icasmthg, EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,Ci), 'plot','off');
            [~,idx_f] = min(abs(freq-60));
            subplot(2,2,2),plot(freq(1:idx_f),spectra(1:idx_f),'LineWidth',2),axis([0 60 min(spectra(1:idx_f))-1 max(spectra(1:idx_f))+1]);
            title('power spectra');

            % 3: single trial values from -500 to 500ms
            dataT = squeeze(EEG.icaact(Ci,idx_t(1):idx_t(2),:));
            dataT = bsxfun(@minus,dataT,mean(dataT(:)));
            dataT = dataT.^2;
            clims = [0 mean(dataT(:))+4*std(dataT(:))];
            xlim = timePeriod;
            ylim = [1 size(EEG.icaact,3)];
            subplot(2,2,3),imagesc(xlim,ylim,flipud(dataT'),clims);
            colorbar; set(gca,'YDir','normal');
            title('single trial power values');

            % 4: Topographie plot 2d
            subplot(2,2,4),topoplot(EEG.icawinv(:,Ci),EEG.chanlocs,'electrodes','on'); colorbar
            title('2d headplot');
            suptitle(['CI ' int2str(Ci)])

            waitforbuttonpress
            close(h)
        else
            % open figure with 3 subplots (1 to 3)
            h = figure;

            % 1: power spectra from 0 to 100 hz
            icasmthg = (EEG.icaweights(Ci,:)*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            [spectra,freq] = spectopo( icasmthg, EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,Ci), 'plot','off');
            [~,idx_f] = min(abs(freq-60));
            subplot(2,2,2),plot(freq(1:idx_f),spectra(1:idx_f),'LineWidth',2),axis([0 60 min(spectra(1:idx_f))-1 max(spectra(1:idx_f))+1]);
            title('power spectra');
            
            % 4: Topographie plot 2d
            subplot(2,2,2),topoplot(EEG.icawinv(:,Ci),EEG.chanlocs,'electrodes','on'); colorbar
            title('2d headplot');
            suptitle(['CI ' int2str(Ci)])
            
            % 1: ongoin time plot
            subplot(2,2,[3,4]), plot(EEG.times,EEG.icaact(Ci,:));
            grid
            title('continous time');
            ylabel('Arbitrary Units'); xlabel('ms');

            waitforbuttonpress
            close(h)
        end
    end
    
elseif strcmp(locFile{end}{1},'ica_cleaned')
    eegplot(EEG.data);
    
    % compute power spectral density
    concat_data = reshape(EEG.data,1,numel(EEG.data));
    transform = fft(concat_data,EEG.srate/2);
    PSD = transform.*conj(transform)/(EEG.srate/2);
    f = 1000/(EEG.srate/2)*(0:127);
    figure;
    plot(f,PSD(1:128))
    title('Power spectral density over all channels')
    xlabel('Frequency (Hz)')
elseif strcmp(locFile{end}{1},'after_inverse_model')
    figure;
    imagesc(EEG.significant_sources)
%     axis([EEG.times(1) EEG.times(end) 1 size(EEG.significant_sources,1)]);
elseif strcmp(locFile{end}{1},'after_pci')  
    SSsum = sum(EEG.significant_sources',1);
    [~,index]=sort(SSsum);
    sorted=EEG.significant_sources(index,:);
    
    % Downsampling to make "PCI" calculation faster.
    sorted = sorted(1:str2double(data_struct.downsample_pci):end,1:end);
    SST = sorted;
    
    % cut time domain to 100ms before and 350ms after stimulus
    [~,ind_zero] = min(abs(EEG.times-0));
    time_ind = 1/EEG.srate;
    SST = SST(:,ind_zero-(0.1/time_ind):ind_zero+(0.35/time_ind));
    
    figure;
    imagesc(SST);
    title(['PCI: ' num2str(EEG.PCI)]);
    colormap([0,0,0;1,1,1]);
end

pause(1.5);

end