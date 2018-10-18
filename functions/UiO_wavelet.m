% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,logFile] = UiO_wavelet(data_struct,subj_name,EEG,logFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the 'ica_cleaned' data (if available)
% subj_name: subject name according to csvfile
% logFile: logFile of previous function. If empty [] this function will
%       load the 'ica_cleaned' logFile (if availeble)
%
% This function decomposes the EEG data into the time-frequency space using 
% complex morelt wavelets. Options (set in csv file) are (1) the wavelet length
% (2) frequency range (3) number of frequencies (4) number of cycles.
% Wavelet decomposition works on both, continuous and epoched data.
% However, please make shure that epoch length are not too short. Otherwise
% edge artifacts may occur.
% 
% by questions:
% b.e.juel@medisin.uio.no
% sevenius.nilsen@gmail.com
% benjamin.thuerer@kit.edu
%

function [EEG,logFile] = UiO_wavelet(data_struct,subj_name,EEG,logFile)


if nargin < 2
    error('provide at least data_struct and subject name. See help UiO_wavelet')
end

% check if EEG structure is provided. If not, load previous data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,'ica_cleaned');   
    else
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end

%% remove EOG channels
if find(cell2mat(strfind({EEG.chanlocs.labels},'EOG')))
    label_cell = strfind({EEG.chanlocs.labels},'EOG');
    for i = 1:length(label_cell)
        if ~isempty(label_cell{i} == 1)
            label_idx(i) = 1;
        end
    end
    
    label_idx = label_idx == 1;
    
%     if ndims(EEG.data)==3
%         EEG.data(label_idx,:,:) = [];
%     elseif ndims(EEG.data)==2
%         EEG.data(label_idx,:) = [];
%     else
%         error('not 2 and not 3 dimensions. What else could it be?')
%     end
    EEG.chanlocs(label_idx) = [];
end


%% if 3D, restrict to mean over trials after wavelet transform
% remove trials and keep only the relevent trials according to the csv-file
% the programming is complicated and maybe week because it is late and I am tired now ;-)

if ndims(EEG.data) > 2
    csv_trials = data_struct.mean_trials;
    delim_idx1 = strfind(csv_trials,',');
    delim_idx2 = strfind(csv_trials,'-');

    delim_add = [delim_idx1,delim_idx2];
    
    if isempty(delim_add)
        good_trials = str2double(csv_trials);
    else
        delim_add = sort(delim_add);

        good_trials = [];

        for i = 1:length(delim_idx2)
            delimDiff_down = find(delim_add < delim_idx2(i));
            delimDiff_up = find(delim_add > delim_idx2(i));
            if isempty(delimDiff_up) && isempty(delimDiff_down)
                good_trials = [good_trials, str2double(csv_trials(1:delim_idx2(i)-1)):str2double(csv_trials(delim_idx2(i)+1:end))];
            elseif isempty(delimDiff_up) && ~isempty(delimDiff_down)
                good_trials = [good_trials, str2double(csv_trials(delim_add(delimDiff_down(end))+1:delim_idx2(i)-1)):str2double(csv_trials(delim_idx2(i)+1:end))];
            elseif ~isempty(delimDiff_up) && isempty(delimDiff_down)
                good_trials = [good_trials, str2double(csv_trials(1:delim_idx2(i)-1)):str2double(csv_trials(delim_idx2(i)+1:delim_add(delimDiff_up(1))-1))];
            elseif ~isempty(delimDiff_up) && ~isempty(delimDiff_down)
                good_trials = [good_trials, str2double(csv_trials(delim_add(delimDiff_down(end))+1:delim_idx2(i)-1)):str2double(csv_trials(delim_idx2(i)+1:delim_add(delimDiff_up(1))-1))];
            end
        end

        numb_idx = 1:length(csv_trials);
        numb_idx(delim_add) = [];

        consecutive_idx = diff(numb_idx);
        consecutive_idx = find(consecutive_idx == 1);

        i = 1;
        k = length(numb_idx);

        while i <= k
            if ~isempty(find(consecutive_idx == i)) && ~isempty(find(consecutive_idx == i+1))
                good_trials = [good_trials,str2double(csv_trials(numb_idx(i):numb_idx(i+2)))];
                i = i+3;
            elseif ~isempty(find(consecutive_idx == i)) && isempty(find(consecutive_idx == i+1))
                good_trials = [good_trials,str2double(csv_trials(numb_idx(i):numb_idx(i+1)))];
                i = i+2;
            else
                good_trials = [good_trials,str2double(csv_trials(numb_idx(i)))];
                i = i+1;
            end
        end
    end

    good_trials = sort(unique(good_trials));

    % now, good_trials are determined but we already deleted bad trials. So
    % we have to add the deleted trials before cutting for good-trials
    for i = 1:length(EEG.accBadEpochs);
        if find(EEG.accBadEpochs(i) == good_trials)
            good_trials(EEG.accBadEpochs(i) == good_trials) = [];
        end
        EEG.data = cat(3,EEG.data(:,:,1:EEG.accBadEpochs(i)-1), zeros(size(EEG.data,1),size(EEG.data,2)), EEG.data(:,:,EEG.accBadEpochs(i):end));
        fake_epoch = EEG.epoch(1);
        first_epochs = EEG.epoch(1:EEG.accBadEpochs(i)-1);
        second_epochs = EEG.epoch(EEG.accBadEpochs(i):end);
        EEG.epoch = [first_epochs,fake_epoch,second_epochs];
    end

    EEG.data = EEG.data(:,:,good_trials);
    EEG.epoch = EEG.epoch(good_trials);
end

%% baseline indices because of fixation cross (only for Study DAVOS)

% for trialI = 1:length(EEG.epoch)
%     bl_idx_cell = [];
%     bl_idx = [];
%     bl_idx_cell = strfind(EEG.epoch(trialI).eventtype,'S  1');
%     for i = 1:length(bl_idx_cell)
%         if ~isempty(bl_idx_cell{i} == 1)
%             bl_idx(i) = 1;
%         end
%     end
%     bl_idx = bl_idx == 1;
%     if isempty(bl_idx)
%         bl_trial(trialI) = {-2000};
%     else
%         bl_trial(trialI) = EEG.epoch(trialI).eventlatency(bl_idx);
%     end
% end

%% start with wavelet transformation for each channel
dataTransformed = zeros(size(EEG.data,1),str2double(data_struct.wavelet_f_numb), ...
        size(EEG.data,2));
    
% perform whole wavelet analysis for each EEG channel
for chan = 1:size(EEG.data,1)
    chan_idx = chan;

    % set lower and higher frequency edge from csv-file (2)
    delim_idx = strfind(data_struct.wavelet_f_range,'-');
    min_f = str2double(data_struct.wavelet_f_range(1:delim_idx-1));
    max_f = str2double(data_struct.wavelet_f_range(delim_idx+1:end));
    
    % min to max f with account of frequencies (3) from csv-file in logarithmic
    % space
    num_freq = str2double(data_struct.wavelet_f_numb);
    frequencies = logspace(log10(min_f),log10(max_f),num_freq); %log
    [EEG.f] = frequencies;
    
    % set wavelet length (1), convolution and number of cycles
    start_length = str2double(data_struct.wavelet_length)/2;
    
    if start_length < 0.5
        error('You do not want to have a wavelet length shorter than 1 sec!')
    end
    
    time = -start_length:1/EEG.srate:start_length;
    half_of_wavelet_size = (length(time)-1)/2;
    n_wavelet = length(time);
    
    if ndims(EEG.data) > 2 %check if epoched or continous data
        n_data = size(EEG.data,2)*size(EEG.data,3);
    else
        n_data = size(EEG.data,2);
    end
    
    n_convolution = n_wavelet+n_data-1;
    n_conv_pow2   = pow2(nextpow2(n_convolution));
    
    % get wavelet cycles (4)
    delim_idx = strfind(data_struct.wavelet_cycles,',');
    low_cycle = str2double(data_struct.wavelet_cycles(1:delim_idx-1));
    high_cycle = str2double(data_struct.wavelet_cycles(delim_idx+1:end));
    wavelet_cycles = logspace(log10(low_cycle),log10(high_cycle),num_freq);
    
    tf_data = zeros(length(frequencies),size(EEG.data,2));
    if ndims(EEG.data) > 2 %check if epoched or continous data
        fft_data = fft(reshape(EEG.data(chan_idx,:,:),1,n_data),n_conv_pow2);
    else
        fft_data = fft(reshape(EEG.data(chan_idx,:),1,n_data),n_conv_pow2);
    end
    

    % now, perform wavelet decomposition for each frequency of each channel
    for fi=1:length(frequencies)

        % create complex wavelet (wavelet * complex exponential * gaussian window) and perform FFT on wavelet
        wavelet = (pi*frequencies(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles(fi) /(2*pi*frequencies(fi)))^2))/frequencies(fi);
        fft_wavelet = fft(wavelet,n_conv_pow2);

        % get convolution (use inverse fft)
        convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
        convolution_result_fft = convolution_result_fft(1:n_convolution);
        convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        
        if ndims(EEG.data) > 2
            convolution_result_fft = reshape(convolution_result_fft,size(EEG.data,2),size(EEG.data,3));
            tf_data(fi,:) = mean(abs(convolution_result_fft).^2,2);
            
            % only for DAVOS study
%             bl_copy = [];
%             sampleBins = 502;
%             for trl = 1:size(EEG.data,3)
%                 step_trial = [];
%                 [~,blIdx(1)] = min(abs(EEG.times - (bl_trial{trl} - (0.5/(1/EEG.srate)))));
%                 [~,blIdx(2)] = min(abs(EEG.times - bl_trial{trl} - 1));
%                 step_trial = abs(convolution_result_fft(blIdx(1):blIdx(2),trl)).^2;
%                 if length(step_trial) < sampleBins
%                     step_trial(end+1:sampleBins) = NaN;
%                 end
%                 bl_copy(trl,:) = step_trial;
%             end
%             bl_mean(fi) = nanmedian(nanmean(bl_copy,1),2);
            
        else
            tf_data(fi,:) = abs(convolution_result_fft).^2;
        end
        
        % tf_data is now frequency * time
    end
    
    dataTransformed(chan,:,:) = tf_data;
%     blTransformed(chan,:) = bl_mean;

    disp(['Wavelet decomposition done for channel ' int2str(chan) ' of ' int2str(size(EEG.data,1)) ])
end

EEG.data = dataTransformed;
% EEG.bl = blTransformed;

%%
% loc file entry
logFile{end+1} = {'after_wavelet',['Wavelet decomposition is done between ' num2str(min_f) ' and ' num2str(max_f) ...
    ' Hz in logarithmic space with ' num2str(num_freq) ' bins and the number of cycles changed according to the frequency from ' ...
    num2str(low_cycle) ' cycles to ' num2str(high_cycle) ' cycles.']};

if str2double(data_struct.plot_always)==1
    UiO_plots(data_struct,subj_name,EEG,logFile);
end

disp('wavelet decomposition is done')


end