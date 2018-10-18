% Convert EEG-recording so that it is compatible with SpiSoP for sleep
% scoring. This function will save a '.edf' and a '.mat' file in the same
% folder as the original file (filepath)

% convertToSpisop(filepath, filename)

% optional inputs:
% filepath: path of the folder where the file is located
% filename: whole filename with ending (e.g. .vhdr) 


function [] = UiO_convertToSpisop(filepath, filename)

if nargin < 1
    [filename, filepath] = uigetfile('*.*','select the EEG file');
end

addpath('Z:\Matlab_Scripts\eeglab14_1_2b\');
eeglab
close

EEG = pop_loadbv(filepath, filename, [], []);

% 0.1 Hz highpass filter + downsample to 256 hz
EEG = pop_resample(EEG, 125, 0.8, 0.4);
EEG = pop_eegfiltnew(EEG, 'locutoff',0.3,'plotfreqz',0);

% Change TP9 and TP10 to M1 and M2
RField = strcmp(lower({EEG.chanlocs.labels}),'tp9');
EEG.chanlocs(RField).labels = 'M1';
RField = strcmp(lower({EEG.chanlocs.labels}),'tp10');
EEG.chanlocs(RField).labels = 'M2';

% centering and detrending
EEG.data = bsxfun(@minus,EEG.data,mean(EEG.data,2));
EEG.data = detrend(EEG.data);

% eog1 eog2 c3-m2 f3-m2 c4-m1 f4-m1 ueog-m1 leog-m2 cz-(m1-m2)
newEEG = [];
newEEG(1,:) = EEG.data(strcmp(lower({EEG.chanlocs.labels}),'veogu'),:);
EEG.chanlocs(1).labels = 'eog1';
newEEG(2,:) = EEG.data(strcmp(lower({EEG.chanlocs.labels}),'veogl'),:);
EEG.chanlocs(2).labels = 'eog2';
newEEG(3,:) = EEG.data(strcmp(lower({EEG.chanlocs.labels}),'c3'),:) - EEG.data(strcmp(lower({EEG.chanlocs.labels}),'m2'),:);
EEG.chanlocs(3).labels = 'c3-m2';
newEEG(4,:) = EEG.data(strcmp(lower({EEG.chanlocs.labels}),'f3'),:) - EEG.data(strcmp(lower({EEG.chanlocs.labels}),'m2'),:);
EEG.chanlocs(4).labels = 'f3-m2';
newEEG(5,:) = EEG.data(strcmp(lower({EEG.chanlocs.labels}),'c4'),:) - EEG.data(strcmp(lower({EEG.chanlocs.labels}),'m1'),:);
EEG.chanlocs(5).labels = 'c4-m1';
newEEG(6,:) = EEG.data(strcmp(lower({EEG.chanlocs.labels}),'f4'),:) - EEG.data(strcmp(lower({EEG.chanlocs.labels}),'m1'),:);
EEG.chanlocs(6).labels = 'f4-m1';
newEEG(7,:) = EEG.data(strcmp(lower({EEG.chanlocs.labels}),'veogu'),:) - EEG.data(strcmp(lower({EEG.chanlocs.labels}),'m1'),:);
EEG.chanlocs(7).labels = 'eog1-m1';
newEEG(8,:) = EEG.data(strcmp(lower({EEG.chanlocs.labels}),'veogl'),:) - EEG.data(strcmp(lower({EEG.chanlocs.labels}),'m2'),:);
EEG.chanlocs(8).labels = 'eog2-m2';
newEEG(9,:) = EEG.data(strcmp(lower({EEG.chanlocs.labels}),'cz'),:) - (EEG.data(strcmp(lower({EEG.chanlocs.labels}),'m2'),:) - EEG.data(strcmp(lower({EEG.chanlocs.labels}),'m1'),:));
EEG.chanlocs(9).labels = 'Cz-lm';

EEG.chanlocs(10:end) = [];
EEG.data = newEEG;
EEG = eeg_checkset(EEG);

%center data
EEG.data = bsxfun(@minus,EEG.data,mean(EEG.data,2));

new_filename = [filename(1:strfind(filename,'.')-1) '_for_Spisop'];

pop_writeeeg(EEG, [filepath '\' new_filename '.edf'], 'TYPE','EDF');

pop_saveset(EEG, 'filename', new_filename, 'filepath', filepath);

save([filepath '\' new_filename '.mat'],'EEG','-v7.3')

disp('files saved as .edf and .set');

end