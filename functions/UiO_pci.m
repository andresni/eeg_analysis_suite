% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,locFile] = UiO_pci(data_struct,subj_name,EEG,locFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the 'after_ica' data (if availeble)
% subj_name: subject name according to csvfile
% locFile: locFile of previous function. If empty [] this function will
%       load the 'after_ica' locFile (if availeble)
%
% A function for calculating the Perturbational Complexity Index, given a
% binary matrix of significant sources of cortical activity, and an
% integer value for the downsampling (1: no downsampling, 2: every other
% point, 3: every third etc)
%
% by questions: b.e.juel@medisin.uio.no or benjamin.thuerer@kit.edu
%
function [EEG,locFile] = UiO_pci(data_struct,subj_name,EEG,locFile)


if nargin < 2
    error('provide at least data_struct and subject name. See help UiO_pca')
end

% check if EEG structure is provided. If not, load previous data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,locFile] = UiO_load_data(data_struct,subj_name,'after_inverse_model');   
    else
        [EEG,locFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end

fprintf('\n ****************************************\n');
fprintf('             CALCULATING PCI');
fprintf('\n ****************************************\n');


% Reorganizing the SS matrix to the most ordered state
%SSsum = sum(SS,2);
SSsum = sum(EEG.significant_sources',1);
[~,index]=sort(SSsum);
sorted=EEG.significant_sources(index,:);


% Downsampling to make "PCI" calculation faster.
sorted = sorted(1:str2double(data_struct.downsample_pci_source):end,1:str2double(data_struct.downsample_pci_time):end);
SST = sorted;

% find 0 time index and restrict SST from 0ms to 300ms after
% stimulus
time_line = EEG.times(1:str2double(data_struct.downsample_pci_time):end);
[~,ind_zero] = min(abs(time_line-0));
[~,ind_350] = min(abs(time_line-300));
SST = SST(:,ind_zero:ind_350);

% Printing out a figure, with the significant source matrix, for
% visualization.
figure;imagesc(SST);
title('Who knows what the PCI will be?!');
colormap([0,0,0;1,1,1]);
hold on
pause(0.2);


disp('The fun begins!');

PCI = UiO_calc_lz_complexity(SST(:),'exhaustive',1);

title(['\fontsize{60} \color{red} PCI = ' num2str(PCI)],'interpreter','tex')

EEG.PCI = PCI;

locFile{end+1} = {'after_pci',['pci is calculated for this dataset and is' num2str(PCI)]};

if str2double(data_struct.plot_always)==1
    UiO_plots(data_struct,subj_name,EEG,locFile);
end

end