% summary Plot function
% 
% by questions:
% b.e.juel@medisin.uio.no
% sevenius.nilsen@gmail.com
% benjamin.thuerer@kit.edu

function [EEG,logFile] = UiO_summaryplot(data_struct,subj_name,EEG,logFile)
% 
% %% load raw-data for ERP:

% check if the filepath is seperated by / or \ and and seperate file-path
% from file-name
% if isempty(strfind(data_struct.vhdrsource,'\'))
%     char_idx = strfind(data_struct.vhdrsource,'/'); 
% else
%     char_idx = strfind(data_struct.vhdrsource,'\');
% end
% 
% data_name = data_struct.vhdrsource(char_idx(end)+1:end);
% data_path = data_struct.vhdrsource(1:char_idx(end));
% 
% % check if header ending is provided, otherwise add .vhdr to the file-name
% if strcmp(data_name(end-4:end),'.vhdr')
%     RAWEEG = pop_loadbv(data_path, data_name, [], []);
% else
%     RAWEEG = pop_loadbv(data_path, [data_name '.vhdr'], [], []);
% end

% 
% % sampling rate
% if str2double(data_struct.downsample_rate) > 1
%     Nsrate = str2double(data_struct.downsample_rate); %which sampling rate for resample?
% elseif str2double(data_struct.downsample_rate) == 0
%     Nsrate = [];
% else
%     warning('sampling rate to low. New sampling rate: 1000 hz')
%     Nsrate = 1000;
% end
% 
% RAWEEG = pop_resample(RAWEEG, Nsrate, 0.8, 0.4);
% 


% RAWEEG = pop_epoch( RAWEEG, {  'R128'  }, [str2double(data_struct.trial_start)/EEG.srate str2double(data_struct.trial_end)/EEG.srate] ...
%     , 'newname', ' resampled epochs', 'epochinfo', 'yes');
% 
% %remove  EOG channels
% RAWEEG.data(end-1:end,:,:) = [];
% 

% % rawERP = mean(RAWEEG.data(:,ERPstart:ERPend-1,:),3); 
% rawERP = squeeze(mean(mean(RAWEEG.data,1),3));
% rawERP = rawERP(ERPstart+1:ERPend);
% rawERP = rawERP-mean(rawERP(1:95));
% % rawERP = abs(rawERP);
% % rawERP = rawERP-mean(rawERP);


%%
% check if EEG structure is provided. If not, load previous data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,'after_pci');   
    else
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end

% necessary variables and parameters
% for title:
PCI = EEG.PCI; %PCI value
chanNum = 62-length(EEG.CHremoved); % number of channels left (non-interpolated)
trialNum = size(EEG.data,3); % remaining trials
comps = size(EEG.icaweights,1);% number of accepted components


% for first plot
% perhaps load proper file?
fs = EEG.srate;
ERPstart = -(EEG.xmin+0.1)*fs;
ERPend = ERPstart+0.4/(1/fs);   
times = EEG.times(ERPstart+1:ERPend);
ERPaxis = [-10,10];

% for second plot
% perhaps load proper file?
cleanERP = mean(EEG.data(:,ERPstart:ERPend-1,:),3); % 

% for third plot
GMFP = mean(abs(cleanERP));

% for topo-plots
% first 50
first50 = abs(GMFP(round(0.115*fs):round(0.135*fs)));
[aux, ind] = max(first50);
ind = ind + round(0.115*fs);
topo50 = cleanERP(:,ind);
ind1 = ind;

% 50-100
next100 = abs(GMFP(round(0.13*fs):round(0.17*fs)));
[aux, ind] = max(next100);
ind = ind + round(0.13*fs);
topo100 = cleanERP(:,ind);
ind2 = ind;

% 100-200
next200 = abs(GMFP(round(0.175*fs):round(0.25*fs)));
[aux, ind] = max(next200);
ind = ind + round(0.175*fs);
topo200 = cleanERP(:,ind);
ind3 = ind;

% 200-300
next300 = abs(GMFP(round(0.3*fs):round(0.4*fs)));
[aux, ind] = max(next300);
ind = ind + round(0.3*fs);
topo300 = cleanERP(:,ind);
ind4 = ind;


% for last plot
Sources = EEG.significant_sources(:,ERPstart:ERPend-1); % 
SSsum = sum(Sources',1);
[~,index]=sort(SSsum);
sorted=Sources(index,:);
sourceProp = sum(sum(Sources(:,0.1*fs:end)))/numel(Sources(:,0.1*fs:end)); %proportion of active sources


% making plots
rows = 6;
cols = 4;

figure;
set(gcf,'units','centimeter','position',[10,1,15,24]);

% subplot(rows,cols,1:4); % Here the "raw" ERP goes
% semilogy(times,rawERP);ylim(ERPaxis);ylabel({'"raw" ERP'; '(\mu V)'});
% mainTitle = {data_struct.session;['PCI score: ' num2str(PCI) ', active source proportion = ' num2str(sourceProp)];...
%     ['channels: ' num2str(chanNum) ', trials: ' num2str(trialNum) ', components: ' num2str(comps)]};
% 
% set(gca,'XTickLabel',[]);

% Here the cleaned/final ERP goes
% first row plot
subplot(rows,cols,1:cols); 
plot(times,cleanERP);ylim(ERPaxis);ylabel({'clean ERP'; '(\mu V)'});
mainTitle = {['Session name: ' data_struct.session];['PCI score: ' num2str(PCI) ', active sources = ' num2str(sourceProp*100) '%'];...
    ['Accepted # of:   channels: ' num2str(chanNum) ', trials: ' num2str(trialNum) ', components: ' num2str(comps)]};
title(mainTitle, 'Interpreter', 'none');
set(gca,'XTickLabel',[]);

% Here the cleaned GMFP goes
% second row
subplot(rows,cols,cols+1:2*cols); 
semilogy(times,GMFP,'LineWidth',2);ylim(ERPaxis);ylabel({'GMFP'; '(\mu V)'});
set(gca,'XTickLabel',[]);


% Topography plots
% third row
subplot(rows,cols,2*cols+1); % Here the topography from the peak in 0-50ms goes
topoplot(topo50,EEG.chanlocs);title(['peak at ' num2str((ind1-0.1*fs)) 'ms'])
subplot(rows,cols,2*cols+2); % Here the topography from the peak in 50-100ms goes
topoplot(topo100,EEG.chanlocs);title(['peak at ' num2str((ind2-0.1*fs)) 'ms'])
subplot(rows,cols,2*cols+3); % Here the topography from the peak in 100-200ms goes
topoplot(topo200,EEG.chanlocs);title(['peak at ' num2str((ind3-0.1*fs)) 'ms'])
subplot(rows,cols,2*cols+4); % Here the topography from the peak in 200-300ms goes
topoplot(topo300,EEG.chanlocs);title(['peak at ' num2str((ind4-0.1*fs)) 'ms'])


subplot(rows,cols,3*cols+1:6*cols); % Here the downsampled significant source matrix goes
imagesc(times,1:size(Sources,1),sorted);
xlabel('time (ms)'); ylabel('Sources');
set(gca,'YTickLabel',[]);


% saving figure

% check if file_save path is provided and check for / or \
if str2double(data_struct.save_folder) == 0
    if isempty(strfind(data_struct.vhdrsource,'\'))
        char_idx = strfind(data_struct.vhdrsource,'/'); 
    else
        char_idx = strfind(data_struct.vhdrsource,'\');
    end
    data_path = data_struct.vhdrsource(1:char_idx(end));
else
    if isempty(strfind(data_struct.save_folder,'\'))
        char_idx = strfind(data_struct.save_folder,'/'); 
    else
        char_idx = strfind(data_struct.save_folder,'\');
    end
    data_path = [data_struct.save_folder(1:char_idx(end))];
end

figName = [data_path subj_name{1} '\' data_struct.session '\' subj_name{1} '_' data_struct.session '_after_pci']
export_fig(figName, '-transparent', '-tiff', '-nofontswap', '-m2')



