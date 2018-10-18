% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,logFile] = UiO_inverse_model(data_struct,subj_name,EEG,logFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the 'after_ica' data (if availeble)
% subj_name: subject name according to csvfile
% logFile: logFile of previous function. If empty [] this function will
%       load the 'after_ica' logFile (if availeble)
%
% This function will do inverse modelling on the EEG data. Make sure that
% you provide the lead field matrix (lft), surface (srf) and location (loc)
% files of the subject (see csvfile).
% 
% by questions:
% b.e.juel@medisin.uio.no
% sevenius.nilsen@gmail.com
% benjamin.thuerer@kit.edu
%
function [EEG,logFile] = UiO_inverse_model(data_struct,subj_name,EEG,logFile)

if nargin < 2
    error('provide at least data_struct and subject name. See help UiO_pca')
end


% check if EEG structure is provided. If not, load previous data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,'ica_cleaned');   
    else
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end

% remove EOG channels if exist
if find(strcmp({EEG.chanlocs.labels},'HEOG'))
    HEOG_idx = strcmp({EEG.chanlocs.labels},'HEOG');
    EEG.chanlocs(HEOG_idx).labels = [];
    if ndims(EEG.data)==3
        EEG.data(HEOG_idx,:,:) = [];
    elseif ndims(EEG.data)==2
        EEG.data(HEOG_idx,:) = [];
    else
        error('not 2 and not 3 dimensions. What do you want?')
    end
end
if find(strcmp({EEG.chanlocs.labels},'VEOG'))
    VEOG_idx = strcmp({EEG.chanlocs.labels},'VEOG');
    EEG.chanlocs(VEOG_idx).labels = [];
    if ndims(EEG.data)==3
        EEG.data(VEOG_idx,:,:) = [];
    elseif ndims(EEG.data)==2
        EEG.data(VEOG_idx,:) = [];
    else
        error('not 2 and not 3 dimensions. What do you want?')
    end
end

clear VEOG_idx HEOG_idx

% set baseline Length in ms before pulse
baselineLength = str2double(data_struct.baseline_sources);

% load lead field matrix, surface and locations from MRI (using BESA)
[~, L]=readBESAlft(data_struct.lftsource);
% srf = readBESAsrf(data_struct.srfsource);
% locs = readBESAloc(data_struct.locsource);

% downsampling rate 
baselineS = EEG.srate*baselineLength/1000;

sources = size(L,2)/3;
channels = size(L,1);

% lead field matrix of BESA MRI is in blocks: x, y, z
% Separating the X, Y, and Z directions in the lead field matrix
Lx = L(:,1:sources);
Ly = L(:,sources+1:sources*2);
Lz = L(:,2*sources+1:sources*3);

% lead field matrix of BESA matlab might be different: x, y, z, x, y, z
% % oprion 2: x1,y1,z1,x2,y2,z2,...,xn,yn,zn
% Lx = L(:,1:3:end);
% Ly = L(:,2:3:end);
% Lz = L(:,3:3:end);

%%
% finding nodes that most closely matches surface locations
% if surface_sources
%     disp('Picking sources that are close to the surface. Downsampling.');
%     tic;
%     for i = 1:size(srf.CoordsVertices,1)
%         aux = repmat(srf.CoordsVertices(i,:),size(locs,1),1);
%         [diff(i), ind(i)] = min(sum(sqrt((aux-locs).^2),2));
%     end
%     %Keepin only the unique elements
%     [ind, keep] = unique(ind);
%     diff = diff(keep);
%     toc
% else
%     ind = 1:sources;
% end
% % Drawing random sources we will use in our reconstruction. 
% % CONSIDER USING AVERAGE OF NEIGHBOURS IN STEAD
% usedSources = sort(datasample(1:sources,floor(sources/downsample),'Replace', false));
% locs = locs(usedSources,:);

% Setting the only sources considered in reconstruction, by drawing random 
% nodes (number based on downsampling factor) within the sources near the surface

% downsampling
% indices = sort(datasample(ind,total_sources,'Replace', false));
% usedSources = indices;
% locDownsampled = locs(usedSources,:);

%% Projection vector values and downsampling
% this creates the length of the projection vector according to the lead field matrix
Projection_orig = sqrt(Lx.^2 + Ly.^2 + Lz.^2);

clear Lx Ly Lz

% this downsamples the new lead field matrix (projections) to the max influential
% sources to the electrodes (these are the sources which affect the EEG
% most)
ds_idx = str2double(data_struct.downsample_source); % how much percentage to keep

Max_proj = max(Projection_orig,[],1);
sorted_max = sort(Max_proj);
Thresh = sorted_max(end-floor((length(sorted_max)/100*ds_idx))+1:end);
source_idx = Max_proj >= min(Thresh);

EEG.source_idx_keep = find(source_idx); % store the indices of keeping sources
L_d = Projection_orig(:,source_idx); % downsampled lead field matrix

clear Projection_orig sorted_max Thresh Max_proj source_idx

%% Generate the inverse model
disp('generating inverse model')

% regularization term (put in lambda)
lambda = 10^(-6);
I = eye(channels); % generating identity matrix for regularization

W = (L_d*L_d' + lambda*I);

disp('calculating source activity')

% for plotting process stage in percentage in while loop
t_i = size(EEG.data,3);
t = 0:floor(t_i/10):t_i;

S = zeros(size(L_d,2),size(EEG.data,2),size(EEG.data,3));

i = 1;
while i <= size(EEG.data,3)
    Ytrial = squeeze(EEG.data(:,:,i));
    S(:,:,i) = L_d'*(W\Ytrial);
    i = i+1;
    if find(i==t)
        fprintf([' ' int2str(100/t_i*i+1) '%% processed \n']);
    end
end

n_size = size(L_d,2);
clear L_d L W Ytrial I locs

%% Renormalizing the source activity matrix (old)
% disp('Renormalizing the Source activity');

% baselineMatrix = S(:,1:baselineS-1,:);
% 
% reBL = reshape(baselineMatrix,size(baselineMatrix,1),size(baselineMatrix,2)*size(baselineMatrix,3));
% meanB = mean(reBL,2);
% 
% tic;
% 
% Si = 1;
% while Si <= size(S,1)
%     S(Si,:,:) = (squeeze(S(Si,:,:))-meanB(Si))./meanB(Si);
%     Si = Si+1;
% end
% 
% toc

%% new renormalizing and thresholding (you don't neet UiO_bootstrap now!

% set the baseline and find the mean and std
disp('computing baseline and standard deviation for renormalizing')
disp('This may take a while...');

baselineMatrix = S(:,1:baselineS-1,:);
i_std = 1;
meanB = [];
stdB = [];
while i_std <= size(baselineMatrix,1)
    sqzd_bl = squeeze(baselineMatrix(i_std,:,:));
    meanB(i_std,:) = 1/numel(sqzd_bl)*sum(sqzd_bl(:));
    stdB(i_std,:) = sqrt(1/(numel(sqzd_bl)-1)*(sum(sum((sqzd_bl-meanB(i_std,:)).^2))));
    i_std = i_std+1;
end

clear baselineMatrix i_std sqzd_bl


% create a loop which finds the best threshold according to 99.9 % of the
% significant sources are after the Pulse and then use the threshold which
% is nearest to 50% of significant sources within -100ms to 300ms compared
% all sources from -100ms to 300ms (0ms=pulse)

k = 2; % where the multiplier for the std should start
i = 100; % how much percent of significant sources should the loop start with
A = []; % stores the percentage of significant sources within -100 and 300ms
m = 0; % iteration

Average = squeeze(mean(S,3)); 
t = 0:2:50; %for plotting progess in significant sources

clear S

% find zero time index and sample bins for time window in while loop
[~,ind_zero] = min(abs(EEG.times-0));
time_ind = 1/EEG.srate;

while i >= 99.99
    m=m+1;
    Thresh = meanB+k*stdB; 
    SA = abs(Average)>abs(repmat(Thresh,1,size(EEG.data,2)));
    i = sum(sum(SA(:,ind_zero:end)))/sum(SA(:))*100;
    if isnan(i)
        i = 100;
    end
    k = k-0.01;    
    A(m) = sum(sum(SA(:,ind_zero:ind_zero+(0.30/time_ind))))/numel(SA(:,ind_zero:ind_zero+(0.30/time_ind)))*100;
    if find(ceil(A(m))==t)
        fprintf([ num2str(A(m)) '%% significant sources \n']);
    end
    if A(m) > 50 %break if 50% (best tradeoff) achieved
        break
    end    
end

Anew = A(end);
disp([ num2str(Anew) '% significant sources between 0 and 350ms detected']);
% figure
% imagesc(SA)

EEG.significant_sources = SA;

% loc file entry
logFile{end+1} = {'after_inverse_model',['data now provides source activity of ' ...
    num2str(n_size) ' sources. The ratio is ' num2str(Anew)]};

if str2double(data_struct.plot_always)==1
    UiO_plots(data_struct,subj_name,EEG,logFile);
end

disp('data inverse modelling is done')

end

