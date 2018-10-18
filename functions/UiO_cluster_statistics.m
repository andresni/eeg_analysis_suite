% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,logFile] = UiO_cluster_statistics(time_window,frequency,perm_s,conditions)
% 
% time_window: time window in ms for which statistics should be performed
%   (e.g. [-500 0])
% frequency: frequency band width for which statistics should be performed
%   (e.g. [8 13])
% perm_s: number of randomizations for permutation test
% sessions: two time-points or sessions which should be computed (e.g.
%   {'practice','savings'})
%
% This function performs cluster based statistics on the interaction
% (group*time) effect. The function asks for a group table. This should be
% a .csv file with the first column containing subject IDs and the second
% column containing group IDs
% 
% by questions:
% b.e.juel@medisin.uio.no
% sevenius.nilsen@gmail.com
% benjamin.thuerer@kit.edu
%

function UiO_cluster_statistics(time_window, frequency, perm_s,sessions)

% set folder
%folderpath = uigetdir('','choose the folder with the EEG data');
folderpath = 'C:\Users\Promotion\03_Studie-III_DAVOS\data\EEG_computed';

c_dir = dir(folderpath);
dots_idx = strfind({c_dir.name},'.');
dots_idx = sum(cell2mat(dots_idx) == 1);

%% import group file
% [groupFileName, groupFilePath] = uigetfile('','choose the file with the group table');
if strcmp(sessions{2},'savings')
    groupFileName = 'group_ID_sleep_savings.csv';
else
    groupFileName = 'group_ID_sleep_transfer.csv';
end
groupFilePath = 'C:\Users\Promotion\03_Studie-III_DAVOS\data\';

groupFile = [groupFilePath groupFileName];

fid = fopen(groupFile,'r');
lineIdx = 1;
nextLine = fgetl(fid);
group_ID = cell(120,2);

while ~isequal(nextLine,-1)
    idx = strfind(nextLine,';');
    group_ID{lineIdx,1} = nextLine(1:idx-1);
    group_ID{lineIdx,2} = nextLine(idx+1:end);
    lineIdx = lineIdx + 1;
    nextLine = fgetl(fid);
end

group_ID(lineIdx:end,:) = [];
gidx = {};

data_session_one = [];

for i = 1:length(group_ID)
    sbj_name = group_ID{i,1};
    gidx(i) = group_ID(find(strcmp({group_ID{:,1}},sbj_name)),2);
    disp(['Load dataset ' sbj_name ' for session ' sessions{1}]);
    load_file = [folderpath '\' sbj_name '\' sessions{1} '\' sbj_name '_' sessions{1} '_after_norm_power.mat'];
    load(load_file)
    [~,t_wndw(1)] = min(abs(EEG.times-time_window(1)));
    [~,t_wndw(2)] = min(abs(EEG.times-time_window(2)));
    [~,f_rng(1)] = min(abs(EEG.f-frequency(1)));
    [~,f_rng(2)] = min(abs(EEG.f-frequency(2)));
    data_session_one(i,:) = mean(mean(EEG.data(:,f_rng(1):f_rng(2),t_wndw(1):t_wndw(2)),2),3);
end

data_session_two = [];

for i = 1:length(group_ID)
    sbj_name = group_ID{i,1};
    disp(['Load dataset ' sbj_name ' for session ' sessions{2}]);
    load_file = [folderpath '\' sbj_name '\' sessions{2} '\' sbj_name '_' sessions{2} '_after_norm_power.mat'];
    load(load_file)
    [~,t_wndw(1)] = min(abs(EEG.times-time_window(1)));
    [~,t_wndw(2)] = min(abs(EEG.times-time_window(2)));
    [~,f_rng(1)] = min(abs(EEG.f-frequency(1)));
    [~,f_rng(2)] = min(abs(EEG.f-frequency(2)));
    data_session_two(i,:) = mean(mean(EEG.data(:,f_rng(1):f_rng(2),t_wndw(1):t_wndw(2)),2),3);
end

chanlocs = EEG.chanlocs;

%% create table for RM ANOVA
Meas = dataset([1 2]','VarNames',{'Measurements'});
    
P_map = zeros(size(EEG.data,1),1);
F_map = zeros(size(EEG.data,1),1);

for i = 1:size(data_session_one,2)
    t = table(gidx',data_session_one(:,i),data_session_two(:,i),'VariableNames',{'group','meas_1','meas_2'});   
    
    % group effect
%     [tbl_p,tbl_m] = anova1(t.meas_2,gidx','off');
%     P_map(i) = tbl_p;
%     F_map(i) = abs(tbl_m{2,5});
    
    % Interaction effect (group*time)
    rm = fitrm(t,'meas_1-meas_2 ~ group','WithinDesign',Meas);
    ranovatbl = ranova(rm);
    P_map(i) = ranovatbl.pValue(2);
    F_map(i) = abs(ranovatbl.F(2));
end

%%
% topoplot grids of P-values
[~,topo_P] = topoplot(P_map,chanlocs,'noplot','on');

% topoplot grids of F-values
[~,topo_F,plotrad,xmesh,ymesh] = topoplot(F_map,chanlocs,'noplot','on');

% building clusters
topo_P = topo_P <= .05;
[obs_P,nblobs] = bwlabeln(topo_P);

clustsum_obs = zeros(1,nblobs);

abs_topo_F = abs(topo_F);
for i=1:nblobs
    clustsum_obs(i) = sum(abs_topo_F(obs_P(:)==i));
end

%% Interaction effect:
if ~isempty(clustsum_obs)
    permP_map = ones(size(EEG.data,1),1);
    permF_map = ones(size(EEG.data,1),1);


    %combine datasets and shuffle data

    % for group difference
    % perm_set = data_session_one; % or data_session_two;

    %for interaction
    perm_set = [data_session_one;data_session_two];


    t_p = t;
    n = size(perm_set,1);
    percentages = perm_s/10:perm_s/10:perm_s;
    percentages2 = 10:10:100;

    perm_i = 1;
    max_perm = zeros(1,perm_s);
    disp('start permutation testing. This may take a while');

    while perm_i <= perm_s
       for i = 1:size(EEG.data,1)
           rpt = randperm(n);
           permData = perm_set(rpt,i);
           t_p.meas_1 = permData(1:length(group_ID));
           t_p.meas_2 = permData(length(group_ID)+1:end);

           % group effect
    %         [perm_tbl_p,perm_tbl_m] = anova1(t_p.meas_2,gidx','off');
    %         permP_map(i) = perm_tbl_p;
    %         permF_map(i) = abs(perm_tbl_m{2,5});

            % Interaction effect (group*time)
            rm_p = fitrm(t_p,'meas_1-meas_2 ~ group','WithinDesign',Meas);
            perm_tbl = ranova(rm_p);
            permP_map(i) = perm_tbl.pValue(2);
            permF_map(i) = abs(perm_tbl.F(2));
       end
           % topoplot grids of P-values
        [~,perm_P] = topoplot(permP_map,chanlocs,'noplot','on');

        % topoplot grids of F-values
        [~,perm_F] = topoplot(permF_map,chanlocs,'noplot','on');

        % building clusters
        perm_P = perm_P <= .05;
        [perm_P,perm_nblobs] = bwlabeln(perm_P);

        clustsum = zeros(1,perm_nblobs);
        perm_F = abs(perm_F);
        if find(perm_nblobs)
            for iii=1:perm_nblobs
                clustsum(iii) = sum(perm_F(perm_P(:)==iii));
            end
            max_perm(1,perm_i) = max(clustsum);
        end

        if ~isempty(find(percentages == perm_i))
            disp([num2str(percentages2(percentages == perm_i)) ' % done']);
        end
        perm_i = perm_i +1;

        clear rpt permData rm_p perm_tbl permP_map permF_map perm_map nblobs clustsum perm_nblobs
    end


    stat_cluster = sum(repmat(abs(max_perm)',1,length(clustsum_obs)) > repmat(clustsum_obs,length(max_perm),1), 1) / perm_s;
else
    stat_cluster = [];
end

% start plotting
if find(stat_cluster)
    combi = [];
    [~,idx] = find(stat_cluster < 0.05);
    sum_map = zeros(length(idx),size(topo_F,1),size(topo_F,2));
    for i = 1:length(idx)
        sum_map(i,:,:) = obs_P == idx(i);
    end
    sum_map = squeeze(sum(sum_map,1));
    sum_map( sum_map >= 1 ) = 100;
    combi(1,:,:) = sum_map;
    combi(2,:,:) = topo_F;
    combi = squeeze(sum(combi,1));
    h = figure;
    toporeplot(combi,'plotrad',plotrad,'xsurface',xmesh,'ysurface',ymesh,'chanlocs',chanlocs,'electrodes','on','headrad',0.54,'style','map');
    set(gca,'clim',[0 12]);
    title(['topoplot for ' num2str(frequency(1)) ' to ' num2str(frequency(2)) ' hz for the time window ' num2str(time_window(1)) ' to ' num2str(time_window(2)) ' ms.' ...
        ' and session ' sessions{2}]);
else
    disp('no significant clusters');
end

end