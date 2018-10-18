%EEGANALYSIS Main analysis script. Calls stuff
%and reads stuff.
% 
% by questions:
% b.e.juel@medisin.uio.no
% sevenius.nilsen@gmail.com
% benjamin.thuerer@kit.edu

function [] = eeganalysis(sbj,csvname,fldnames1,param)

clear

csvname = uigetfile('*.csv','choose the csv data in private folder (change file format to "All Files")');
csvname(end-3:end) = [];

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath('functions\')

addpath('Z:\Matlab_Scripts\eeglab14_1_2b\');

eeglab

close

%%Calling the read csv file function

param=readcsv(csvname);

% distcomp.feature( 'LocalUseMpiexec', false );
% pool = parpool();


% subject names
fldnames1 = fieldnames(param);

% runs functions for every subject and session provided in the csvfile
for sbj = 1:size(fldnames1,1)
    for ses = 1:size(fieldnames(eval(['param.' fldnames1{sbj}])),1)
        % finding session names and find indices before and after function
        % list
        fldnames2 = fieldnames(eval(['param.' fldnames1{sbj}]));
        fldnames3 = fieldnames(eval(['param.' fldnames1{sbj} '.' fldnames2{ses}]));
        idx_src1 = find(strcmp(fldnames3,'vhdrsource'));
        idx_src2 = find(strcmp(fldnames3,'downsample_rate'));
        
        % define data_struct: specific structure for a specific subject and
        % session which should be processed
        data_struct = eval(['param.' fldnames1{sbj} '.' fldnames2{ses}]);
        
        % finding the right order (n) of the functions. This is complicated
        % because numbers are strings and can have one or two characters.
        % In addition numbers can be delimited by ',' if the function
        % should be run several times
        si_edge = idx_src2-idx_src1-1;
        si = 1;
        n = zeros(si_edge-1,1);
        while si < si_edge
            numbers_string = eval(['param.' fldnames1{sbj} '.' fldnames2{ses} '.' fldnames3{si+idx_src1}]);
            x = strfind(numbers_string, ',');
            if isempty(x)
                n(si,1) = str2double(numbers_string);
            elseif ~isempty(x)
                x_it = 1:length(eval(['param.' fldnames1{sbj} '.' fldnames2{ses} '.' fldnames3{si+idx_src1}]));
                x_diff = setdiff(x_it,x);
                diff_x_diff = diff(x_diff);
                for n_comma = 1:length(diff_x_diff)
                    if n_comma == 1
                        n(si,n_comma) = str2double(numbers_string(x_diff(1:diff_x_diff(n_comma)-1)));
                    else
                        n(si,n_comma) = str2double(numbers_string(x_diff(diff_x_diff(n_comma-1):diff_x_diff(n_comma)-1)));
                    end
                    if length(diff_x_diff) == n_comma
                        n(si,n_comma+1) = str2double(numbers_string(x_diff(diff_x_diff(n_comma):end)));
                    end
                end
            end
            si = si+1;
        end
        
        % now n defines the order of the functions and feval will start
        % each function in the right order
        t = 1;
        logFile = [];
        EEG = [];
        while 1
            if isempty(find(n==t))
                break
            end
            [c,~] = find(n==t);
            [EEG,logFile] = feval(fldnames3{c(1)+idx_src1},data_struct,fldnames1(sbj),EEG,logFile);
            t = t+1;
        end
    end
end

% delete(pool)

end

