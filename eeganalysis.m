function [] = eeganalysis( csvname )
%EEGANALYSIS Main analysis script. Calls stuff
%and reads stuff.
%By André Sevenius Nilsen & Benjamin Thürer
% sevenius.nilsen@gmail.com
% benjamin.thuerer@kit.edu

%%Calling the read csv file function
param=readcsv(csvname);

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
        
        % finding the right order (n) of the functions. This is complex
        % because numbers are strings and can have one or two characters.
        % In addition numbers can be delimited by ',' if the function
        % should be run several times
        si_edge = idx_src2-idx_src1-1;
        si = 1;
        n = zeros(si_edge,1);
        while si < si_edge+1
            numbers_string = eval(['param.' fldnames1{sbj} '.' fldnames2{ses} '.' fldnames3{si+idx_src1}]);
            x = strfind(numbers_string, ',');
            if isempty(x)
                n(si,1) = str2double(numbers_string);
            elseif ~isempty(x)
                x_it = 1:length(eval(['param.' fldnames1{sbj} '.' fldnames2{ses} '.' fldnames3{si+idx_src1}]));
                x_diff = setdiff(x_it,x);
                diff_x_diff = diff(x_diff);
                for n_comma = 1:length(diff_x_diff)
                    if diff_x_diff(n_comma) > 1
                        n(si,n_comma) = str2double(numbers_string(x_diff(1)));
                        x_diff(1) = [];
                    else
                        n(si,n_comma) = str2double(numbers_string(x_diff(1:2)));
                        x_diff(1:2) = [];
                    end
                end
            end
            si = si+1;
        end
        
        % now n defines the order of the functions and feval will start
        % each function in the right order
        t = 1;
        locFile = [];
        EEG = [];
        while 1
            if isempty(find(n==t))
                break
            end
            [c,~] = find(n==t);
            [EEG,locFile] = feval(fldnames3{c(1)+idx_src1},data_struct,fldnames1(sbj),EEG,locFile);
            t = t+1;
        end
    end
end

end

