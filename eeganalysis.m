function [] = eeganalysis( csvname )
%EEGANALYSIS Main analysis script. Calls stuff
%and reads stuff.
%By Andr� Sevenius Nilsen & Benjamin Th�rer
% sevenius.nilsen@gmail.com
% benjamin.thuerer@kit.edu

%%Calling the read csv file function
param=readcsv(csvname);

fldnames1 = fieldnames(param);

for sbj = 1:size(fldnames1,1)
    for ses = 1:size(fieldnames(eval(['param.' fldnames1{sbj}])),1)
        fldnames2 = fieldnames(eval(['param.' fldnames1{sbj}]));
        fldnames3 = fieldnames(eval(['param.' fldnames1{sbj} '.' fldnames2{ses}]));
        idx_src1 = find(strcmp(fldnames3,'start_of_function_list'));
        idx_src2 = find(strcmp(fldnames3,'end_of_function_list'));
        
        data_struct = eval(['param.' fldnames1{sbj} '.' fldnames2{ses}]);
        
        si_edge = idx_src2-idx_src1-1;
        si = 1;
        n = zeros(si_edge,1);
        while si < si_edge+1
            n(si) = str2num(eval(['param.' fldnames1{sbj} '.' fldnames2{ses} '.' fldnames3{si+idx_src1}]));
            si = si+1;
        end
        
        t = 1;
        locFile = [];
        EEG = [];
        while 1
            if isempty(find(n==t))
                break
            end
            [EEG,locFile] = feval(fldnames3{find(n==t)+idx_src1},data_struct,fldnames1(sbj),EEG,locFile);
            t = t+1;
        end
    end
end

end

