function [param] = readcsv(csvname)
%READCSV This function reads a csvfile 
%and puts its data into a struct, dynamically
% 'C:\Users\Promotion\Matlab_Scripts\Consciousness\EEG_processing\UiO_eeganalysis\csvs\'
%%Read file
fileName=[csvname '.csv']; %Specifies the csvs folder as storage
fid = fopen(fileName,'r');   
  lineArray = cell(200,1);     % Preallocate a cell matrix
  lineIndex = 1;               % Iterator, cause while loops are quicker
  nextLine = fgetl(fid);       % Get line
while ~isequal(nextLine,-1)         
    lineArray{lineIndex} = nextLine;  % Add the line to the cell matrix
    lineIndex = lineIndex+1;          
    nextLine = fgetl(fid);            %Get line
end
fclose(fid);        

%%Format the read file  
lineArray = lineArray(1:lineIndex-1);  % Remove empty cells, but this doesn't work properly
for iLine = 1:lineIndex-1              
    lineData = textscan(lineArray{iLine},'%s',...  %# Read strings
                        'Delimiter',';');
    lineData = lineData{1};              % Remove cell encapsulation
    if strcmp(lineArray{iLine}(end),';')  % Account for when the line ends with a delimiter
      lineData{end+1} = '';                     
    end
    lineArray(iLine,1:numel(lineData)) = lineData;  % Overwrite line data
end
la=lineArray;

%%Create a struct out of read data
for i=2:length(la(:,1))    % 1st row is sacred, and there shouldn't be more than 200 csv lines. 
     for x=2:length(la(1,:)) % 1st col is sacred, and there shouldn't be more than 50 columns (48 participants). Can be increased later
         if isempty(la{1,x})==0 && isempty(la{2,x})==0 && isempty(la{i,2}) == 0 && isempty(la{i,x}) == 0 % If col and row header, and cell, is not empty
            param.(la{1,x}).(la{2,x}).(la{i,1})=la{i,x}; % And now i realise this is wrong and should crash. Old version.
         end
     end
end
  

end

