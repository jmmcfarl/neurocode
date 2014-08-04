function filenames = SearchFileNames(directory, keywords)
% Search for files containing 'word' in a given directory
% directory: string
% keywords
% Yuwei Cui, Mar 7 2012


datafiles = dir(directory);

filenames = {};
for i=1:length(datafiles)
    if isempty(strfind(datafiles(i).name, keywords)) 
        
        continue;
    else
        if ~exist('filenames')
            filenames = datafiles(i).name;
        else
            filenames = [filenames; datafiles(i).name;];
        end
    end
end
fprintf(' %d files have been found\n',length(filenames));