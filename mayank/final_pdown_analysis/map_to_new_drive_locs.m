function new_dir = map_to_new_drive_locs(old_dir,is_file,force_second)
%% helper function for mapping old data locs to locations spread out over two external HDs)

if nargin < 2 || isempty(is_file)
    is_file = false;
end
if nargin < 3
    force_second = false;
end

%first try to find the directory on this drive
new_dir = strrep(old_dir,'C:','/Volumes/Expansion Drive/James_laptop_backup');
new_dir = strrep(new_dir,'\','/');

if is_file
    if ~strcat(new_dir(end-3:end),'.mat')
    new_file = [new_dir '.mat'];
    else
       new_file = new_dir; 
    end
    if ~exist(new_file,'file') || force_second
        new_dir = strrep(old_dir,'C:','/Volumes/anne is awesome');
        new_dir = strrep(new_dir,'\','/');
        if ~strcat(new_dir(end-3:end),'.mat')
            new_file = [new_dir '.mat'];
        else
            new_file = new_dir;
        end
    end
    new_dir = new_file;
else
    %for a few of them, have to look on this other drive
    if ~exist(new_dir,'dir') || force_second
        new_dir = strrep(old_dir,'C:','/Volumes/anne is awesome');
        new_dir = strrep(new_dir,'\','/');
    end
    
    if ~exist(new_dir,'dir')
        fprintf('Couldnt find %s!\n',new_dir)
    end
end