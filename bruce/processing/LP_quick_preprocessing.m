clear all
close all

addpath(genpath('~/Bruce_matlab'));

Expt_name = 'M013';
monName = 'jbe';
data_dir = '/media/NTlab_data3/Data/bruce/';
stim_dir = strcat(data_dir,Expt_name,'/stims');

cd([data_dir Expt_name]);

matfiles = what; matfiles = matfiles.mat;
[startinds,~,token,expt_chunk,tokenname] = regexp(matfiles,sprintf('%s%s.([0-9]{1,3}).mat',monName,Expt_name));
is_match = cellfun(@(x) ~isempty(x),startinds);
block_nums = cellfun(@(x) str2num(cell2mat((x{1}))),tokenname(is_match));
block_nums = sort(block_nums);

fprintf('Found %d blocks ranging from %d to %d\n',length(unique(block_nums)),min(block_nums),max(block_nums));

%%
ExptFileName = strcat(monName,Expt_name,'Expts.mat');
if ~exist(ExptFileName);
    [Expts,Ex] = ReadExptDir_jmm([data_dir Expt_name],monName);
    fprintf('Saving file %s\n',ExptFileName);
    save(ExptFileName,Expts);
else
   fprintf('Loading existing file %s\n',ExptFileName);
   load(ExptFileName);
end

%%
force_rls_process = true;

cd(stim_dir);

cur_files = dir('*.rc*'); 
cur_file_names = arrayfun(@(x) x.name,cur_files,'uniformoutput',0);
[startinds,~,token,expt_chunk,tokenname] = regexp(cur_file_names,'rls.rc([0-9]{1,3})');
rls_list = sort(cellfun(@(x) str2num(cell2mat((x{1}))),tokenname));

disparity_flag = false;
process_rls_files(stim_dir,rls_list,disparity_flag);

%%
cd(stim_dir);
force_rls_align = true;

align_rls_data(stim_dir,Expts,block_nums);

%%
rmpath(genpath('~/Bruce_matlab'));

%%


