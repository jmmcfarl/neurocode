clear all
close all

Expt_name = 'M005';
monName = 'jbe';
data_dir = '/media/NTlab_data3/Data/bruce/';

cd([data_dir Expt_name]);

matfiles = what; matfiles = matfiles.mat;
[startinds,~,token,expt_chunk,tokenname] = regexp(matfiles,sprintf('%s%s.([0-9]{1,3}).mat',monName,Expt_name));
is_match = cellfun(@(x) ~isempty(x),startinds);
expt_nums = cellfun(@(x) str2num(cell2mat((x{1}))),tokenname(is_match));

fprintf('Found %d expts ranging from %d to %d\n',length(unique(expt_nums)),min(expt_nums),max(expt_nums));
%%



block_nums = [3 4 5];
force_rls_process = true;
rls_list = 1:12;
force_rls_align = true;

n_probes = 24;

Expts = {};
for bb = 1:length(block_nums)
    dat_name = [pwd sprintf('/%s%s.%d.mat',monName,exp_name,block_nums(bb))];
    [a,cur_Expt] = APlaySpkFile(dat_name,'nospikes','noerrs');
    Expts = {Expts{:} cur_Expt{:}};
end

if ~exist([pwd '/stim_data.mat'],'file') || force_rls_process
    process_rls_files(pwd,rls_list);
else
    fprintf('Load stim_data.mat\n');
    load('./stim_data.mat');
end

if ~exist([pwd '/expt_data.mat'],'file') || force_rls_align
    align_rls_data(pwd,Expts,block_nums);
end
fprintf('Loading expt_data.mat\n');
load('./expt_data.mat');

