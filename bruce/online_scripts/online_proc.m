clear all

data_dir = '~/James_scripts/bruce/online_processing/';
cd(data_dir);

addpath('~/James_scripts/bruce/bruce_code/AllV_backup/AllV_old/');
addpath('~/James_scripts/bruce/processing/');
addpath(genpath('~/James_scripts/iCSD/'))
monName = 'lem';
exp_name = 'M266';
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

