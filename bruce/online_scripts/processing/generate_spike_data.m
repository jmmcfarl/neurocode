close all
clear all
Expt_name = 'M270';

dir_prefix = '~';
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
cd(data_dir)
%%

FullV_data = dir('*FullV.mat');
n_expts = length(FullV_data);
for ee = 1:n_expts
    fprintf('Running AllVPcs for file %d of %d\n',ee,n_expts);
    AllVPcs([data_dir '/' FullV_data(ee).name],'tchan',1:24,'reapply','savespikes','noninteract');
end