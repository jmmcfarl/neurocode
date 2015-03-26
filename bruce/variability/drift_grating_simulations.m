clear all

global Expt_name monk_name bar_ori rec_type

Expt_name = 'G086';
bar_ori = 0; %bar orientation to use (only for UA recs)
monk_name = 'jbe';

mod_data_name = 'corrected_models2';
compact_data_name = 'packaged_data';
et_anal_name = 'full_eyetrack_Rinit';

et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
mod_data_dir = ['~/Analysis/bruce/' Expt_name '/models'];

%%
data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    error('Couldnt find data directory');
end
Edata_file = strcat(data_dir,'/',monk_name,Expt_name,'Expts');
fprintf('Loading %s\n',Edata_file);
load(Edata_file);

%is this a laminar probe or utah array rec?
if strcmp(Expts{1}.Header.DataType,'GridData 96')
    rec_type = 'UA';
elseif strcmp(Expts{1}.Header.DataType,'Spike2')
    rec_type = 'LP';
end

compact_file = [data_dir sprintf('/%s_ori%d',compact_data_name,bar_ori)];
fprintf('Loading %s\n',compact_file);
load(compact_file);
%%
%if using coil info
if any(params.use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end
et_file = [et_dir sprintf('/%s_ori%d',et_anal_name,bar_ori)];
fprintf('Loading %s\n',et_file);
load(et_file);

mod_file = [mod_data_dir sprintf('/%s_ori%d',mod_data_name,bar_ori)];
fprintf('Loading %s\n',mod_file);
load(mod_file);

%%

