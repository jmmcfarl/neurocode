
clear all
fig_dir = '/home/james/Analysis/bruce/ET_final/';
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
% Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
bar_ori = 0;

if bar_ori == 0
    ov_SU_data = load('~/Analysis/bruce/ET_final/SUA_sum_data_hor.mat');
    full_expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
else
    full_expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
    ov_SU_data = load('~/Analysis/bruce/ET_final/SUA_sum_data_ver.mat');
end
ov_SU_data = ov_SU_data.Expt_sua_data;

all_su_data = [];
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee}
    cur_expt_num = str2num(Expt_name(2:end));
    data_dir = ['~/Data/bruce/' Expt_name];
    anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
    cd(anal_dir);
    
    if bar_ori == 0
        %for horizontal
        mod_data_name = 'monoc_eyecorr_hbar_mods.mat';
        anal_name = 'monoc_eyecorr_hbar2.mat';
        mod_data_name_EC = 'monoc_eyecorr_hbar_mods_Cinit.mat';
        mod_data_name_ECproc = 'monoc_eyecorr_hbar_mods_CPinit.mat';
    elseif bar_ori == 90
        %for vertical
        mod_data_name = 'monoc_eyecorr_vbar_mods.mat';
        anal_name = 'monoc_eyecorr_vbar2.mat';
        mod_data_name_EC = 'monoc_eyecorr_vbar_mods_Cinit.mat';
        mod_data_name_ECproc = 'monoc_eyecorr_vbar_mods_CPinit.mat';
    end
    
    data = load(anal_name,'dit_R2_LOO','it_R2','it_R2_LOO','et_tr_set');
    mod_data = load(mod_data_name,'all_mod_SU*');
    mod_data_EC = load(mod_data_name_EC,'all_mod_R2');
    mod_data_ECproc = load(mod_data_name_ECproc,'all_mod_R2');
    
    cur_SU_sumdata = ov_SU_data{cellfun(@(x) x.Expt_num,ov_SU_data) == cur_expt_num};
    cur_all_su_nums = [cur_SU_sumdata.sua_data(:).su_num];
    
    clear cur_su_data
    su_numbers = mod_data.all_mod_SUnum(data.et_tr_set(mod_data.all_mod_SUnum(data.et_tr_set) > 0));
    su_inds = find(ismember(mod_data.all_mod_SUnum,su_numbers));
    for ii = 1:length(su_numbers)
        cur_imp = data.dit_R2_LOO(ii,end,su_inds(ii));
        
        cur_su_data(ii).usable = true;
        cur_su_data(ii).before_R2 = data.it_R2(1,su_inds(ii));
        cur_su_data(ii).after_R2 = cur_imp;
        cur_su_data(ii).R2_EC = mod_data_EC.all_mod_R2(su_inds(ii));
        cur_su_data(ii).R2_ECproc = mod_data_ECproc.all_mod_R2(su_inds(ii));
        cur_su_data(ii).probe_num = mod_data.all_mod_SU(su_inds(ii));
        cur_su_data(ii).expt_num = Expt_name;
        
        cur_spkloc = find(cur_all_su_nums == su_numbers(ii));
        cur_su_data(ii).num_blocks = cur_SU_sumdata.sua_data(cur_spkloc).num_blocks;
        cur_su_data(ii).n_spikes = cur_SU_sumdata.sua_data(cur_spkloc).n_spikes;
        cur_su_data(ii).mean_rate = cur_SU_sumdata.sua_data(cur_spkloc).mean_rate;
    end
    
    all_su_data = cat(2,all_su_data,cur_su_data);
end
%%
min_num_blocks = 3;
num_blocks = [all_su_data(:).num_blocks];
n_spikes = [all_su_data(:).n_spikes];

used_sus = find(num_blocks >= min_num_blocks);

%%
all_R2 = [[all_su_data(used_sus).before_R2]' [all_su_data(used_sus).R2_EC]' [all_su_data(used_sus).R2_ECproc]' [all_su_data(used_sus).after_R2]'];
all_R2 = bsxfun(@rdivide,all_R2(:,2:end),all_R2(:,1));

h = figure;
boxplot(all_R2,'labels',{'Coil corrected','Processed coil','Inferred'});
xl = xlim();
line(xl,[1 1],'color','k','linestyle','--');
ylabel('Relative R2');
ylim([0 4])
%%
fig_width = 3.27; rel_height = 0.8;
figufy(h);
fname = [fig_dir 'hor_modR2_init_compare_SU.pdf'];
exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h);

