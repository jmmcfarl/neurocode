clear all
% close all

data_dir_base = '~/Analysis/bruce';
Expt_nums = [86 87 88 89 91 93];
%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat
mahal_thresh = su_data.mah_thresh;

bin_mua_data = [];
bin_mu_probes = [];
bin_mu_exnums = [];
bin_sua_data = [];
bin_su_probes = [];
bin_su_exnums = [];

all_mua_data = [];
all_mu_probes = [];
all_mu_exnums = [];
all_sua_data = [];
all_su_probes = [];
all_su_exnums = [];
for ex = 1:length(Expt_nums)
    fprintf('Loading data from expt %d\n',Expt_nums(ex));
    fname = [data_dir_base sprintf('/G0%d',Expt_nums(ex)) '/sac_mod/full_binoc_sacmod_data.mat'];
    load(fname);
    
    bin_mua_data = cat(1,bin_mua_data,mua_data');
    bin_mu_probes = cat(1,bin_mu_probes,(1:96)');
    bin_mu_exnums = cat(1,bin_mu_exnums,Expt_nums(ex)*ones(96,1));

    bin_sua_data = cat(1,bin_sua_data,sua_data');
    bin_su_probes = cat(1,bin_su_probes,su_probes');
    bin_su_exnums = cat(1,bin_su_exnums,Expt_nums(ex)*ones(length(su_probes),1));
   
    fname = [data_dir_base sprintf('/G0%d',Expt_nums(ex)) '/sac_mod/full_sacmod_data.mat'];
    load(fname);
    
    all_mua_data = cat(1,all_mua_data,mua_data');
    all_mu_probes = cat(1,all_mu_probes,(1:96)');
    all_mu_exnums = cat(1,all_mu_exnums,Expt_nums(ex)*ones(96,1));

    all_sua_data = cat(1,all_sua_data,sua_data');
    all_su_probes = cat(1,all_su_probes,su_probes');
    all_su_exnums = cat(1,all_su_exnums,Expt_nums(ex)*ones(length(su_probes),1));

end
lags = anal_params.lags*anal_params.dt;
dt = anal_params.dt;

%% view all expt-avg mu sacmod
bad_probes = 16;
close all
f1 = figure();
f2 = figure();
for ex = 1:length(Expt_nums)
    fprintf('Expt %d\n',Expt_nums(ex));
    cur_mu_set = find(bin_mu_exnums == Expt_nums(ex));
    cur_su_set = find(bin_su_exnums == Expt_nums(ex));
    if length(cur_mu_set) ~= 96
        error('Wrong mu Num');
    end
    exclude = [bad_probes; bin_su_probes(cur_su_set)];
    cur_mu_set(exclude) = [];
    
    mu_dense_gsac_avgs = get_struct_data(bin_mua_data,cur_mu_set,'dense_binoc_gsac_avg');
    mu_sparse_gsac_avgs = get_struct_data(bin_mua_data,cur_mu_set,'sparse_binoc_gsac_avg');
    mu_dense_msac_avgs = get_struct_data(bin_mua_data,cur_mu_set,'dense_binoc_msac_avg');
    mu_sparse_msac_avgs = get_struct_data(bin_mua_data,cur_mu_set,'sparse_binoc_msac_avg');
    
    mu_gsac = get_struct_data(all_mua_data,cur_mu_set,'gsac_avg');
    mu_msac = get_struct_data(all_mua_data,cur_mu_set,'msac_avg');
    
    figure(f1);hold on; grid on
    shadedErrorBar(lags,nanmean(mu_dense_gsac_avgs),nanstd(mu_dense_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','r'});
    shadedErrorBar(lags,nanmean(mu_sparse_gsac_avgs),nanstd(mu_sparse_gsac_avgs)/sqrt(length(cur_mu_set)),{'color','b'});
    shadedErrorBar(lags,nanmean(mu_gsac),nanstd(mu_gsac)/sqrt(length(cur_mu_set)),{'color','k'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    ylim([0.7 1.4]);
    
    figure(f2);hold on; grid on
    shadedErrorBar(lags,nanmean(mu_dense_msac_avgs),nanstd(mu_dense_msac_avgs)/sqrt(length(cur_mu_set)),{'color','r'});
    shadedErrorBar(lags,nanmean(mu_sparse_msac_avgs),nanstd(mu_sparse_msac_avgs)/sqrt(length(cur_mu_set)),{'color','b'});
    shadedErrorBar(lags,nanmean(mu_msac),nanstd(mu_msac)/sqrt(length(cur_mu_set)),{'color','k'});
    xlim([-0.2 0.5]);
    xl = xlim(); yl = ylim();
    line(xl,[1 1],'color','k','linestyle','--');
    line([0 0],yl,'color','k','linestyle','--');
    ylim([0.7 1.4]);

    pause
    figure(f1); clf;
    figure(f2); clf;
end

