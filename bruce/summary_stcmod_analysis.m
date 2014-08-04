clear all
% close all

data_dir_base = '~/Analysis/bruce';
Expt_nums = [81 85 86 87 88 89 91 92 93];
% Expt_nums = [85 86];

bad_probes = 16;
%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat
mahal_thresh = su_data.mah_thresh;

all_sua_data = [];
all_su_probes = [];
all_su_exnums = [];
for ex = 1:length(Expt_nums)
    fprintf('Loading data from expt %d\n',Expt_nums(ex));
    fname = [data_dir_base sprintf('/G0%d',Expt_nums(ex)) '/sac_mod/full_bar_stcmod_data2.mat'];
    load(fname);
            
    all_sua_data = cat(1,all_sua_data,sua_data');
    all_su_probes = cat(1,all_su_probes,su_probes');
    all_su_exnums = cat(1,all_su_exnums,Expt_nums(ex)*ones(length(sua_data),1));
        
end
% lags = anal_params.lags*anal_params.dt;
dt = anal_params.dt;
sac_bin_cents = anal_params.sac_bin_cents;

%%
close all
n_sus = length(all_sua_data);
use_sua_hor = get_struct_data(all_sua_data,1:n_sus,'stc_hor_use');
use_sua_ver = get_struct_data(all_sua_data,1:n_sus,'stc_ver_use');
xl = [-0.2 0.5];
f1 = figure();
f2 = figure();
f3 = figure();
for ss = 1:n_sus
% for ss = [23 33 46]
if all_su_exnums(ss) == 81
    nPix = 17;
else
    nPix = 24;
end
    fprintf('Expt %d su probe %d\n',all_su_exnums(ss),all_su_probes(ss));
    if use_sua_hor(ss)
        figure(f1);
        subplot(6,3,1)
        cur_sta = all_sua_data(ss).sta_hor; ca = max(abs(cur_sta));
        imagesc(reshape(cur_sta,12,nPix)); caxis([-ca ca]);
        cur_stcs = all_sua_data(ss).stcs_hor; ca = max(abs(cur_stcs(:)));
        for ii = 1:3
        subplot(6,3,3+ii)
        imagesc(reshape(cur_stcs(:,ii),12,nPix)); caxis([-ca ca]);
        end
        for ii = 1:3
        subplot(6,3,6+ii)
        imagesc(reshape(cur_stcs(:,end-3+ii),12,nPix)); caxis([-ca ca]);
        end
        
        figure(f2);
        subplot(2,4,1);hold on
        plot(sac_bin_cents,all_sua_data(ss).gsac_ta_hor_sm);xlim(xl);
        title(sprintf('Avg hor rate: %.2f',all_sua_data(ss).avg_rate_hor/dt));
        subplot(2,4,2); hold on
        plot(sac_bin_cents,all_sua_data(ss).gsac_tkern_hor);xlim(xl);
        subplot(2,4,3); hold on
        plot(sac_bin_cents,all_sua_data(ss).gsac_skern_hor);xlim(xl);
        subplot(2,4,4); hold on
        plot(sac_bin_cents,all_sua_data(ss).gsacdep_info_hor);xlim(xl);hold on
        line(xl,all_sua_data(ss).ov_info_hor([1 1]),'color','k','linestyle','--');
        title(sprintf('Nspikes: %d  NSacs: %d',all_sua_data(ss).nspks_hor,all_sua_data(ss).nsacs_hor(1)));
        
        figure(f3);
        subplot(2,4,1);hold on
        plot(sac_bin_cents,all_sua_data(ss).msac_ta_hor_sm,'r');xlim(xl);
        title(sprintf('Avg hor rate: %.2f',all_sua_data(ss).avg_rate_hor/dt));
        subplot(2,4,2); hold on
        plot(sac_bin_cents,all_sua_data(ss).msac_tkern_hor,'r');xlim(xl);
        subplot(2,4,3); hold on
        plot(sac_bin_cents,all_sua_data(ss).msac_skern_hor,'r');xlim(xl);
        subplot(2,4,4); hold on
        plot(sac_bin_cents,all_sua_data(ss).msacdep_info_hor,'r');xlim(xl);
        line(xl,all_sua_data(ss).ov_info_hor([1 1]),'color','k','linestyle','--');
    end
    if use_sua_ver(ss)
        figure(f1);
        subplot(6,3,9+1)
        cur_sta = all_sua_data(ss).sta_ver; ca = max(abs(cur_sta));
        imagesc(reshape(cur_sta,12,nPix)); caxis([-ca ca]);
        cur_stcs = all_sua_data(ss).stcs_ver; ca = max(abs(cur_stcs(:)));
        for ii = 1:3
        subplot(6,3,9+3+ii)
        imagesc(reshape(cur_stcs(:,ii),12,nPix)); caxis([-ca ca]);
        end
        for ii = 1:3
        subplot(6,3,9+6+ii)
        imagesc(reshape(cur_stcs(:,end-3+ii),12,nPix)); caxis([-ca ca]);
        end
        
        figure(f2);
        subplot(2,4,5);hold on
        plot(sac_bin_cents,all_sua_data(ss).gsac_ta_ver_sm);xlim(xl);
        title(sprintf('Avg ver rate: %.2f',all_sua_data(ss).avg_rate_ver/dt));
        subplot(2,4,6);hold on
        plot(sac_bin_cents,all_sua_data(ss).gsac_tkern_ver);xlim(xl);
        subplot(2,4,7);hold on
        plot(sac_bin_cents,all_sua_data(ss).gsac_skern_ver);xlim(xl);
        subplot(2,4,8);hold on
        plot(sac_bin_cents,all_sua_data(ss).gsacdep_info_ver);xlim(xl);
        line(xl,all_sua_data(ss).ov_info_ver([1 1]),'color','k','linestyle','--');
        title(sprintf('Nspikes: %d  NSacs: %d',all_sua_data(ss).nspks_ver,all_sua_data(ss).nsacs_ver(1)));

        figure(f3);
        subplot(2,4,5);hold on
        plot(sac_bin_cents,all_sua_data(ss).msac_ta_ver_sm,'r');xlim(xl);
        title(sprintf('Avg ver rate: %.2f',all_sua_data(ss).avg_rate_ver/dt));
        subplot(2,4,6);hold on
        plot(sac_bin_cents,all_sua_data(ss).msac_tkern_ver,'r');xlim(xl);
        subplot(2,4,7);hold on
        plot(sac_bin_cents,all_sua_data(ss).msac_skern_ver,'r');xlim(xl);
        subplot(2,4,8);hold on
        plot(sac_bin_cents,all_sua_data(ss).msacdep_info_ver,'r');xlim(xl);
        line(xl,all_sua_data(ss).ov_info_ver([1 1]),'color','k','linestyle','--');
    end
     
    pause
    figure(f1);clf;
    figure(f2);clf;
    figure(f3);clf;
end