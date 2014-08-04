clear all
cd ~/Data/bruce/G075
%%
use_lfps = [1:8:96];
load ./fortyhz_lfp_avgs_v2.mat
forty_all_phase = all_block_stim_phase;
forty_spk_avg = avg_smspks;
forty_spk_sem = sem_smspks;
forty_spk_std = std_smspks;
% forty_spk_avg = avg_spks;
% forty_spk_sem = sem_spks;
% forty_spk_std = std_spks;
forty_lfp_avg = avg_lfpsf;
forty_lfp_sem = sem_lfpsf;
forty_lfp_std = std_lfpsf;



load ./loc_glob_lfp_avgs.mat
movie_ids = [1e3 2e3 4e3 5e3 6e3 7e3 8e3 9e3 10e3 11e3 12e3 13e3];
loc_spk_avg = avg_smspk;
loc_spk_sem = sem_smspk;
loc_spk_std = std_smspk;
% for ss =1:length(movie_ids)
%     n_trials(ss) = sum(~isnan(all_block_stim_spks(ss,:,1,1)),2);
%     avg_spk(ss,:,:) = squeeze(nanmean(all_block_stim_spks(ss,:,:,:)));
%     sem_spk(ss,:,:) = squeeze(nanstd(all_block_stim_spks(ss,:,:,:)))/sqrt(n_trials(ss));
%     std_spk(ss,:,:) = squeeze(nanstd(all_block_stim_spks(ss,:,:,:)));
% end
% loc_spk_avg = avg_spk;
% loc_spk_sem = sem_spk;
% loc_spk_std = std_spk;
loc_lfp_avg = avg_lfpsf;
loc_lfp_sem = sem_lfpsf;
loc_lfp_std = std_lfpsf;

%%
stim_fs = 1e4/117.5;
stim_times = (0:40:320)/stim_fs;
for ii = 1:96
    ii
    cur_unit = ii;
    cur_el = 1;
    forty_stim = 4;
    % forty_stim2 = 5;
    loc_stim = 10;
    % loc_stim2 = 11;
    figure(1)
    shadedErrorBar(lags/Fsd,squeeze(forty_spk_avg(forty_stim,cur_unit,:)),squeeze(forty_spk_sem(forty_stim,cur_unit,:)),{'color','b'});
    hold on
    shadedErrorBar(lags/Fsd,squeeze(loc_spk_avg(loc_stim,cur_unit,:)),squeeze(loc_spk_sem(loc_stim,cur_unit,:)),{'color','r'});
    
    stim_fs = 1e4/117.5;
    set(gcf,'position',[00 1000 1300 500])
    
    yl = ylim();
    for i = 1:length(stim_times)
        line(stim_times([i i]),yl,'color','k')
    end
    ylim(yl)
    xlim([1 2.5])
    
    
    % loc_stim2 = 11;
    figure(2)
    shadedErrorBar(lags/Fsd,squeeze(forty_lfp_avg(forty_stim,cur_el,:)),squeeze(forty_lfp_sem(forty_stim,cur_el,:)),{'color','b'});
    hold on
    shadedErrorBar(lags/Fsd,squeeze(loc_lfp_avg(loc_stim,cur_el,:)),squeeze(loc_lfp_sem(loc_stim,cur_el,:)),{'color','r'});
    temp_forty_avg = squeeze(forty_lfp_avg(forty_stim,cur_el,:));
    [tempb,tempa] = butter(2,15/(Fsd/2),'low');
    temp_forty_avg = filtfilt(tempb,tempa,temp_forty_avg);
    plot(lags/Fsd,temp_forty_avg,'g','linewidth',2)
    set(gcf,'position',[00 00 1300 500])
    
    xlim([1 2.5])
    yl = ylim();
    for i = 1:length(stim_times)
        line(stim_times([i i]),yl,'color','k')
    end
    ylim(yl)
    
    pause
    figure(1);clf;figure(2);clf;
end

%%
load ./noise_40_phase_mod_data


load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
nearest_lfps = nan(length(use_lfps),1);
for ll = 1:96
   all_dists = sqrt((X_pos-X_pos(ll)).^2 + (Y_pos-Y_pos(ll)).^2); 
   all_dists(ll) = inf;
   [~,best_loc] = min(all_dists(use_lfps));
   nearest_lfps(ll) = best_loc;
end


stim_fs = 1e4/117.5;
stim_times = (0:40:320)/stim_fs;
cur_unit = 19;
cur_el = nearest_lfps(cur_unit);
% cur_el2 = 10;
forty_stim = 6;
% forty_stim2 = 5;
loc_stim = 12;
% loc_stim2 = 11;


forty_avg_phase = squeeze(forty_all_phase{forty_stim}(:,cur_el,:,:));
n_trials = size(forty_avg_phase,1);
forty_phase_mod_out = nan(n_trials,length(lags));
for i = 1:n_trials
    cur_alpha_phase = squeeze(forty_avg_phase(i,:,:));
    Pmat = nan(length(lags),length(wfreqs)*nbins);
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat(:,cur_set) = cur_tb;
    end
    forty_phase_mod_out(i,:) = Pmat*stim_ind_phase_pfilt(cur_unit,:)';
end

loc_avg_phase = squeeze(all_block_stim_phase(loc_stim,:,cur_el,:,:));
n_trials = size(loc_avg_phase,1);
loc_phase_mod_out = nan(n_trials,length(lags));
for i = 1:n_trials
    cur_alpha_phase = squeeze(loc_avg_phase(i,:,:));
    Pmat = nan(length(lags),length(wfreqs)*nbins);
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat(:,cur_set) = cur_tb;
    end
    loc_phase_mod_out(i,:) = Pmat*stim_ind_phase_pfilt(cur_unit,:)';
end

%%
figure(1)
shadedErrorBar(lags/Fsd,squeeze(forty_spk_avg(forty_stim,cur_unit,:)),squeeze(forty_spk_sem(forty_stim,cur_unit,:)),{'color','b'});
hold on
shadedErrorBar(lags/Fsd,squeeze(loc_spk_avg(loc_stim,cur_unit,:)),squeeze(loc_spk_sem(loc_stim,cur_unit,:)),{'color','r'});
stim_fs = 1e4/117.5;
% set(gcf,'position',[00 1000 1300 500])
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
ylim(yl)
% xlim([1 2.5])



% z1 = (squeeze(forty_spk_avg(forty_stim,cur_unit,:)));
% z2 = (squeeze(loc_spk_avg(loc_stim,cur_unit,:)));
% ov_avg = 0.5*mean(z1)+0.5*mean(z2);
figure(2)
hold on
plot(lags/Fsd,zscore(mean(forty_phase_mod_out)),'b')
plot(lags/Fsd,zscore(mean(loc_phase_mod_out)),'r')
% plot(lags/Fsd,z1/ov_avg,'b');
% plot(lags/Fsd,z2/ov_avg,'r');
% plot(lags/Fsd,(z1-z2)/ov_avg,'k');
% 
z1 = squeeze(forty_lfp_avg(forty_stim,cur_el,:));
% [tempb,tempa] = butter(2,40/(Fsd/2),'low');
% z1 = filtfilt(tempb,tempa,z1);
% z2 = squeeze(forty_lfp_avg(forty_stim,cur_el2,:));
z2 = squeeze(loc_lfp_avg(loc_stim,cur_el,:));
norm = 0.5*std(z1)+0.5*std(z2);
plot(lags/Fsd,z1/norm-3,'c','linewidth',2);
plot(lags/Fsd,z2/norm-3,'g','linewidth',2);

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
ylim(yl)


