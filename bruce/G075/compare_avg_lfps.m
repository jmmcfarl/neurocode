clear all
cd ~/Data/bruce/G075
%%
load ./fortyhz_lfp_avgs
forty_avg = avg_lfpsf;
forty_sem = sem_lfpsf;
forty_phase_locking_time = phase_locking_time;

load ./sim_sac_lfp_avgs.mat
sac_avg = avg_lfpsf;
sac_sem = sem_lfpsf;
sac_phase_locking_time = phase_locking_time;

load ./loc_glob_lfp_avgs.mat
loc_avg = avg_lfpsf;
loc_sem = sem_lfpsf;
loc_phase_locking_time = phase_locking_time;

%%
stim_fs = 1e4/117.5;
stim_times = (0:40:320)/stim_fs;
cur_el = 2;
forty_stim = 4;
forty_stim2 = 5;
loc_stim = 10;
loc_stim2 = 11;
figure
shadedErrorBar(lags/Fsd,squeeze(forty_avg(forty_stim,cur_el,:)),squeeze(forty_sem(forty_stim,cur_el,:)),{'color','b'});
hold on
shadedErrorBar(lags/Fsd,squeeze(forty_avg(forty_stim2,cur_el,:)),squeeze(forty_sem(forty_stim2,cur_el,:)),{'color','g'});
shadedErrorBar(lags/Fsd,squeeze(loc_avg(loc_stim,cur_el,:)),squeeze(loc_sem(loc_stim,cur_el,:)),{'color','r'});
shadedErrorBar(lags/Fsd,squeeze(loc_avg(loc_stim2,cur_el,:)),squeeze(loc_sem(loc_stim2,cur_el,:)),{'color','k'});

stim_fs = 1e4/117.5;

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
ylim(yl)

%%
stim_fs = 1e4/117.5;
stim_times = (0:40:320)/stim_fs;
cur_el = 2;
forty_stim = 2;
loc_stim = 8;
sac_stim = 8;
figure
shadedErrorBar(lags/Fsd,squeeze(forty_avg(forty_stim,cur_el,:)),squeeze(forty_sem(forty_stim,cur_el,:)),{'color','b'});
hold on
shadedErrorBar(lags/Fsd,squeeze(loc_avg(loc_stim,cur_el,:)),squeeze(loc_sem(loc_stim,cur_el,:)),{'color','r'});
shadedErrorBar(lags/Fsd,squeeze(sac_avg(sac_stim,cur_el,:)),squeeze(sac_sem(sac_stim,cur_el,:)),{'color','k'});
temp_forty_avg = squeeze(forty_avg(forty_stim,cur_el,:));
% temp_sac_avg = squeeze(sac_avg(sac_stim,cur_el,:));
[tempb,tempa] = butter(2,15/(Fsd/2),'low');
temp_forty_avg = filtfilt(tempb,tempa,temp_forty_avg);
% plot(lags/Fsd,temp_forty_avg,'g','linewidth',2)
% plot(lags/Fsd,temp_forty_avg+temp_sac_avg,'c','linewidth',2)
stim_fs = 1e4/117.5;

yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
ylim(yl)

%%
stim_fs = 1e4/117.5;
stim_times = (0:40:320)/stim_fs;
cur_el = 2;
forty_stim = 1;
loc_stim = 7;
sac_stim = 7;
cl = [0.1 0.7];

figure
subplot(3,1,1)
pcolor(lags/Fsd,wfreqs,squeeze(forty_phase_locking_time(forty_stim,:,:))');shading interp
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
ylim(yl);
caxis(cl);
set(gca,'yscale','log');

subplot(3,1,2)
pcolor(lags/Fsd,wfreqs,squeeze(sac_phase_locking_time(sac_stim,:,:))');shading interp
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
ylim(yl);
caxis(cl);
set(gca,'yscale','log');

subplot(3,1,3)
pcolor(lags/Fsd,wfreqs,squeeze(loc_phase_locking_time(loc_stim,:,:))');shading interp
yl = ylim();
for i = 1:length(stim_times)
    line(stim_times([i i]),yl,'color','k')
end
ylim(yl)
caxis(cl);
set(gca,'yscale','log');
