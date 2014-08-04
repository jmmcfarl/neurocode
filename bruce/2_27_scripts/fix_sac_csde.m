clear all
% close all
addpath(genpath('~/James_scripts'));

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-5 5;-5 5];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

%%
Fs = 1000;
dsf = 4;
Fsd = Fs/dsf;
niqf = Fs/2;
% [b,a] = butter(4,[45 60]/niqf);
[b,a] = butter(2,[1 110]/niqf);
[b_hf,a_hf] = butter(2,[100]/niqf,'high');
[b_gam,a_gam] = butter(2,[30 60]/niqf);

%%
all_stim_inds = [];
all_sac_inds = [];
all_peak_inds = [];
all_fix_inds = [];
all_fix_end_inds = [];
all_lfp_data = [];
all_lfp_data_hf = [];
all_lfp_data_gam = [];
cur_length = 0;
for blockid = 1:3;
    fprintf('Block %d of %d\n',blockid,3);
    
    %%
    cd /Users/James/Data/bruce/2_27_12
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
    %get start times of each LFP trial
    n_lfp_trials = length(LFP.Trials);
    lfp_trial_start = nan(n_lfp_trials,1);
    lfp_trial_stop = nan(n_lfp_trials,1);
    for i = 1:n_lfp_trials
        lfp_trial_start(i) = LFP.Trials(i).Start/1e4; %originally in tenths of ms
        lfp_trial_stop(i) = LFP.Trials(i).End/1e4;
        lfp_dur(i) = size(LFP.Trials(i).LFP,1)/1000;
    end
    
    lfp_time = [];
    lfp_samps = [];
    for i = 1:n_lfp_trials
        lfp_time = [lfp_time linspace(lfp_trial_start(i),lfp_trial_start(i)+lfp_dur(i),size(LFP.Trials(i).LFP,1))];
        lfp_samps = [lfp_samps; LFP.Trials(i).LFP];
    end
    
    lfp_samps = filtfilt(b,a,lfp_samps);
    lfp_samps_hf = filtfilt(b_hf,a_hf,lfp_samps);
    lfp_samps_hf = abs(hilbert(lfp_samps_hf));
    lfp_sampsd_hf = downsample(lfp_samps_hf,dsf);
    lfp_samps_gam = filtfilt(b_gam,a_gam,lfp_samps);
    lfp_samps_gam = abs(hilbert(lfp_samps_gam));
    lfp_sampsd_gam = downsample(lfp_samps_gam,dsf);
    lfp_sampsd = downsample(lfp_samps,dsf);
    lfp_timed = downsample(lfp_time,dsf);
       
    cur_set = find(blockids==blockid);
    cur_fix_times = all_fix_start_times(cur_set);
    cur_sac_times = all_fix_start_times(cur_set)-all_sac_durs(cur_set);
    cur_peak_times = all_sac_peaklocs(cur_set);
    
    cur_stim_inds = round(interp1(lfp_timed,1:length(lfp_timed),Blocks{blockid}.stimtime));
    cur_fix_inds = round(interp1(lfp_timed,1:length(lfp_timed),cur_fix_times));
    cur_fix_end_inds = round(interp1(lfp_timed,1:length(lfp_timed),all_fix_stop_times(cur_set)));
    cur_sac_inds = round(interp1(lfp_timed,1:length(lfp_timed),cur_sac_times));
    cur_peak_inds = round(interp1(lfp_timed,1:length(lfp_timed),cur_peak_times));
%     cur_fix_inds(isnan(cur_fix_inds)) = [];
%     cur_sac_inds(isnan(cur_sac_inds)) = [];
    
    all_stim_inds = [all_stim_inds; cur_length + cur_stim_inds(:)];
    all_sac_inds = [all_sac_inds; cur_length + cur_sac_inds(:)];
    all_fix_inds = [all_fix_inds; cur_length + cur_fix_inds(:)];
    all_fix_end_inds = [all_fix_end_inds; cur_length + cur_fix_end_inds(:)];
    all_peak_inds = [all_peak_inds; cur_length + cur_peak_inds(:)];
%     all_sac_ids = [all_sac_ids; cur_set(2:end)];
%     all_fix_ids = [all_fix_ids; cur_set];
    all_lfp_data = [all_lfp_data; lfp_sampsd];
    all_lfp_data_hf = [all_lfp_data_hf; lfp_sampsd_hf];
    all_lfp_data_gam = [all_lfp_data_gam; lfp_sampsd_gam];
    cur_length = size(all_lfp_data,1);
    
end

%%
cur_set = find(blockids < 4);
bad = find(all_sac_durs(cur_set) > 0.1);
all_sac_inds(bad) = [];
all_fix_inds(bad) = [];
all_fix_end_inds(bad) = [];
all_peak_inds(bad) = [];
blockids(cur_set(bad)) = [];
all_sac_durs(cur_set(bad)) = [];
all_sac_amps(cur_set(bad)) = [];

all_lfp_data_hf_norm = bsxfun(@minus,all_lfp_data_hf,mean(all_lfp_data_hf));
all_lfp_data_hf_norm = bsxfun(@rdivide,all_lfp_data_hf_norm,std(all_lfp_data_hf));

all_lfp_data_gam_norm = bsxfun(@minus,all_lfp_data_gam,mean(all_lfp_data_gam));
all_lfp_data_gam_norm = bsxfun(@rdivide,all_lfp_data_gam_norm,std(all_lfp_data_gam));


%%
% used_hf = all_lfp_data_hf;
% used_gam = all_lfp_data_gam;
used_hf = all_lfp_data_hf_norm;
used_gam = all_lfp_data_gam_norm;


backlag = round(Fsd*0.1);
forlag = round(Fsd*0.5);
lags = -backlag:forlag;
% all_fix_inds(all_fix_inds < backlag | all_fix_inds > size(all_lfp_data,1)-forlag) = [];
n_fixs = length(all_fix_inds);
fix_trg_avg = zeros(24,length(lags));
fix_trg_mat = zeros(24,n_fixs,length(lags));
fix_trg_avg_hf = zeros(24,length(lags));
fix_trg_avg_gam = zeros(24,length(lags));
% fix_trg_mat_hf = zeros(24,n_fixs,length(lags));
cur_len = zeros(length(lags),1);
for i = 1:n_fixs
    cur_set = (all_fix_inds(i)-backlag):(all_fix_inds(i)+forlag);
    cur_end = all_fix_end_inds(i);
    cur_set(cur_set > cur_end) = [];
    fix_trg_mat(:,i,1:length(cur_set)) = all_lfp_data(cur_set,:)';
    fix_trg_avg(:,1:length(cur_set)) = fix_trg_avg(:,1:length(cur_set)) + all_lfp_data(cur_set,:)';
%     fix_trg_mat_hf(:,i,:) = used_hf(cur_set,:)';
    fix_trg_avg_hf(:,1:length(cur_set)) = fix_trg_avg_hf(:,1:length(cur_set)) + used_hf(cur_set,:)';
    fix_trg_avg_gam(:,1:length(cur_set)) = fix_trg_avg_gam(:,1:length(cur_set)) + used_gam(cur_set,:)';
    cur_len(1:length(cur_set)) = cur_len(1:length(cur_set)) + 1;
end
fix_trg_avg = bsxfun(@rdivide,fix_trg_avg,cur_len');
fix_trg_avg_hf = bsxfun(@rdivide,fix_trg_avg_hf,cur_len');
fix_trg_avg_gam = bsxfun(@rdivide,fix_trg_avg_gam,cur_len');

% all_sac_inds(all_sac_inds < backlag | all_sac_inds > size(all_lfp_data,1)-forlag) = [];
n_sacs = length(all_sac_inds);
sac_trg_avg = zeros(24,length(lags));
sac_trg_mat = zeros(24,n_fixs,length(lags));
sac_trg_avg_hf = zeros(24,length(lags));
sac_trg_avg_gam = zeros(24,length(lags));
% sac_trg_mat_hf = zeros(24,n_fixs,length(lags));
for i = 1:n_sacs
    cur_set = (all_sac_inds(i)-backlag):(all_sac_inds(i)+forlag);
    sac_trg_avg = sac_trg_avg + all_lfp_data(cur_set,:)';
    sac_trg_mat(:,i,:) = all_lfp_data(cur_set,:)';
    sac_trg_avg_hf = sac_trg_avg_hf + used_hf(cur_set,:)';
    sac_trg_avg_gam = sac_trg_avg_gam + used_gam(cur_set,:)';
%     sac_trg_mat_hf(:,i,:) = used_hf(cur_set,:)';
end
sac_trg_avg = sac_trg_avg/n_sacs;
sac_trg_avg_hf = sac_trg_avg_hf/n_sacs;
sac_trg_avg_gam = sac_trg_avg_gam/n_sacs;

n_sacs = length(all_peak_inds);
mid_trg_avg = zeros(24,length(lags));
mid_trg_mat = zeros(24,n_fixs,length(lags));
mid_trg_avg_hf = zeros(24,length(lags));
mid_trg_avg_gam = zeros(24,length(lags));
% mid_trg_mat_hf = zeros(24,n_fixs,length(lags));
for i = 1:n_sacs
    cur_set = (all_peak_inds(i)-backlag):(all_peak_inds(i)+forlag);
    mid_trg_avg = mid_trg_avg + all_lfp_data(cur_set,:)';
    mid_trg_mat(:,i,:) = all_lfp_data(cur_set,:)';
    mid_trg_avg_hf = mid_trg_avg_hf + used_hf(cur_set,:)';
    mid_trg_avg_gam = mid_trg_avg_gam + used_gam(cur_set,:)';
%     mid_trg_mat_hf(:,i,:) = used_hf(cur_set,:)';
end
mid_trg_avg = mid_trg_avg/n_sacs;
mid_trg_avg_hf = mid_trg_avg_hf/n_sacs;
mid_trg_avg_gam = mid_trg_avg_gam/n_sacs;

n_stims = length(all_stim_inds);
stim_trg_avg = zeros(24,length(lags));
stim_trg_mat = zeros(24,n_stims,length(lags));
stim_trg_avg_hf = zeros(24,length(lags));
stim_trg_avg_gam = zeros(24,length(lags));
% stim_trg_mat_hf = zeros(24,n_stims,length(lags));
for i = 1:n_stims
    cur_set = (all_stim_inds(i)-backlag):(all_stim_inds(i)+forlag);
    stim_trg_avg = stim_trg_avg + all_lfp_data(cur_set,:)';
    stim_trg_mat(:,i,:) = all_lfp_data(cur_set,:)';
    stim_trg_avg_hf = stim_trg_avg_hf + used_hf(cur_set,:)';
    stim_trg_avg_gam = stim_trg_avg_gam + used_gam(cur_set,:)';
%     stim_trg_mat_hf(:,i,:) = used_hf(cur_set,:)';
end
stim_trg_avg = stim_trg_avg/n_stims;
stim_trg_avg_hf = stim_trg_avg_hf/n_stims;
stim_trg_avg_gam = stim_trg_avg_gam/n_stims;

% use_fta = fix_trg_avg(1:1:end,:);
% % use_fta2 = fix_trg_avg(2:2:end,:);
% fix_trg_csd = use_fta(1:end-2,:) + use_fta(3:end,:) - 2*use_fta(2:end-1,:);
% % fix_trg_csd2 = use_fta2(1:end-2,:) + use_fta2(3:end,:) - 2*use_fta2(2:end-1,:);
%% COMPARE TRIGGERED AVG LFP AMPS
close all
subplot(2,1,1)
plot(lags/Fsd,fix_trg_avg(1,:),'linewidth',1)
hold on
plot(lags/Fsd,fix_trg_avg(14,:),'r','linewidth',1)
plot(lags/Fsd,fix_trg_avg(22,:),'k','linewidth',1)
xlim([-0.05 0.25])
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time since fixation onset (s)','fontsize',20)
ylabel('Average amplitude','fontsize',20)
subplot(2,1,2)
plot(lags/Fsd,stim_trg_avg(1,:),'linewidth',1)
hold on
plot(lags/Fsd,stim_trg_avg(14,:),'r','linewidth',1)
plot(lags/Fsd,stim_trg_avg(22,:),'k','linewidth',1)
xlim([-0.05 0.25])
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time since stimulus onset (s)','fontsize',20)
ylabel('Average amplitude','fontsize',20)

figure
subplot(2,1,1)
plot(lags/Fsd,fix_trg_avg(1,:)/max(abs(fix_trg_avg(1,:))),'linewidth',1)
hold on
plot(lags/Fsd,fix_trg_avg(14,:)/max(abs(fix_trg_avg(14,:))),'r','linewidth',1)
plot(lags/Fsd,fix_trg_avg(22,:)/max(abs(fix_trg_avg(22,:))),'k','linewidth',1)
xlim([-0.05 0.25])
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time since fixation onset (s)','fontsize',16)
ylabel('Average amplitude','fontsize',16)
subplot(2,1,2)
plot(lags/Fsd,stim_trg_avg(1,:)/max(abs(stim_trg_avg(1,:))),'linewidth',1)
hold on
plot(lags/Fsd,stim_trg_avg(14,:)/max(abs(stim_trg_avg(14,:))),'r','linewidth',1)
plot(lags/Fsd,stim_trg_avg(22,:)/max(abs(stim_trg_avg(22,:))),'k','linewidth',1)
xlim([-0.05 0.25])
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time since stimulus onset (s)','fontsize',16)
ylabel('Average amplitude','fontsize',16)

close all
subplot(2,1,1)
plot(lags/Fsd,mid_trg_avg(1,:),'linewidth',1)
hold on
plot(lags/Fsd,mid_trg_avg(14,:),'r','linewidth',1)
plot(lags/Fsd,mid_trg_avg(22,:),'k','linewidth',1)
xlim([-0.05 0.25])
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time since fixation onset (s)','fontsize',16)
ylabel('Average amplitude','fontsize',16)
subplot(2,1,2)
plot(lags/Fsd,stim_trg_avg(1,:),'linewidth',1)
hold on
plot(lags/Fsd,stim_trg_avg(14,:),'r','linewidth',1)
plot(lags/Fsd,stim_trg_avg(22,:),'k','linewidth',1)
xlim([-0.05 0.25])
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time since stimulus onset (s)','fontsize',16)
ylabel('Average amplitude','fontsize',16)

%% COMPARE TRG AVG HF AMPS
close all
subplot(2,1,1)
plot(lags/Fsd,fix_trg_avg_hf(1,:),'linewidth',1)
hold on
plot(lags/Fsd,fix_trg_avg_hf(14,:),'r','linewidth',1)
plot(lags/Fsd,fix_trg_avg_hf(22,:),'k','linewidth',1)
xlim([-0.05 0.25])
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time since fixation onset (s)','fontsize',16)
ylabel('Average amplitude','fontsize',16)
subplot(2,1,2)
plot(lags/Fsd,stim_trg_avg_hf(1,:),'linewidth',1)
hold on
plot(lags/Fsd,stim_trg_avg_hf(14,:),'r','linewidth',1)
plot(lags/Fsd,stim_trg_avg_hf(22,:),'k','linewidth',1)
xlim([-0.05 0.25])
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time since stimulus onset (s)','fontsize',16)
ylabel('Average amplitude','fontsize',16)

%%


use_ch = 2:2:24;
cmap = jet(length(use_ch));
for i = 1:length(use_ch)
%     plot(lags/Fsd,fix_trg_avg(use_ch(i),:)+i*1,'color',cmap(i,:));
    plot(lags/Fsd,fix_trg_avg(use_ch(i),:)+i*0.5,'k');
    hold on
end

figure
use_ch = 2:2:24;
cmap = jet(length(use_ch));
for i = 1:length(use_ch)
%     plot(lags/Fsd,fix_trg_avg(use_ch(i),:)+i*1,'color',cmap(i,:));
    plot(lags/Fsd,sac_trg_avg(use_ch(i),:)+i*0.5,'k');
    hold on
end

%%
use_ch = 2:2:24;
cmap = jet(length(use_ch));
for i = 1:length(use_ch)
%     plot(lags/Fsd,fix_trg_avg(use_ch(i),:)+i*1,'color',cmap(i,:));
    plot(lags/Fsd,fix_trg_avg_hf(use_ch(i),:)+i*0.1,'k');
    hold on
end

figure
use_ch = 2:2:24;
cmap = jet(length(use_ch));
for i = 1:length(use_ch)
%     plot(lags/Fsd,fix_trg_avg(use_ch(i),:)+i*1,'color',cmap(i,:));
    plot(lags/Fsd,sac_trg_avg_hf(use_ch(i),:)+i*0.1,'k');
    hold on
end

%%
cur_sac_durs = all_sac_durs(blockids < 4);
% cur_sac_durs = all_sac_amps(blockids < 4);
n_bins = 10;
bin_edges = prctile(cur_sac_durs,linspace(0,100,n_bins+1));
bin_cents = 0.5*bin_edges(1:end-1)+0.5*bin_edges(2:end);
for i = 1:n_bins
    cur = find(cur_sac_durs >= bin_edges(i) & cur_sac_durs < bin_edges(i+1));
%     cur_fix_trg_avg(:,i,:) = mean(fix_trg_mat(:,cur,:),2);
%     cur_sac_trg_avg(:,i,:) = mean(sac_trg_mat(:,cur,:),2);
%     cur_mid_trg_avg(:,i,:) = mean(mid_trg_mat(:,cur,:),2);
    cur_fix_trg_avg(:,i,:) = mean(fix_trg_mat_hf(:,cur,:),2);
    cur_sac_trg_avg(:,i,:) = mean(sac_trg_mat_hf(:,cur,:),2);
    cur_mid_trg_avg(:,i,:) = mean(mid_trg_mat_hf(:,cur,:),2);
end

%%
close all
ch = 13;
cmap = jet(n_bins);
for i = 1:n_bins
    subplot(3,1,1)
    plot(lags/Fsd,squeeze(cur_fix_trg_avg(ch,i,:)),'color',cmap(i,:))
    hold on
    xlim([0 0.3])
    xlabel('Time since fixation onset (s)','fontsize',14)
    ylabel('Average Amplitude','fontsize',14)
    subplot(3,1,2)
    plot(lags/Fsd,squeeze(cur_sac_trg_avg(ch,i,:)),'color',cmap(i,:))
    hold on
    xlim([0 0.3])
    xlabel('Time since saccade onset (s)','fontsize',14)
    ylabel('Average Amplitude','fontsize',14)
    subplot(3,1,3)
    plot(lags/Fsd,squeeze(cur_mid_trg_avg(ch,i,:)),'color',cmap(i,:))
    hold on
    xlim([0 0.3])
    xlabel('Time since saccade peak (s)','fontsize',14)
    ylabel('Average Amplitude','fontsize',14)
end

% figure
% plot(lags/Fsd,squeeze(mean(fix_trg_mat(ch,:,:),2)),'k')
% hold on
% plot(lags/Fsd,squeeze(mean(sac_trg_mat(ch,:,:),2)),'b')
% plot(lags/Fsd,squeeze(mean(mid_trg_mat(ch,:,:),2)),'r')
%     xlim([0 0.3])
% xlabel('Time since onset (s)','fontsize',14)
% ylabel('Average amplitude','fontsize',14)
% legend('Fixation onset','Saccade onset','Saccade Peak')

figure
plot(lags/Fsd,squeeze(mean(fix_trg_mat_hf(ch,:,:),2)),'k')
hold on
plot(lags/Fsd,squeeze(mean(sac_trg_mat_hf(ch,:,:),2)),'b')
plot(lags/Fsd,squeeze(mean(mid_trg_mat_hf(ch,:,:),2)),'r')
    xlim([0 0.3])
xlabel('Time since onset (s)','fontsize',14)
ylabel('Average amplitude','fontsize',14)
legend('Fixation onset','Saccade onset','Saccade Peak')

%%
subplot(3,1,[1 2])
temp = abs(fix_trg_avg);
temp = bsxfun(@rdivide,temp,sum(temp));
imagesc(lags/Fsd,1:24,temp)
shg
xlim([0 0.3])
xlabel('Time (s)','fontsize',14)
ylabel('Channel Number','fontsize',14)
set(gca,'ydir','normal')
subplot(3,1,3)
plot(lags/Fsd,squeeze(fix_trg_avg(22,:)),'k')
hold on
plot(lags/Fsd,squeeze(fix_trg_avg(14,:)),'b')
plot(lags/Fsd,squeeze(fix_trg_avg(2,:)),'r')
xlim([0 0.3])
xlabel('Time since fixation onset (s)','fontsize',14)
ylabel('Amplitude','fontsize',14)
legend('Deep','Middle','Superficial')