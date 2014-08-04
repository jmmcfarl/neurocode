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

Pix2Deg = 0.018837;

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

cd /Users/James/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

%%
Fs = 1000;
dsf = 4;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,120/niqf,'low');
[b_hf,a_hf] = butter(2,120/niqf,'high');
[b_gam,a_gam] = butter(2,[40 60]/niqf);
[b_lf,a_lf] = butter(2,60/niqf,'low');
[b_alpha,a_alpha] = butter(2,[5 12]/niqf);
gam_sm = round(Fsd*0.015);

backlag = round(0.1*Fsd);
forwardlag = round(0.5*Fsd);
lags = -backlag:forwardlag;

scales = logspace(log10(2.5),log10(50),30);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);


%%
all_lfp_mat = [];
all_lfphf_mat = [];
all_lfpgam_mat = [];
all_used_fixs = [];
all_lfpampgram_mat = [];
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);

    %%
    cd /Users/James/Data/bruce/2_27_12
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
    
    if blockid == 4
        LFP.Trials = LFP.Trials(1:5);
    end

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
        cur_len(i) = size(LFP.Trials(i).LFP,1);
        if i < n_lfp_trials
            next_start = lfp_trial_start(i+1);
            start_len = length(lfp_trial_start(i):1/Fs:next_start);
        else
            next_start = Inf;
            start_len = Inf;
        end
        cur_end(i) = min(cur_len(i),start_len);
        cur_t = lfp_trial_start(i):1/Fs:(lfp_trial_start(i)+cur_end(i)/Fs);
        cur_t(cur_end(i)+1:end) = [];
        lfp_time = [lfp_time cur_t];
        lfp_samps = [lfp_samps; LFP.Trials(i).LFP(1:cur_end(i),:)];
    end
    
    lfp_samps = zscore(lfp_samps);
    
    lfp_timed = downsample(lfp_time,dsf);
    
    lfp_samps_hf = filtfilt(b_hf,a_hf,lfp_samps);
    lfp_samps_hf = abs(hilbert(lfp_samps_hf));
    
    lfp_samps_lf = filtfilt(b_lf,a_lf,lfp_samps_hf);
    lfp_samps_lf = downsample(lfp_samps_lf,dsf);
    %     lfp_samps_lf = zscore(lfp_samps_lf);
    
    lfp_samps_gam = filtfilt(b_gam,a_gam,lfp_samps);
    %     lfp_samps_gam = abs(hilbert(lfp_samps_gam));
    
    lfp_samps_gam = downsample(lfp_samps_gam,dsf);
    for i = 1:24
        lfp_samps_gam(:,i) = sqrt(jmm_smooth_1d_cor(lfp_samps_gam(:,i).^2,gam_sm));
    end
    lfp_samps_gam = zscore(lfp_samps_gam);
    
    lfp_sampsd = filtfilt(b,a,lfp_samps);
    lfp_sampsd = downsample(lfp_sampsd,dsf);
%     lfp_sampsd = zscore(lfp_sampsd);
    lfp_timed = downsample(lfp_time,dsf);
    
    
    all_ampgram = zeros(24,length(lfp_timed),length(wfreqs));
    for cc = 1:24
        %     cc = 17;
        fprintf('Channel %d of %d\n',cc,24);
        temp = cwt(lfp_sampsd(:,cc),scales,'cmor1-1');
        tampgram = abs(temp);
        all_ampgram(cc,:,:) = tampgram';
    end
    clear temp tampgram
    all_ampgram = bsxfun(@minus,all_ampgram,mean(all_ampgram,2));
    all_ampgram = bsxfun(@rdivide,all_ampgram,std(all_ampgram,[],2));
    
    %%
    cur_fix_set = find(blockids == blockid);
    cur_fix_start_inds = round(interp1(lfp_timed,1:length(lfp_timed),all_fix_start_times(cur_fix_set)));
    cur_fix_stop_inds = round(interp1(lfp_timed,1:length(lfp_timed),all_fix_stop_times(cur_fix_set)));
   
    bad = find(cur_fix_start_inds < backlag | cur_fix_stop_inds > length(lfp_timed) - forwardlag);
    cur_fix_start_inds(bad) = [];
    cur_fix_stop_inds(bad) = [];
    cur_fix_set(bad) = [];

    cur_fix_durs = cur_fix_stop_inds - cur_fix_start_inds;
    cur_n_fixs = length(cur_fix_start_inds);
    
    fix_trig_lfp_mat = nan(24,cur_n_fixs,length(lags));
    fix_trig_lfphf_mat = nan(24,cur_n_fixs,length(lags));
    fix_trig_lfpgam_mat = nan(24,cur_n_fixs,length(lags));
    fix_trig_ampgram_mat = nan(24,cur_n_fixs,length(lags),length(wfreqs));
    for i = 1:cur_n_fixs
        cur_set = (cur_fix_start_inds(i)-backlag):(cur_fix_start_inds(i)+forwardlag);
        if cur_fix_durs(i) < forwardlag
            cur_set(backlag+cur_fix_durs(i):end) = [];
        end
        fix_trig_lfp_mat(:,i,1:length(cur_set)) = lfp_sampsd(cur_set,:)';
        fix_trig_lfphf_mat(:,i,1:length(cur_set)) = lfp_samps_hf(cur_set,:)';
        fix_trig_lfpgam_mat(:,i,1:length(cur_set)) = lfp_samps_gam(cur_set,:)';
        fix_trig_ampgram_mat(:,i,1:length(cur_set),:) = squeeze(all_ampgram(:,cur_set,:));
    end
    
    cur_sac_amps = all_sac_amps(cur_fix_set);
    cur_sac_durs = all_sac_durs(cur_fix_set);
    cur_sac_dirs = atan2(all_sac_dy(cur_fix_set),all_sac_dx(cur_fix_set));
    
    %%
    all_lfp_mat = cat(2,all_lfp_mat,fix_trig_lfp_mat);
    all_lfphf_mat = cat(2,all_lfphf_mat,fix_trig_lfphf_mat);
    all_lfpgam_mat = cat(2,all_lfpgam_mat,fix_trig_lfpgam_mat);
    all_lfpampgram_mat = cat(2,all_lfpampgram_mat,fix_trig_ampgram_mat);
    all_used_fixs = [all_used_fixs; cur_fix_set];
    
end


%%
used_sac_amps = all_sac_amps(all_used_fixs);
used_sac_durs = all_sac_durs(all_used_fixs);
used_sac_dirs = atan2(all_sac_dy(all_used_fixs),all_sac_dx(all_used_fixs));

%%
sac_trg_avg = squeeze(nanmean(all_lfp_mat,2));
sac_trg_avg_gam = squeeze(nanmean(all_lfpgam_mat,2));
sac_trg_avg_ampmat = squeeze(nanmean(all_lfpampgram_mat,2));

% sm_sacs = find(used_sac_amps < 1);
% big_sacs = find(used_sac_amps > 1);
% 
% smsac_trg_avg = squeeze(nanmean(all_lfp_mat(:,sm_sacs,:),2));
% bisac_trg_avg = squeeze(nanmean(all_lfp_mat(:,big_sacs,:),2));
% 
% %%
% for ch = 1:24;
% plot(lags/Fsd,smsac_trg_avg(ch,:))
% hold on
% plot(lags/Fsd,bisac_trg_avg(ch,:),'r')
% pause
% clf
% end
%%
% n_bins = 15;
% bin_edges = prctile(used_sac_durs,linspace(0,100,n_bins+1));
% for i = 1:n_bins
%     cur_set = find(used_sac_durs >= bin_edges(i) & used_sac_durs < bin_edges(i+1));
%     avg_lfp(i,:,:) = nanmean(all_lfp_mat(:,cur_set,:),2);
%     avg_lfp_hf(i,:,:) = nanmean(all_lfphf_mat(:,cur_set,:),2);
%     avg_lfp_gam(i,:,:) = nanmean(all_lfpgam_mat(:,cur_set,:),2);
%     
%     avg_lfp_p(i,:,:) = nanmean(all_lfp_mat_p(:,cur_set,:),2);
%     avg_lfp_hf_p(i,:,:) = nanmean(all_lfphf_mat_p(:,cur_set,:),2);
%     avg_lfp_gam_p(i,:,:) = nanmean(all_lfpgam_mat_p(:,cur_set,:),2);
%     
%     avg_durs(i) = mean(used_sac_durs(cur_set));
% end
% 

%%
figure
imagesc(lags/Fsd,1:24,sac_trg_avg_gam);
set(gca,'ydir','normal','fontsize',16,'fontname','arial')
xlabel('Time (s)','fontname','arial')
ylabel('Channel','fontname','arial')
colorbar

figure
imagesc(lags/Fsd,1:24,sac_trg_avg);
set(gca,'ydir','normal','fontsize',16,'fontname','arial')
xlabel('Time (s)','fontname','arial')
ylabel('Channel','fontname','arial')
colorbar

%%
% figure
% pcolor(lags/Fsd,wfreqs,squeeze(sac_trg_avg_ampmat(1,:,:))');shading flat
% set(gca,'ydir','normal','fontsize',16,'fontname','arial')
% xlabel('Time (s)','fontname','arial')
% ylabel('Frequency','fontname','arial')
% colorbar

figure
pcolor(lags/Fsd,wfreqs,squeeze(sac_trg_avg_ampmat(14,:,:))');shading flat
set(gca,'fontsize',16,'fontname','arial')
xlabel('Time (s)','fontname','arial')
ylabel('Frequency','fontname','arial')
colorbar