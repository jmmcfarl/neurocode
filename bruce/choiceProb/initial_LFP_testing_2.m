%% LOAD PROCESSED DATA
clear all
% cd('/home/james/Data/bruce/ChoiceProb/')
cd('~/Data/bruce/ChoiceProb/')

Expt_name = 'M230';

if strcmp(Expt_name,'M239')
    load('M239/lemM239.image.ORBW.LFP.mat')
elseif strcmp(Expt_name,'M230')
    load('M230/lemM230.image.ORBW.LFP.mat')
end

%% SET BLOCK TRIAL RANGES
if strcmp(Expt_name,'M239')
    block_trial_boundaries = [1 157;
        158 265;
        266 380;
        381 499;
        500 627;
        628 747;
        748 870;
        871 988;
        989 1104];
elseif strcmp(Expt_name,'M230')
    block_trial_boundaries = [1 140;
        141 275;
        276 409;
        410 533;
        534 668;
        669 801;
        802 936;
        937 1058;
        1059 1116;
        1117 1163;
        1164 1240;
        1241 1375;
        1376 1415;
        1416 1480;
        1481 1612;
        1613 1740];
end
n_blocks = size(block_trial_boundaries,1);

%% Assemble stim, Robs, LFP into trial structure
Ntrials = length(AllExpt.Expt.Trials);
LFP_Fs = 1.0 / AllExpt.Expt.Header.LFPsamplerate;  % this is 1 kHz
dt = 0.005;  % assuming 10 ms frame time exactly
trialDur = 2; %in sec.
Nframes = 200; %number of video frames (@100Hz) per trial
NT = ceil(trialDur/dt); %number of time bins per trial

%% SU and MU data (combined)
Nunits = length(AllExpt.Header);
Ucellnumbers = [AllExpt.Header(:).cellnumber];
Uchans = round([AllExpt.Header(:).probe]);

%% relevant trial variables
trialNums = [AllExpt.Expt.Trials(:).Trial];
trialOR = [AllExpt.Expt.Trials(:).or];
trialOB = [AllExpt.Expt.Trials(:).ob];
trialrwDir = [AllExpt.Expt.Trials(:).rwdir];
trialrespDir = [AllExpt.Expt.Trials(:).RespDir];
trialSe = [AllExpt.Expt.Trials(:).se];

%% get total trial-by-trial spike counts for each unit
tot_spks_per_trial = zeros(Ntrials,Nunits);
for tr = 1:Ntrials
    for nn = 1:Nunits
        trindx = find( AllExpt.Spikes{nn}.Trial == trialNums(tr));
        if ~isempty(trindx)
            tot_spks_per_trial(tr,nn) = sum(AllExpt.Spikes{nn}.Spikes{trindx} <= trialDur/1e-4 & AllExpt.Spikes{nn}.Spikes{trindx} >= 0);
        end
    end
end

%find any blocks where the unit has no spikes and set these values to nans
spikes_per_block = nan(n_blocks,Nunits);
for bb = 1:n_blocks
    spikes_per_block(bb,:) = sum(tot_spks_per_trial(block_trial_boundaries(bb,1):block_trial_boundaries(bb,2),:));
end

tot_spks_per_trial_norm = tot_spks_per_trial;
for cc = 1:Nunits
    bad_blocks = find(spikes_per_block(:,cc) == 0);
    for ii = bad_blocks'
        tot_spks_per_trial_norm(block_trial_boundaries(ii,1):block_trial_boundaries(ii,2),cc) = nan;
    end
end

%% Get Binned Spikes
fullRobs = zeros(NT,Ntrials,Nunits);
for tr = 1:Ntrials
    for nn = 1:Nunits
        trindx = find( AllExpt.Spikes{nn}.Trial == trialNums(tr));
        if ~isempty(trindx)
            spks = double( AllExpt.Spikes{nn}.Spikes{trindx} ) * 1e-4;  % no adjustment from start of trial
            cur_Robs = histc( spks, (0:(NT))*dt);
            fullRobs(:,tr,nn) = cur_Robs(1:end-1);
        end
    end
end
for cc = 1:Nunits
    fullRobs(:,isnan(tot_spks_per_trial_norm(:,cc)),cc) = nan;
end

trialSpkCnts = squeeze(nansum(fullRobs,1));
totSpkCnts = squeeze(nansum(trialSpkCnts,1));
avg_spk_rates = nanmean(reshape(fullRobs,[],Nunits));

%% CALCULATE CHOICE PROBS
tot_spks_per_trial_norm = tot_spks_per_trial;
tot_spks_per_trial_norm(tot_spks_per_trial == 0) = nan;

sig_trials = find(trialOB < 130);
stim_up = sig_trials(trialrwDir(sig_trials)==1);
stim_down = sig_trials(trialrwDir(sig_trials)==-1);
sig_prob = nan(Nunits,1);
for cc = 1:Nunits
    cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(sig_trials,cc)));
    up_hist = histc(tot_spks_per_trial_norm(stim_up,cc),cur_resp_ax);
    down_hist = histc(tot_spks_per_trial_norm(stim_down,cc),cur_resp_ax);
    fract_cor = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_up,cc)));
    fract_incor = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_down,cc)));
    sig_prob(cc) = trapz(fract_incor,fract_cor);
end

test_trials = find(trialOB == 130 & trialrespDir ~= 0);
% test_trials = find(trialrespDir ~= 0);
resp_up = test_trials(trialrespDir(test_trials) == 1);
resp_down = test_trials(trialrespDir(test_trials) == -1);
choice_prob = nan(Nunits,1);
cp_pval = nan(Nunits,1);
for cc = 1:Nunits
    cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(test_trials,cc)));
    up_hist = histc(tot_spks_per_trial_norm(resp_up,cc),cur_resp_ax);
    down_hist = histc(tot_spks_per_trial_norm(resp_down,cc),cur_resp_ax);
    true_pos = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(resp_up,cc)));
    false_pos = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(resp_down,cc)));
    choice_prob(cc) = trapz(false_pos,true_pos);
    [~,cp_pval(cc)] = ttest2(tot_spks_per_trial_norm(resp_up,cc),tot_spks_per_trial_norm(resp_down,cc));
end

nboot = 500;
boot_sig_prob = nan(Nunits,nboot);
boot_choice_prob = nan(Nunits,nboot);
for nn = 1:nboot
    randrwdir = trialrwDir(randperm(Ntrials));
    stim_up = sig_trials(randrwdir(sig_trials)==1);
    stim_down = sig_trials(randrwdir(sig_trials)==-1);
    
    randrespdir = trialrespDir(randperm(Ntrials));
    resp_up = test_trials(randrespdir(test_trials) == 1);
    resp_down = test_trials(randrespdir(test_trials) == -1);
    for cc = 1:Nunits
        cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(sig_trials,cc)));
        up_hist = histc(tot_spks_per_trial_norm(stim_up,cc),cur_resp_ax);
        down_hist = histc(tot_spks_per_trial_norm(stim_down,cc),cur_resp_ax);
        fract_cor = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_up,cc)));
        fract_incor = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_down,cc)));
        boot_sig_prob(cc,nn) = trapz(fract_incor,fract_cor);
        
        cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(test_trials,cc)));
        up_hist = histc(tot_spks_per_trial_norm(resp_up,cc),cur_resp_ax);
        down_hist = histc(tot_spks_per_trial_norm(resp_down,cc),cur_resp_ax);
        fract_cor = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(resp_up,cc)));
        fract_incor = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(resp_down,cc)));
        boot_choice_prob(cc,nn) = trapz(fract_incor,fract_cor);
    end
end
choice_prob_ci = prctile(boot_choice_prob',[5 95]);
sig_prob_ci = prctile(boot_sig_prob',[5 95]);
choice_prob_z = (choice_prob - nanmean(boot_choice_prob,2))./nanstd(boot_choice_prob,[],2);
sig_prob_z = (sig_prob - nanmean(boot_sig_prob,2))./nanstd(boot_sig_prob,[],2);

%% USE WITHIN-CHOICE AVGS
test_trials = find(trialOB == 130 & trialrespDir ~= 0);

tot_spks_per_trial_ms = bsxfun(@minus,tot_spks_per_trial_norm,nanmean(tot_spks_per_trial_norm(test_trials,:)));

cur_resp = trialrespDir(test_trials);
resp_up = test_trials(cur_resp == 1);
resp_down = test_trials(cur_resp == -1);

avg_up_resp = nanmean(tot_spks_per_trial_ms(resp_up,:));
avg_down_resp = nanmean(tot_spks_per_trial_ms(resp_down,:));
choice_cond_avgs = [avg_up_resp; avg_down_resp];
choice_cond_covmat = squeeze(nanmean(bsxfun(@times,choice_cond_avgs,reshape(choice_cond_avgs,2,1,Nunits))));
% choice_cond_covmat = nancov(choice_cond_avgs,1);

normfac = sqrt(nanvar(tot_spks_per_trial_ms(test_trials,:))'*nanvar(tot_spks_per_trial_ms(test_trials,:)));

choice_cond_corrmat = choice_cond_covmat./normfac;
%% AVG OF WITHIN-CHOICE COVMATS
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
% test_trials = find(trialrespDir~= 0);

tot_spks_per_trial_ms = bsxfun(@minus,tot_spks_per_trial_norm,nanmean(tot_spks_per_trial_norm(test_trials,:)));
tot_spks2 = reshape(tot_spks_per_trial_ms,[],1,Nunits);

cur_resp = trialrespDir(test_trials);
resp_up = test_trials(cur_resp == 1);
resp_down = test_trials(cur_resp == -1);

%covariance of responses on unequal trials with up choice
[II,JJ] = meshgrid(1:length(resp_up));
uset = find(II > JJ);
avg_up_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(resp_up(II(uset)),:),tot_spks2(resp_up(JJ(uset)),:,:)),1));

%covariance of responses on unequal trials with up choice
[II,JJ] = meshgrid(1:length(resp_down));
uset = find(II ~= JJ);
avg_down_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(resp_down(II(uset)),:),tot_spks2(resp_down(JJ(uset)),:,:)),1));

%avg chioce-conditional covariances
avg_choice_covmat = 0.5*avg_down_covmat + 0.5*avg_up_covmat;

[II,JJ] = meshgrid(1:length(test_trials));
cur_Dmat = abs(squareform(pdist(cur_resp')));
cur_Dmat(JJ == II) = nan;
uset = find(cur_Dmat == 0);
avg_same_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(test_trials(II(uset)),:),tot_spks2(test_trials(JJ(uset)),:,:)),1));

noise_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(test_trials,:),tot_spks2(test_trials,:,:)),1));

%normalization matrix of variance products
normfac = sqrt(nanvar(tot_spks_per_trial_ms(test_trials,:))'*nanvar(tot_spks_per_trial_ms(test_trials,:)));

avg_choice_corrmat = avg_choice_covmat./normfac;
% opp_corrmat = avg_opp_covmat./normfac;
avg_same_corrmat = avg_same_covmat./normfac;

noise_corrmat = noise_covmat./normfac;


%% USE WAVELET ANALYSIS TO COMPUTE PHASE-LOCKING SPECTRA FOR EACH UNIT
LFP_trial_taxis = AllExpt.Expt.Header.LFPtimes*1e-4;

LFP_dsf = 2;
LFP_Fsd = LFP_Fs/LFP_dsf;
%anti-aliasing filter and high-pass filter
aa_hcf = LFP_Fsd/2*0.8;
% [b_aa,a_aa] = butter(4,aa_hcf/(LFP_Fs/2),'low');
aa_lcf = 0.5;
[b_aa,a_aa] = butter(2,[aa_lcf aa_hcf]/(LFP_Fs/2));

LFP_trial_taxis_ds = downsample(LFP_trial_taxis,LFP_dsf);

nprobes = 24;
uprobes = 1:1:nprobes;

%wavelet parameters
nwfreqs = 30;
min_freq = 1.; max_freq = 100;
% nwfreqs = 5;
% min_freq = 1; max_freq = 6;
min_scale = 1/max_freq*LFP_Fsd;
max_scale = 1/min_freq*LFP_Fsd;
wavetype = 'cmor1-1';
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,wavetype,1/LFP_Fsd);

% beg_buffer = round(0.15/dt);
beg_buffer = round(0/dt);
end_buffer = round(0/dt);
trial_dur = round(2/dt);
R_trial_taxis = (beg_buffer:(trial_dur - end_buffer))*dt;
TLEN = length(R_trial_taxis);

all_spk_id = [];
all_spk_phases = nan(sum(tot_spks_per_trial(:)),length(wfreqs),length(uprobes));
% trial_LFP_cwt = nan(Ntrials,TLEN,length(wfreqs),length(uprobes));
trial_LFP_real = nan(Ntrials,TLEN,length(wfreqs),length(uprobes));
trial_LFP_imag = nan(Ntrials,TLEN,length(wfreqs),length(uprobes));
trial_LFPs = nan(Ntrials,TLEN,length(uprobes));
for tr = 1:Ntrials
    fprintf('Trial %d of %d\n',tr,Ntrials);
    cur_LFPs = double(AllExpt.Expt.Trials(tr).LFP(:,uprobes));
    bad_LFPs = isnan(cur_LFPs(:,1));
    cur_LFPs(isnan(cur_LFPs)) = 0;
    cur_LFPs = filtfilt(b_aa,a_aa,cur_LFPs);
    cur_LFPs = downsample(cur_LFPs,LFP_dsf);
    bad_LFPs = downsample(bad_LFPs,LFP_dsf);
    
    cur_cwt = nan(size(cur_LFPs,1),length(wfreqs),length(uprobes));
    for cc = 1:length(uprobes)
        cur_cwt(:,:,cc) = cwt(cur_LFPs(:,cc),scales,wavetype)';
    end
    
    cur_LFPs(bad_LFPs,:) = nan;
    cur_cwt(bad_LFPs,:,:) = nan;
%     trial_LFP_cwt(tr,:,:,:) = interp1(LFP_trial_taxis_ds,cur_cwt,R_trial_taxis);
    trial_LFP_real(tr,:,:,:) = interp1(LFP_trial_taxis_ds,real(cur_cwt),R_trial_taxis);
    trial_LFP_imag(tr,:,:,:) = interp1(LFP_trial_taxis_ds,imag(cur_cwt),R_trial_taxis);
    trial_LFPs(tr,:,:) = interp1(LFP_trial_taxis_ds,cur_LFPs,R_trial_taxis);
end
% trial_LFP_cwt = permute(trial_LFP_cwt,[2 1 3 4]);
trial_LFP_real = permute(trial_LFP_real,[2 1 3 4]);
trial_LFP_imag = permute(trial_LFP_imag,[2 1 3 4]);
trial_LFPs = permute(trial_LFPs,[2 1 3]);

% trial_LFP_cwt = reshape(trial_LFP_cwt,[],length(wfreqs)*length(uprobes));
trial_LFP_real = reshape(trial_LFP_real,[],length(wfreqs)*length(uprobes));
trial_LFP_imag = reshape(trial_LFP_imag,[],length(wfreqs)*length(uprobes));
trial_LFPs = reshape(trial_LFPs,[],length(uprobes));

%%
trial_LFPs = nanzscore(trial_LFPs);
trial_LFP_mag = sqrt(trial_LFP_real.^2 + trial_LFP_imag.^2);

% trial_LFP_cwt = bsxfun(@rdivide,trial_LFP_cwt,nanstd(abs(trial_LFP_cwt)));
trial_LFP_real = bsxfun(@rdivide,trial_LFP_real,nanstd(trial_LFP_mag));
trial_LFP_imag = bsxfun(@rdivide,trial_LFP_imag,nanstd(trial_LFP_mag));


%%
tbt_LFPs = reshape(trial_LFPs,[TLEN Ntrials length(uprobes)]);
% tbt_LFP_cwt = reshape(trial_LFP_cwt,[TLEN Ntrials length(wfreqs) length(uprobes)]);
tbt_LFP_real = reshape(trial_LFP_real,[TLEN Ntrials length(wfreqs) length(uprobes)]);
tbt_LFP_imag = reshape(trial_LFP_imag,[TLEN Ntrials length(wfreqs) length(uprobes)]);

% %%
% test_trials = find(trialOB == 130 & trialrespDir ~= 0);
% % test_trials = find(trialrespDir ~= 0);
% un_stim_seeds = unique(trialSe(test_trials));
% 
% beg_buff = 0; %number of bins from beginning of trial to exclude
% end_buff = 0;
% 
% % maxlag_ED = 4;
% % ED_space = 0.2;
% % ED_bin_edges = 0:ED_space:maxlag_ED;
% % ED_bin_centers = (ED_bin_edges(1:end-1)+ED_bin_edges(2:end))/2;
% 
% %subtract off avg rates for times within this set of trials
% [RR,TT] = meshgrid(1:Ntrials,1:NT);
% istrain = (ismember(RR(:),test_trials)) & (TT(:) > beg_buff);
% fullRobs_resh = reshape(fullRobs,[],Nunits);
% fullRobs_ms = bsxfun(@minus,fullRobs,reshape(nanmean(fullRobs_resh(istrain,:)),[1 1 Nunits]));
% % fullRobs_ms = bsxfun(@minus,fullRobs,nanmean(fullRobs(:,test_trials,:),2));
% tot_var = nanvar(fullRobs_resh(istrain,:));
% fullRobs_ms = fullRobs_ms(:,test_trials,:);
% 
% u_LFP_chs = [1:8];
% for ww = 1:length(wfreqs)
% ww
% u_lfp_times = find(R_trial_taxis > 0 & R_trial_taxis <= (trial_dur)*dt);
% % tbt_LFP_data = squeeze(tbt_LFP_cwt(u_lfp_times,test_trials,ww,:));
% LFP_real = squeeze(tbt_LFP_real(u_lfp_times,test_trials,ww,u_LFP_chs));
% LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,test_trials,ww,u_LFP_chs));
% LFP_phase = atan2(LFP_imag,LFP_real);
% LFP_amp = sqrt(LFP_imag.^2 + LFP_real.^2);
% 
% LFP_real = LFP_real./LFP_amp;
% LFP_imag = LFP_imag./LFP_amp;
% 
% %compute bin edge locations for LFP metric
% all_LFP_dists = [];
% for tr = 1:length(un_stim_seeds)
%     cur_tr_set = find(trialSe(test_trials) == un_stim_seeds(tr));
% %     fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
%     [II,JJ] = meshgrid(1:length(cur_tr_set));
%     if length(cur_tr_set) >= 2
%         for tt = (beg_buff+1):(NT-end_buff)
% %             cur_LFP_state = squeeze(tbt_LFP_data(tt,cur_tr_set,:));
% %             cur_LFP_state = cat(2,real(cur_LFP_state),imag(cur_LFP_state));
%             cur_LFP_state = cat(2,squeeze(LFP_real(tt,cur_tr_set,:)),squeeze(LFP_imag(tt,cur_tr_set,:)));
% %             cur_LFP_state = squeeze(LFP_phase(tt,cur_tr_set,:));
% %             cur_LFP_state = squeeze(LFP_amp(tt,cur_tr_set,:));
% 
%             cur_Dmat = squareform(pdist(cur_LFP_state))/sqrt(length(uprobes));
% %             cur_Dmat = abs(circ_dist2(cur_LFP_state));
% %             cur_Dmat = squareform(pdist(cur_LFP_state));
%             cur_Dmat(JJ >= II) = nan;
%             all_LFP_dists = cat(1,all_LFP_dists,cur_Dmat(~isnan(cur_Dmat)));
%         end
%     end
% end
% n_LFP_bins = 6;
% ED_bin_edges = prctile(all_LFP_dists,0:100/n_LFP_bins:100);
% ED_bin_centers = (ED_bin_edges(1:end-1)+ED_bin_edges(2:end))/2;
% 
% LFP_metric_XC = zeros(length(ED_bin_centers),Nunits,Nunits);
% LFP_metric_cnt = zeros(length(ED_bin_centers),Nunits,Nunits);
% rand_XC = zeros(Nunits,Nunits);
% rand_cnt = zeros(Nunits,Nunits);
% for tr = 1:length(un_stim_seeds)
%     
%     cur_tr_set = find(trialSe(test_trials) == un_stim_seeds(tr));
% %     fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
%     [II,JJ] = meshgrid(1:length(cur_tr_set));
%     
%     if length(cur_tr_set) >= 2
%     
%     for tt = (beg_buff+1):(NT-end_buff)
%         cur_Robs = squeeze(fullRobs_ms(tt,cur_tr_set,:));
%         cur_Robs2 = reshape(cur_Robs,[],1,Nunits);
%         
% %         cur_LFP_state = squeeze(tbt_LFP_data(tt,cur_tr_set,:));
% %         cur_LFP_state = cat(2,real(cur_LFP_state),imag(cur_LFP_state));
%             cur_LFP_state = cat(2,squeeze(LFP_real(tt,cur_tr_set,:)),squeeze(LFP_imag(tt,cur_tr_set,:)));
% %             cur_LFP_state = squeeze(LFP_phase(tt,cur_tr_set,:));
% %             cur_LFP_state = squeeze(LFP_amp(tt,cur_tr_set,:));
%          
%         cur_Dmat = squareform(pdist(cur_LFP_state))/sqrt(length(uprobes));
% %             cur_Dmat = abs(circ_dist2(cur_LFP_state));
%         cur_Dmat(JJ >= II) = nan;
%         for jj = 1:length(ED_bin_centers)
%             curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
%             temp = bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:));
%             LFP_metric_XC(jj,:,:) = LFP_metric_XC(jj,:,:) + nansum(temp,1);
%             LFP_metric_cnt(jj,:,:) = LFP_metric_cnt(jj,:,:) + sum(~isnan(temp),1);
%             %         cur_XC(tt,jj,:,:) = (nanmean(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
%         end
%         curset = ~isnan(cur_Dmat);
%         temp = bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:));
%         rand_XC = rand_XC + squeeze(nansum(temp,1));
%         rand_cnt = rand_cnt + squeeze(sum(~isnan(temp),1));
%         %     rand_XC(tt,:,:) = squeeze(nanmean(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
%     end
%     end
% end
% 
% LFP_metric_XC = LFP_metric_XC./LFP_metric_cnt;
% rand_XC = rand_XC./rand_cnt;
% 
% normfac = sqrt(tot_var'*tot_var);
% 
% LFP_expl_cov = squeeze(LFP_metric_XC(1,:,:)) - rand_XC;
% LFP_expl_corr = LFP_expl_cov./normfac;
% 
% all_test(ww,:,:) = LFP_expl_corr;
% end
%%
% close all
% for ii = 1:Nunits
% plot(squeeze(LFP_metric_XC(:,ii,ii)),'o-')
% xl = xlim();
% line(xl,rand_XC(ii,ii)+[0 0],'color','r')
% pause
% clf
% end

%%
% test_trials = find(trialOB == 130 & trialrespDir ~= 0);
% ut = (beg_buff+1):(NT-end_buff);
% 
% resh_LFP_real = reshape(tbt_LFP_real(ut,test_trials,:,:),length(ut)*length(test_trials),length(wfreqs),length(uprobes));
% resh_LFP_imag = reshape(tbt_LFP_imag(ut,test_trials,:,:),length(ut)*length(test_trials),length(wfreqs),length(uprobes));
% resh_LFP_phase = atan2(resh_LFP_imag,resh_LFP_real);
% 
% for cc = 1:Nunits;
% cur_Robs = reshape(squeeze(fullRobs(ut,test_trials,cc)),[],1);
% uset = find(~isnan(cur_Robs) & ~isnan(resh_LFP_phase(:,1,1)));
% 
% cur_spk_inds = convert_to_spikebins(cur_Robs(uset));
% cur_spk_r = squeeze(circ_r(resh_LFP_phase(uset(cur_spk_inds),:,:)));
% all_spk_r(cc,:,:) = cur_spk_r;
% end

%%
test_trials = find(trialOB == 130 & trialrespDir ~= 0);

beg_buff = 0; %number of bins from beginning of trial to exclude
end_buff = 0;

u_lfp_times = find(R_trial_taxis > 0 & R_trial_taxis <= (trial_dur)*dt);
LFP_real = squeeze(tbt_LFP_real(u_lfp_times,test_trials,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,test_trials,:,:));
LFP_amp = sqrt(LFP_imag.^2 + LFP_real.^2);

trial_pow = squeeze(nanmean(LFP_amp));
trial_pow = bsxfun(@minus,trial_pow,nanmean(trial_pow));

Robs = fullRobs((1+beg_buff):(NT-end_buff),test_trials,:);
Robs_tot = squeeze(nanmean(Robs));
Robs_tot = bsxfun(@minus,Robs_tot,nanmean(Robs_tot));

clear unit_tot_avg_amp
for cc = 1:Nunits
   
    cur_Robs = reshape(Robs_tot(:,cc),[],1);
     
    trig_avg_amp = squeeze(nanmean(bsxfun(@times,trial_pow,cur_Robs)));

    unit_tot_avg_amp(cc,:,:) = trig_avg_amp;
end

%%
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
% test_trials = find(trialrespDir ~= 0);

beg_buff = 0; %number of bins from beginning of trial to exclude
end_buff = 0;

u_lfp_times = find(R_trial_taxis > 0 & R_trial_taxis <= (trial_dur)*dt);
LFP_real = squeeze(tbt_LFP_real(u_lfp_times,test_trials,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,test_trials,:,:));
LFP_amp = sqrt(LFP_imag.^2 + LFP_real.^2);

LFP_real = bsxfun(@minus,LFP_real,nanmean(LFP_real,2));
LFP_imag = bsxfun(@minus,LFP_imag,nanmean(LFP_imag,2));
LFP_amp = bsxfun(@minus,LFP_amp,nanmean(LFP_amp,2));

LFP_amp = reshape(LFP_amp,[],length(uprobes)*length(wfreqs));
LFP_real = reshape(LFP_real,[],length(uprobes)*length(wfreqs));
LFP_imag = reshape(LFP_imag,[],length(uprobes)*length(wfreqs));

LFP_real = bsxfun(@rdivide,LFP_real,nanstd(LFP_amp));
LFP_imag = bsxfun(@rdivide,LFP_imag,nanstd(LFP_amp));

Robs = fullRobs((1+beg_buff):(NT-end_buff),test_trials,:);
Robs = reshape(Robs,[],Nunits);
Robs = bsxfun(@minus,Robs,nanmean(Robs));
Robs = reshape(Robs,[],length(test_trials),Nunits);

for cc = 1:Nunits
   
    cur_Robs = reshape(Robs(:,:,cc),[],1);
    
%     uset = find(~isnan(cur_Robs));
%     cur_spkbins = convert_to_spikebins(cur_Robs(uset));
%     trig_avg_real = nanmean(LFP_real(uset(cur_spkbins),:));
%     trig_avg_imag = nanmean(LFP_imag(uset(cur_spkbins),:));
%     trig_avg_amp = nanmean(LFP_amp(uset(cur_spkbins),:));
 
    trig_avg_real = nanmean(bsxfun(@times,LFP_real,cur_Robs));
    trig_avg_imag = nanmean(bsxfun(@times,LFP_imag,cur_Robs));
    trig_avg_amp = nanmean(bsxfun(@times,LFP_amp,cur_Robs));

    unit_trig_avg_real(cc,:) = trig_avg_real;
    unit_trig_avg_imag(cc,:) = trig_avg_imag;
    unit_trig_avg_amp(cc,:) = trig_avg_amp;
end


LFP_real = reshape(LFP_real,[],length(test_trials),length(wfreqs),length(uprobes));
LFP_imag = reshape(LFP_imag,[],length(test_trials),length(wfreqs),length(uprobes));
LFP_amp = reshape(LFP_amp,[],length(test_trials),length(wfreqs),length(uprobes));

un_stim_seeds = unique(trialSe(test_trials));
all_LFP_Aavgs = zeros(Nunits,length(wfreqs),length(uprobes));
all_LFP_ravgs = zeros(Nunits,length(wfreqs),length(uprobes));
all_LFP_iavgs = zeros(Nunits,length(wfreqs),length(uprobes));
all_LFP_cnts = zeros(Nunits,length(wfreqs),length(uprobes));
for tr = 1:length(un_stim_seeds)
    cur_tr_set = find(trialSe(test_trials) == un_stim_seeds(tr));
    
    if length(cur_tr_set) >= 2
        [II,JJ] = meshgrid(1:length(cur_tr_set));
        uset = find(II ~= JJ);
        
        cur_LFP_real = LFP_real(:,cur_tr_set,:,:);
        cur_LFP_imag = LFP_imag(:,cur_tr_set,:,:);
        cur_LFP_amp = LFP_amp(:,cur_tr_set,:,:);
        
        for cc = 1:Nunits
            cur_Robs = Robs(:,cur_tr_set,cc);
            
            cur_real = squeeze(nanmean(bsxfun(@times,cur_LFP_real(:,II(uset),:,:),cur_Robs(:,JJ(uset)))));
            cur_imag = squeeze(nanmean(bsxfun(@times,cur_LFP_imag(:,II(uset),:,:),cur_Robs(:,JJ(uset)))));
            cur_amp = squeeze(nanmean(bsxfun(@times,cur_LFP_amp(:,II(uset),:,:),cur_Robs(:,JJ(uset)))));
            all_LFP_ravgs(cc,:,:) = all_LFP_ravgs(cc,:,:) + nansum(cur_real);
            all_LFP_iavgs(cc,:,:) = all_LFP_iavgs(cc,:,:) + nansum(cur_imag);
            all_LFP_Aavgs(cc,:,:) = all_LFP_Aavgs(cc,:,:) + nansum(cur_amp);
            all_LFP_cnts(cc,:,:) = all_LFP_cnts(cc,:,:) + sum(~isnan(cur_real));
        end
    end
end

all_LFP_ravgs = all_LFP_ravgs./all_LFP_cnts;
all_LFP_iavgs = all_LFP_iavgs./all_LFP_cnts;
all_LFP_Aavgs = all_LFP_Aavgs./all_LFP_cnts;

controlled_ravgs = reshape(unit_trig_avg_real,Nunits,length(wfreqs),length(uprobes)) - all_LFP_ravgs;
controlled_iavgs = reshape(unit_trig_avg_imag,Nunits,length(wfreqs),length(uprobes)) - all_LFP_iavgs;
controlled_Aavgs = reshape(unit_trig_avg_amp,Nunits,length(wfreqs),length(uprobes)) - all_LFP_Aavgs;
uncontrolled_ravgs = reshape(unit_trig_avg_real,Nunits,length(wfreqs),length(uprobes));
uncontrolled_iavgs = reshape(unit_trig_avg_imag,Nunits,length(wfreqs),length(uprobes));
uncontrolled_Aavgs = reshape(unit_trig_avg_amp,Nunits,length(wfreqs),length(uprobes));

%%

%%
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
up_trials = trialrespDir(test_trials) == 1;
down_trials = trialrespDir(test_trials) == -1;

beg_buff = 0; %number of bins from beginning of trial to exclude
end_buff = 0;

u_lfp_times = find(R_trial_taxis > 0 & R_trial_taxis <= (trial_dur)*dt);
LFP_real = squeeze(tbt_LFP_real(u_lfp_times,test_trials,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,test_trials,:,:));
LFP_amp = sqrt(LFP_imag.^2 + LFP_real.^2);

LFP_real = bsxfun(@rdivide,LFP_real,nanstd(LFP_amp));
LFP_imag = bsxfun(@rdivide,LFP_imag,nanstd(LFP_amp));

LFP_amp_upavg = squeeze(nanmean(LFP_amp(:,up_trials,:,:),2));
LFP_amp_downavg = squeeze(nanmean(LFP_amp(:,down_trials,:,:),2));

LFP_real_upavg = squeeze(nanmean(LFP_real(:,up_trials,:,:),2));
LFP_real_downavg = squeeze(nanmean(LFP_real(:,down_trials,:,:),2));
LFP_imag_upavg = squeeze(nanmean(LFP_imag(:,up_trials,:,:),2));
LFP_imag_downavg = squeeze(nanmean(LFP_imag(:,down_trials,:,:),2));

LFP_PC_up = sqrt(LFP_real_upavg.^2 + LFP_imag_upavg.^2);
LFP_PC_down = sqrt(LFP_real_downavg.^2 + LFP_imag_downavg.^2);

%%
uu = find(R_trial_taxis(u_lfp_times) > 0.4);

mean_LFP_amp = squeeze(nanmean(LFP_amp,4));
mean_trial_pow = squeeze(nanmean(mean_LFP_amp(uu,:,:)));

trial_pow = squeeze(nanmean(LFP_amp(uu,:,:,:)));

resp_up = find(trialrespDir(test_trials) == 1);
resp_down = find(trialrespDir(test_trials) == -1);
for cc = 1:length(uprobes)
    for ww = 1:length(wfreqs)
        cur_resp_ax = prctile(trial_pow(:,ww,cc),[0:5:100]);
        up_hist = histc(trial_pow(resp_up,ww,cc),cur_resp_ax);
        down_hist = histc(trial_pow(resp_down,ww,cc),cur_resp_ax);
        true_pos = cumsum(up_hist)/sum(~isnan(trial_pow(resp_up,ww,cc)));
        false_pos = cumsum(down_hist)/sum(~isnan(trial_pow(resp_down,ww,cc)));
        all_LFPpow_CP(ww,cc) = trapz(false_pos,true_pos);
    end
end

resp_up = find(trialrespDir(test_trials) == 1);
resp_down = find(trialrespDir(test_trials) == -1);
for ww = 1:length(wfreqs)
    cur_resp_ax = prctile(mean_trial_pow(:,ww),[0:5:100]);
    up_hist = histc(mean_trial_pow(resp_up,ww),cur_resp_ax);
    down_hist = histc(mean_trial_pow(resp_down,ww),cur_resp_ax);
    true_pos = cumsum(up_hist)/sum(~isnan(mean_trial_pow(resp_up,ww)));
    false_pos = cumsum(down_hist)/sum(~isnan(mean_trial_pow(resp_down,ww)));
    LFPpow_CP(ww) = trapz(false_pos,true_pos);
end

n_boot_samps = 1e3;
for bb = 1:n_boot_samps
    fprintf('Boot %d of %d\n',bb,n_boot_samps);
    randrespdir = round(rand(length(test_trials),1));
    resp_up = find(randrespdir == 1);
    resp_down = find(randrespdir == 0);
    for ww = 1:length(wfreqs)
        cur_resp_ax = prctile(mean_trial_pow(:,ww),[0:5:100]);
        up_hist = histc(mean_trial_pow(resp_up,ww),cur_resp_ax);
        down_hist = histc(mean_trial_pow(resp_down,ww),cur_resp_ax);
        true_pos = cumsum(up_hist)/sum(~isnan(mean_trial_pow(resp_up,ww)));
        false_pos = cumsum(down_hist)/sum(~isnan(mean_trial_pow(resp_down,ww)));
        rand_LFPpow_CP(ww,bb) = trapz(false_pos,true_pos);
    end
end

LFPpow_CI = prctile(rand_LFPpow_CP',[5 95]);

%%
fig_dir = '/home/james/Desktop/CPfigs/';
close all
fig_width = 4;
rel_height = 0.8;

wfreqs_ls = linspace(min(wfreqs),max(wfreqs),1000);
all_LFPpow_CP_ls = interp1(fliplr(wfreqs),flipud(all_LFPpow_CP),wfreqs_ls);

f1 = figure();
% imagesc(wfreqs_ls,1:24,1-all_LFPpow_CP_ls');
% imagesc(wfreqs,1:24,1-all_LFPpow_CP');
pcolor(wfreqs,1:24,1-all_LFPpow_CP');shading interp
xlim([1.2 100]);
colorbar
% set(gca,'xticklabel',fliplr(wfreqs));
caxis([0.35 0.65]);
xlabel('Frequency (Hz)');
ylabel('Depth');
set(gca,'xscale','log');
figufy(f1);
fname = [fig_dir sprintf('LFPpow_CP_%s.pdf',Expt_name)];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

ucells = find(totSpkCnts > 1e4);
avg_controlled = squeeze(nanmean(controlled_Aavgs(ucells,:,:)));

f2 = figure();
% imagescnan(wfreqs,1:24,avg_controlled');shading flat
pcolor(wfreqs,1:24,avg_controlled');shading interp
xlim([1.2 100]);
caxis([-0.005 0.005])
set(gca,'xscale','log');
colorbar
xlabel('Frequency (Hz)');
ylabel('Depth');
figufy(f2);
fname = [fig_dir 'Gsac_GRIM_rates.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

%%
poss_trials = find(trialOB == 130 & trialrespDir ~= 0);
un_stim_seeds = unique(trialSe(poss_trials));

beg_buff = 0; %number of bins from beginning of trial to exclude
end_buff = 0;

u_lfp_times = find(R_trial_taxis > 0 & R_trial_taxis <= (trial_dur)*dt);
LFP_real = squeeze(tbt_LFP_real(u_lfp_times,test_trials,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,test_trials,:,:));
LFP_comp = LFP_real + sqrt(-1)*LFP_imag;
LFP_comp = bsxfun(@minus,LFP_comp,nanmean(LFP_comp,2));
LFP_Acomp = abs(LFP_comp);
LFP_Acomp = bsxfun(@minus,LFP_Acomp,nanmean(LFP_Acomp,2));

for cc = 1:8;
cc
all_cwt_same = zeros(length(u_lfp_times),length(wfreqs));
all_cwt_diff = all_cwt_same;
all_same_cnt = all_cwt_same;
all_diff_cnt = all_cwt_same;
all_Acwt_same = zeros(length(u_lfp_times),length(wfreqs));
all_Acwt_diff = all_Acwt_same;
all_A2cwt_same = zeros(length(u_lfp_times),length(wfreqs));
all_A2cwt_diff = all_A2cwt_same;
for ss = 1:length(un_stim_seeds)
    cur_trials = find(trialSe(poss_trials) == un_stim_seeds(ss));
    if length(cur_trials) >= 2
        cur_resp = trialrespDir(poss_trials(cur_trials));
        [II,JJ] = meshgrid(1:length(cur_trials));
        cur_Y = squeeze(LFP_comp(:,cur_trials,:,cc));
        cur_AY = squeeze(LFP_Acomp(:,cur_trials,:,cc));
        
        cur_Dmat = abs(squareform(pdist(cur_resp')));
        cur_Dmat(JJ == II) = nan;
        
        uset = find(cur_Dmat == 0);
        temp = cur_Y(:,II(uset),:).*conj(cur_Y(:,JJ(uset),:));
%         temp = cur_Y(:,II(uset),:).*(cur_Y(:,JJ(uset),:));
        tempm = squeeze(nanmean(temp,2));
        all_cwt_same(~isnan(tempm)) = all_cwt_same(~isnan(tempm))  + tempm(~isnan(tempm)) ;
        temp = abs(cur_Y(:,II(uset),:)).*abs(cur_Y(:,JJ(uset),:));
         tempm = squeeze(nanmean(temp,2));
        all_Acwt_same(~isnan(tempm)) = all_Acwt_same(~isnan(tempm))  + tempm(~isnan(tempm)) ;
        temp = cur_AY(:,II(uset),:).*cur_AY(:,JJ(uset),:);
         tempm = squeeze(nanmean(temp,2));
        all_A2cwt_same(~isnan(tempm)) = all_A2cwt_same(~isnan(tempm))  + tempm(~isnan(tempm)) ;
        all_same_cnt = all_same_cnt + squeeze(sum(~isnan(temp),2));
        
%         uset = find(cur_Dmat ~= 0);
        uset = find(~isnan(cur_Dmat));
        temp = cur_Y(:,II(uset),:).*conj(cur_Y(:,JJ(uset),:));
%         temp = cur_Y(:,II(uset),:).*(cur_Y(:,JJ(uset),:));
        tempm = squeeze(nanmean(temp,2));
        all_cwt_diff(~isnan(tempm)) = all_cwt_diff(~isnan(tempm)) + tempm(~isnan(tempm));
        temp = abs(cur_Y(:,II(uset),:)).*abs(cur_Y(:,JJ(uset),:));
        tempm = squeeze(nanmean(temp,2));
        all_Acwt_diff(~isnan(tempm)) = all_Acwt_diff(~isnan(tempm)) + tempm(~isnan(tempm));
        temp = cur_AY(:,II(uset),:).*cur_AY(:,JJ(uset),:);
        tempm = squeeze(nanmean(temp,2));
        all_A2cwt_diff(~isnan(tempm)) = all_A2cwt_diff(~isnan(tempm)) + tempm(~isnan(tempm));
        all_diff_cnt = all_diff_cnt + squeeze(sum(~isnan(temp),2));
    end
    
end

all_cwt_same = abs(all_cwt_same)./all_same_cnt;
all_cwt_diff = abs(all_cwt_diff)./all_diff_cnt;
all_Acwt_same = all_Acwt_same./all_same_cnt;
all_Acwt_diff = all_Acwt_diff./all_diff_cnt;
all_A2cwt_same = all_A2cwt_same./all_same_cnt;
all_A2cwt_diff = all_A2cwt_diff./all_diff_cnt;

% all_cwt_choice = all_cwt_same - all_cwt_diff;
% norm_fac = squeeze(nanvar(abs(tbt_LFP_cwt(:,poss_trials,:,cc)),[],2));
% all_cwt_corr = all_cwt_choice./norm_fac;
all_cwt_same_norm = all_cwt_same./all_Acwt_same;
all_cwt_choice = all_cwt_same - all_cwt_diff;
all_cwt_choice_norm(cc,:,:) = all_cwt_choice./(0.5*all_Acwt_same + 0.5*all_Acwt_diff);
all_A2cwt_choice(cc,:,:) = (all_A2cwt_same - all_A2cwt_diff)./(0.5*all_Acwt_same + 0.5*all_Acwt_diff);
end

%%
LFP_dsf = 4;
LFP_Fsd = LFP_Fs/LFP_dsf;
%anti-aliasing filter and high-pass filter
aa_hcf = LFP_Fsd/2*0.8;
% [b_aa,a_aa] = butter(4,aa_hcf/(LFP_Fs/2),'low');
aa_lcf = 1;
[b_aa,a_aa] = butter(2,[aa_lcf aa_hcf]/(LFP_Fs/2));

LFP_trial_taxis_ds = downsample(LFP_trial_taxis,LFP_dsf);

nprobes = 24;
uprobes = 1:1:nprobes;


beg_buffer = -round(0.2/dt);
end_buffer = round(0/dt);
trial_dur = round(2/dt);
R_trial_taxis = (beg_buffer:(trial_dur - end_buffer))*dt;
TLEN = length(R_trial_taxis);

all_spk_id = [];
all_spk_phases = nan(sum(tot_spks_per_trial(:)),length(wfreqs),length(uprobes));
trial_LFPs = nan(Ntrials,TLEN,length(uprobes));
for tr = 1:Ntrials
    fprintf('Trial %d of %d\n',tr,Ntrials);
    cur_LFPs = double(AllExpt.Expt.Trials(tr).LFP(:,uprobes));
    bad_LFPs = isnan(cur_LFPs(:,1));
    cur_LFPs(isnan(cur_LFPs)) = 0;
    cur_LFPs = filtfilt(b_aa,a_aa,cur_LFPs);
    cur_LFPs = downsample(cur_LFPs,LFP_dsf);
    bad_LFPs = downsample(bad_LFPs,LFP_dsf);

    cur_LFPs(bad_LFPs,:) = nan;
    trial_LFPs(tr,:,:) = interp1(LFP_trial_taxis_ds,cur_LFPs,R_trial_taxis);
end
trial_LFPs = permute(trial_LFPs,[2 1 3]);


%%
addpath(genpath('~/James_scripts/iCSD/'))

cur_LFPs = permute(trial_LFPs,[3 1 2]);
cur_LFPs = cur_LFPs(:,1:150,:);
csd_params.Fs = LFP_Fsd; %sample freq
csd_params.BrainBound = 1; %first channel that is in the brain
csd_params.ChanSep = 0.05; %channel sep in mm
csd_params.diam = 2; %current disc diameter (in mm)
csd_method = 'spline';

csd_mat = PettersenCSD(cur_LFPs,csd_method,csd_params);
avg_csd = squeeze(nanmean(csd_mat,3));
csd_tax = R_trial_taxis(1:150);


