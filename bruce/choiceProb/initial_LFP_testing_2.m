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
dt = 0.01;  % assuming 10 ms frame time exactly
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

%wavelet parameters
nwfreqs = 20;
min_freq = 1.; max_freq = 50;
% nwfreqs = 5;
% min_freq = 1; max_freq = 6;
min_scale = 1/max_freq*LFP_Fsd;
max_scale = 1/min_freq*LFP_Fsd;
wavetype = 'cmor1-1';
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,wavetype,1/LFP_Fsd);

% beg_buffer = round(0.15/dt);
beg_buffer = round(-0.2/dt);
end_buffer = round(-0.2/dt);
trial_dur = round(2/dt);
R_trial_taxis = (beg_buffer:(trial_dur - end_buffer))*dt;
TLEN = length(R_trial_taxis);

nprobes = 24;
uprobes = 1:3:nprobes;
all_spk_id = [];
all_spk_phases = nan(sum(tot_spks_per_trial(:)),length(wfreqs),length(uprobes));
trial_LFP_cwt = nan(Ntrials,TLEN,length(wfreqs),length(uprobes));
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
    trial_LFP_cwt(tr,:,:,:) = interp1(LFP_trial_taxis_ds,cur_cwt,R_trial_taxis);
    trial_LFPs(tr,:,:) = interp1(LFP_trial_taxis_ds,cur_LFPs,R_trial_taxis);
end
trial_LFP_cwt = permute(trial_LFP_cwt,[2 1 3 4]);
trial_LFPs = permute(trial_LFPs,[2 1 3]);

trial_LFP_cwt = reshape(trial_LFP_cwt,[],length(wfreqs)*length(uprobes));
trial_LFPs = reshape(trial_LFPs,[],length(uprobes));

%%
trial_LFPs = nanzscore(trial_LFPs);
trial_LFP_cwt = bsxfun(@rdivide,trial_LFP_cwt,nanstd(abs(trial_LFP_cwt)));


%%
tbt_LFPs = reshape(trial_LFPs,[TLEN Ntrials length(uprobes)]);
tbt_LFP_cwt = reshape(trial_LFP_cwt,[TLEN Ntrials length(wfreqs) length(uprobes)]);

%%
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
un_stim_seeds = unique(trialSe(test_trials));

beg_buff = 0; %number of bins from beginning of trial to exclude
end_buff = 0;

maxlag_ED = 5;
ED_space = 0.5;
ED_bin_edges = 0:ED_space:maxlag_ED;
ED_bin_centers = (ED_bin_edges(1:end-1)+ED_bin_edges(2:end))/2;

%subtract off avg rates for times within this set of trials
[RR,TT] = meshgrid(1:Ntrials,1:NT);
istrain = (ismember(RR(:),test_trials)) & (TT(:) > beg_buff);
fullRobs_resh = reshape(fullRobs,[],Nunits);
fullRobs_ms = bsxfun(@minus,fullRobs,reshape(nanmean(fullRobs_resh(istrain,:)),[1 1 Nunits]));
tot_var = nanvar(fullRobs_resh(istrain,:));
fullRobs_ms = fullRobs_ms(:,test_trials,:);

ww = 16;

u_lfp_times = find(R_trial_taxis > beg_buff & R_trial_taxis <= dt*(trial_dur - end_buff));
tbt_LFP_data = squeeze(tbt_LFP_cwt(u_lfp_times,test_trials,ww,:));

LFP_metric_XC = zeros(length(ED_bin_centers),Nunits,Nunits);
LFP_metric_cnt = zeros(length(ED_bin_centers),Nunits,Nunits);
rand_XC = zeros(Nunits,Nunits);
rand_cnt = zeros(Nunits,Nunits);
for tr = 1:length(un_stim_seeds)
    
    cur_tr_set = find(trialSe(test_trials) == un_stim_seeds(tr));
    fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
    [II,JJ] = meshgrid(1:length(cur_tr_set));
    
    
    for tt = (beg_buff+1):NT
        cur_Robs = squeeze(fullRobs_ms(tt,cur_tr_set,:));
        cur_Robs2 = reshape(cur_Robs,[],1,Nunits);
        
        cur_LFP_state = squeeze(tbt_LFP_data(tt,cur_tr_set,:));
        cur_LFP_state = cat(2,real(cur_LFP_state),imag(cur_LFP_state));
        
        cur_Dmat = squareform(pdist(cur_LFP_state))/sqrt(length(uprobes));
        cur_Dmat(JJ >= II) = nan;
        for jj = 1:length(ED_bin_centers)
            curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
            LFP_metric_XC(jj,:,:) = LFP_metric_XC(jj,:,:) + nansum(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1);
            %         cur_XC(tt,jj,:,:) = (nanmean(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
        end
        curset = ~isnan(cur_Dmat);
        rand_XC = rand_XC + squeeze(nansum(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
        %     rand_XC(tt,:,:) = squeeze(nanmean(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
    end
end

