%% LOAD PROCESSED DATA
clear all
% cd('/home/james/Data/bruce/ChoiceProb/')
cd('~/Data/bruce/ChoiceProb/')

Expt_name = 'G037';

% if strcmp(Expt_name,'M239')
%     load('M239/lemM239.image.ORBW.LFP.mat')
% elseif strcmp(Expt_name,'M230')
%     load('M230/lemM230.image.ORBW.LFP.mat')
% end
load('G037/jbeG037.Cells.image.ORBW.LFP.mat');

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
elseif strcmp(Expt_name,'G037')
    block_trial_boundaries = [1 148;
        149 454;
        455 1361];
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
spk_delay_buff = 0.05;
for tr = 1:Ntrials
    for nn = 1:Nunits
        trindx = find( AllExpt.Spikes{nn}.Trial == trialNums(tr));
        if ~isempty(trindx)
            tot_spks_per_trial(tr,nn) = sum(AllExpt.Spikes{nn}.Spikes{trindx}*1e-4 <= (trialDur + spk_delay_buff) & ...
                (AllExpt.Spikes{nn}.Spikes{trindx}*1e-4 >= spk_delay_buff));
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

%%
un_OB = unique(trialOB);
psths = nan(NT,Nunits,length(un_OB));
for ii = 1:length(un_OB)
   cur_trials = find(trialOB == un_OB(ii));
   psths(:,:,ii) = squeeze(nanmean(fullRobs(:,cur_trials,:),2));
end

close all
sm_win = 2;
t_axis = (1:NT)*dt;
for ii = 1:Nunits
   plot(t_axis,jmm_smooth_1d_cor(squeeze(psths(:,ii,3)),sm_win));
   hold on
   plot(t_axis,jmm_smooth_1d_cor(squeeze(psths(:,ii,4)),sm_win),'r');
   ii
   pause
   clf
end
%% CALCULATE CHOICE PROBS
tot_spks_per_trial_norm = tot_spks_per_trial;
tot_spks_per_trial_norm(tot_spks_per_trial == 0) = nan;

sig_trials = find(abs(trialOB) == 60);
stim_up = sig_trials(trialrwDir(sig_trials)==1);
stim_down = sig_trials(trialrwDir(sig_trials)==-1);
sig_prob = nan(Nunits,1);
sp_pval = nan(Nunits,1);
for cc = 1:Nunits
    cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(sig_trials,cc)));
    up_hist = histc(tot_spks_per_trial_norm(stim_up,cc),cur_resp_ax);
    down_hist = histc(tot_spks_per_trial_norm(stim_down,cc),cur_resp_ax);
    fract_cor = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_up,cc)));
    fract_incor = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_down,cc)));
    sig_prob(cc) = trapz(fract_incor,fract_cor);
    [~,sp_pval(cc)] = ttest2(tot_spks_per_trial_norm(stim_up,cc),tot_spks_per_trial_norm(stim_down,cc));
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
    [~,~,~,CP_Neuron(cc)]=perfcurve(trialrespDir(test_trials), tot_spks_per_trial_norm(test_trials,cc),1);
    [~,cp_pval(cc)] = ttest2(tot_spks_per_trial_norm(resp_up,cc),tot_spks_per_trial_norm(resp_down,cc));
end

nboot = 1000;
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
choice_prob_ci = prctile(boot_choice_prob',[2.5 97.5]);
sig_prob_ci = prctile(boot_sig_prob',[2.5 97.5]);
choice_prob_z = (choice_prob - nanmean(boot_choice_prob,2))./nanstd(boot_choice_prob,[],2);
sig_prob_z = (sig_prob - nanmean(boot_sig_prob,2))./nanstd(boot_sig_prob,[],2);

boot_sig_sig = find(sig_prob > sig_prob_ci(2,:)' | sig_prob < sig_prob_ci(1,:)');
boot_choice_sig = find(choice_prob > choice_prob_ci(2,:)' | choice_prob < choice_prob_ci(1,:)');

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
% LFP_trial_taxis = AllExpt.Expt.Header.LFPtimes*1e-4;
LFP_Fs = 1/AllExpt.Expt.Header.LFPsamplerate;
LFP_offset = AllExpt.Expt.Header.preperiod/1e4;
LFP_trial_taxis = (1:length(AllExpt.Expt.Header.LFPtimes))/LFP_Fs - LFP_offset;

LFP_dsf = 2;
LFP_Fsd = LFP_Fs/LFP_dsf;
%anti-aliasing filter and high-pass filter
aa_hcf = LFP_Fsd/2*0.8;
% [b_aa,a_aa] = butter(4,aa_hcf/(LFP_Fs/2),'low');
aa_lcf = 0.5;
[b_aa,a_aa] = butter(2,[aa_lcf aa_hcf]/(LFP_Fs/2));

LFP_trial_taxis_ds = downsample(LFP_trial_taxis,LFP_dsf);

nprobes = 96;
uprobes = 1:2:nprobes;

%wavelet parameters
nwfreqs = 25;
min_freq = 1.5; max_freq = 100;
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
% all_spk_phases = nan(sum(tot_spks_per_trial(:)),length(wfreqs),length(uprobes));
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

%% compute power-based spike-lfp coupling directly (simultaneous trials)
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

%% COMPUTE SPIKE-LFP COUPLING WITH AND WITHOUT CONTROLLING FOR STIM_DEPENDENCE
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
test_trials = find(trialOB == 130 & trialrespDir ~= 0);

up_trials = find(trialrespDir(test_trials) == 1);
down_trials = find(trialrespDir(test_trials) == -1);

% train_trials = find(trialOB ~= 130);
train_trials = find(trialOB == 70 | trialOB==80);

beg_buff = 10; %number of bins from beginning of trial to exclude
end_buff = 0;
% beg_buff = 0; %number of bins from beginning of trial to exclude
% end_buff = 0;

%range of times over which to compute lfp power
lfp_bt = 0.4;
lfp_et = 0.2;
u_lfp_times = find(R_trial_taxis > 0 & R_trial_taxis <= (trial_dur)*dt);
u_lfp_times = u_lfp_times(R_trial_taxis(u_lfp_times) >= lfp_bt & R_trial_taxis(u_lfp_times) <= (trial_dur-lfp_et));
LFP_real = squeeze(tbt_LFP_real(u_lfp_times,:,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,:,:,:));
LFP_amp = sqrt(LFP_imag.^2 + LFP_real.^2);
trial_pow = squeeze(nanmean(LFP_amp));
% trial_pow = bsxfun(@minus,trial_pow,nanmean(trial_pow));
trial_pow = bsxfun(@rdivide,trial_pow,nanstd(trial_pow));

Robs = fullRobs((1+beg_buff):(NT-end_buff),:,:);
trial_Robs = squeeze(sum(Robs));
% trial_Robs = bsxfun(@minus,trial_Robs,nanmean(trial_Robs));

%% DETERMINE THE CP PREDICTION FROM POWER AT EACH INDIVIDUAL FREQ
uchs = 1:1:length(uprobes);
reg_params = NMMcreate_reg_params('lambda_L2',5);
stim_params = NMMcreate_stim_params(length(uchs));
silent = 1;

trial_pow_train = trial_pow(train_trials,:,uchs);
trial_pow_test = trial_pow(test_trials,:,uchs);
trial_Robs_train = trial_Robs(train_trials,:);
trial_Robs_test = trial_Robs(test_trials,:);
clear mod_params
for cc = 1:Nunits
    cc
    for ww = 1:length(wfreqs)
        cur_X = squeeze(trial_pow_train(:,ww,:));
        cur_Y = trial_Robs_train(:,cc);
        uindx = find(~isnan(cur_Y) & ~any(isnan(cur_X),2));
        
        init_mod = NMMinitialize_model(stim_params,1,{'lin'},reg_params);
        init_mod = NMMfit_filters(init_mod,cur_Y,cur_X,[],uindx,silent);
%         B = glmfit(cur_X,trial_Robs_train(:,cc),'poisson');
        
        mod_params(cc,ww,:) = init_mod.mods(1).filtK;

        cur_X = squeeze(trial_pow_test(:,ww,:));
        cur_Y = trial_Robs_test(:,cc);
        
        [~,~,Yhat] = NMMeval_model(init_mod,cur_Y,cur_X);
%         Yhat = glmval(B,cur_X,'log');
        Yhat(isnan(cur_Y)) = nan;
        
        cur_resp_ax = prctile(Yhat,[0:5:100]);
        up_hist = histc(Yhat(up_trials),cur_resp_ax);
        down_hist = histc(Yhat(down_trials),cur_resp_ax);
        true_pos = cumsum(up_hist)/sum(~isnan(Yhat(up_trials)));
        false_pos = cumsum(down_hist)/sum(~isnan(Yhat(down_trials)));
        all_LFPpow_pred_CP(cc,ww) = trapz(false_pos,true_pos);

    end
    
    test_R = trial_Robs_test(:,cc);
    cur_resp_ax = 0:(nanmax(test_R));
    up_hist = histc(test_R(up_trials),cur_resp_ax);
    down_hist = histc(test_R(down_trials),cur_resp_ax);
    true_pos = cumsum(up_hist)/sum(~isnan(test_R(up_trials)));
    false_pos = cumsum(down_hist)/sum(~isnan(test_R(down_trials)));
    new_choice_prob(cc) = trapz(false_pos,true_pos);
   
end

% for ww = 1:length(wfreqs)
% [tempc(ww),tempp(ww)] = corr(choice_prob(ucells),all_LFPpow_pred_CP(ucells,ww),'type','spearman');
% end

%% COMPUTE CP PREDICTIONS USING POWER ACROSS A RANGE OF (LOWER) FREQS
ufreqs = find(wfreqs >= 1.5 & wfreqs <= 10);
uchs = 1:1:length(uprobes);
reg_params = NMMcreate_reg_params('lambda_L2',1,'lambda_d2XT',1);
stim_params = NMMcreate_stim_params([length(ufreqs) length(uchs)]);
silent = 1;

for cc = 1:Nunits
        cur_X = squeeze(trial_pow_train(:,ufreqs,:));
        cur_Y = trial_Robs_train(:,cc);
        uindx = find(~isnan(cur_Y) & ~any(isnan(reshape(cur_X,length(train_trials),[])),2));
    
        init_mod = NMMinitialize_model(stim_params,1,{'lin'},reg_params);
        init_mod = NMMfit_filters(init_mod,cur_Y,cur_X,[],uindx,silent);
        
        cur_X = squeeze(trial_pow_test(:,ufreqs,:));
        cur_X = reshape(cur_X,length(test_trials),[]);
        cur_Y = trial_Robs_test(:,cc);
        [~,~,Yhat] = NMMeval_model(init_mod,cur_Y,cur_X);
        Yhat(isnan(cur_Y)) = nan;
        
        Y_resid = cur_Y - Yhat;
                
        cur_resp_ax = prctile(Yhat,[0:5:100]);
        up_hist = histc(Yhat(up_trials),cur_resp_ax);
        down_hist = histc(Yhat(down_trials),cur_resp_ax);
        true_pos = cumsum(up_hist)/sum(~isnan(Yhat(up_trials)));
        false_pos = cumsum(down_hist)/sum(~isnan(Yhat(down_trials)));
        all_LFPpow_allw_pred_CP(cc) = trapz(false_pos,true_pos);

        cur_resp_ax = prctile(Y_resid,[0:5:100]);
        up_hist = histc(Y_resid(up_trials),cur_resp_ax);
        down_hist = histc(Y_resid(down_trials),cur_resp_ax);
        true_pos = cumsum(up_hist)/sum(~isnan(Y_resid(up_trials)));
        false_pos = cumsum(down_hist)/sum(~isnan(Y_resid(down_trials)));
        all_LFPpow_allw_resid_CP(cc) = trapz(false_pos,true_pos);

%         cur_resp_ax = prctile(cur_Y,[0:5:100]);
        cur_resp_ax = min(cur_Y):max(cur_Y);
        up_hist = histc(cur_Y(up_trials),cur_resp_ax);
        down_hist = histc(cur_Y(down_trials),cur_resp_ax);
        true_pos = cumsum(up_hist)/sum(~isnan(cur_Y(up_trials)));
        false_pos = cumsum(down_hist)/sum(~isnan(cur_Y(down_trials)));
        all_LFPpow_allw_actual_CP(cc) = trapz(false_pos,true_pos);
        
        cur_Y = cur_Y - nanmean(cur_Y);
        Yhat = Yhat - nanmean(Yhat);
        Y_resid = Y_resid - nanmean(Y_resid);
        
        act_choice_related_var(cc) = 0.5*nanmean(cur_Y(up_trials)).^2 + 0.5*nanmean(cur_Y(down_trials)).^2;
        pred_choice_related_var(cc) = 0.5*nanmean(Yhat(up_trials)).^2 + 0.5*nanmean(Yhat(down_trials)).^2;
        resid_choice_related_var(cc) = 0.5*nanmean(Y_resid(up_trials)).^2 + 0.5*nanmean(Y_resid(down_trials)).^2;
        
        total_var(cc) = nanvar(cur_Y);

end

% [a,b] = corr(all_LFPpow_allw_pred_CP(ucells)',new_choice_prob(ucells)','type','spearman')
%% COMPUTE CHOICE_CONDITIONAL TRIALAVG POW SPECTRA

test_trials = find(trialOB == 130 & trialrespDir ~= 0);
sig_trials = find(abs(trialOB) == 60 & trialrespDir ~= 0);
cup_trials = trialrespDir(test_trials) == 1;
cdown_trials = trialrespDir(test_trials) == -1;
sup_trials = trialrwDir(sig_trials) == 1;
sdown_trials = trialrwDir(sig_trials) == -1;

bt = 0;
et = 0;

u_lfp_times = find(R_trial_taxis > bt & R_trial_taxis <= (trial_dur)*dt-et);
LFP_real = squeeze(tbt_LFP_real(u_lfp_times,:,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,:,:,:));
LFP_amp = sqrt(LFP_imag.^2 + LFP_real.^2);
% amp_SDs = nanstd(reshape(LFP_amp,[],length(wfreqs),length(uprobes)));
% LFP_amp = bsxfun(@rdivide,LFP_amp,reshape(amp_SDs,1,1,length(wfreqs),length(uprobes)));

LFP_amp_upavg = squeeze(nanmean(LFP_amp(:,test_trials(cup_trials),:,:),2));
LFP_amp_downavg = squeeze(nanmean(LFP_amp(:,test_trials(cdown_trials),:,:),2));
LFP_choice_powdiff = LFP_amp_upavg - LFP_amp_downavg;

LFP_amp_upavg = squeeze(nanmean(LFP_amp(:,sig_trials(sup_trials),:,:),2));
LFP_amp_downavg = squeeze(nanmean(LFP_amp(:,sig_trials(sdown_trials),:,:),2));
LFP_sig_powdiff = LFP_amp_upavg - LFP_amp_downavg;

avg_choice_specdiff = squeeze(nanmean(LFP_choice_powdiff));
avg_sig_specdiff = squeeze(nanmean(LFP_sig_powdiff));

raw_LFP_upavg = squeeze(nanmean(tbt_LFPs(u_lfp_times,test_trials(cup_trials),:),2));
raw_LFP_downavg = squeeze(nanmean(tbt_LFPs(u_lfp_times,test_trials(cdown_trials),:),2));
raw_LFP_choicediff = raw_LFP_upavg - raw_LFP_downavg;

raw_LFP_upavg = squeeze(nanmean(tbt_LFPs(u_lfp_times,sig_trials(sup_trials),:),2));
raw_LFP_downavg = squeeze(nanmean(tbt_LFPs(u_lfp_times,sig_trials(sdown_trials),:),2));
raw_LFP_sigdiff = raw_LFP_upavg - raw_LFP_downavg;

nboots = 500;
boot_samps = nan(nboots,length(wfreqs));
boot_amp_samps = nan(nboots,length(wfreqs));
for ii = 1:nboots
    ii
rr = rand(length(test_trials),1) > 0.5;
% raw_LFP_rupavg = squeeze(nanmean(tbt_LFPs(u_lfp_times,test_trials(rr),:),2));
% raw_LFP_rdownavg = squeeze(nanmean(tbt_LFPs(u_lfp_times,test_trials(~rr),:),2));
% raw_LFP_rchoicediff = raw_LFP_rupavg - raw_LFP_rdownavg;
LFP_amp_upavg = squeeze(nanmean(LFP_amp(:,test_trials(rr),:,:),2));
LFP_amp_downavg = squeeze(nanmean(LFP_amp(:,test_trials(~rr),:,:),2));
LFP_rchoice_powdiff = LFP_amp_upavg - LFP_amp_downavg;
avg_rchoice_specdiff = squeeze(nanmean(LFP_rchoice_powdiff));
boot_samps(ii,:) = nanmean(avg_rchoice_specdiff,2);
boot_amp_samps(ii,:) = nanmean(abs(avg_rchoice_specdiff),2);
end
boot_CI = prctile(boot_samps,[2.5 97.5]);
amp_boot_CI = prctile(boot_amp_samps,[2.5 97.5]);

nboots = 500;
rboot_samps = nan(nboots,length(wfreqs));
rboot_amp_samps = nan(nboots,length(wfreqs));
for ii = 1:nboots
    ii
rr = randperm(length(sig_trials));
LFP_amp_upavg = squeeze(nanmean(LFP_amp(:,sig_trials(rr(1:sum(sup_trials))),:,:),2));
LFP_amp_downavg = squeeze(nanmean(LFP_amp(:,sig_trials(rr((sum(sup_trials)+1):end)),:,:),2));
LFP_rchoice_powdiff = LFP_amp_upavg - LFP_amp_downavg;
avg_rchoice_specdiff = squeeze(nanmean(LFP_rchoice_powdiff));
rboot_samps(ii,:) = nanmean(avg_rchoice_specdiff,2);
rboot_amp_samps(ii,:) = nanmean(abs(avg_rchoice_specdiff),2);
end
sig_boot_CI = prctile(rboot_samps,[2.5 97.5]);
sig_amp_boot_CI = prctile(rboot_amp_samps,[2.5 97.5]);

%%
close all
for cc = 1:length(uprobes)
   cc
%    subplot(2,1,1);
%    pcolor(R_trial_taxis(u_lfp_times),wfreqs,squeeze(LFP_choice_powdiff(:,:,cc))');shading flat
%     caxis([-0.4 0.4])
%     set(gca,'yscale','log')
%    subplot(2,1,2);
%    pcolor(R_trial_taxis(u_lfp_times),wfreqs,squeeze(LFP_sig_powdiff(:,:,cc))');shading flat
%     caxis([-0.4 0.4])
%     set(gca,'yscale','log')
    plot(wfreqs,avg_choice_specdiff(:,cc),wfreqs,avg_sig_specdiff(:,cc),'r')
    set(gca,'xscale','log'); axis tight
    ylim([-0.2 0.2]);
    xl = xlim(); line(xl,[0 0],'color','k');
     pause
   clf
end

%% COMPUTE CP OF TRIAL_AVG POW SPECTRA
lfp_bt = 0.3;
lfp_et = 0.3;
u_lfp_times = find(R_trial_taxis > 0 & R_trial_taxis <= (trial_dur)*dt);
u_lfp_times = u_lfp_times(R_trial_taxis(u_lfp_times) >= lfp_bt & R_trial_taxis(u_lfp_times) <= (trial_dur-lfp_et));
test_trials = find(trialOB == 130 & trialrespDir ~= 0);

LFP_real = squeeze(tbt_LFP_real(u_lfp_times,test_trials,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,test_trials,:,:));
LFP_amp = sqrt(LFP_imag.^2 + LFP_real.^2);

mean_LFP_amp = squeeze(nanmean(LFP_amp,4));
mean_trial_pow = squeeze(nanmean(mean_LFP_amp));

trial_pow = squeeze(nanmean(LFP_amp));

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

% n_boot_samps = 1e3;
% for bb = 1:n_boot_samps
%     fprintf('Boot %d of %d\n',bb,n_boot_samps);
%     randrespdir = round(rand(length(test_trials),1));
%     resp_up = find(randrespdir == 1);
%     resp_down = find(randrespdir == 0);
%     for ww = 1:length(wfreqs)
%         cur_resp_ax = prctile(mean_trial_pow(:,ww),[0:5:100]);
%         up_hist = histc(mean_trial_pow(resp_up,ww),cur_resp_ax);
%         down_hist = histc(mean_trial_pow(resp_down,ww),cur_resp_ax);
%         true_pos = cumsum(up_hist)/sum(~isnan(mean_trial_pow(resp_up,ww)));
%         false_pos = cumsum(down_hist)/sum(~isnan(mean_trial_pow(resp_down,ww)));
%         rand_LFPpow_CP(ww,bb) = trapz(false_pos,true_pos);
%     end
% end
% 
% LFPpow_CI = prctile(rand_LFPpow_CP',[5 95]);

%%


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

for cc = 1:length(uprobes);
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
            tempm = squeeze(nanmean(temp,2));
            all_cwt_same(~isnan(tempm)) = all_cwt_same(~isnan(tempm))  + tempm(~isnan(tempm)) ;
            
            temp = abs(cur_Y(:,II(uset),:)).*abs(cur_Y(:,JJ(uset),:));
            tempm = squeeze(nanmean(temp,2));
            all_Acwt_same(~isnan(tempm)) = all_Acwt_same(~isnan(tempm))  + tempm(~isnan(tempm)) ;
            
            temp = cur_AY(:,II(uset),:).*cur_AY(:,JJ(uset),:);
            tempm = squeeze(nanmean(temp,2));
            all_A2cwt_same(~isnan(tempm)) = all_A2cwt_same(~isnan(tempm))  + tempm(~isnan(tempm)) ;
            all_same_cnt = all_same_cnt + squeeze(sum(~isnan(temp),2));
            
            uset = find(~isnan(cur_Dmat));
            temp = cur_Y(:,II(uset),:).*conj(cur_Y(:,JJ(uset),:));
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
    all_cwt_choice_norm(cc,:,:) = all_cwt_choice./(all_Acwt_diff);
    all_A2cwt_choice(cc,:,:) = (all_A2cwt_same - all_A2cwt_diff)./(all_Acwt_diff);
end

%%
fig_dir = '/home/james/Desktop/CPfigs/';
close all
fig_width = 5;
rel_height = 0.8;
to_print = false;

ucells = find(nansum(trial_Robs_test) > 1e3);
% ucells = find(totSpkCnts > 1e4);
uSUs = ucells(Ucellnumbers(ucells) > 0);
uMUs = setdiff(ucells,uSUs);

weval = logspace(log10(2),log10(80.001),500);
dinterp = 1:0.25:24;
[Xo,Yo] = meshgrid(wfreqs,uprobes);
[Xq,Yq] = meshgrid(weval,dinterp);
interp_CPmap = interp2(Xo,Yo,1-all_LFPpow_CP',Xq,Yq);

freq_markers = [5 10 20 40 80];
freq_inds = interp1(weval,1:length(weval),freq_markers);

f1 = figure();
imagesc(1:length(weval),(dinterp-1)*0.05,interp_CPmap);
% xlim([2 80]);
colorbar
set(gca,'xtick',freq_inds,'xticklabel',freq_markers);
caxis([0.35 0.65]);
xlabel('Frequency (Hz)');
ylabel('Depth');
if to_print
figufy(f1);
fname = [fig_dir sprintf('%s_LFPpow_CP.pdf',Expt_name)];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);
end

avg_controlled = squeeze(nanmean(controlled_Aavgs(ucells,:,:)));
avg_uncontrolled = squeeze(nanmean(uncontrolled_Aavgs(ucells,:,:)));

interp_map = interp2(Xo,Yo,avg_controlled',Xq,Yq);

f2 = figure();
imagescnan(1:length(weval),(dinterp-1)*0.05,interp_map);
% xlim([1.5 80]);
caxis([-0.005 0.005])
colorbar
set(gca,'xtick',freq_inds,'xticklabel',freq_markers);
xlabel('Frequency (Hz)');
ylabel('Depth');
if to_print
figufy(f2);
fname = [fig_dir sprintf('%s_rate_LFP_dep.pdf',Expt_name)];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);
end

f3 = figure();
shadedErrorBar(wfreqs,1-nanmean(all_LFPpow_pred_CP(ucells,:)),nanstd(all_LFPpow_pred_CP(ucells,:))/sqrt(length(ucells)),{'color','k'})
hold on
set(gca,'xscale','log');
xlim([2 80])
ylim([0.35 0.65]);
line([1.5 100],[0.5 0.5],'color','k');
xl = xlim();
xlabel('Frequency (Hz)');
ylabel('Predicted CP');
if to_print
figufy(f3);
fname = [fig_dir sprintf('%s_Freq_pred_CP.pdf',Expt_name)];
exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f3);
end

msize = 4;
f4 = figure();
plot(1-new_choice_prob(uMUs),1-all_LFPpow_allw_pred_CP(uMUs),'o','markersize',msize)
hold on
plot(1-new_choice_prob(uSUs),1-all_LFPpow_allw_pred_CP(uSUs),'ro','markersize',msize)
line([0 1],[0 1],'color','k')
xlim([0.25 0.75]); ylim([0.25 0.75])
line([0 1],[0.5 0.5],'color','k','linestyle','--')
line([0.5 0.5],[0 1],'color','k','linestyle','--')
xlabel('Measured CP');
ylabel('Predicted CP');
if to_print
figufy(f4);
fname = [fig_dir sprintf('%s_Pred_meas_CP_scatter.pdf',Expt_name)];
exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f4);
end

tavg_choice_norm = squeeze(nanmean(all_cwt_choice_norm,2));
f5 = figure();
shadedErrorBar(wfreqs,nanmean(tavg_choice_norm),nanstd(tavg_choice_norm));
set(gca,'xscale','log');
xlim([2 80]);
if to_print
figufy(f5);
fname = [fig_dir sprintf('%s_Phasedep_choiceLFPmod.pdf',Expt_name)];
exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f5);o
end
%%
% data.Expt = Expt_name;
% data.ucells = ucells;
% data.uSUs = uSUs;
% data.uMUs = uMUs;
% data.LFP_pred_CP = 1-all_LFPpow_allw_pred_CP;
% data.actual_CP = 1-new_choice_prob;
% sname = sprintf('%s_tempdata.mat',Expt_name);
% cd ~/Desktop
% save(sname,'data');

%%
fig_dir = '/home/james/Desktop/CPfigs/';
close all
E1_name = 'M230_tempdata.mat';
E2_name = 'M239_tempdata.mat';

cd ~/Desktop
data1 = load(E1_name);
data2 = load(E2_name);
all_pred_CP = [data1.data.LFP_pred_CP data2.data.LFP_pred_CP];
all_act_CP = [data1.data.actual_CP data2.data.actual_CP];
type1 = zeros(size(data1.data.actual_CP));
type1(data1.data.uSUs) = 2; type1(data1.data.uMUs) = 1;
type2 = zeros(size(data2.data.actual_CP));
type2(data2.data.uSUs) = 2; type2(data2.data.uMUs) = 1;

all_types = [type1 type2];

[a,b] = corr(all_pred_CP(all_types > 0)',all_act_CP(all_types > 0)','type','spearman')

msize = 6;
f1 = figure(); hold on
plot(all_act_CP(all_types == 1),all_pred_CP(all_types == 1),'o','markersize',msize);
plot(all_act_CP(all_types == 2),all_pred_CP(all_types == 2),'ro','markersize',msize);
line([0 1],[0 1],'color','k')
xlim([0.25 0.75]); ylim([0.25 0.75])
line([0 1],[0.5 0.5],'color','k','linestyle','--')
line([0.5 0.5],[0 1],'color','k','linestyle','--')
xlabel('Measured CP');
ylabel('Predicted CP');

% fig_width = 4; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'Combined_predact_CP_scatter.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);


%%
close all

fig_dir = '/home/james/Desktop/CPfigs/';
fig_width = 5;
rel_height = 0.8;

bt = 0.2;
et = 1.8;
uu = find(R_trial_taxis(u_lfp_times) > bt & R_trial_taxis(u_lfp_times) < et);
ch_avg_coherence = squeeze(nanmean(all_cwt_choice_norm(:,uu,:),2));

weval = logspace(log10(2),log10(80.001),500);
dinterp = 1:0.25:24;
[Xo,Yo] = meshgrid(wfreqs,uprobes);
[Xq,Yq] = meshgrid(weval,dinterp);
interp_map = interp2(Xo,Yo,ch_avg_coherence,Xq,Yq);

freq_markers = [5 10 20 40 80];
freq_inds = interp1(weval,1:length(weval),freq_markers);

f1 = figure();
imagesc(1:length(weval),(dinterp-1)*0.05,interp_map);
colorbar
set(gca,'xtick',freq_inds,'xticklabel',freq_markers);
cam = max(abs(caxis())); caxis([0 cam]);
xlabel('Frequency (Hz)');
ylabel('Depth (mm)');

fig_width = 5; rel_height = 0.8;
figufy(f1);
fname = [fig_dir sprintf('%s_LFPchoice_coherence2.pdf',Expt_name)];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

