%% LOAD PROCESSED DATA
clear all
% cd('/home/james/Data/bruce/ChoiceProb/')
cd('~/Data/bruce/ChoiceProb/')

Expt_name = 'M239';

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
min_freq = 1.; max_freq = 60;
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
trial_LFP_phases = nan(Ntrials,TLEN,length(wfreqs),length(uprobes));
trial_LFP_amps = nan(Ntrials,TLEN,length(wfreqs),length(uprobes));
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
    cur_phasegram = unwrap(angle(cur_cwt));
    cur_ampgram = abs(cur_cwt);
    
    cur_phasegram(bad_LFPs,:,:) = nan;
    cur_ampgram(bad_LFPs,:,:) = nan;
    cur_LFPs(bad_LFPs,:) = nan;
    cur_cwt(bad_LFPs,:,:) = nan;
    
    trial_LFP_phases(tr,:,:,:) = mod(interp1(LFP_trial_taxis_ds,cur_phasegram,R_trial_taxis),2*pi)-pi;
    trial_LFP_amps(tr,:,:,:) = interp1(LFP_trial_taxis_ds,cur_ampgram,R_trial_taxis);
    trial_LFP_cwt(tr,:,:,:) = interp1(LFP_trial_taxis_ds,cur_cwt,R_trial_taxis);
    trial_LFPs(tr,:,:) = interp1(LFP_trial_taxis_ds,cur_LFPs,R_trial_taxis);
end
trial_LFP_phases = permute(trial_LFP_phases,[2 1 3 4]);
trial_LFP_amps = permute(trial_LFP_amps,[2 1 3 4]);
trial_LFP_cwt = permute(trial_LFP_cwt,[2 1 3 4]);
trial_LFPs = permute(trial_LFPs,[2 1 3]);

trial_LFP_phases = reshape(trial_LFP_phases,[],length(wfreqs)*length(uprobes));
trial_LFP_amps = reshape(trial_LFP_amps,[],length(wfreqs)*length(uprobes));
trial_LFP_cwt = reshape(trial_LFP_cwt,[],length(wfreqs)*length(uprobes));
trial_LFPs = reshape(trial_LFPs,[],length(uprobes));

trial_LFP_amps = nanzscore(trial_LFP_amps);
trial_LFPs = nanzscore(trial_LFPs);

LFP_tbt = reshape(trial_LFP_amps(:,end),[TLEN Ntrials]);

%%
tbt_LFP_amp = reshape(trial_LFP_amps,[TLEN Ntrials length(wfreqs) length(uprobes)]);
tbt_LFP_phases = reshape(trial_LFP_phases,[TLEN Ntrials length(wfreqs) length(uprobes)]);
tbt_LFPs = reshape(trial_LFPs,[TLEN Ntrials length(uprobes)]);
tbt_LFP_cwt = reshape(trial_LFP_cwt,[TLEN Ntrials length(wfreqs) length(uprobes)]);

%%
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
train_trials = find(trialOB ~= 130);

resp_up = trialrespDir(test_trials) == 1;
resp_down = trialrespDir(test_trials) == -1;
up_LFP_amp = squeeze(nanmean(tbt_LFP_amp(:,test_trials(resp_up),:,:),2));
down_LFP_amp = squeeze(nanmean(tbt_LFP_amp(:,test_trials(resp_down),:,:),2));
choice_LFP_ampdiff = up_LFP_amp - down_LFP_amp;

cur_phase = permute(tbt_LFP_phases(:,test_trials(resp_up),:,:),[2 1 3 4]);
bad_trials = find(any(isnan(cur_phase(:,:,1,1)),2));
cur_phase(bad_trials,:,:,:) = [];
up_LFP_r = squeeze(circ_r(cur_phase));
% up_LFP_r = squeeze(circ_mean(cur_phase));

cur_phase = permute(tbt_LFP_phases(:,test_trials(resp_down),:,:),[2 1 3 4]);
bad_trials = find(any(isnan(cur_phase(:,:,1,1)),2));
cur_phase(bad_trials,:,:,:) = [];
% down_LFP_r(162:end,:,:,:) = [];
down_LFP_r = squeeze(circ_r(cur_phase));
% down_LFP_r = squeeze(circ_mean(cur_phase));
choice_LFP_rdiff = up_LFP_r - down_LFP_r;

sresp_up = trialrespDir(train_trials) == 1;
sresp_down = trialrespDir(train_trials) == -1;
sup_LFP_amp = squeeze(nanmean(tbt_LFP_amp(:,train_trials(resp_up),:,:),2));
sdown_LFP_amp = squeeze(nanmean(tbt_LFP_amp(:,train_trials(resp_down),:,:),2));
schoice_LFP_ampdiff = sup_LFP_amp - sdown_LFP_amp;

up_LFPs = squeeze(nanmean(tbt_LFPs(:,test_trials(resp_up),:,:),2));
down_LFPs = squeeze(nanmean(tbt_LFPs(:,test_trials(resp_down),:,:),2));
choice_LFP_diff = up_LFPs - down_LFPs;

norm = squeeze(nanstd(tbt_LFP_amp(:,test_trials,:,:),[],2));
choice_LFP_ampdiff_z = choice_LFP_ampdiff./norm;

norm = squeeze(nanstd(tbt_LFPs(:,test_trials,:),[],2));
choice_LFP_diff_z = choice_LFP_diff./norm;

beg_buff = 0.15; end_buff = 0.05;
utimes = find(R_trial_taxis > beg_buff & 2 - R_trial_taxis > end_buff);
tot_trial_pow = squeeze(nanmean(tbt_LFP_amp(utimes,:,:,:)));
tot_trial_powX = reshape(tot_trial_pow(test_trials,:,:),length(test_trials),[]);
tot_trial_powX = bsxfun(@minus,tot_trial_powX,nanmean(tot_trial_powX));
tot_trial_powX = bsxfun(@rdivide,tot_trial_powX,nanstd(tot_trial_powX));

tot_spkY = squeeze(nansum(fullRobs(utimes,test_trials,:)));
% tot_spkY = tot_spks_per_trial_norm(test_trials,:);
tot_spkY = bsxfun(@minus,tot_spkY,nanmean(tot_spkY));
pow_spk_cov = squeeze(nanmean(bsxfun(@times,tot_spkY,reshape(tot_trial_powX,length(test_trials),1,[]))));
normfac = sqrt(squeeze(bsxfun(@times,nanvar(tot_spkY),nanvar(tot_trial_powX)')))';
pow_spk_corr = pow_spk_cov./normfac;
pow_spk_corr = reshape(pow_spk_corr,[Nunits length(wfreqs) length(uprobes)]);

trial_pow_choice_diff = nanmean(tot_trial_powX(resp_up,:)) - nanmean(tot_trial_powX(resp_down,:));
% trial_pow_sepvec = inv(nancov(tot_trial_powX))*trial_pow_choice_diff';
%%


% temp_Robs = fullRobs(utimes,test_trials,:);
% temp_Robs = bsxfun(@minus,temp_Robs,nanmean(temp_Robs));
% 
% cc = 4;
% LFP_ampX = reshape(tbt_LFP_amp(utimes,test_trials,:,cc),[],length(wfreqs));
% spk_Y = reshape(temp_Robs,[],1,Nunits);
% LFP_ampX = bsxfun(@minus,LFP_ampX,nanmean(LFP_ampX));
% spk_Y = bsxfun(@minus,spk_Y,nanmean(spk_Y));
% 
% LFP_spk_cov = squeeze(nanmean(bsxfun(@times,LFP_ampX,spk_Y)));
% normfac = sqrt(bsxfun(@times,nanvar(LFP_ampX),squeeze(nanvar(spk_Y)))');
% LFP_spk_corr = LFP_spk_cov./normfac;
% 
% temp_Robs = fullRobs(utimes,test_trials,:);
% spk_Y = reshape(temp_Robs,[],1,Nunits);

%%
% ch = 6;
% LFPmod_choice_prob = nan(Nunits,length(wfreqs));
% for ww = 1:length(wfreqs)
% ww
% LFP_ampXtr = reshape(tbt_LFP_amp(utimes,train_trials,ww,:),[],length(uprobes));
% LFP_ampX = reshape(tbt_LFP_amp(utimes,test_trials,ww,:),[],length(uprobes));
% % LFP_ampX = reshape(tbt_LFP_amp(utimes,test_trials,ww,ch),[],1);
% clear B pred_rate
% for cc = 1:Nunits
%     cur_Y = squeeze(spk_Y(:,:,cc));
%     cur_uset = find(~isnan(cur_Y));
%     B(cc,:) = glmfit(LFP_ampX(cur_uset,:),cur_Y(cur_uset),'poisson');
% %     B(cc,:) = glmfit(LFP_ampX(cur_uset,:),cur_Y(cur_uset),'poisson');
%     pred_rate(:,cc) = glmval(B(cc,:)',LFP_ampX,'log');
% end
% tbt_LFP_pred_rate = reshape(pred_rate,[length(utimes) length(test_trials) Nunits]);
% tbt_LFP_sumrate = squeeze(nansum(tbt_LFP_pred_rate));
% 
% 
% tbt_LFP_sumrate_ms = bsxfun(@minus,tbt_LFP_sumrate,nanmean(tbt_LFP_sumrate));
% tot_spks2 = reshape(tbt_LFP_sumrate_ms,[],1,Nunits);
% 
% [II,JJ] = meshgrid(1:length(test_trials));
% cur_Dmat = abs(squareform(pdist(trialrespDir(test_trials)')));
% cur_Dmat(JJ == II) = nan;
% uset = find(cur_Dmat == 0);
% avg_LFPpred_choice_covmat = squeeze(nanmean(bsxfun(@times,tbt_LFP_sumrate_ms(II(uset),:),tot_spks2(JJ(uset),:,:)),1));
% normfac = sqrt(nanvar(tot_spks_per_trial_ms(test_trials,:))'*nanvar(tot_spks_per_trial_ms(test_trials,:)));
% avgLFPpred_choice_corrmat = avg_LFPpred_choice_covmat./normfac;
% 
% %%
% resp_up = find(trialrespDir(test_trials) == 1);
% resp_down = find(trialrespDir(test_trials) == -1);
% for cc = 1:Nunits
%     cur_resp_ax = prctile(tbt_LFP_sumrate(:,cc),[0:5:100]);
%     up_hist = histc(tbt_LFP_sumrate(resp_up,cc),cur_resp_ax);
%     down_hist = histc(tbt_LFP_sumrate(resp_down,cc),cur_resp_ax);
%     true_pos = cumsum(up_hist)/sum(~isnan(tbt_LFP_sumrate(resp_up,cc)));
%     false_pos = cumsum(down_hist)/sum(~isnan(tbt_LFP_sumrate(resp_down,cc)));
%     LFPmod_choice_prob(cc,ww) = trapz(false_pos,true_pos);
% end
% end

%%
beg_buff = 0.5; end_buff = 0.5;
utimes = find(R_trial_taxis > beg_buff & 2 - R_trial_taxis > end_buff);
tot_trial_pow = squeeze(nanmean(tbt_LFP_amp(utimes,:,:,:)));

tot_trial_powX = reshape(tot_trial_pow(train_trials,:,:),length(train_trials),[]);
tot_trial_powX = bsxfun(@minus,tot_trial_powX,nanmean(tot_trial_powX));
tot_trial_powX = bsxfun(@rdivide,tot_trial_powX,nanstd(tot_trial_powX));
all_powX_tr = reshape(tot_trial_powX,length(train_trials),length(wfreqs),[]);

tot_trial_powX = reshape(tot_trial_pow(test_trials,:,:),length(test_trials),[]);
tot_trial_powX = bsxfun(@minus,tot_trial_powX,nanmean(tot_trial_powX));
tot_trial_powX = bsxfun(@rdivide,tot_trial_powX,nanstd(tot_trial_powX));
all_powX = reshape(tot_trial_powX,length(test_trials),length(wfreqs),[]);

tot_spkY = squeeze(nansum(fullRobs(utimes,test_trials,:)));
tot_spkY_tr = squeeze(nansum(fullRobs(utimes,train_trials,:)));

% ch = 4;

LFPmod_choice_prob = nan(Nunits,length(wfreqs));
% tot_spkY = tot_spks_per_trial_norm(test_trials,:);
for ww = 1:length(wfreqs)
ww
% LFP_ampX = reshape(all_powX(:,ww,ch),[],1);
LFP_ampX = reshape(all_powX(:,ww,:),[],length(uprobes));
LFP_ampX_tr = reshape(all_powX_tr(:,ww,:),[],length(uprobes));
% LFP_ampX = reshape(tbt_LFP_amp(utimes,test_trials,ww,ch),[],1);
clear B 
pred_rate = nan(length(test_trials),Nunits);
for cc = 1:Nunits
    cur_Y = squeeze(tot_spkY_tr(:,cc));
    cur_uset = find(~isnan(cur_Y));
    B(cc,:) = glmfit(LFP_ampX_tr(cur_uset,:),cur_Y(cur_uset),'poisson');
%     cur_Y = squeeze(tot_spkY(:,cc));
%     cur_uset = find(~isnan(cur_Y));
%     B(cc,:) = glmfit(LFP_ampX(cur_uset,:),cur_Y(cur_uset),'poisson');
    cur_Y = squeeze(tot_spkY(:,cc));
    cur_uset = find(~isnan(cur_Y));
    pred_rate(cur_uset,cc) = glmval(B(cc,:)',LFP_ampX(cur_uset,:),'log');
end
tbt_LFP_sumrate = reshape(pred_rate,[length(test_trials) Nunits]);


tbt_LFP_sumrate_ms = bsxfun(@minus,tbt_LFP_sumrate,nanmean(tbt_LFP_sumrate));
tot_spks2 = reshape(tbt_LFP_sumrate_ms,[],1,Nunits);

[II,JJ] = meshgrid(1:length(test_trials));
cur_Dmat = abs(squareform(pdist(trialrespDir(test_trials)')));
cur_Dmat(JJ == II) = nan;
uset = find(cur_Dmat == 0);
avg_LFPpred_choice_covmat = squeeze(nanmean(bsxfun(@times,tbt_LFP_sumrate_ms(II(uset),:),tot_spks2(JJ(uset),:,:)),1));
normfac = sqrt(nanvar(tot_spks_per_trial_ms(test_trials,:))'*nanvar(tot_spks_per_trial_ms(test_trials,:)));
avgLFPpred_choice_corrmat = avg_LFPpred_choice_covmat./normfac;

%%
resp_up = find(trialrespDir(test_trials) == 1);
resp_down = find(trialrespDir(test_trials) == -1);
for cc = 1:Nunits
    cur_resp_ax = prctile(tbt_LFP_sumrate(:,cc),[0:5:100]);
    up_hist = histc(tbt_LFP_sumrate(resp_up,cc),cur_resp_ax);
    down_hist = histc(tbt_LFP_sumrate(resp_down,cc),cur_resp_ax);
    true_pos = cumsum(up_hist)/sum(~isnan(tbt_LFP_sumrate(resp_up,cc)));
    false_pos = cumsum(down_hist)/sum(~isnan(tbt_LFP_sumrate(resp_down,cc)));
    LFPmod_choice_prob(cc,ww) = trapz(false_pos,true_pos);
end
end


%%
poss_trials = find(trialOB == 130 & trialrespDir ~= 0);
un_stim_seeds = unique(trialSe(poss_trials));
for cc = 1:8;
cc
all_cwt_same = zeros(length(R_trial_taxis),length(wfreqs));
all_cwt_diff = all_cwt_same;
all_same_cnt = all_cwt_same;
all_diff_cnt = all_cwt_same;
all_Acwt_same = zeros(length(R_trial_taxis),length(wfreqs));
all_Acwt_diff = all_Acwt_same;
for ss = 1:length(un_stim_seeds)
    cur_trials = poss_trials(trialSe(poss_trials) == un_stim_seeds(ss));
    if length(cur_trials) >= 2
        cur_resp = trialrespDir(cur_trials);
        [II,JJ] = meshgrid(1:length(cur_trials));
        cur_Y = squeeze(tbt_LFP_cwt(:,cur_trials,:,cc));
        
        cur_Dmat = abs(squareform(pdist(cur_resp')));
        cur_Dmat(JJ == II) = nan;
        
        uset = find(cur_Dmat == 0);
        temp = cur_Y(:,II(uset),:).*cur_Y(:,JJ(uset),:);
        tempm = squeeze(nanmean(temp,2));
        all_cwt_same(~isnan(tempm)) = all_cwt_same(~isnan(tempm))  + tempm(~isnan(tempm)) ;
        temp = abs(cur_Y(:,II(uset),:)).*abs(cur_Y(:,JJ(uset),:));
         tempm = squeeze(nanmean(temp,2));
        all_Acwt_same(~isnan(tempm)) = all_Acwt_same(~isnan(tempm))  + tempm(~isnan(tempm)) ;
        all_same_cnt = all_same_cnt + squeeze(sum(~isnan(temp),2));
        
%         uset = find(cur_Dmat ~= 0);
        uset = find(~isnan(cur_Dmat));
        temp = cur_Y(:,II(uset),:).*cur_Y(:,JJ(uset),:);
        tempm = squeeze(nanmean(temp,2));
        all_cwt_diff(~isnan(tempm)) = all_cwt_diff(~isnan(tempm)) + tempm(~isnan(tempm));
        temp = abs(cur_Y(:,II(uset),:)).*abs(cur_Y(:,JJ(uset),:));
        tempm = squeeze(nanmean(temp,2));
        all_Acwt_diff(~isnan(tempm)) = all_Acwt_diff(~isnan(tempm)) + tempm(~isnan(tempm));
        all_diff_cnt = all_diff_cnt + squeeze(sum(~isnan(temp),2));
    end
    
end

all_cwt_same = abs(all_cwt_same)./all_same_cnt;
all_cwt_diff = abs(all_cwt_diff)./all_diff_cnt;
all_Acwt_same = all_Acwt_same./all_same_cnt;
all_Acwt_diff = all_Acwt_diff./all_diff_cnt;

% all_cwt_choice = all_cwt_same - all_cwt_diff;
% norm_fac = squeeze(nanvar(abs(tbt_LFP_cwt(:,poss_trials,:,cc)),[],2));
% all_cwt_corr = all_cwt_choice./norm_fac;
all_cwt_same_norm = all_cwt_same./all_Acwt_same;
all_cwt_choice = abs(all_cwt_same) - abs(all_cwt_diff);
all_cwt_choice_norm(cc,:,:) = all_cwt_choice./(0.5*all_Acwt_same + 0.5*all_Acwt_diff);
end

%%
poss_trials = find(trialOB == 130 & trialrespDir ~= 0);
un_stim_seeds = unique(trialSe(poss_trials));
for cc = 1:8;
    cc
    all_cwt_same = zeros(length(R_trial_taxis),length(wfreqs));
    all_cwt_diff = all_cwt_same;
    all_same_cnt = all_cwt_same;
    all_diff_cnt = all_cwt_same;
    all_Acwt_same = zeros(length(R_trial_taxis),length(wfreqs));
    all_Acwt_diff = all_Acwt_same;
    cur_trials = poss_trials;
    cur_resp = trialrespDir(cur_trials);
    [II,JJ] = meshgrid(1:length(cur_trials));
    cur_Y = squeeze(tbt_LFP_cwt(:,cur_trials,:,cc));
    
    cur_Dmat = abs(squareform(pdist(cur_resp')));
    cur_Dmat(JJ == II) = nan;
    
    same_uset = find(cur_Dmat == 0);
    
    chunk_size = 1e3;
    n_chunks = floor(length(same_uset)/chunk_size);
    for ii = 1:n_chunks
        fprintf('%d of %d\n',ii,n_chunks);
        cur_set = (ii-1)*chunk_size + (1:chunk_size);
        
        temp = cur_Y(:,II(same_uset(cur_set)),:).*cur_Y(:,JJ(same_uset(cur_set)),:);
        tempm = squeeze(nanmean(temp,2));
        all_cwt_same(~isnan(tempm)) = all_cwt_same(~isnan(tempm))  + tempm(~isnan(tempm)) ;
        temp = abs(cur_Y(:,II(same_uset(cur_set)),:)).*abs(cur_Y(:,JJ(same_uset(cur_set)),:));
        tempm = squeeze(nanmean(temp,2));
        all_Acwt_same(~isnan(tempm)) = all_Acwt_same(~isnan(tempm))  + tempm(~isnan(tempm)) ;
        all_same_cnt = all_same_cnt + squeeze(sum(~isnan(temp),2));
    end
    
    rand_uset = find(~isnan(cur_Dmat));
    n_chunks = floor(length(rand_uset)/chunk_size);
    for ii = 1:n_chunks
        fprintf('%d of %d\n',ii,n_chunks);
        cur_set = (ii-1)*chunk_size + (1:chunk_size);
        
        temp = cur_Y(:,II(rand_uset(cur_set)),:).*cur_Y(:,JJ(rand_uset(cur_set)),:);
        tempm = squeeze(nanmean(temp,2));
        all_cwt_diff(~isnan(tempm)) = all_cwt_diff(~isnan(tempm))  + tempm(~isnan(tempm)) ;
        temp = abs(cur_Y(:,II(rand_uset(cur_set)),:)).*abs(cur_Y(:,JJ(rand_uset(cur_set)),:));
        tempm = squeeze(nanmean(temp,2));
        all_Acwt_diff(~isnan(tempm)) = all_Acwt_diff(~isnan(tempm))  + tempm(~isnan(tempm)) ;
        all_diff_cnt = all_diff_cnt + squeeze(sum(~isnan(temp),2));
    end
    
    all_cwt_same = abs(all_cwt_same)./all_same_cnt;
    all_cwt_diff = abs(all_cwt_diff)./all_diff_cnt;
    all_Acwt_same = all_Acwt_same./all_same_cnt;
    all_Acwt_diff = all_Acwt_diff./all_diff_cnt;
    
    % all_cwt_choice = all_cwt_same - all_cwt_diff;
    % norm_fac = squeeze(nanvar(abs(tbt_LFP_cwt(:,poss_trials,:,cc)),[],2));
    % all_cwt_corr = all_cwt_choice./norm_fac;
    all_cwt_same_norm = all_cwt_same./all_Acwt_same;
    all_cwt_choice = abs(all_cwt_same) - abs(all_cwt_diff);
    all_cwt_choice_norm(cc,:,:) = all_cwt_choice./(0.5*all_Acwt_same + 0.5*all_Acwt_diff);
end

