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
           tot_spks_per_trial(tr,nn) = length(AllExpt.Spikes{nn}.Spikes{trindx}); 
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

trialSpkCnts = squeeze(nansum(fullRobs));
totSpkCnts = squeeze(nansum(trialSpkCnts));
avg_spk_rates = nanmean(reshape(fullRobs,[],Nunits));

fullRobs_ms = bsxfun(@minus,fullRobs,reshape(avg_spk_rates,[1 1 Nunits]));

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

%%
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
uset = find(II > JJ);
avg_down_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(resp_down(II(uset)),:),tot_spks2(resp_down(JJ(uset)),:,:)),1));

%avg chioce-conditional covariances
avg_choice_covmat = 0.5*avg_down_covmat + 0.5*avg_up_covmat;

%normalization matrix of variance products
normfac = sqrt(nanvar(tot_spks_per_trial_ms(test_trials,:))'*nanvar(tot_spks_per_trial_ms(test_trials,:)));

avg_choice_corrmat = avg_choice_covmat./normfac;

noise_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(test_trials,:),tot_spks2(test_trials,:,:)),1));
noise_corrmat = noise_covmat./normfac;

%% USE WAVELET ANALYSIS TO COMPUTE PHASE-LOCKING SPECTRA FOR EACH UNIT
LFP_trial_taxis = AllExpt.Expt.Header.LFPtimes*1e-4;

LFP_dsf = 2;
LFP_Fsd = LFP_Fs/LFP_dsf;
%anti-aliasing filter and high-pass filter
aa_hcf = LFP_Fsd/2*0.8;
% [b_aa,a_aa] = butter(4,aa_hcf/(LFP_Fs/2),'low');
aa_lcf = 0.5;
[b_aa,a_aa] = butter(4,[aa_lcf aa_hcf]/(LFP_Fs/2));

LFP_trial_taxis_ds = downsample(LFP_trial_taxis,LFP_dsf);

%wavelet parameters
nwfreqs = 8;
min_freq = 1; max_freq = 20;
% nwfreqs = 5;
% min_freq = 1; max_freq = 6;
min_scale = 1/max_freq*LFP_Fsd;
max_scale = 1/min_freq*LFP_Fsd;
wavetype = 'cmor1-1';
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,wavetype,1/LFP_Fsd);

beg_buffer = round(0.15/dt);
R_trial_taxis = (beg_buffer:(NT-1))*dt + dt/2;
TLEN = length(R_trial_taxis);

nprobes = 24;
uprobes = 1:3:nprobes;
all_spk_id = [];
all_spk_phases = nan(sum(tot_spks_per_trial(:)),length(wfreqs),length(uprobes));
trial_LFP_phases = nan(Ntrials,TLEN,length(wfreqs),length(uprobes));
trial_LFP_amps = nan(Ntrials,TLEN,length(wfreqs),length(uprobes));
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
    
    trial_LFP_phases(tr,:,:,:) = mod(interp1(LFP_trial_taxis_ds,cur_phasegram,R_trial_taxis),2*pi)-pi;
    trial_LFP_amps(tr,:,:,:) = interp1(LFP_trial_taxis_ds,cur_ampgram,R_trial_taxis);
    trial_LFPs(tr,:,:) = interp1(LFP_trial_taxis_ds,cur_LFPs,R_trial_taxis);
end
trial_LFP_phases = permute(trial_LFP_phases,[2 1 3 4]);
trial_LFP_amps = permute(trial_LFP_amps,[2 1 3 4]);
trial_LFPs = permute(trial_LFPs,[2 1 3]);

trial_LFP_phases = reshape(trial_LFP_phases,[],length(wfreqs)*length(uprobes));
trial_LFP_amps = reshape(trial_LFP_amps,[],length(wfreqs)*length(uprobes));
trial_LFPs = reshape(trial_LFPs,[],length(uprobes));

trial_LFP_amps = nanzscore(trial_LFP_amps);
trial_LFPs = nanzscore(trial_LFPs);

LFP_tbt = reshape(trial_LFP_amps(:,end),[TLEN Ntrials]);


%%
good_LFP_trials = find(all(~isnan(LFP_tbt)));
utrials = good_LFP_trials(trialOB(good_LFP_trials) == 130);
% utrials = 1:length(trialOB);
[RR,TT] = meshgrid(1:Ntrials,1:TLEN);
RR = RR(:);

tdata = find(~ismember(RR,utrials) & ismember(RR,good_LFP_trials));
% tdata = find(ismember(RR,utrials));
udata = find(ismember(RR,utrials));
% xv_frac = 0.;
% xv_trials = randperm(length(utrials));
% xv_trials = xv_trials(1:round(xv_frac*length(utrials)));
% tr_trials = setdiff(1:length(utrials),xv_trials);

% trdata = find(ismember(RR(udata),utrials(tr_trials)));
% xvdata = find(ismember(RR(udata),utrials(xv_trials)));

lambda_L2 = 5e3;
mod_reg_params = NMMcreate_reg_params('lambda_L2',lambda_L2);
silent = 1;

uchs = 1:length(uprobes);
[CC,WW] = meshgrid(1:length(uprobes),1:length(wfreqs));
ch_use = find(ismember(CC(:),uchs));
mod_stim_params(1:2) = NMMcreate_stim_params(length(wfreqs)*length(uchs));

X{1} = trial_LFP_amps(:,ch_use).*cos(trial_LFP_phases(:,ch_use));
X{2} = trial_LFP_amps(:,ch_use).*sin(trial_LFP_phases(:,ch_use));

all_prate_outs = nan(length(udata),Nunits);
for cc = 1:Nunits
    cc
    cur_Robs = reshape(fullRobs((beg_buffer+1):end,:,cc),[],1);
%     uinds = find(~isnan(cur_Robs));
    truinds = tdata(~isnan(cur_Robs(tdata)));
    xvuinds = udata(~isnan(cur_Robs(udata)));
    if ~isempty(truinds)
        lfp_mod(cc) = NMMinitialize_model(mod_stim_params,[1 1],{'lin','lin'},mod_reg_params,[1 2]);
        lfp_mod(cc) = NMMfit_filters(lfp_mod(cc),cur_Robs,X,[],truinds,silent);
        [LL, nullLL, pred_rate] = NMMeval_model( lfp_mod(cc), cur_Robs, X, [], truinds);
        [xvLL, nullxvLL, pred_rate] = NMMeval_model( lfp_mod(cc), cur_Robs, X, [], xvuinds);
        lfp_mod_LLimp(cc) = (LL - nullLL)/log(2);
        lfp_mod_xvLLimp(cc) = (xvLL - nullxvLL)/log(2);
        [~, ~, all_prate_outs(:,cc)] = NMMeval_model( lfp_mod(cc), cur_Robs, X,[],udata);
    else
        lfp_mod(cc) = [];
    end
end

%%
all_prate_outs_ms = bsxfun(@minus,all_prate_outs,nanmean(all_prate_outs));
all_prate_outs_ms = reshape(all_prate_outs_ms,[TLEN length(utrials) Nunits]);

all_sum_outs = squeeze(nansum(all_prate_outs_ms));
zscore_sumouts = nanmean(nanzscore(all_sum_outs),2);
bad_trials = find(zscore_sumouts >4);
all_prate_outs_ms(:,bad_trials,:) = nan;
all_sum_outs(bad_trials,:) = nan;

all_sum_outs_ms = bsxfun(@minus,all_sum_outs,nanmean(all_sum_outs));

unique_ses = unique(trialSe(utrials));
samestim_covmat = zeros(Nunits,Nunits);
samestim_sumcovmat = zeros(Nunits,Nunits);
cnt = 0;
for ss = 1:length(unique_ses)
   curset = find(trialSe(utrials) == unique_ses(ss)); 
   if length(curset) >= 2
   temp_prates = squeeze(all_prate_outs_ms(:,curset,:));
   temp_prates2 = reshape(temp_prates,TLEN,[],1,Nunits);
   
   [II,JJ] = meshgrid(1:length(curset));
    uset = find(II > JJ);
    temp = squeeze(nansum(nanmean(bsxfun(@times,temp_prates(:,II(uset),:),temp_prates2(:,JJ(uset),:,:)),1),2));
    if ~any(isnan(temp(:)))
    samestim_covmat = samestim_covmat + temp;
    end
    
   temp_prates = squeeze(all_sum_outs_ms(curset,:));
   temp_prates2 = reshape(temp_prates,[],1,Nunits);
    temp = squeeze(nansum(bsxfun(@times,temp_prates(II(uset),:),temp_prates2(JJ(uset),:,:)),1));
    if ~any(isnan(temp(:)))
    samestim_sumcovmat = samestim_sumcovmat + temp;
    cnt = cnt + length(uset);
    end
    
   end
end
samestim_covmat = samestim_covmat/cnt;
samestim_sumcovmat = samestim_sumcovmat/cnt;

cur_Robs = reshape(fullRobs((beg_buffer+1:end),utrials,:),[],Nunits);
ov_LFP_covmat = nancov(all_prate_outs);
LFP_noise_covmat = ov_LFP_covmat - samestim_covmat;

ov_LFP_sumcovmat = nancov(all_sum_outs);
LFP_noise_sumcovmat = ov_LFP_sumcovmat - samestim_sumcovmat;

% LFP_normfac = sqrt(var(all_prate_outs)'*var(all_prate_outs));
LFP_normfac = sqrt(nanvar(cur_Robs)'*nanvar(cur_Robs));
LFP_noise_corrmat = LFP_noise_covmat./LFP_normfac;

LFP_sumnormfac = sqrt(nanvar(tot_spks_per_trial_norm)'*nanvar(tot_spks_per_trial_norm));
LFP_noise_sumcorrmat = LFP_noise_sumcovmat./LFP_sumnormfac;

%%
all_prate_outs_resh = reshape(all_prate_outs,[TLEN length(utrials) Nunits]);
LFPmod_pred_spkcnts = squeeze(nansum(all_prate_outs_resh));
LFPmod_pred_spkcnts(bad_trials,:) = nan;

resp_up = find(trialrespDir(utrials) == 1);
resp_down = find(trialrespDir(utrials) == -1);
LFPmod_choice_prob = nan(Nunits,1);
for cc = 1:Nunits
   cur_resp_ax = prctile(LFPmod_pred_spkcnts(:,cc),[0:5:100]);
   up_hist = histc(LFPmod_pred_spkcnts(resp_up,cc),cur_resp_ax);
   down_hist = histc(LFPmod_pred_spkcnts(resp_down,cc),cur_resp_ax);
   true_pos = cumsum(up_hist)/sum(~isnan(LFPmod_pred_spkcnts(resp_up,cc)));
   false_pos = cumsum(down_hist)/sum(~isnan(LFPmod_pred_spkcnts(resp_down,cc)));
   LFPmod_choice_prob(cc) = trapz(false_pos,true_pos);
end

LFPmod_choice_diff = nanmean(LFPmod_pred_spkcnts(resp_up,:)) - nanmean(LFPmod_pred_spkcnts(resp_down,:));
act_choice_diff = nanmean(tot_spks_per_trial_norm(utrials(resp_up),:)) - nanmean(tot_spks_per_trial_norm(utrials(resp_down),:));

%% CALCULATE CHOICE-predictable covariance
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
un_stim_seeds = unique(trialSe(test_trials));

beg_buff = 15; %number of bins from beginning of trial to exclude

%subtract off avg rates for times within this set of trials
[RR,TT] = meshgrid(1:Ntrials,1:NT);
istrain = (ismember(RR(:),test_trials)) & (TT(:) > beg_buff);
fullRobs_resh = reshape(fullRobs,[],Nunits);
fullRobs_ms = bsxfun(@minus,fullRobs,reshape(nanmean(fullRobs_resh(istrain,:)),[1 1 Nunits]));

cur_XC = zeros(Nunits,Nunits);
cur_cnt = zeros(Nunits,1);
cur_XC2 = zeros(Nunits,Nunits);
cur_cnt2 = zeros(Nunits,1);
for tr = 1:length(un_stim_seeds)
    cur_tr_set = test_trials(trialSe(test_trials) == un_stim_seeds(tr));
    fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
    if length(cur_tr_set) >= 2
        cur_resp = trialrespDir(cur_tr_set);
        [II,JJ] = meshgrid(1:length(cur_tr_set));
        cur_Robs = squeeze(fullRobs((beg_buff+1):end,cur_tr_set,:));
        cur_Robs2 = reshape(cur_Robs,[],length(cur_tr_set),1,Nunits);
        
        cur_Dmat = abs(squareform(pdist(cur_resp')));
        cur_Dmat(logical(eye(length(cur_tr_set)))) = nan;
        cur_Dmat(JJ > II) = nan;
        
        curset = find(cur_Dmat == 0);
        cur_XC = cur_XC + squeeze(nansum(nanmean(bsxfun(@times,cur_Robs2(:,II(curset),:,:),cur_Robs(:,JJ(curset),:))),2));
        cur_cnt = cur_cnt + squeeze(sum(~(isnan(cur_Robs(1,II(curset),:))) & ~(isnan(cur_Robs(1,JJ(curset),:))),2));
        
        curset = find(~isnan(cur_Dmat));
        cur_XC2 = cur_XC2 + squeeze(nansum(nanmean(bsxfun(@times,cur_Robs2(:,II(curset),:,:),cur_Robs(:,JJ(curset),:))),2));
        cur_cnt2 = cur_cnt2 + squeeze(sum(~(isnan(cur_Robs(1,II(curset),:))) & ~(isnan(cur_Robs(1,JJ(curset),:))),2));
    end
end

sameVar = bsxfun(@rdivide,cur_XC,cur_cnt);
randVar = bsxfun(@rdivide,cur_XC2,cur_cnt2);

choice_covmat = sameVar - randVar;
SU_choice_varfrac = diag(choice_covmat)./diag(randVar);
normfac = sqrt(diag(randVar)*diag(randVar)');
choice_corrmat = choice_covmat./normfac;

%% LFP state-explained covariance
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
un_stim_seeds = unique(trialSe(test_trials));
targChan = 10;

R_trial_taxis = (beg_buffer:(NT-1))*dt + dt/2;
TLEN = length(R_trial_taxis);
trial_LFP_state = nan(TLEN,Ntrials,nprobes);

%5-15 looks good, see similar structure in choice- and LFP-explained cov
%mats
lcf = 1.5; hcf = 4; %cut-off freqs for LFP filter
[filt_bb,filt_aa] = butter(2,[lcf hcf]/(LFP_Fs/2));
for tr = 1:Ntrials
    cur_LFPs = double(AllExpt.Expt.Trials(tr).LFP);
    cur_LFPs(isnan(cur_LFPs)) = 0;
    cur_LFPs = filtfilt(filt_bb,filt_aa,cur_LFPs);
%     cur_phase = unwrap(angle(hilbert(cur_LFPs)));
%     cur_phase_interp = mod(interp1(LFP_trial_taxis,cur_phase,R_trial_taxis),2*pi);
    cur_lfp_interp = interp1(LFP_trial_taxis,cur_LFPs,R_trial_taxis);
    trial_LFP_state(:,tr,:) = cur_lfp_interp;
end

beg_buff = 15; %number of bins from beginning of trial to exclude

maxlag_ED = 2;
ED_space = 0.25;
ED_bin_edges = 0:ED_space:maxlag_ED;
ED_bin_centers = (ED_bin_edges(1:end-1)+ED_bin_edges(2:end))/2;

%subtract off avg rates for times within this set of trials
[RR,TT] = meshgrid(1:Ntrials,1:NT);
istrain = (ismember(RR(:),test_trials)) & (TT(:) > beg_buff);
fullRobs_resh = reshape(fullRobs,[],Nunits);
fullRobs_ms = bsxfun(@minus,fullRobs,reshape(nanmean(fullRobs_resh(istrain,:)),[1 1 Nunits]));
tot_var = nanvar(fullRobs_resh(istrain,:));

cur_tr_set = test_trials;
fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
[II,JJ] = meshgrid(1:length(cur_tr_set));

cur_LFP_data = squeeze(trial_LFP_state(:,cur_tr_set,targChan));
cur_XC = nan(NT,length(ED_bin_centers),Nunits,Nunits);
rand_XC = nan(NT,Nunits,Nunits);
for tt = (beg_buff+1):NT
    tt
    cur_Robs = squeeze(fullRobs_ms(tt,cur_tr_set,:));
    cur_Robs2 = reshape(cur_Robs,[],1,Nunits);
    cur_Dmat = abs(circ_dist2(cur_LFP_data(tt,:),cur_LFP_data(tt,:)));
    cur_Dmat(JJ >= II) = nan;
    for jj = 1:length(ED_bin_centers)
        curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
        %                 cur_XC(jj,:) = cur_XC(jj,:) + squeeze(nansum(bsxfun(@times,cur_Robs(II(curset),:),cur_Robs(JJ(curset),:)),1));
        cur_XC(tt,jj,:,:) = (nanmean(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
    end
    curset = ~isnan(cur_Dmat);
    %             rand_XC = rand_XC + squeeze(nansum(bsxfun(@times,cur_Robs(II(curset),:),cur_Robs(JJ(curset),:)),1));
    rand_XC(tt,:,:) = squeeze(nanmean(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
end
cur_XC = squeeze(nanmean(cur_XC));
rand_XC = squeeze(nanmean(rand_XC));

LFP_explained(ppp,bbb,:,:) = squeeze(stateDepVar(1,:,:)) - randVar; %dont worry about spline interpolation for now, just take the bin with most similar LFP states
