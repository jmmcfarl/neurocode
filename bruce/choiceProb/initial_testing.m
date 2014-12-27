%% LOAD PROCESSED DATA
clear all
% cd('/home/james/Data/bruce/ChoiceProb/')
cd('~/Data/bruce/ChoiceProb/')

load('M239/lemM239.image.ORBW.LFP.mat')

%% FOR M239
block_trial_boundaries = [1 157;
 158 265;
 266 380;
 381 499;
 500 627;
 628 747;
 748 870;
 871 988;
 989 1104];
n_blocks = size(block_trial_boundaries,1);

%% Assemble stim, Robs, LFP into trial structure
Ntrials = length(AllExpt.Expt.Trials);
LFP_Fs = 1.0 / AllExpt.Expt.Header.LFPsamplerate;  % this is 1 kHz 
dt = 0.005;  % assuming 10 ms frame time exactly
trialDur = 2; %in sec.
Nframes = 200; %number of video frames (@100Hz) per trial
NT = ceil(trialDur/dt); %number of time bins per trial

% SU and MU data
Nunits = length(AllExpt.Header);
SUindx = []; MUindx = [];  SUchans = []; MUchans = [];  SUcellnumbers = [];
for nn = 1:Nunits
  if AllExpt.Header(nn).cellnumber > 0
    SUcellnumbers(end+1) = AllExpt.Header(nn).cellnumber;
    SUindx(end+1) = nn;
    SUchans(end+1) = AllExpt.Header(nn).probe;
  else
    MUindx(end+1) = nn;
    MUchans(end+1) = AllExpt.Header(nn).probe;
  end    
end
SUchans = round(SUchans);

%%
trialNums = [AllExpt.Expt.Trials(:).Trial];
trialOR = [AllExpt.Expt.Trials(:).or];
trialOB = [AllExpt.Expt.Trials(:).ob];
trialrwDir = [AllExpt.Expt.Trials(:).rwdir];
trialrespDir = [AllExpt.Expt.Trials(:).RespDir];
trialSe = [AllExpt.Expt.Trials(:).se];

%% get total trial-by-trial spike counts for each SU
tot_spks_per_trial = zeros(Ntrials,length(SUindx));
for tr = 1:Ntrials
    for nn = 1:length(SUindx)
        trindx = find( AllExpt.Spikes{SUindx(nn)}.Trial == trialNums(tr));
        if ~isempty(trindx)
           tot_spks_per_trial(tr,nn) = length(AllExpt.Spikes{SUindx(nn)}.Spikes{trindx}); 
        end
    end
end

tot_spks_per_trial_MU = zeros(Ntrials,length(MUindx));
for tr = 1:Ntrials
    for nn = 1:length(MUindx)
        trindx = find( AllExpt.Spikes{MUindx(nn)}.Trial == trialNums(tr));
        if ~isempty(trindx)
           tot_spks_per_trial_MU(tr,nn) = length(AllExpt.Spikes{MUindx(nn)}.Spikes{trindx}); 
        end
    end
end

%find any blocks where the unit has no spikes and set these values to nans
spikes_per_block = nan(n_blocks,length(SUindx));
spikes_per_block_MU = nan(n_blocks,length(MUindx));
for bb = 1:n_blocks
   spikes_per_block(bb,:) = sum(tot_spks_per_trial(block_trial_boundaries(bb,1):block_trial_boundaries(bb,2),:)); 
   spikes_per_block_MU(bb,:) = sum(tot_spks_per_trial_MU(block_trial_boundaries(bb,1):block_trial_boundaries(bb,2),:)); 
end

tot_spks_per_trial_norm = tot_spks_per_trial;
for cc = 1:length(SUindx)
   bad_blocks = find(spikes_per_block(:,cc) == 0);
   for ii = bad_blocks'
       tot_spks_per_trial_norm(block_trial_boundaries(ii,1):block_trial_boundaries(ii,2),cc) = nan;
   end
end
 
tot_spks_per_trial_norm_MU = tot_spks_per_trial_MU;
for cc = 1:length(MUindx)
   bad_blocks = find(spikes_per_block_MU(:,cc) == 0);
   for ii = bad_blocks'
       tot_spks_per_trial_norm_MU(block_trial_boundaries(ii,1):block_trial_boundaries(ii,2),cc) = nan;
   end
end

%% USE WAVELET ANALYSIS TO COMPUTE PHASE-LOCKING SPECTRA FOR EACH UNIT
LFP_trial_taxis = AllExpt.Expt.Header.LFPtimes*1e-4;
R_trial_taxis = (0:(NT-1))*dt + dt/2;

LFP_dsf = 2;
LFP_Fsd = LFP_Fs/LFP_dsf;
%anti-aliasing filter
aa_hcf = LFP_Fsd/2*0.8;
[b_aa,a_aa] = butter(4,aa_hcf/(LFP_Fs/2),'low');

LFP_trial_taxis_ds = downsample(LFP_trial_taxis,LFP_dsf);

%wavelet parameters
nwfreqs = 15;
min_freq = 5; max_freq = 80;
min_scale = 1/max_freq*LFP_Fsd;
max_scale = 1/min_freq*LFP_Fsd;
wavetype = 'cmor1-1';
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,wavetype,1/LFP_Fsd);

nprobes = 24;
uprobes = 1:3:nprobes;
all_spk_id = [];
all_spk_phases = nan(sum(tot_spks_per_trial(:)),length(wfreqs),length(uprobes));
% trial_LFP_phases = nan(NT,Ntrials,length(wfreqs),length(uprobes));
cnt = 0;
for tr = 1:Ntrials
    fprintf('Trial %d of %d\n',tr,Ntrials);
    cur_LFPs = double(AllExpt.Expt.Trials(tr).LFP(:,uprobes));
    cur_LFPs(isnan(cur_LFPs)) = 0;
    cur_LFPs = filtfilt(b_aa,a_aa,cur_LFPs);
    cur_LFPs = downsample(cur_LFPs,LFP_dsf);
    
    cur_cwt = nan(size(cur_LFPs,1),length(wfreqs),length(uprobes));
    for cc = 1:length(uprobes)
        cur_cwt(:,:,cc) = cwt(cur_LFPs(:,cc),scales,wavetype)';
    end
    cur_phasegram = angle(cur_cwt);
    
    cur_spks = [];
    for nn = 1:length(SUindx)
        trindx = find( AllExpt.Spikes{SUindx(nn)}.Trial == trialNums(tr));
        if ~isempty(trindx)
            spks = double( AllExpt.Spikes{SUindx(nn)}.Spikes{trindx} ) * 1e-4;  % no adjustment from start of trial
            cur_spks = cat(1,cur_spks,[spks ones(length(spks),1)*nn]);
        end
    end
    cur_spk_inds = round(interp1(LFP_trial_taxis_ds,1:length(LFP_trial_taxis_ds),cur_spks(:,1)));
    all_spk_phases(cnt + (1:length(cur_spk_inds)),:,:) = cur_phasegram(cur_spk_inds,:,:);
    %     all_spk_phases = cat(1,all_spk_phases,cur_phasegram(cur_spk_inds,:,:));
    all_spk_id = cat(1,all_spk_id,cur_spks(:,2));
    cnt = cnt + length(cur_spk_inds);
end

%% USE WAVELET ANALYSIS TO COMPUTE PHASE-LOCKING SPECTRA FOR EACH UNIT
LFP_trial_taxis = AllExpt.Expt.Header.LFPtimes*1e-4;

LFP_dsf = 2;
LFP_Fsd = LFP_Fs/LFP_dsf;
%anti-aliasing filter and high-pass filter
aa_hcf = LFP_Fsd/2*0.8;
% [b_aa,a_aa] = butter(4,aa_hcf/(LFP_Fs/2),'low');
aa_lcf = 1;
[b_aa,a_aa] = butter(4,[aa_lcf aa_hcf]/(LFP_Fs/2));

LFP_trial_taxis_ds = downsample(LFP_trial_taxis,LFP_dsf);

%wavelet parameters
nwfreqs = 12;
min_freq = 2; max_freq = 70;
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
udata = find(ismember(RR,utrials));
xv_frac = 0.;
xv_trials = randperm(length(utrials));
xv_trials = xv_trials(1:round(xv_frac*length(utrials)));
tr_trials = setdiff(1:length(utrials),xv_trials);

trdata = find(ismember(RR(udata),utrials(tr_trials)));
xvdata = find(ismember(RR(udata),utrials(xv_trials)));

lambda_L2 = 5e3;
mod_reg_params = NMMcreate_reg_params('lambda_L2',lambda_L2);
silent = 1;

uchs = 1:length(uprobes);
[CC,WW] = meshgrid(1:length(uprobes),1:length(wfreqs));
ch_use = find(ismember(CC(:),uchs));
mod_stim_params(1:2) = NMMcreate_stim_params(length(wfreqs)*length(uchs));

X{1} = trial_LFP_amps(udata,ch_use).*cos(trial_LFP_phases(udata,ch_use));
X{2} = trial_LFP_amps(udata,ch_use).*sin(trial_LFP_phases(udata,ch_use));

all_prate_outs = nan(length(udata),length(SUindx));
for cc = 1:length(SUindx)
    cc
    cur_Robs = reshape(fullRobs((beg_buffer+1):end,:,cc),[],1);
    cur_Robs = cur_Robs(udata);
    uinds = find(~isnan(cur_Robs));
    truinds = trdata(~isnan(cur_Robs(trdata)));
    xvuinds = xvdata(~isnan(cur_Robs(xvdata)));
    if ~isempty(uinds)
        lfp_mod(cc) = NMMinitialize_model(mod_stim_params,[1 1],{'lin','lin'},mod_reg_params,[1 2]);
%         lfp_mod(cc) = NMMfit_filters(lfp_mod(cc),cur_Robs,X,[],uinds,silent);
        lfp_mod(cc) = NMMfit_filters(lfp_mod(cc),cur_Robs,X,[],truinds,silent);
        [LL, nullLL, pred_rate] = NMMeval_model( lfp_mod(cc), cur_Robs, X, [], truinds);
        [xvLL, nullxvLL, pred_rate] = NMMeval_model( lfp_mod(cc), cur_Robs, X, [], xvuinds);
%         [LL,~,~,~,~,~,nullLL] = NMMmodel_eval( lfp_mod(cc), cur_Robs(truinds),get_Xcell_tInds(X,truinds));
%         [xvLL,~,~,~,~,~,nullxvLL] = NMMmodel_eval( lfp_mod(cc), cur_Robs(xvuinds),get_Xcell_tInds(X,xvuinds));
        lfp_mod_LLimp(cc) = (LL - nullLL)/log(2);
        lfp_mod_xvLLimp(cc) = (xvLL - nullxvLL)/log(2);
%         cur_filts = reshape([lfp_mod(cc).mods(:).filtK],[],1);
%         cur_filts = Xcov*cur_filts;
        [~, ~, all_prate_outs(:,cc)] = NMMeval_model( lfp_mod(cc), cur_Robs, X);
    else
        lfp_mod(cc) = [];
    end
end

%%
all_prate_outs_ms = bsxfun(@minus,all_prate_outs,nanmean(all_prate_outs));
all_prate_outs_ms = reshape(all_prate_outs_ms,[TLEN length(utrials) length(SUindx)]);

all_sum_outs = squeeze(nansum(all_prate_outs_ms/dt));
zscore_sumouts = nanmean(nanzscore(all_sum_outs),2);
bad_trials = find(zscore_sumouts >4);
utrials(bad_trials) = [];
all_prate_outs_ms(:,bad_trials,:) = [];
all_sum_outs(bad_trials,:) = [];

unique_ses = unique(trialSe(utrials));
samestim_covmat = zeros(length(SUindx),length(SUindx));
samestim_sumcovmat = zeros(length(SUindx),length(SUindx));
cnt = 0;
for ss = 1:length(unique_ses)
   curset = find(trialSe(utrials) == unique_ses(ss)); 
   if length(curset) >= 2
   temp_prates = squeeze(all_prate_outs_ms(:,curset,:));
   temp_prates2 = reshape(temp_prates,TLEN,[],1,length(SUindx));
   
   [II,JJ] = meshgrid(1:length(curset));
    uset = find(II > JJ);
    temp = squeeze(nansum(nanmean(bsxfun(@times,temp_prates(:,II(uset),:),temp_prates2(:,JJ(uset),:,:)),1),2));
    samestim_covmat = samestim_covmat + temp;

   temp_prates = squeeze(all_sum_outs(curset,:));
   temp_prates2 = reshape(temp_prates,[],1,length(SUindx));
    temp = squeeze(nansum(bsxfun(@times,temp_prates(II(uset),:),temp_prates2(JJ(uset),:,:)),1));
    samestim_sumcovmat = samestim_sumcovmat + temp;

    cnt = cnt + length(uset);
   end
end
samestim_covmat = samestim_covmat/cnt;
samestim_sumcovmat = samestim_sumcovmat/cnt;

cur_Robs = reshape(fullRobs((beg_buffer+1:end),utrials,:),[],length(SUindx));
ov_LFP_covmat = cov(all_prate_outs);
LFP_noise_covmat = ov_LFP_covmat - samestim_covmat;

ov_LFP_sumcovmat = cov(all_sum_outs);
LFP_noise_sumcovmat = ov_LFP_sumcovmat - samestim_sumcovmat;

% LFP_normfac = sqrt(var(all_prate_outs)'*var(all_prate_outs));
LFP_normfac = sqrt(nanvar(cur_Robs)'*nanvar(cur_Robs));
LFP_noise_corrmat = LFP_noise_covmat./LFP_normfac;

LFP_sumnormfac = sqrt(nanvar(all_sum_outs)'*nanvar(all_sum_outs));
LFP_noise_sumcorrmat = LFP_noise_sumcovmat./LFP_sumnormfac;
%%
XX = [trial_LFP_amps.*cos(trial_LFP_phases) trial_LFP_amps.*sin(trial_LFP_phases)];
YY = reshape(fullRobs(beg_buffer:end,:,:),[],length(SUindx));
for ii = 1:length(SUindx)
    YY(isnan(YY(:,ii)),ii) = nanmean(YY(:,ii));
end
bad_units = [13 14];
YY(:,bad_units) = [];
YY = nanzscore(sqrt(YY));

%%
[A,B,R] = canoncorr(XX,YY);
%%
close all
nmodes = length(R);
for ii = 1:nmodes
%    temp = reshape(XXcov*A(:,ii),length(wfreqs)*length(uprobes),2);
   temp = reshape(A(:,ii),length(wfreqs)*length(uprobes),2);
   temp = sqrt(temp(:,1).^2 + temp(:,2).^2);
   amp_spec = reshape(temp,length(wfreqs),length(uprobes));
   pcolor(wfreqs,uprobes,amp_spec');shading flat
   set(gca,'xscale','log');
   pause
   clf
end
%%
all_phase_lock = nan(length(SUindx),length(wfreqs),length(uprobes));
for nn = 1:length(SUindx)
   cur_spk_set = find(all_spk_id == nn);
   for cc = 1:length(uprobes)
       for ww = 1:length(wfreqs)
           all_phase_lock(nn,ww,cc) = circ_kappa(all_spk_phases(cur_spk_set,ww,cc));
       end
   end    
end

%%
close all
for nn = 1:length(SUindx)
    fprintf('SU from probe %d\n',SUchans(nn));
    pcolor(wfreqs,uprobes,squeeze(all_phase_lock(nn,:,:))'); shading flat; colorbar
    pause
    clf
end
%% Get Binned Spikes
fullRobs = zeros(NT,Ntrials,length(SUindx));
for tr = 1:Ntrials
    for nn = 1:length(SUindx)
        trindx = find( AllExpt.Spikes{SUindx(nn)}.Trial == trialNums(tr));
        if ~isempty(trindx)
            spks = double( AllExpt.Spikes{SUindx(nn)}.Spikes{trindx} ) * 1e-4;  % no adjustment from start of trial
            cur_Robs = histc( spks, (0:(NT))*dt);
            fullRobs(:,tr,nn) = cur_Robs(1:end-1);
        end
    end
end
for cc = 1:length(SUindx)
    fullRobs(:,isnan(tot_spks_per_trial_norm(:,cc)),cc) = nan;
end

trialSpkCnts = squeeze(nansum(fullRobs));
totSpkCnts = squeeze(nansum(trialSpkCnts));
avg_spk_rates = nanmean(reshape(fullRobs,[],length(SUindx)));
% bad_trials = (tot_spks_per_trial == 0);
% bad_trials = permute(repmat(bad_trials,[1 1 NT]),[3 1 2]);
% fullRobs(bad_trials) = nan;
fullRobs_ms = bsxfun(@minus,fullRobs,reshape(avg_spk_rates,[1 1 length(SUindx)]));

%%
clear LFP_explained
% poss_tchans = [1:2:24];
% LFP_bands = [2 4; 5 10;10 20; 20 40;40 80];
poss_tchans = [20];
LFP_bands = [40 80];
for bbb = 1:size(LFP_bands,1)
    fprintf('Band %d of %d\n',bbb,size(LFP_bands,1));
trial_LFP_state = nan(NT,Ntrials,nprobes);

%5-15 looks good, see similar structure in choice- and LFP-explained cov
%mats
% lcf = 6; hcf = 15; %cut-off freqs for LFP filter
lcf = LFP_bands(bbb,1); hcf = LFP_bands(bbb,2); %cut-off freqs for LFP filter
[filt_bb,filt_aa] = butter(2,[lcf hcf]/(LFP_Fs/2));
for tr = 1:Ntrials
    cur_LFPs = double(AllExpt.Expt.Trials(tr).LFP);
    cur_LFPs(isnan(cur_LFPs)) = 0;
    cur_LFPs = filtfilt(filt_bb,filt_aa,cur_LFPs);
    cur_phase = unwrap(angle(hilbert(cur_LFPs)));
    cur_phase_interp = mod(interp1(LFP_trial_taxis,cur_phase,R_trial_taxis),2*pi);
    trial_LFP_state(:,tr,:) = cur_phase_interp;
end

% rand_tperm = randperm(Ntrials);
% trial_LFP_state = trial_LFP_state(:,rand_tperm,:);
%% LFP state-explained covariance
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
% test_trials = find(trialrespDir ~= 0);
un_stim_seeds = unique(trialSe(test_trials));
for ppp = 1:length(poss_tchans)
    fprintf('Chan %d of %d\n',ppp,length(poss_tchans));
%targChan = 1;
targChan = poss_tchans(ppp);

beg_buff = 15; %number of bins from beginning of trial to exclude

maxlag_ED = pi;
ED_space = pi/10;
ED_bin_edges = 0:ED_space:maxlag_ED;
ED_bin_centers = (ED_bin_edges(1:end-1)+ED_bin_edges(2:end))/2;

%subtract off avg rates for times within this set of trials
[RR,TT] = meshgrid(1:Ntrials,1:NT);
istrain = (ismember(RR(:),test_trials)) & (TT(:) > beg_buff);
fullRobs_resh = reshape(fullRobs,[],length(SUindx));
fullRobs_ms = bsxfun(@minus,fullRobs,reshape(nanmean(fullRobs_resh(istrain,:)),[1 1 length(SUindx)]));
tot_var = nanvar(fullRobs_resh(istrain,:));

% % cur_XC = zeros(length(ED_bin_centers),length(SUindx));
% cur_XC = zeros(length(ED_bin_centers),length(SUindx),length(SUindx));
% cur_cnt = zeros(length(ED_bin_centers),length(SUindx));
% % rand_XC = zeros(1,length(SUindx));
% rand_XC = zeros(length(SUindx),length(SUindx));
% rand_cnt = zeros(1,length(SUindx));
% for tr = 1:length(un_stim_seeds)
%     cur_tr_set = test_trials(trialSe(test_trials) == un_stim_seeds(tr));
%     fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
%     if length(cur_tr_set) >= 2
%         [II,JJ] = meshgrid(1:length(cur_tr_set));
%         
%         cur_LFP_data = squeeze(trial_LFP_state(:,cur_tr_set,targChan));
%         for tt = (beg_buff+1):NT
%             cur_Robs = squeeze(fullRobs_ms(tt,cur_tr_set,:));
%             cur_Robs2 = reshape(cur_Robs,[],1,length(SUindx));
% %             cur_Robs = squeeze(fullRobs(tt,cur_tr_set,:));
%             cur_Dmat = abs(circ_dist2(cur_LFP_data(tt,:),cur_LFP_data(tt,:)));
% %             cur_Dmat(logical(eye(length(cur_tr_set)))) = nan;
%             cur_Dmat(JJ >= II) = nan;
%             for jj = 1:length(ED_bin_centers)
%                 curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
% %                 cur_XC(jj,:) = cur_XC(jj,:) + squeeze(nansum(bsxfun(@times,cur_Robs(II(curset),:),cur_Robs(JJ(curset),:)),1));
%                 cur_XC(jj,:,:) = cur_XC(jj,:,:) + (nansum(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
% %                 cur_XC(jj,:) = cur_XC(jj,:) + squeeze(nansum(cur_Robs(II(curset),:).*cur_Robs(JJ(curset),:),1));
%                 cur_cnt(jj,:) = cur_cnt(jj,:) + sum(~(isnan(cur_Robs(II(curset),:))) & ~(isnan(cur_Robs(JJ(curset),:))),1);
%             end
%             curset = ~isnan(cur_Dmat);
% %             rand_XC = rand_XC + squeeze(nansum(bsxfun(@times,cur_Robs(II(curset),:),cur_Robs(JJ(curset),:)),1));
%             rand_XC = rand_XC + squeeze(nansum(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
%             rand_cnt = rand_cnt + sum(~(isnan(cur_Robs(II(curset),:))) & ~(isnan(cur_Robs(JJ(curset),:))),1);
%         end
%     end
% end
% % stateDepVar = cur_XC./cur_cnt;
% stateDepVar = bsxfun(@rdivide,cur_XC,cur_cnt);
% % randVar = rand_XC./rand_cnt;
% randVar = bsxfun(@rdivide,rand_XC,rand_cnt);


cur_tr_set = test_trials;
fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
[II,JJ] = meshgrid(1:length(cur_tr_set));

cur_LFP_data = squeeze(trial_LFP_state(:,cur_tr_set,targChan));
cur_XC = nan(NT,length(ED_bin_centers),length(SUindx),length(SUindx));
rand_XC = nan(NT,length(SUindx),length(SUindx));
for tt = (beg_buff+1):NT
    tt
    cur_Robs = squeeze(fullRobs_ms(tt,cur_tr_set,:));
    cur_Robs2 = reshape(cur_Robs,[],1,length(SUindx));
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

end
end
normfac = sqrt(tot_var'*tot_var);
LFP_covmat = squeeze(LFP_explained)./normfac;
%% CALCULATE CHOICE PROBS
tot_spks_per_trial_norm = tot_spks_per_trial;
tot_spks_per_trial_norm(tot_spks_per_trial == 0) = nan;

sig_trials = find(trialOB < 130);
stim_up = sig_trials(trialrwDir(sig_trials)==1);
stim_down = sig_trials(trialrwDir(sig_trials)==-1);
sig_prob = nan(length(SUindx),1);
for cc = 1:length(SUindx)
   cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(sig_trials,cc)));
   up_hist = histc(tot_spks_per_trial_norm(stim_up,cc),cur_resp_ax);
   down_hist = histc(tot_spks_per_trial_norm(stim_down,cc),cur_resp_ax);
   fract_cor = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_up,cc)));
   fract_incor = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_down,cc)));
   sig_prob(cc) = trapz(fract_incor,fract_cor);
end

% test_trials = find(trialOB == 130 & trialrespDir ~= 0);
test_trials = find(trialrespDir ~= 0);
resp_up = test_trials(trialrespDir(test_trials) == 1);
resp_down = test_trials(trialrespDir(test_trials) == -1);
choice_prob = nan(length(SUindx),1);
cp_pval = nan(length(SUindx),1);
for cc = 1:length(SUindx)
   cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(test_trials,cc)));
   up_hist = histc(tot_spks_per_trial_norm(resp_up,cc),cur_resp_ax);
   down_hist = histc(tot_spks_per_trial_norm(resp_down,cc),cur_resp_ax);
   true_pos = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(resp_up,cc)));
   false_pos = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(resp_down,cc)));
   choice_prob(cc) = trapz(false_pos,true_pos);
   [~,cp_pval(cc)] = ttest2(tot_spks_per_trial_norm(resp_up,cc),tot_spks_per_trial_norm(resp_down,cc));
end

nboot = 500;
boot_sig_prob = nan(length(SUindx),nboot);
boot_choice_prob = nan(length(SUindx),nboot);
for nn = 1:nboot
    randrwdir = trialrwDir(randperm(Ntrials));
    stim_up = sig_trials(randrwdir(sig_trials)==1);
    stim_down = sig_trials(randrwdir(sig_trials)==-1);
    
    randrespdir = trialrespDir(randperm(Ntrials));
    resp_up = test_trials(randrespdir(test_trials) == 1);
    resp_down = test_trials(randrespdir(test_trials) == -1);
    for cc = 1:length(SUindx)
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
%% CALCULATE CHOICE-predictable covariance
test_trials = find(trialOB == 130);
% test_trials = 1:Ntrials;
un_stim_seeds = unique(trialSe(test_trials));

beg_buff = 15; %number of bins from beginning of trial to exclude
randrespDir = trialrespDir(randperm(Ntrials));

%subtract off avg rates for times within this set of trials
[RR,TT] = meshgrid(1:Ntrials,1:NT);
istrain = (ismember(RR(:),test_trials)) & (TT(:) > beg_buff);
fullRobs_resh = reshape(fullRobs,[],length(SUindx));
fullRobs_ms = bsxfun(@minus,fullRobs,reshape(nanmean(fullRobs_resh(istrain,:)),[1 1 length(SUindx)]));

cur_XC = zeros(length(SUindx),length(SUindx));
cur_cnt = zeros(length(SUindx),1);
cur_XC2 = zeros(length(SUindx),length(SUindx));
cur_cnt2 = zeros(length(SUindx),1);
rand_XC = zeros(length(SUindx),length(SUindx));
rand_cnt = zeros(1,length(SUindx));
for tr = 1:length(un_stim_seeds)
    cur_tr_set = find(trialSe == un_stim_seeds(tr) & ismember(1:length(trialOB),test_trials));
    fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
    if length(cur_tr_set) >= 2
        cur_resp = trialrespDir(cur_tr_set);
%         cur_resp = randrespDir(cur_tr_set);
        
        [II,JJ] = meshgrid(1:length(cur_tr_set));
        
        for tt = (beg_buff+1):NT
            cur_Robs = squeeze(fullRobs_ms(tt,cur_tr_set,:));
            cur_Robs2 = reshape(cur_Robs,[],1,length(SUindx));
            
            cur_Dmat = abs(squareform(pdist(cur_resp')));
            cur_Dmat(logical(eye(length(cur_tr_set)))) = nan;
            cur_Dmat(JJ > II) = nan;
            
            curset = find(cur_Dmat == 0);
            cur_XC = cur_XC + squeeze(nansum(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
            cur_cnt = cur_cnt + sum(~(isnan(cur_Robs(II(curset),:))) & ~(isnan(cur_Robs(JJ(curset),:))),1)';
            
            curset = find(cur_Dmat == 2);
            cur_XC2 = cur_XC2 + squeeze(nansum(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
            cur_cnt2 = cur_cnt2 + sum(~(isnan(cur_Robs(II(curset),:))) & ~(isnan(cur_Robs(JJ(curset),:))),1)';

            curset = ~isnan(cur_Dmat);
            rand_XC = rand_XC + squeeze(nansum(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
            rand_cnt = rand_cnt + sum(~(isnan(cur_Robs(II(curset),:))) & ~(isnan(cur_Robs(JJ(curset),:))),1);
        end
    end
end
sameVar = bsxfun(@rdivide,cur_XC,cur_cnt);
oppVar = bsxfun(@rdivide,cur_XC2,cur_cnt2);
randVar = bsxfun(@rdivide,rand_XC,rand_cnt);

choice_explained = sameVar - randVar;
SU_choice_varfrac = diag(choice_explained)./diag(randVar);

%%
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
% test_trials = find(trialrespDir~= 0);
un_stim_seeds = unique(trialSe(test_trials));

tot_spks_per_trial_ms = bsxfun(@minus,tot_spks_per_trial_norm,nanmean(tot_spks_per_trial_norm(test_trials,:)));
tot_spks2 = reshape(tot_spks_per_trial_ms,[],1,length(SUindx));

all_same_covmat = nan(length(un_stim_seeds),length(SUindx),length(SUindx));
all_diff_covmat = nan(length(un_stim_seeds),length(SUindx),length(SUindx));
all_rand_covmat = nan(length(un_stim_seeds),length(SUindx),length(SUindx));
for tr = 1:length(un_stim_seeds)
    cur_tr_set = test_trials(trialSe(test_trials) == un_stim_seeds(tr));
    fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
    if length(cur_tr_set) >= 2
        cur_resp = trialrespDir(cur_tr_set);
        %         cur_resp = randrespDir(cur_tr_set);
        
        [II,JJ] = meshgrid(1:length(cur_tr_set));
        cur_Dmat = abs(squareform(pdist(cur_resp')));
        cur_Dmat(logical(eye(length(cur_tr_set)))) = nan;
        cur_Dmat(JJ > II) = nan;
        curset = find(cur_Dmat == 0);
        all_same_covmat(tr,:,:) = nanmean(bsxfun(@times,tot_spks_per_trial_ms(cur_tr_set(II(curset)),:),tot_spks2(cur_tr_set(JJ(curset)),:,:)),1);

        curset = find(cur_Dmat == 2);
        all_diff_covmat(tr,:,:) = nanmean(bsxfun(@times,tot_spks_per_trial_ms(cur_tr_set(II(curset)),:),tot_spks2(cur_tr_set(JJ(curset)),:,:)),1);

        curset = ~isnan(cur_Dmat);
        all_rand_covmat(tr,:,:) = nanmean(bsxfun(@times,tot_spks_per_trial_ms(cur_tr_set(II(curset)),:),tot_spks2(cur_tr_set(JJ(curset)),:,:)),1);
    end
end
avg_same_covmat = squeeze(nanmean(all_same_covmat));
avg_diff_covmat = squeeze(nanmean(all_diff_covmat));
avg_rand_covmat = squeeze(nanmean(all_rand_covmat));

%%
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
% test_trials = find(trialrespDir~= 0);

tot_spks_per_trial_ms = bsxfun(@minus,tot_spks_per_trial_norm,nanmean(tot_spks_per_trial_norm(test_trials,:)));
tot_spks2 = reshape(tot_spks_per_trial_ms,[],1,length(SUindx));

cur_resp = trialrespDir(test_trials);
[II,JJ] = meshgrid(1:length(test_trials));
cur_Dmat = abs(squareform(pdist(cur_resp')));
cur_Dmat(logical(eye(length(test_trials)))) = nan;
cur_Dmat(JJ > II) = nan;

curset = find(cur_Dmat == 0);
avg_same_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(test_trials(II(curset)),:),tot_spks2(test_trials(JJ(curset)),:,:)),1));

% curset = find(cur_Dmat == 2);
curset = ~isnan(cur_Dmat);
avg_rand_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(test_trials(II(curset)),:),tot_spks2(test_trials(JJ(curset)),:,:)),1));

curset = find(cur_Dmat == 2);
avg_diff_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(test_trials(II(curset)),:),tot_spks2(test_trials(JJ(curset)),:,:)),1));

avg_choice_covmat = avg_same_covmat - avg_rand_covmat;
ov_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(test_trials,:),tot_spks2(test_trials,:,:)),1));
normfac = sqrt(nanvar(tot_spks_per_trial_ms(test_trials,:))'*nanvar(tot_spks_per_trial_ms(test_trials,:)));
ov_covmat_norm = avg_choice_covmat./normfac;

%%
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
% test_trials = find(trialrespDir~= 0);

tot_spks_per_trial_ms = bsxfun(@minus,tot_spks_per_trial_norm,nanmean(tot_spks_per_trial_norm(test_trials,:)));
tot_spks2 = reshape(tot_spks_per_trial_ms,[],1,length(SUindx));

cur_resp = trialrespDir(test_trials);
resp_up = test_trials(cur_resp == 1);
resp_down = test_trials(cur_resp == -1);

[II,JJ] = meshgrid(1:length(resp_up));
uset = find(II > JJ);
avg_up_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(resp_up(II(uset)),:),tot_spks2(resp_up(JJ(uset)),:,:)),1));

[II,JJ] = meshgrid(1:length(resp_down));
uset = find(II > JJ);
avg_down_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(resp_down(II(uset)),:),tot_spks2(resp_down(JJ(uset)),:,:)),1));

avg_choice_covmat = 0.5*avg_down_covmat + 0.5*avg_up_covmat;
ov_covmat = squeeze(nanmean(bsxfun(@times,tot_spks_per_trial_ms(test_trials,:),tot_spks2(test_trials,:,:)),1));
normfac = sqrt(nanvar(tot_spks_per_trial_ms(test_trials,:))'*nanvar(tot_spks_per_trial_ms(test_trials,:)));
ov_covmat_norm = avg_choice_covmat./normfac;

%%
% test_trials = find(trialOB == 130);
test_trials = 1:Ntrials;
cur_Robs = fullRobs_ms(:,test_trials,:);
cur_Robs = reshape(cur_Robs,[],length(SUindx));
emp_xcovs = squeeze(nanmean(bsxfun(@times,cur_Robs,reshape(cur_Robs,[],1,length(SUindx)))));
noise_covs = emp_xcovs - randVar;
noise_covs(logical(eye(length(SUindx)))) = nan;
%%
for nn = 1:length(SUindx)
%     pcolor(wfreqs,1:length(uprobes),squeeze(all_phase_lock(nn,:,:))'); shading flat; colorbar
    figure
    plot(ED_bin_centers,stateDepVar(:,nn,nn))
    
    xl = xlim();
    line(xl,randVar(nn,nn) + [0 0],'color','r')
    pause
    close all
end
