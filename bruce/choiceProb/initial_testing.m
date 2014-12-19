%% LOAD PROCESSED DATA
clear all
% cd('/home/james/Data/bruce/ChoiceProb/')
cd('~/Data/bruce/ChoiceProb/')

load('M239/lemM239.image.ORBW.LFP.mat')

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

%%
tot_spks_per_trial = zeros(Ntrials,length(SUindx));
for tr = 1:Ntrials
    for nn = 1:length(SUindx)
        trindx = find( AllExpt.Spikes{SUindx(nn)}.Trial == trialNums(tr));
        if ~isempty(trindx)
           tot_spks_per_trial(tr,nn) = length(AllExpt.Spikes{SUindx(nn)}.Spikes{trindx}); 
        end
    end
end
%%
LFP_trial_taxis = AllExpt.Expt.Header.LFPtimes*1e-4;
R_trial_taxis = (0:(NT-1))*dt + dt/2;

LFP_dsf = 2;
LFP_Fsd = LFP_Fs/LFP_dsf;
aa_hcf = LFP_Fsd/2*0.8;
[b_aa,a_aa] = butter(4,aa_hcf/(LFP_Fs/2),'low');

LFP_trial_taxis_ds = downsample(LFP_trial_taxis,LFP_dsf);

nwfreqs = 10;
min_freq = 2; max_freq = 60;
min_scale = 1/max_freq*LFP_Fsd;
max_scale = 1/min_freq*LFP_Fsd;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/LFP_Fsd);

nprobes = 24;
uprobes = 1:2:nprobes;
all_spk_phases = [];
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
        cur_cwt(:,:,cc) = cwt(cur_LFPs(:,cc),scales,'cmor1-1')';
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
for nn = 1:length(SUindx)
    pcolor(wfreqs,1:length(uprobes),squeeze(all_phase_lock(nn,:,:))'); shading flat; colorbar
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
totSpkCnts = squeeze(nansum(fullRobs));
avg_spk_rates = nanmean(reshape(fullRobs,[],length(SUindx)));
fullRobs_ms = bsxfun(@minus,fullRobs,reshape(avg_spk_rates,[1 1 length(SUindx)]));

%%
bad_trials = (tot_spks_per_trial == 0);
bad_trials = permute(repmat(bad_trials,[1 1 NT]),[3 1 2]);
fullRobs_ms(bad_trials) = nan;
%%
trial_LFP_state = nan(NT,Ntrials,nprobes);

lcf = 5; hcf = 15; %cut-off freqs for LFP filter
[filt_bb,filt_aa] = butter(2,[lcf hcf]/(LFP_Fs/2));
for tr = 1:Ntrials
    cur_LFPs = double(AllExpt.Expt.Trials(tr).LFP);
    cur_LFPs(isnan(cur_LFPs)) = 0;
    cur_LFPs = filtfilt(filt_bb,filt_aa,cur_LFPs);
    cur_phase = unwrap(angle(hilbert(cur_LFPs)));
    cur_phase_interp = mod(interp1(LFP_trial_taxis,cur_phase,R_trial_taxis),2*pi);
    trial_LFP_state(:,tr,:) = cur_phase_interp;
end

%%
test_trials = find(trialOB == 130);
% test_trials = 1:Ntrials;
un_stim_seeds = unique(trialSe(test_trials));

targChan = 1;

beg_buff = 15; %number of bins from beginning of trial to exclude

maxlag_ED = pi;
ED_space = pi/10;
ED_bin_edges = 0:ED_space:maxlag_ED;
ED_bin_centers = (ED_bin_edges(1:end-1)+ED_bin_edges(2:end))/2;

% cur_XC = zeros(length(ED_bin_centers),length(SUindx));
cur_XC = zeros(length(ED_bin_centers),length(SUindx),length(SUindx));
cur_cnt = zeros(length(ED_bin_centers),length(SUindx));
% rand_XC = zeros(1,length(SUindx));
rand_XC = zeros(length(SUindx),length(SUindx));
rand_cnt = zeros(1,length(SUindx));
for tr = 1:length(un_stim_seeds)
    cur_tr_set = find(trialSe == un_stim_seeds(tr));
    cur_tr_set(trialOB(cur_tr_set) ~= 130) = [];
    fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
    if length(cur_tr_set) >= 2
        [II,JJ] = meshgrid(1:length(cur_tr_set));
        
        cur_LFP_data = squeeze(trial_LFP_state(:,cur_tr_set,targChan));
        for tt = (beg_buff+1):NT
            cur_Robs = squeeze(fullRobs_ms(tt,cur_tr_set,:));
            cur_Robs2 = reshape(cur_Robs,[],1,length(SUindx));
%             cur_Robs = squeeze(fullRobs(tt,cur_tr_set,:));
            cur_Dmat = abs(circ_dist2(cur_LFP_data(tt,:),cur_LFP_data(tt,:)));
            cur_Dmat(logical(eye(length(cur_tr_set)))) = nan;
            cur_Dmat(JJ > II) = nan;
            for jj = 1:length(ED_bin_centers)
                curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
%                 cur_XC(jj,:) = cur_XC(jj,:) + squeeze(nansum(bsxfun(@times,cur_Robs(II(curset),:),cur_Robs(JJ(curset),:)),1));
                cur_XC(jj,:,:) = cur_XC(jj,:,:) + (nansum(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
%                 cur_XC(jj,:) = cur_XC(jj,:) + squeeze(nansum(cur_Robs(II(curset),:).*cur_Robs(JJ(curset),:),1));
                cur_cnt(jj,:) = cur_cnt(jj,:) + sum(~(isnan(cur_Robs(II(curset),:))) & ~(isnan(cur_Robs(JJ(curset),:))),1);
            end
            curset = ~isnan(cur_Dmat);
%             rand_XC = rand_XC + squeeze(nansum(bsxfun(@times,cur_Robs(II(curset),:),cur_Robs(JJ(curset),:)),1));
            rand_XC = rand_XC + squeeze(nansum(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
            rand_cnt = rand_cnt + sum(~(isnan(cur_Robs(II(curset),:))) & ~(isnan(cur_Robs(JJ(curset),:))),1);
        end
    end
end
% stateDepVar = cur_XC./cur_cnt;
stateDepVar = bsxfun(@rdivide,cur_XC,cur_cnt);
% randVar = rand_XC./rand_cnt;
randVar = bsxfun(@rdivide,rand_XC,rand_cnt);

LFP_explained = squeeze(stateDepVar(1,:,:)) - randVar; %dont worry about spline interpolation for now, just take the bin with most similar LFP states
%% CALCULATE CHOICE PROBS
tot_spks_per_trial_norm = tot_spks_per_trial;
tot_spks_per_trial_norm(tot_spks_per_trial == 0) = nan;

test_trials = find(trialOB == 130);
resp_up = test_trials(trialrespDir(test_trials) == 1);
resp_down = test_trials(trialrespDir(test_trials) == -1);

choice_prob = nan(length(SUindx),1);
for cc = 1:length(SUindx)
   cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(test_trials,cc)));
   up_hist = histc(tot_spks_per_trial_norm(resp_up,cc),cur_resp_ax);
   down_hist = histc(tot_spks_per_trial_norm(resp_down,cc),cur_resp_ax);
   fract_cor = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(resp_up,cc)));
   fract_incor = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(resp_down,cc)));
   choice_prob(cc) = trapz(fract_incor,fract_cor);
end

%% CALCULATE CHOICE-predictable covariance
test_trials = find(trialOB == 130);
% test_trials = 1:Ntrials;
% un_stim_seeds = unique(trialSe(test_trials));

targChan = 1;

beg_buff = 15; %number of bins from beginning of trial to exclude

cur_XC = zeros(length(SUindx),length(SUindx));
cur_cnt = zeros(length(SUindx),1);
cur_XC2 = zeros(length(SUindx),length(SUindx));
cur_cnt2 = zeros(length(SUindx),1);
rand_XC = zeros(length(SUindx),length(SUindx));
rand_cnt = zeros(1,length(SUindx));
for tr = 1:length(un_stim_seeds)
    cur_tr_set = find(trialSe == un_stim_seeds(tr));
    cur_tr_set(trialOB(cur_tr_set) ~= 130) = [];
    fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
    if length(cur_tr_set) >= 2
        cur_resp = trialrespDir(cur_tr_set);
        
        [II,JJ] = meshgrid(1:length(cur_tr_set));
        
        cur_LFP_data = squeeze(trial_LFP_state(:,cur_tr_set,targChan));
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
