%% LOAD PROCESSED DATA
cd('/home/james/Data/bruce/ChoiceProb/')

load('M239/lemM239.image.ORBW.LFP.mat')

%% Assemble stim, Robs, LFP into trial structure
Ntrials = length(AllExpt.Expt.Trials);
LFP_Fs = 1.0 / AllExpt.Expt.Header.LFPsamplerate;  % this is 1 kHz 
dt = 0.01;  % assuming 10 ms frame time exactly
trialDur = 2; %in sec.
NT = 200; %number of video frames (@100Hz) per trial

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

%% Get Binned Spikes
fullRobs = nan(NT,Ntrials,length(SUindx));
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
LFP_trial_taxis = AllExpt.Expt.Header.LFPtimes*1e-4;
R_trial_taxis = (0:(NT-1))*dt + dt/2;

nprobes = 24;

trial_LFP_state = nan(NT,Ntrials,nprobes);

lcf = 10; hcf = 20; %cut-off freqs for LFP filter
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
% test_trials = find(trialOB == 130);
test_trials = 1:Ntrials;
un_stim_seeds = unique(trialSe(test_trials));

targChan = 22;

beg_buff = 15; %number of bins from beginning of trial to exclude

maxlag_ED = pi/2;
ED_space = pi/20;
ED_bin_edges = 0:ED_space:maxlag_ED;
ED_bin_centers = (ED_bin_edges(1:end-1)+ED_bin_edges(2:end))/2;

cur_XC = zeros(length(ED_bin_centers),length(SUindx));
cur_cnt = zeros(length(ED_bin_centers),length(SUindx));
rand_XC = zeros(1,length(SUindx));
rand_cnt = zeros(1,length(SUindx));
for tr = 1:length(un_stim_seeds)
    cur_tr_set = find(trialSe == un_stim_seeds(tr));
%     cur_tr_set(trialOB(cur_tr_set) ~= 130) = [];
    fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
    if length(cur_tr_set) >= 2
        [II,JJ] = meshgrid(1:length(cur_tr_set));
        
        cur_LFP_data = squeeze(trial_LFP_state(:,cur_tr_set,targChan));
        for tt = (beg_buff+1):NT
            cur_Robs = squeeze(fullRobs_ms(tt,cur_tr_set,:));
%             cur_Robs = squeeze(fullRobs(tt,cur_tr_set,:));
            cur_Dmat = abs(circ_dist2(cur_LFP_data(tt,:),cur_LFP_data(tt,:)));
            cur_Dmat(logical(eye(length(cur_tr_set)))) = nan;
            
            for jj = 1:length(ED_bin_centers)
                curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
                cur_XC(jj,:) = cur_XC(jj,:) + squeeze(nansum(bsxfun(@times,cur_Robs(II(curset),:),cur_Robs(JJ(curset),:)),1));
                cur_cnt(jj,:) = cur_cnt(jj,:) + sum(~(isnan(cur_Robs(II(curset),:)) & ~(isnan(cur_Robs(JJ(curset),:)))));
            end
            curset = ~isnan(cur_Dmat);
            rand_XC = rand_XC + squeeze(nansum(bsxfun(@times,cur_Robs(II(curset),:),cur_Robs(JJ(curset),:)),1));
            rand_cnt = rand_cnt + nansum(~(isnan(cur_Robs(II(curset),:))) & ~(isnan(cur_Robs(JJ(curset),:))));
        end
    end
end
stateDepVar = cur_XC./cur_cnt;
randVar = rand_XC./rand_cnt;

