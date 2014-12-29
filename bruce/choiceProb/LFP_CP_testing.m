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

%% FOR M239
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
TLEN = length(LFP_trial_taxis_ds);

nprobes = 24;
uprobes = 1:1:nprobes;
all_spk_id = [];
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
    trial_LFPs(tr,:,:) = cur_LFPs;
end
trial_LFPs = permute(trial_LFPs,[2 1 3]);

%%
beg_buff = 0.15; end_buff = 0.05; trial_dur = 2;
upts = find(LFP_trial_taxis_ds >= beg_buff & LFP_trial_taxis_ds <= trial_dur-end_buff);

poss_trials = find(trialOB == 130 & trialSe > 0);
% poss_trials = find(trialOB > 50 & trialSe > 0);
% un_seeds = unique(trialSe(poss_trials));
un_seeds = unique([trialSe(poss_trials)' trialOB(poss_trials)'],'rows');
same_pairs = []; diff_pairs = [];
for ii = 1:length(un_seeds)
    cur_trials = find(trialSe(poss_trials) == un_seeds(ii,1) & trialrespDir(poss_trials) ~= 0 ...
    & trialOB(poss_trials) == un_seeds(ii,2));
    cur_resp = trialrespDir(poss_trials(cur_trials));
    
    if length(cur_resp) >= 2
        poss_pairs = nchoosek(1:length(cur_resp),2);
        resp_pairs = cur_resp(poss_pairs);
        same_set = find(diff(resp_pairs,[],2) == 0);
        opp_set = find(diff(resp_pairs,[],2) ~= 0);
        same_pairs = cat(1,same_pairs,cur_trials(poss_pairs(same_set,:)));
        diff_pairs = cat(1,diff_pairs,cur_trials(poss_pairs(opp_set,:)));
    end
end
rand_pairs = randperm(length(poss_trials))';
rand_pairs = [rand_pairs circshift(rand_pairs,[1 0])];

same_OB = trialOB(poss_trials(same_pairs));
diff_OB = trialOB(poss_trials(diff_pairs));
%%
poss_OB = unique(same_OB(:));

params.Fs = LFP_Fsd;
params.tapers = [5 9];
params.trialave = 1;
params.err = [2 0.05];
clear Csame Cdiff
for cc = 1:length(uprobes)
    cc
    for ob = 1:length(poss_OB)
        curset = find(same_OB(:,1) == poss_OB(ob));
        
    data1 = squeeze(trial_LFPs(upts,poss_trials(same_pairs(curset,1)),cc));
    data2 = squeeze(trial_LFPs(upts,poss_trials(same_pairs(curset,2)),cc));
    bad_trials = find(any(isnan(data1),1) | any(isnan(data2),1));
    data1(:,bad_trials) = []; data2(:,bad_trials) = [];
    [Csame(ob,cc,:),~,S12same,S1same,S2same,f]=coherencyc(data1,data2,params);
    
         curset = find(diff_OB(:,1) == poss_OB(ob));
   data1 = squeeze(trial_LFPs(upts,poss_trials(diff_pairs(curset,1)),cc));
    data2 = squeeze(trial_LFPs(upts,poss_trials(diff_pairs(curset,2)),cc));
    bad_trials = find(any(isnan(data1),1) | any(isnan(data2),1));
    data1(:,bad_trials) = []; data2(:,bad_trials) = [];
    [Cdiff(ob,cc,:),~,S12diff,S1diff,S2diff,f]=coherencyc(data1,data2,params);
    
    end
    
    data1 = squeeze(trial_LFPs(upts,poss_trials(rand_pairs(:,1)),cc));
    data2 = squeeze(trial_LFPs(upts,poss_trials(rand_pairs(:,2)),cc));
    bad_trials = find(any(isnan(data1),1) | any(isnan(data2),1));
    data1(:,bad_trials) = []; data2(:,bad_trials) = [];
    [Crand(cc,:),~,S12rand,S1rand,S2rand,f]=coherencyc(data1,data2,params);
end
Csame = squeeze(Csame);
Cdiff = squeeze(Cdiff);
%%
poss_lws = [1 2 1 3];
for cc = 1:length(uprobes)
    plot(f,Csame(cc,:),f,Cdiff(cc,:),'r',f,Crand(cc,:),'k');
% for ob = [1 4]
%     hold on
%     plot(f,squeeze(Csame(ob,cc,:)),f,squeeze(Cdiff(ob,cc,:)),'r','linewidth',poss_lws(ob));
% end
xlim([0 100])
    pause
    clf
end
%%
up_choices = poss_trials(trialrespDir(poss_trials) == 1);
down_choices = poss_trials(trialrespDir(poss_trials) == -1);

params.tapers = [2 3];
for cc = 1:length(uprobes)
    data = squeeze(trial_LFPs(upts,up_choices,cc));
    bad_trials = find(any(isnan(data),1));
    data(:,bad_trials) = [];
    [S_up(cc,:),f]=mtspectrumc(data,params);
    
    data = squeeze(trial_LFPs(upts,down_choices,cc));
    bad_trials = find(any(isnan(data),1));
    data(:,bad_trials) = [];
    [S_down(cc,:),f]=mtspectrumc(data,params);
end
%%
for cc = 1:length(uprobes)
    plot(f,log10(abs(S_up(cc,:))),f,log10(S_down(cc,:)),'r');
    xlim([0 100]);
    pause
    clf
    
end

%%
cc = 8;
        
    data1 = squeeze(trial_LFPs(upts,poss_trials(same_pairs(:,1)),cc));
    data2 = squeeze(trial_LFPs(upts,poss_trials(same_pairs(:,2)),cc));
    bad_trials = find(any(isnan(data1),1) | any(isnan(data2),1));
    data1(:,bad_trials) = []; data2(:,bad_trials) = [];
    [Csame(cc,:),~,S12same,S1same,S2same,f,~,~,Csameerr]=coherencyc(data1,data2,params);
    
   data1 = squeeze(trial_LFPs(upts,poss_trials(diff_pairs(:,1)),cc));
    data2 = squeeze(trial_LFPs(upts,poss_trials(diff_pairs(:,2)),cc));
    bad_trials = find(any(isnan(data1),1) | any(isnan(data2),1));
    data1(:,bad_trials) = []; data2(:,bad_trials) = [];
    [Cdiff(cc,:),~,S12diff,S1diff,S2diff,f,~,~,Cdifferr]=coherencyc(data1,data2,params);

    %%
    close all
    hold on
    plot(f,Csame(cc,:),f,Csameerr,'b--')
    plot(f,Cdiff(cc,:),'r',f,Cdifferr,'r--')