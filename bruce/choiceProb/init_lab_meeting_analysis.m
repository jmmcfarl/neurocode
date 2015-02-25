%% LOAD PROCESSED DATA
clear all
% cd('/home/james/Data/bruce/ChoiceProb/')
cd('~/Data/bruce/ChoiceProb/')

Expt_name = 'M230';
% Expt_name = 'G037';

if strcmp(Expt_name,'M239')
    load('M239/lemM239.image.ORBW.LFP.mat')
elseif strcmp(Expt_name,'M230')
    load('M230/lemM230.image.ORBW.LFP.mat')
elseif strcmp(Expt_name,'G037')
    load('G037/jbeG037.Cells.image.ORBW.LFP.mat');
end

fig_dir = '/home/james/Desktop/lab_meeting_CP/';

%% SET BLOCK TRIAL RANGES
flip = false;
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
        1241 1374;
        1376 1415;
        1416 1480;
        1481 1612;
        1613 1740];
elseif strcmp(Expt_name,'G037')
    block_trial_boundaries = [1 148;
        149 454;
        455 1361];
    flip = true;
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
spk_delay_buff = 0.05;
tot_spks_per_trial = zeros(Ntrials,Nunits);
for tr = 1:Ntrials
    for nn = 1:Nunits
        trindx = find( AllExpt.Spikes{nn}.Trial == trialNums(tr));
        if ~isempty(trindx)
            cur_st = double(AllExpt.Spikes{nn}.Spikes{trindx})*1e-4;
            tot_spks_per_trial(tr,nn) = sum(cur_st < (trialDur) & cur_st > spk_delay_buff);
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

%% PSYCHOPHYSICAL PERFORMANC
un_OBs = [60 70 80 130 -80 -70 -60];
sig_strength = [3 2 1 0 -1 -2 -3];
fract_one = nan(length(un_OBs),1);
neuromet_avgs = nan(Nunits,length(un_OBs));
for ss = 1:length(un_OBs)
    if flip
    cur_OB_neg = find(trialOB == un_OBs(ss) & trialrespDir ~= 0);
    else
    cur_OB_neg = find(trialOB == abs(un_OBs(ss)) & trialrwDir == sign(un_OBs(ss)) & trialrespDir ~= 0);
    end
    fract_one(ss) = sum(trialrespDir(cur_OB_neg) == 1)/length(cur_OB_neg);
    neuromet_avgs(:,ss) = nanmean(tot_spks_per_trial_norm(cur_OB_neg,:));
end
neuromet_avgs = bsxfun(@rdivide,neuromet_avgs,nanmean(tot_spks_per_trial_norm)');

cmet_avgs = nan(Nunits,length(un_OBs));
for ss = 1:length(un_OBs)
    cur_OB_neg = find(trialOB == abs(un_OBs(ss)) & trialrespDir == sign(un_OBs(ss)) & trialrespDir ~= 0);
    cmet_avgs(:,ss) = nanmean(tot_spks_per_trial_norm(cur_OB_neg,:));
end
cmet_avgs = bsxfun(@rdivide,cmet_avgs,nanmean(tot_spks_per_trial_norm)');

if strcmp(Expt_name,'G037')
    fract_one = flipud(fract_one);
    neuromet_avgs = fliplr(neuromet_avgs);    
end
f1 = figure();
plot(sig_strength,fract_one,'o-');
xlabel('Signal strength');
ylabel('Fraction chose pos');
ylim([0 1]);
xl = xlim();
line(xl,[0.5 0.5],'color','k');

f2 = figure();
errorbar(sig_strength,nanmean(neuromet_avgs),nanstd(neuromet_avgs)/sqrt(Nunits));
hold on
errorbar(sig_strength,nanmean(cmet_avgs),nanstd(cmet_avgs)/sqrt(Nunits),'r');
legend('Stim-based','Choice-based');
% 
% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir sprintf('Psych_perf_%s.pdf',Expt_name)];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir sprintf('Neur_perf_%s.pdf',Expt_name)];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% CALCULATE SIGNAL PROBS
poss_sig_strengths = unique(abs(trialOB));
poss_sig_strengths(poss_sig_strengths == 30) = [];
% poss_sig_strengths(poss_sig_strengths == 130) = [];
all_sig_prob = nan(length(poss_sig_strengths),Nunits);
for ss = 1:length(poss_sig_strengths)
    sig_trials = find(abs(trialOB) == poss_sig_strengths(ss));
    stim_up = sig_trials(trialrwDir(sig_trials)==1);
    stim_down = sig_trials(trialrwDir(sig_trials)==-1);
    for cc = 1:Nunits
        cur_resp_ax = 0:(nanmax(tot_spks_per_trial_norm(sig_trials,cc)));
        if max(cur_resp_ax) > 0
        up_hist = histc(tot_spks_per_trial_norm(stim_up,cc),cur_resp_ax);
        down_hist = histc(tot_spks_per_trial_norm(stim_down,cc),cur_resp_ax);
        fract_cor = cumsum(up_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_up,cc)));
        fract_incor = cumsum(down_hist)/sum(~isnan(tot_spks_per_trial_norm(stim_down,cc)));
        all_sig_prob(ss,cc) = trapz(fract_incor,fract_cor);
        end
    end
end
sig_sign = sign(all_sig_prob(1,:)-0.5);

sig_sign(isnan(sig_sign)) = sign(all_sig_prob(2,isnan(sig_sign))-0.5);

corrected_sig_prob = bsxfun(@times,all_sig_prob-0.5,sig_sign)+0.5;
sig_prob = all_sig_prob(1,:);
%% CALCULATE CHOICE PROBS
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

%%
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
sig_trials = find(abs(trialOB) == 60 & trialrespDir ~= 0);

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
choice_prob_ci = prctile(boot_choice_prob',[2.5 97.5]);
choice_prob_z = (choice_prob - nanmean(boot_choice_prob,2))./nanstd(boot_choice_prob,[],2);
sig_prob_ci = prctile(boot_sig_prob',[2.5 97.5]);
sig_prob_z = (sig_prob' - nanmean(boot_sig_prob,2))./nanstd(boot_sig_prob,[],2);

%%
f1 = figure();
plot(choice_prob,all_sig_prob(1,:),'o')
hold on
% plot(choice_prob,sig_prob(2,:),'r.')
plot(choice_prob,all_sig_prob(3,:),'r.');
line([0 1],[0 1],'color','k');
line([0 1],[0.5 0.5],'color','k','linestyle','--');
line([0.5 0.5],[0 1],'color','k','linestyle','--');
xlabel('Choice prob');
ylabel('Signal prob');

xax = linspace(0,1,100);
sig_dist = histc(sig_prob,xax)/Nunits;
choice_dist = histc(choice_prob,xax)/Nunits;
f2 = figure();
stairs(xax,sig_dist,'r');
hold on
stairs(xax,choice_dist,'b');
yl = ylim();
line([0.5 0.5],yl,'color','k');

% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir sprintf('CP_vs_SP_%s.pdf',Expt_name)];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir sprintf('CP_SP_dists_%s.pdf',Expt_name)];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%%
sp_sign = sign(0.5-sig_prob);
cp_sign = sign(0.5-choice_prob);

test_trials = find(trialOB == 130 & trialrespDir ~= 0);
resp_up = test_trials(trialrespDir(test_trials) == 1);
resp_down = test_trials(trialrespDir(test_trials) == -1);

cup_psths = squeeze(nanmean(fullRobs(:,resp_up,:),2));
cdown_psths = squeeze(nanmean(fullRobs(:,resp_down,:),2));
choice_psth_diff = cup_psths - cdown_psths;

test_trials = find(abs(trialOB) == 60 & trialrespDir ~= 0);
resp_up = test_trials(trialrwDir(test_trials) == 1);
resp_down = test_trials(trialrwDir(test_trials) == -1);

up_psths = squeeze(nanmean(fullRobs(:,resp_up,:),2));
down_psths = squeeze(nanmean(fullRobs(:,resp_down,:),2));
sig_psth_diff = up_psths - down_psths;

avg_rates = nanmean(reshape(fullRobs(:,test_trials,:),[],Nunits));
avg_rates_sp = avg_rates.*sp_sign;
avg_rates_cp = avg_rates.*cp_sign';

choice_psth_diff = bsxfun(@rdivide,choice_psth_diff,avg_rates_cp);
cup_psths = bsxfun(@rdivide,cup_psths,avg_rates);
cdown_psths = bsxfun(@rdivide,cdown_psths,avg_rates);
sig_psth_diff = bsxfun(@rdivide,sig_psth_diff,avg_rates_sp);
up_psths = bsxfun(@rdivide,up_psths,avg_rates);
down_psths = bsxfun(@rdivide,down_psths,avg_rates);


sm_win = round(0.03/dt);
for cc = 1:Nunits
    choice_psth_diff(:,cc) = jmm_smooth_1d_cor(choice_psth_diff(:,cc),sm_win);
    sig_psth_diff(:,cc) = jmm_smooth_1d_cor(sig_psth_diff(:,cc),sm_win);
    cup_psths(:,cc) = jmm_smooth_1d_cor(cup_psths(:,cc),sm_win);
    cdown_psths(:,cc) = jmm_smooth_1d_cor(cdown_psths(:,cc),sm_win);
    up_psths(:,cc) = jmm_smooth_1d_cor(up_psths(:,cc),sm_win);
    down_psths(:,cc) = jmm_smooth_1d_cor(down_psths(:,cc),sm_win);
end

uset = find(avg_rates/dt > 5);
f1 = figure(); hold on
shadedErrorBar((1:NT)*dt,nanmean(choice_psth_diff(:,uset)'),nanstd(choice_psth_diff(:,uset)')/sqrt(length(uset)),{'color','r'});
shadedErrorBar((1:NT)*dt,nanmean(sig_psth_diff(:,uset)'),nanstd(sig_psth_diff(:,uset)')/sqrt(length(uset)),{'color','b'});
%% CALCULATE CHOICE-predictable covariance at dt resolution (controlling for stim-dep)
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
un_stim_seeds = unique(trialSe(test_trials));

% beg_buff = round(0.15/dt); %number of bins from beginning of trial to exclude
beg_buff = round(0.05/dt); %number of bins from beginning of trial to exclude

%subtract off avg rates for times within this set of trials
[RR,TT] = meshgrid(1:Ntrials,1:NT);
istrain = (ismember(RR(:),test_trials)) & (TT(:) > beg_buff);
fullRobs_resh = reshape(fullRobs,[],Nunits);
fullRobs_ms = bsxfun(@minus,fullRobs,reshape(nanmean(fullRobs_resh(istrain,:)),[1 1 Nunits]));

cur_XC = zeros(Nunits,Nunits);
cur_cnt = zeros(Nunits,Nunits);
cur_XC2 = zeros(Nunits,Nunits);
cur_cnt2 = zeros(Nunits,Nunits);
for tr = 1:length(un_stim_seeds)
    cur_tr_set = test_trials(trialSe(test_trials) == un_stim_seeds(tr));
    fprintf('Tr %d, %d rpts\n',tr,length(cur_tr_set));
    if length(cur_tr_set) >= 2
        cur_resp = trialrespDir(cur_tr_set);
        [II,JJ] = meshgrid(1:length(cur_tr_set));
        cur_Robs = (fullRobs_ms((beg_buff+1):end,cur_tr_set,:));
        cur_Robs2 = reshape(cur_Robs,[],length(cur_tr_set),1,Nunits);
        
        cur_Dmat = abs(squareform(pdist(cur_resp')));
        cur_Dmat(JJ == II) = nan;
        
        %same response
        curset = find(cur_Dmat == 0);
        temp = nanmean(bsxfun(@times,cur_Robs2(:,II(curset),:,:),cur_Robs(:,JJ(curset),:)),1);
        cur_XC = cur_XC + squeeze(nansum(temp,2));
        cur_cnt = cur_cnt + squeeze(sum(~isnan(temp),2));
%         cur_XC = cur_XC + squeeze(nansum(nanmean(bsxfun(@times,cur_Robs2(:,II(curset),:,:),cur_Robs(:,JJ(curset),:)),1),2));
%         cur_cnt = cur_cnt + squeeze(sum(~(isnan(cur_Robs(1,II(curset),:))) & ~(isnan(cur_Robs(1,JJ(curset),:))),2));
        
        %marginal
        curset = find(~isnan(cur_Dmat));
        temp = nanmean(bsxfun(@times,cur_Robs2(:,II(curset),:,:),cur_Robs(:,JJ(curset),:)),1);
        cur_XC2 = cur_XC2 + squeeze(nansum(temp,2));
        cur_cnt2 = cur_cnt2 + squeeze(sum(~isnan(temp),2));
%         cur_XC2 = cur_XC2 + squeeze(nansum(nanmean(bsxfun(@times,cur_Robs2(:,II(curset),:,:),cur_Robs(:,JJ(curset),:)),1),2));
%         cur_cnt2 = cur_cnt2 + squeeze(sum(~(isnan(cur_Robs(1,II(curset),:))) & ~(isnan(cur_Robs(1,JJ(curset),:))),2));
    end
end

% sameVar = bsxfun(@rdivide,cur_XC,cur_cnt);
% randVar = bsxfun(@rdivide,cur_XC2,cur_cnt2);
sameVar = cur_XC./cur_cnt;
randVar = cur_XC2./cur_cnt2;

% ovVar = fullRobs_resh(istrain,:));

choice_covmat = sameVar - randVar;
SU_choice_varfrac = diag(choice_covmat)./diag(randVar);
temp = reshape(fullRobs_ms(:,test_trials,:),[],Nunits);
normfac = sqrt(nanvar(temp)'*nanvar(temp));
choice_corrmat = choice_covmat./normfac;
stim_corrmat = randVar./normfac;

ov_cov = squeeze(nanmean(bsxfun(@times,fullRobs_resh(istrain,:),reshape(fullRobs_resh(istrain,:),sum(istrain),1,[]))));
ov_corr = ov_cov./normfac;

%% CHOICE CONDITIONAL COV (NO STIM CORRECTION) AT DT RESOLUTION
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
poss_choices = [-1 1];
up_resp = find(trialOB == 130 & trialrespDir == 1);
down_resp = find(trialOB == 130 & trialrespDir == -1);

up_avg = squeeze(nanmean(nanmean(fullRobs_ms(:,up_resp,:),2),1));
down_avg = squeeze(nanmean(nanmean(fullRobs_ms(:,down_resp,:),2),1));
choice_cond_avgs = [up_avg'; down_avg'];
dt_choice_cond_covmat = squeeze(nanmean(bsxfun(@times,choice_cond_avgs,reshape(choice_cond_avgs,2,1,Nunits))));

temp = reshape(fullRobs_ms(:,test_trials,:),[],Nunits);
normfac = sqrt(nanvar(temp)'*nanvar(temp));

dt_choice_cond_corrmat = dt_choice_cond_covmat./normfac;

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

%% REMOVE AVG WITHIN SE AND COMPUTE CHOICE COND COV MATS
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
un_stim_seeds = unique(trialSe(test_trials));

tot_spks_seavg = nan(size(tot_spks_per_trial_ms));
tot_spks_senorm = nan(size(tot_spks_per_trial_ms));
for ii = 1:length(un_stim_seeds)
    cur_tr_set = test_trials(trialSe(test_trials) == un_stim_seeds(ii));
    if length(cur_tr_set) >= 2 & length(unique(trialrespDir(cur_tr_set))) > 1
    tr_avgs = nanmean(tot_spks_per_trial_ms(cur_tr_set,:));
    tot_spks_seavg(cur_tr_set,:) = repmat(tr_avgs,length(cur_tr_set),1);
    tot_spks_senorm(cur_tr_set,:) = tot_spks_per_trial_ms(cur_tr_set,:) - repmat(tr_avgs,length(cur_tr_set),1);
    end
end

cur_resp = trialrespDir(test_trials);
resp_up = test_trials(cur_resp == 1);
resp_down = test_trials(cur_resp == -1);

avg_up_resp = nanmean(tot_spks_senorm(resp_up,:));
avg_down_resp = nanmean(tot_spks_senorm(resp_down,:));
choice_cond_avgs = [avg_up_resp; avg_down_resp];
secor_choice_cond_covmat = squeeze(nanmean(bsxfun(@times,choice_cond_avgs,reshape(choice_cond_avgs,2,1,Nunits))));

avg_up_resp = nanmean(tot_spks_seavg(resp_up,:));
avg_down_resp = nanmean(tot_spks_seavg(resp_down,:));
choice_cond_avgs = [avg_up_resp; avg_down_resp];
seavg_choice_cond_covmat = squeeze(nanmean(bsxfun(@times,choice_cond_avgs,reshape(choice_cond_avgs,2,1,Nunits))));

normfac = sqrt(nanvar(tot_spks_senorm(test_trials,:))'*nanvar(tot_spks_senorm(test_trials,:)));

secor_choice_cond_corrmat = secor_choice_cond_covmat./normfac;
seavg_choice_cond_corrmat = seavg_choice_cond_covmat./normfac;


%% USE WAVELET ANALYSIS TO COMPUTE PHASE-LOCKING SPECTRA FOR EACH UNIT
if strcmp(Expt_name,'G037')
    LFP_Fs = 1/AllExpt.Expt.Header.LFPsamplerate;
LFP_offset = AllExpt.Expt.Header.preperiod/1e4;
LFP_trial_taxis = (1:length(AllExpt.Expt.Header.LFPtimes))/LFP_Fs - LFP_offset;
else
LFP_trial_taxis = AllExpt.Expt.Header.LFPtimes*1e-4;
end

LFP_dsf = 2;
LFP_Fsd = LFP_Fs/LFP_dsf;
%anti-aliasing filter and high-pass filter
aa_hcf = LFP_Fsd/2*0.8;
% [b_aa,a_aa] = butter(4,aa_hcf/(LFP_Fs/2),'low');
aa_lcf = 0.5;
[b_aa,a_aa] = butter(2,[aa_lcf aa_hcf]/(LFP_Fs/2));

LFP_trial_taxis_ds = downsample(LFP_trial_taxis,LFP_dsf);

if Expt_name(1) == 'G'
nprobes = 96;
uprobes = 1:4:nprobes;
else
nprobes = 24;
uprobes = 1:2:nprobes;
end

%wavelet parameters
nwfreqs = 25;
min_freq = 1.; max_freq = 80;
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
R_trial_taxis = (1+beg_buffer:(trial_dur - end_buffer))*dt;
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

tbt_LFPs = reshape(trial_LFPs,[TLEN Ntrials length(uprobes)]);
% tbt_LFP_cwt = reshape(trial_LFP_cwt,[TLEN Ntrials length(wfreqs) length(uprobes)]);
tbt_LFP_real = reshape(trial_LFP_real,[TLEN Ntrials length(wfreqs) length(uprobes)]);
tbt_LFP_imag = reshape(trial_LFP_imag,[TLEN Ntrials length(wfreqs) length(uprobes)]);


%% compute power-based spike-lfp coupling directly (on trial-avg basis) 
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
unit_tot_avg_amp = nan(Nunits,length(wfreqs),length(uprobes));
for cc = 1:Nunits
   
    cur_Robs = reshape(Robs_tot(:,cc),[],1);
     
    trig_avg_amp = squeeze(nanmean(bsxfun(@times,trial_pow,cur_Robs)));

    unit_tot_avg_amp(cc,:,:) = trig_avg_amp;
end

%% COMPUTE SPIKE-LFP COUPLING WITH AND WITHOUT CONTROLLING FOR STIM_DEPENDENCE
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
% test_trials = find(trialrespDir ~= 0);

beg_buff = round(0.05/dt); %number of bins from beginning of trial to exclude
end_buff = 0;

u_lfp_times = find(R_trial_taxis > beg_buff*dt & R_trial_taxis <= (trial_dur-end_buff)*dt);
LFP_real = squeeze(tbt_LFP_real(u_lfp_times,test_trials,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,test_trials,:,:));
LFP_amp = sqrt(LFP_imag.^2 + LFP_real.^2);

LFP_real = bsxfun(@minus,LFP_real,nanmean(LFP_real,2));
LFP_imag = bsxfun(@minus,LFP_imag,nanmean(LFP_imag,2));
LFP_amp = bsxfun(@minus,LFP_amp,nanmean(LFP_amp,2));

LFP_amp = reshape(LFP_amp,[],length(uprobes)*length(wfreqs));
% LFP_real = reshape(LFP_real,[],length(uprobes)*length(wfreqs));
% LFP_imag = reshape(LFP_imag,[],length(uprobes)*length(wfreqs));

% LFP_real = bsxfun(@rdivide,LFP_real,nanstd(LFP_amp));
% LFP_imag = bsxfun(@rdivide,LFP_imag,nanstd(LFP_amp));

Robs = fullRobs((1+beg_buff):(NT-end_buff),test_trials,:);
Robs = reshape(Robs,[],Nunits);
Robs = bsxfun(@minus,Robs,nanmean(Robs));
Robs = reshape(Robs,[],length(test_trials),Nunits);

unit_trig_avg_amp = nan(Nunits,length(wfreqs)*length(uprobes));
for cc = 1:Nunits
   cc
    cur_Robs = reshape(Robs(:,:,cc),[],1);
    
%     uset = find(~isnan(cur_Robs));
%     cur_spkbins = convert_to_spikebins(cur_Robs(uset));
%     trig_avg_real = nanmean(LFP_real(uset(cur_spkbins),:));
%     trig_avg_imag = nanmean(LFP_imag(uset(cur_spkbins),:));
%     trig_avg_amp = nanmean(LFP_amp(uset(cur_spkbins),:));
 
%     trig_avg_real = nanmean(bsxfun(@times,LFP_real,cur_Robs));
%     trig_avg_imag = nanmean(bsxfun(@times,LFP_imag,cur_Robs));
    trig_avg_amp = nanmean(bsxfun(@times,LFP_amp,cur_Robs));

%     unit_trig_avg_real(cc,:) = trig_avg_real;
%     unit_trig_avg_imag(cc,:) = trig_avg_imag;
    unit_trig_avg_amp(cc,:) = trig_avg_amp;
end
%%

% LFP_real = reshape(LFP_real,[],length(test_trials),length(wfreqs),length(uprobes));
% LFP_imag = reshape(LFP_imag,[],length(test_trials),length(wfreqs),length(uprobes));
LFP_amp = reshape(LFP_amp,[],length(test_trials),length(wfreqs),length(uprobes));

un_stim_seeds = unique(trialSe(test_trials));
all_LFP_Aavgs = zeros(Nunits,length(wfreqs),length(uprobes));
% all_LFP_ravgs = zeros(Nunits,length(wfreqs),length(uprobes));
% all_LFP_iavgs = zeros(Nunits,length(wfreqs),length(uprobes));
all_LFP_cnts = zeros(Nunits,length(wfreqs),length(uprobes));
for tr = 1:length(un_stim_seeds)
    tr
    cur_tr_set = find(trialSe(test_trials) == un_stim_seeds(tr));
    
    if length(cur_tr_set) >= 2
        [II,JJ] = meshgrid(1:length(cur_tr_set));
        uset = find(II ~= JJ);
        
%         cur_LFP_real = LFP_real(:,cur_tr_set,:,:);
%         cur_LFP_imag = LFP_imag(:,cur_tr_set,:,:);
        cur_LFP_amp = LFP_amp(:,cur_tr_set,:,:);
        
        for cc = 1:Nunits
            cur_Robs = Robs(:,cur_tr_set,cc);
            
%             cur_real = squeeze(nanmean(bsxfun(@times,cur_LFP_real(:,II(uset),:,:),cur_Robs(:,JJ(uset)))));
%             cur_imag = squeeze(nanmean(bsxfun(@times,cur_LFP_imag(:,II(uset),:,:),cur_Robs(:,JJ(uset)))));
            cur_amp = squeeze(nanmean(bsxfun(@times,cur_LFP_amp(:,II(uset),:,:),cur_Robs(:,JJ(uset)))));
%             all_LFP_ravgs(cc,:,:) = all_LFP_ravgs(cc,:,:) + nansum(cur_real);
%             all_LFP_iavgs(cc,:,:) = all_LFP_iavgs(cc,:,:) + nansum(cur_imag);
            all_LFP_Aavgs(cc,:,:) = all_LFP_Aavgs(cc,:,:) + nansum(cur_amp);
            all_LFP_cnts(cc,:,:) = all_LFP_cnts(cc,:,:) + sum(~isnan(cur_amp));
        end
    end
end

% all_LFP_ravgs = all_LFP_ravgs./all_LFP_cnts;
% all_LFP_iavgs = all_LFP_iavgs./all_LFP_cnts;
all_LFP_Aavgs = all_LFP_Aavgs./all_LFP_cnts;

% controlled_ravgs = reshape(unit_trig_avg_real,Nunits,length(wfreqs),length(uprobes)) - all_LFP_ravgs;
% controlled_iavgs = reshape(unit_trig_avg_imag,Nunits,length(wfreqs),length(uprobes)) - all_LFP_iavgs;
controlled_Aavgs = reshape(unit_trig_avg_amp,Nunits,length(wfreqs),length(uprobes)) - all_LFP_Aavgs;
% uncontrolled_ravgs = reshape(unit_trig_avg_real,Nunits,length(wfreqs),length(uprobes));
% uncontrolled_iavgs = reshape(unit_trig_avg_imag,Nunits,length(wfreqs),length(uprobes));
uncontrolled_Aavgs = reshape(unit_trig_avg_amp,Nunits,length(wfreqs),length(uprobes));

%%
% close all
% for cc = 1:Nunits
%    subplot(2,1,1)
%    pcolor(wfreqs,1:length(uprobes),squeeze(uncontrolled_Aavgs(cc,:,:))');shading flat
%    caxis([-5 5]*1e-3);
%    subplot(2,1,2)
%    pcolor(wfreqs,1:length(uprobes),squeeze(controlled_Aavgs(cc,:,:))');shading flat
%    caxis([-5 5]*1e-3);
%    pause
%    clf
%     
% end

f1 = figure();
avg_controlled_Aavgs = squeeze(nanmean(controlled_Aavgs));
avg_uncontrolled_Aavgs = squeeze(nanmean(uncontrolled_Aavgs));
subplot(2,1,1)
pcolor(wfreqs,1:length(uprobes),avg_controlled_Aavgs');shading flat
caxis([-5 5]*1e-3);colorbar
set(gca,'xscale','log');
set(gca,'ydir','reverse');
xlim([2 80]);
title('Stim controlled');
xlabel('Frequency (Hz)');
ylabel('Channel');

subplot(2,1,2)
pcolor(wfreqs,1:length(uprobes),avg_uncontrolled_Aavgs');shading flat
caxis([-5 5]*1e-3);colorbar
set(gca,'xscale','log');
set(gca,'ydir','reverse');
xlim([2 80]);
title('Uncontrolled');
xlabel('Frequency (Hz)');
ylabel('Channel');

% fig_width = 4; rel_height = 1.6;
% figufy(f1);
% fname = [fig_dir sprintf('Avg_Spk_pow_coupling_%s.pdf',Expt_name)];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% LFP power CP predictions
test_trials = find(trialOB == 130 & trialrespDir ~= 0);
train_trials = find(trialOB < 130 & trialrespDir ~= 0);

up_trials = find(trialrespDir(test_trials) == 1);
down_trials = find(trialrespDir(test_trials) == -1);

beg_buff = round(0.05/dt); %number of bins from beginning of trial to exclude
end_buff = 0;
% beg_buff = 0; %number of bins from beginning of trial to exclude
% end_buff = 0;

%range of times over which to compute lfp power
lfp_bt = 0.2;
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
%     for ww = 1:length(wfreqs)
%         cur_X = squeeze(trial_pow_train(:,ww,:));
%         cur_Y = trial_Robs_train(:,cc);
%         uindx = find(~isnan(cur_Y) & ~any(isnan(cur_X),2));
%         
%         init_mod = NMMinitialize_model(stim_params,1,{'lin'},reg_params);
%         init_mod = NMMfit_filters(init_mod,cur_Y,cur_X,[],uindx,silent);
% %         B = glmfit(cur_X,trial_Robs_train(:,cc),'poisson');
%         
%         mod_params(cc,ww,:) = init_mod.mods(1).filtK;
% 
%         cur_X = squeeze(trial_pow_test(:,ww,:));
%         cur_Y = trial_Robs_test(:,cc);
%         
%         [~,~,Yhat] = NMMeval_model(init_mod,cur_Y,cur_X);
% %         Yhat = glmval(B,cur_X,'log');
%         Yhat(isnan(cur_Y)) = nan;
%         
%         cur_resp_ax = prctile(Yhat,[0:5:100]);
%         up_hist = histc(Yhat(up_trials),cur_resp_ax);
%         down_hist = histc(Yhat(down_trials),cur_resp_ax);
%         true_pos = cumsum(up_hist)/sum(~isnan(Yhat(up_trials)));
%         false_pos = cumsum(down_hist)/sum(~isnan(Yhat(down_trials)));
%         all_LFPpow_pred_CP(cc,ww) = trapz(false_pos,true_pos);
% 
%     end
    
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
ufreqs = find(wfreqs >= 2 & wfreqs <= 10); %set of frequencies to use in model
uchs = 1:1:length(uprobes);
reg_params = NMMcreate_reg_params('lambda_L2',1,'lambda_d2XT',1); %very slight reg
stim_params = NMMcreate_stim_params([length(ufreqs) length(uchs)]);
silent = 1;

for cc = 1:Nunits
        %training data (OB ~= 130)
        cur_X = squeeze(trial_pow_train(:,ufreqs,:));
        cur_Y = trial_Robs_train(:,cc);
        uindx = find(~isnan(cur_Y) & ~any(isnan(reshape(cur_X,length(train_trials),[])),2));
    
        %fit model to training data 
        init_mod = NMMinitialize_model(stim_params,1,{'lin'},reg_params);
        init_mod = NMMfit_filters(init_mod,cur_Y,cur_X,[],uindx,silent);
        
        %test data (OB==130)
        cur_X = squeeze(trial_pow_test(:,ufreqs,:));
        cur_X = reshape(cur_X,length(test_trials),[]);
        cur_Y = trial_Robs_test(:,cc);
        
        %get predicted rate on each trial
        [~,~,Yhat] = NMMeval_model(init_mod,cur_Y,cur_X);
        Yhat(isnan(cur_Y)) = nan;
        
        Y_resid = cur_Y - Yhat;
                
        %CP using predicted trial spike count based on power model
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

f1 = figure();
plot(all_LFPpow_allw_pred_CP,new_choice_prob,'o');
line([0 1],[0 1],'color','k');
xlim([0.2 0.8]); ylim([0.2 0.8]);
xlabel('Predicted CP');
ylabel('Measured CP');


% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir sprintf('Prec_Meas_CP_compare_%s.pdf',Expt_name)];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% [a,b] = corr(all_LFPpow_allw_pred_CP(ucells)',new_choice_prob(ucells)','type','spearman')
%% COMPUTE CHOICE_CONDITIONAL TRIALAVG POW SPECTRA

test_trials = find(trialOB == 130 & trialrespDir ~= 0);
up_trials = trialrespDir(test_trials) == 1;
down_trials = trialrespDir(test_trials) == -1;

bt = 0;
et = 0;

u_lfp_times = find(R_trial_taxis > bt & R_trial_taxis <= (trial_dur)*dt-et);
LFP_real = squeeze(tbt_LFP_real(u_lfp_times,test_trials,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,test_trials,:,:));
LFP_amp = sqrt(LFP_imag.^2 + LFP_real.^2);

% LFP_vals = abs(tbt_LFPs(u_lfp_times,test_trials,:));

LFP_amp_upavg = squeeze(nanmean(LFP_amp(:,up_trials,:,:),2));
LFP_amp_downavg = squeeze(nanmean(LFP_amp(:,down_trials,:,:),2));
LFP_choice_powdiff = squeeze(nanmean(LFP_amp_upavg,3) - nanmean(LFP_amp_downavg,3));

% LFP_val_upavg = squeeze(nanmean(LFP_vals(:,up_trials,:),2));
% LFP_val_downavg = squeeze(nanmean(LFP_vals(:,down_trials,:),2));
% LFP_choice_valdiff = squeeze(nanmean(LFP_val_upavg,3) - nanmean(LFP_val_downavg,3));

%%
test_trials = find(abs(trialOB) == 60 & trialrespDir ~= 0);
up_trials = trialrwDir(test_trials) == 1;
down_trials = trialrwDir(test_trials) == -1;

LFP_real = squeeze(tbt_LFP_real(u_lfp_times,test_trials,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,test_trials,:,:));
LFP_amp = sqrt(LFP_imag.^2 + LFP_real.^2);

LFP_samp_upavg = squeeze(nanmean(LFP_amp(:,up_trials,:,:),2));
LFP_samp_downavg = squeeze(nanmean(LFP_amp(:,down_trials,:,:),2));
LFP_sig_powdiff = squeeze(nanmean(LFP_samp_upavg,3) - nanmean(LFP_samp_downavg,3));
% sm_sig = 5;
% for ii = 1:length(wfreqs)
%     LFP_sig_powdiff(:,ii) = jmm_smooth_1d_cor(LFP_sig_powdiff(:,ii),sm_sig);
% end

% ca = [-0.4 0.4];
% close all
% for cc = 1:length(uprobes)
%     pcolor(R_trial_taxis(u_lfp_times),wfreqs,squeeze(LFP_amp_upavg(:,:,cc)-LFP_amp_downavg(:,:,cc))');shading flat
%     caxis(ca);
%     pause
%     clf
% end
% 

f1 = figure();
subplot(2,1,1);
pcolor(R_trial_taxis(u_lfp_times),wfreqs,LFP_choice_powdiff');shading flat
caxis([-0.3 0.3]);
set(gca,'yscale','log');
ylim([2 80]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
subplot(2,1,2);
pcolor(R_trial_taxis(u_lfp_times),wfreqs,LFP_sig_powdiff');shading flat
caxis([-0.3 0.3]);
set(gca,'yscale','log');
ylim([2 80]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');


beg_buff = round(0.2/dt);
end_buff = round(0.2/dt);
choice_specdiff = squeeze(nanmean(LFP_amp_upavg(1+beg_buff:(end-end_buff),:,:)) - nanmean(LFP_amp_downavg(1+beg_buff:(end-end_buff),:,:)));
sig_specdiff = squeeze(nanmean(LFP_samp_upavg(1+beg_buff:(end-end_buff),:,:)) - nanmean(LFP_samp_downavg(1+beg_buff:(end-end_buff),:,:)));
f2 = figure();
subplot(2,1,1);
pcolor(wfreqs,1:24,choice_specdiff');shading flat
caxis([-0.2 0.2]);
set(gca,'xscale','log');
set(gca,'ydir','reverse');
xlim([2 80]);
xlabel('Time (s)');
ylabel('Channel');
subplot(2,1,2);
pcolor(wfreqs,1:24,sig_specdiff');shading flat
caxis([-0.2 0.2]);
set(gca,'xscale','log');
xlim([2 80]);
xlabel('Time (s)');
ylabel('Channel');
set(gca,'ydir','reverse');

% fig_width = 4; rel_height =1.6;
% figufy(f1);
% fname = [fig_dir sprintf('Specgram_powdiff_%s.pdf',Expt_name)];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir sprintf('Specdepth_powdiff_%s.pdf',Expt_name)];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% COMPUTE CP OF TRIAL_AVG POW SPECTRA
lfp_bt = 0.2; %buffers at beginning and end of trial to exclude for power calculation
lfp_et = 0.2;
u_lfp_times = find(R_trial_taxis > 0 & R_trial_taxis <= (trial_dur)*dt);
u_lfp_times = u_lfp_times(R_trial_taxis(u_lfp_times) >= lfp_bt & R_trial_taxis(u_lfp_times) <= (trial_dur-lfp_et));

test_trials = find(trialOB == 130 & trialrespDir ~= 0);

LFP_real = squeeze(tbt_LFP_real(u_lfp_times,test_trials,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,test_trials,:,:));
LFP_amp = sqrt(LFP_imag.^2 + LFP_real.^2);

mean_LFP_amp = squeeze(nanmean(LFP_amp,4));

trial_pow = squeeze(nanmean(LFP_amp));

%compute CP from LFP power at each frequency and channel
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
%         [~,~,~,all_LFPpow_CP(ww,cc)]=perfcurve(trialrespDir(test_trials), trial_pow(:,ww,cc),1);
end
end


%% SIG PROB OF LFPs
lfp_bt = 0.2; %buffers at beginning and end of trial to exclude for power calculation
lfp_et = 0.2;
u_lfp_times = find(R_trial_taxis > 0 & R_trial_taxis <= (trial_dur)*dt);
u_lfp_times = u_lfp_times(R_trial_taxis(u_lfp_times) >= lfp_bt & R_trial_taxis(u_lfp_times) <= (trial_dur-lfp_et));

test_trials = find(abs(trialOB) == 60 & trialrespDir ~= 0);

LFP_real = squeeze(tbt_LFP_real(u_lfp_times,test_trials,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,test_trials,:,:));
LFP_amp = sqrt(LFP_imag.^2 + LFP_real.^2);

mean_LFP_amp = squeeze(nanmean(LFP_amp,4));

trial_pow = squeeze(nanmean(LFP_amp));

%compute CP from LFP power at each frequency and channel
resp_up = find(trialrwDir(test_trials) == 1);
resp_down = find(trialrwDir(test_trials) == -1);
all_LFPpow_SP = nan(length(wfreqs),length(uprobes));
for cc = 1:length(uprobes)
    for ww = 1:length(wfreqs)
        cur_resp_ax = prctile(trial_pow(:,ww,cc),[0:5:100]);
        up_hist = histc(trial_pow(resp_up,ww,cc),cur_resp_ax);
        down_hist = histc(trial_pow(resp_down,ww,cc),cur_resp_ax);
        true_pos = cumsum(up_hist)/sum(~isnan(trial_pow(resp_up,ww,cc)));
        false_pos = cumsum(down_hist)/sum(~isnan(trial_pow(resp_down,ww,cc)));
        all_LFPpow_SP(ww,cc) = trapz(false_pos,true_pos);
    end
end


f1 = figure();
subplot(2,1,1);
pcolor(wfreqs,1:length(uprobes),1-all_LFPpow_CP');shading flat
caxis([0.35 0.65]);
title('Choice prob');
set(gca,'xscale','log','ydir','reverse');
xlabel('Frequency (Hz)');
ylabel('Depth');
xlim([2 80]);

subplot(2,1,2);
pcolor(wfreqs,1:length(uprobes),1-all_LFPpow_SP');shading flat
caxis([0.35 0.65]);
title('Sig prob');
set(gca,'xscale','log','ydir','reverse');
xlim([2 80]);
xlabel('Frequency (Hz)');
ylabel('Depth');

% figufy(f1);
% fname = [fig_dir sprintf('LFPpow_CP_SP_%s.pdf',Expt_name)];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% COMPUTE LFP VARIANCE ATTRIBUTABLE TO CHOICE
poss_trials = find(trialOB == 130 & trialrespDir ~= 0);
un_stim_seeds = unique(trialSe(poss_trials));

beg_buff = 0; %number of bins from beginning of trial to exclude
end_buff = 0;

u_lfp_times = find(R_trial_taxis > 0 & R_trial_taxis <= (trial_dur)*dt);
LFP_real = squeeze(tbt_LFP_real(u_lfp_times,poss_trials,:,:));
LFP_imag = squeeze(tbt_LFP_imag(u_lfp_times,poss_trials,:,:));
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
    for ss = 1:length(un_stim_seeds)
        cur_trials = find(trialSe(poss_trials) == un_stim_seeds(ss)); %trials (of test trials) with current stim seed
        if length(cur_trials) >= 2 %need at least 2
            
            %find which trial pairs had same vs opposite choice dirs
            cur_resp = trialrespDir(poss_trials(cur_trials));
            [II,JJ] = meshgrid(1:length(cur_trials));
            cur_Y = squeeze(LFP_comp(:,cur_trials,:,cc));
            cur_AY = squeeze(LFP_Acomp(:,cur_trials,:,cc));
            
            %difference in choice dir (either 1 or 0)
            cur_Dmat = abs(squareform(pdist(cur_resp')));
            cur_Dmat(JJ == II) = nan; %unequal trial pairs
            
            %trial pairs with same choice
            uset = find(cur_Dmat == 0);
            temp = cur_Y(:,II(uset),:).*conj(cur_Y(:,JJ(uset),:));
            tempm = squeeze(nanmean(temp,2));
            all_cwt_same(~isnan(tempm)) = all_cwt_same(~isnan(tempm))  + tempm(~isnan(tempm)) ;
            
            %compute avg product of amplitudes for normalization
            temp = abs(cur_Y(:,II(uset),:)).*abs(cur_Y(:,JJ(uset),:));
            tempm = squeeze(nanmean(temp,2));
            all_Acwt_same(~isnan(tempm)) = all_Acwt_same(~isnan(tempm))  + tempm(~isnan(tempm)) ;
            
            temp = cur_AY(:,II(uset),:).*cur_AY(:,JJ(uset),:);
            tempm = squeeze(nanmean(temp,2));
            all_same_cnt = all_same_cnt + squeeze(sum(~isnan(temp),2));
            
            %this is the marginal (ignoring choice)
            uset = find(~isnan(cur_Dmat));
            temp = cur_Y(:,II(uset),:).*conj(cur_Y(:,JJ(uset),:));
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
    all_cwt_choice = all_cwt_same - all_cwt_diff;
    all_cwt_choice_norm(cc,:,:) = all_cwt_choice./(all_Acwt_diff);
    all_cwt_seq_norm(cc,:,:) = all_cwt_diff./all_Acwt_diff;
end
%%
close all

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

% fig_width = 5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir sprintf('LFPchoice_coherence_%s.pdf',Expt_name)];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);


%%
