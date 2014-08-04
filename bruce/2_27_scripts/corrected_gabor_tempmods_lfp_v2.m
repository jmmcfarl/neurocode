clear all
% close all
addpath(genpath('~/James_scripts'));

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-5 5;-5 5];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;

Pix2Deg = 0.018837;

% down-sampling fraction for image
dsfrac = 2;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

cd /Users/James/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

SDIM = length(xpatch_inds);
cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_image_patches_corrected_raw

NT = length(used_fixs);
X_resh = reshape(X,NT,SDIM^2);
dt = .005;

cd /Users/James/James_scripts/bruce/modelfits
load pref_oris
cd /Users/James/James_scripts/bruce/modelfits/
% load fixation_gabor_models
load gabor_tempmodfits_en_lin

%%
Fs = 1000;
dsf = 4;
Fsd = Fs/dsf;
niqf = Fs/2;
% [b,a] = butter(4,[45 60]/niqf);
[b,a] = butter(2,120/niqf,'low');
% scales = logspace(log10(2.5),log10(35),20);
scales = logspace(log10(2.5),log10(50),30);
% scales = [5.15 5.7 6.3];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
[b_alpha,a_alpha] = butter(2,[5 15]/niqf);
[b_gamma,a_gamma] = butter(2,[25 60]/niqf);
%%
all_lfp_pow = [];
within_fix_avgs = [];
fix_nums = [];
all_used_inds = [];
all_ampgrams = [];
all_alpha_freq = [];
all_gamma_freq = [];
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    
    %%
    cd /Users/James/Data/bruce/2_27_12
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
    %get start times of each LFP trial
    if blockid == 4
        LFP.Trials = LFP.Trials(1:5);
    end
    n_lfp_trials = length(LFP.Trials);
    lfp_trial_start = nan(n_lfp_trials,1);
    lfp_trial_stop = nan(n_lfp_trials,1);
    for i = 1:n_lfp_trials
        lfp_trial_start(i) = LFP.Trials(i).Start/1e4; %originally in tenths of ms
        lfp_trial_stop(i) = LFP.Trials(i).End/1e4;
        lfp_dur(i) = size(LFP.Trials(i).LFP,1)/1000;
    end
    
    lfp_time = [];
    lfp_samps = [];
    for i = 1:n_lfp_trials
        cur_len(i) = size(LFP.Trials(i).LFP,1);
        if i < n_lfp_trials
            next_start = lfp_trial_start(i+1);
            start_len = length(lfp_trial_start(i):1/Fs:next_start);
        else
            next_start = Inf;
            start_len = Inf;
        end
        cur_end(i) = min(cur_len(i),start_len);
        cur_t = lfp_trial_start(i):1/Fs:(lfp_trial_start(i)+cur_end(i)/Fs);
        cur_t(cur_end(i)+1:end) = [];
        lfp_time = [lfp_time cur_t];
        lfp_samps = [lfp_samps; LFP.Trials(i).LFP(1:cur_end(i),:)];
    end
    
    lfp_samps = filtfilt(b,a,lfp_samps);
    lfp_sampsd = downsample(lfp_samps,dsf);
    lfp_timed = downsample(lfp_time,dsf);
    lfp_alpha = filtfilt(b_alpha,a_alpha,lfp_samps);
    lfp_alpha = zscore(downsample(lfp_alpha,dsf));
    lfp_gamma = filtfilt(b_gamma,a_gamma,lfp_samps);
    lfp_gamma = zscore(downsample(lfp_gamma,dsf));
    
    cur_all_model = find(all_model_blockids==blockid);
    interp_ampgram = zeros(24,length(cur_all_model),length(wfreqs));
    interp_gammafreq = zeros(24,length(cur_all_model));
    interp_alphafreq = zeros(24,length(cur_all_model));
    for cc = 1:24
        %     cc = 17;
        fprintf('Channel %d of %d\n',cc,24);
        temp = cwt(lfp_sampsd(:,cc),scales,'cmor1-1');
        tampgram = abs(temp);
        %         tampgram = bsxfun(@minus,tampgram,mean(tampgram,2));
        %         tampgram = bsxfun(@rdivide,tampgram,std(tampgram,[],2));
        tampgram = tampgram';
        interp_ampgram(cc,:,:) = interp1(lfp_timed,tampgram,all_model_time_axis(cur_all_model));
        
        gamma_freq = dnf_hilbert_instfreq(lfp_gamma(:,cc),Fsd,5);
        alpha_freq = dnf_hilbert_instfreq(lfp_alpha(:,cc),Fsd,5);
        interp_gammafreq(cc,:) = interp1(lfp_timed,gamma_freq.f_smooth,all_model_time_axis(cur_all_model));
        interp_alphafreq(cc,:) = interp1(lfp_timed,alpha_freq.f_smooth,all_model_time_axis(cur_all_model));
    end
    clear temp tampgram
    
    %%
    used_tpoints = find(~isnan(interp_ampgram(1,:,1)));
    all_ampgrams = cat(2,all_ampgrams,interp_ampgram(:,used_tpoints,:));
    all_gamma_freq = cat(2,all_gamma_freq,interp_gammafreq(:,used_tpoints,:));
    all_alpha_freq = cat(2,all_alpha_freq,interp_alphafreq(:,used_tpoints,:));
   all_used_inds = [all_used_inds; cur_all_model(used_tpoints)'];
    
end
% all_ampgrams = sqrt(all_ampgrams);
all_ampgrams = bsxfun(@minus,all_ampgrams,mean(all_ampgrams,2));
all_ampgrams = bsxfun(@rdivide,all_ampgrams,std(all_ampgrams,[],2));

%%
cur_used_fixs = unique(all_model_fixids(all_used_inds));
n_fixs = length(cur_used_fixs);
fix_avg_amps = zeros(24,n_fixs,length(wfreqs));
late_avg_amps = zeros(24,n_fixs,length(wfreqs));
early_avg_amps = zeros(24,n_fixs,length(wfreqs));
fix_avg_gfreq = zeros(24,n_fixs);
late_avg_gfreq = zeros(24,n_fixs);
early_avg_gfreq = zeros(24,n_fixs);
fix_avg_afreq = zeros(24,n_fixs);
late_avg_afreq = zeros(24,n_fixs);
early_avg_afreq = zeros(24,n_fixs);
for i = 1:n_fixs
    cur_set = find(all_model_fixids(all_used_inds) == cur_used_fixs(i));
    fix_avg_amps(:,i,:) = squeeze(mean(all_ampgrams(:,cur_set,:),2));
    fix_avg_gfreq(:,i) = squeeze(median(all_gamma_freq(:,cur_set),2));
    fix_avg_afreq(:,i) = squeeze(median(all_alpha_freq(:,cur_set),2));
    cur_set2 = cur_set(time_since_fix(all_used_inds(cur_set)) >= 0.15 & time_since_fix(all_used_inds(cur_set)) < 0.3);
    late_avg_amps(:,i,:) = squeeze(mean(all_ampgrams(:,cur_set2,:),2));
    late_avg_gfreq(:,i) = squeeze(median(all_gamma_freq(:,cur_set2),2));
    late_avg_afreq(:,i) = squeeze(median(all_alpha_freq(:,cur_set2),2));
    cur_set2 = cur_set(time_since_fix(all_used_inds(cur_set)) < 0.15);
    early_avg_amps(:,i,:) = squeeze(mean(all_ampgrams(:,cur_set2,:),2));
    early_avg_gfreq(:,i) = squeeze(median(all_gamma_freq(:,cur_set2),2));
    early_avg_afreq(:,i) = squeeze(median(all_alpha_freq(:,cur_set2),2));
end
fix_avg_amps = bsxfun(@minus,fix_avg_amps,mean(fix_avg_amps,2));
fix_avg_amps = bsxfun(@rdivide,fix_avg_amps,std(fix_avg_amps,[],2));

late_avg_amps = bsxfun(@minus,late_avg_amps,nanmean(late_avg_amps,2));
late_avg_amps = bsxfun(@rdivide,late_avg_amps,nanstd(late_avg_amps,[],2));
late_avg_amps(isnan(late_avg_amps)) = 0;

early_avg_amps = bsxfun(@minus,early_avg_amps,nanmean(early_avg_amps,2));
early_avg_amps = bsxfun(@rdivide,early_avg_amps,nanstd(early_avg_amps,[],2));
early_avg_amps(isnan(early_avg_amps)) = 0;

%%
X_resh = reshape(X,NT,SDIM^2);
orientations = linspace(0,pi-pi/12,12);
init_params = [5 -5 0 15 6 1];
for i = 1:12
    i
    init_params(3) = orientations(i);
    cur_mask1 = get_pgabor_mask(init_params,0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(init_params,pi/2,[SDIM SDIM]);
    
    mask1_out = X_resh(cur_used_fixs,:)*cur_mask1(:);
    mask2_out = X_resh(cur_used_fixs,:)*cur_mask2(:);
    
    en_out(i,:) = sqrt(mask1_out.^2+mask2_out.^2);
end
tot_cont = mean(en_out);
tot_cont = zscore(tot_cont);
en_out = bsxfun(@minus,en_out,mean(en_out));
en_out = bsxfun(@rdivide,en_out,std(en_out,[],2));

rsq = zeros(24,length(wfreqs),12);
rsq_late = zeros(24,length(wfreqs),12);
rsq_early = zeros(24,length(wfreqs),12);
R_vals = zeros(24,length(wfreqs),12,3);
R_vals_early = zeros(24,length(wfreqs),12,3);
R_vals_late = zeros(24,length(wfreqs),12,3);

for cc = 1:24
    fprintf('Channel %d of %d\n',cc,24);
    for ww = 1:length(wfreqs)
        %     for ww = 9
        for i = 1:12
            temp = [en_out(i,:)' tot_cont' ones(length(cur_used_fixs),1)];
            [B,BINT,R,RINT,STATS] = regress(squeeze(fix_avg_amps(cc,:,ww))',temp);
            rsq(cc,ww,i) = STATS(1);
            R_vals(cc,ww,i,:) = B;
            [B,BINT,R,RINT,STATS] = regress(squeeze(late_avg_amps(cc,:,ww))',temp);
            rsq_late(cc,ww,i) = STATS(1);
            R_vals_late(cc,ww,i,:) = B;
            [B,BINT,R,RINT,STATS] = regress(squeeze(early_avg_amps(cc,:,ww))',temp);
            rsq_early(cc,ww,i) = STATS(1);
            R_vals_early(cc,ww,i,:) = B;
            
        end
    end
end

%%
% rsq_gf = zeros(24,12);
% rsq_late_gf = zeros(24,12);
% rsq_early_gf = zeros(24,12);
% R_vals_gf = zeros(24,12,3);
% R_vals_early_gf = zeros(24,12,3);
% R_vals_late_gf = zeros(24,12,3);
% 
% rsq_af = zeros(24,12);
% rsq_late_af = zeros(24,12);
% rsq_early_af = zeros(24,12);
% R_vals_af = zeros(24,12,3);
% R_vals_early_af = zeros(24,12,3);
% R_vals_late_af = zeros(24,12,3);
% for cc = 1:24
%     for i = 1:12
%     temp = [en_out(i,:)' tot_cont' ones(length(cur_used_fixs),1)];
%     [B,BINT,R,RINT,STATS] = regress(squeeze(fix_avg_gfreq(cc,:))',temp);
%     rsq_gf(cc,i) = STATS(1);
%     R_vals_gf(cc,i,:) = B;
%     [B,BINT,R,RINT,STATS] = regress(squeeze(late_avg_gfreq(cc,:))',temp);
%     rsq_late_gf(cc,i) = STATS(1);
%     R_vals_late_gf(cc,i,:) = B;
%     [B,BINT,R,RINT,STATS] = regress(squeeze(early_avg_gfreq(cc,:))',temp);
%     rsq_early_gf(cc,i) = STATS(1);
%     R_vals_early_gf(cc,i,:) = B;
%     temp = [en_out(i,:)' tot_cont' ones(length(cur_used_fixs),1)];
%     [B,BINT,R,RINT,STATS] = regress(squeeze(fix_avg_afreq(cc,:))',temp);
%     rsq_af(cc,i) = STATS(1);
%     R_vals_af(cc,i,:) = B;
%     [B,BINT,R,RINT,STATS] = regress(squeeze(late_avg_afreq(cc,:))',temp);
%     rsq_late_af(cc,i) = STATS(1);
%     R_vals_late_af(cc,i,:) = B;
%     [B,BINT,R,RINT,STATS] = regress(squeeze(early_avg_afreq(cc,:))',temp);
%     rsq_early_af(cc,i) = STATS(1);
%     R_vals_early_af(cc,i,:) = B;
%     end
% end
%%
cd ~/James_scripts/bruce/modelfits/
load ./lfp_model_fits

% [max_rsq,maxloc] = max(rsq,[],3);
% [max_rsq_late,maxloc_late] = max(rsq_late,[],3);
% [max_rsq_early,maxloc_early] = max(rsq_early,[],3);
cur_oris = linspace(0,pi-pi/8,8);
en_out = nan(8,n_fixs);
X_used = X_resh(cur_used_fixs,:);

used_freq = 8;
hold_const = [0 0 0 0 0 1 0 0];
LB = [-5 -5 0 6 2 0.5 -Inf -Inf];
UB = [10 10 pi SDIM(1)/2 SDIM(1)/6 2 Inf Inf Inf];
% lfp_gabor_params = zeros(24,length(wfreqs),8);
% LL_fin = zeros(24,length(wfreqs));
for cc = 1:24
    fprintf('CC %d of %d\n',cc,24);
    for ww = 1:length(wfreqs)
        %     for ww = 9
        fprintf('freq %d of %d\n',ww,length(wfreqs));
%         init_params(3) = orientations(maxloc(cc,ww));
%         cur_init_params = [init_params zeros(1,2)];
%         [lfp_gabor_params(cc,ww,:),LL_fin(cc,ww)] = fit_gabor_params_lfp(cur_init_params,X_resh(cur_used_fixs,:),...
%             squeeze(fix_avg_amps(cc,:,ww)),[SDIM SDIM],hold_const,LB,UB);
        
        init_params = squeeze(lfp_gabor_params(cc,ww,1:6));
        for i = 1:length(cur_oris)
            init_params(3) = cur_oris(i);
            cur_mask1 = get_pgabor_mask(init_params,0,[SDIM SDIM]);
            cur_mask2 = get_pgabor_mask(init_params,pi/2,[SDIM SDIM]);
            
            mask1_out = X_used*cur_mask1(:);
            mask2_out = X_used*cur_mask2(:);
            
            en_out(i,:) = sqrt(mask1_out.^2+mask2_out.^2);
        end
        tot_cont = mean(en_out);
        tot_cont = zscore(tot_cont);
        en_out = bsxfun(@minus,en_out,mean(en_out));
        en_out = bsxfun(@rdivide,en_out,std(en_out,[],2));
        init_params(3) = lfp_gabor_params(cc,ww,3);
        cur_mask1 = get_pgabor_mask(init_params,0,[SDIM SDIM]);
        cur_mask2 = get_pgabor_mask(init_params,pi/2,[SDIM SDIM]);
        mask1_out = X_used*cur_mask1(:);
        mask2_out = X_used*cur_mask2(:);
        pref_out = zscore(sqrt(mask1_out.^2+mask2_out.^2));
        temp = [pref_out tot_cont' ones(length(cur_used_fixs),1)];
        
        [B_both(cc,ww,:),BINT,R,RINT,STATS] = regress(squeeze(fix_avg_amps(cc,:,ww))',temp);
        rsq_both(cc,ww) = STATS(1);
        
        [B_both_early(cc,ww,:),BINT,R,RINT,STATS] = regress(squeeze(early_avg_amps(cc,:,ww))',temp);
        rsq_both_early(cc,ww) = STATS(1);
        
        [B_both_late(cc,ww,:),BINT,R,RINT,STATS] = regress(squeeze(late_avg_amps(cc,:,ww))',temp);
        rsq_both_late(cc,ww) = STATS(1);
    end
end

%%
[max_rsq,maxloc] = max(rsq,[],3);
[max_rsq_late,maxloc_late] = max(rsq_late,[],3);
[max_rsq_early,maxloc_early] = max(rsq_early,[],3);
cur_oris = linspace(0,pi-pi/8,8);
en_out = nan(8,n_fixs);
X_used = X_resh(cur_used_fixs,:);

used_freq = 8;
hold_const = [0 0 0 0 0 1 0 0];
LB = [-10 -10 0 6 2 0.5 -Inf -Inf];
UB = [10 10 pi SDIM(1)/2 SDIM(1)/6 2 Inf Inf Inf];
lfp_gabor_params = zeros(24,length(wfreqs),8);
LL = zeros(24,length(wfreqs));
lfp_gabor_params_early = zeros(24,length(wfreqs),8);
LL_early = zeros(24,length(wfreqs));
lfp_gabor_params_late = zeros(24,length(wfreqs),8);
LL_late = zeros(24,length(wfreqs));

init_params = zeros(1,6);
for cc = 1:24
    fprintf('CC %d of %d\n',cc,24);
    for ww = 1:length(wfreqs)
        fprintf('freq %d of %d\n',ww,length(wfreqs));
        init_params(3) = orientations(maxloc(cc,ww));
        cur_init_params = [init_params zeros(1,2)];
        [lfp_gabor_params(cc,ww,:),LL(cc,ww)] = fit_gabor_params_lfp(cur_init_params,X_resh(cur_used_fixs,:),...
            squeeze(fix_avg_amps(cc,:,ww)),[SDIM SDIM],hold_const,LB,UB);
        
        init_params(3) = orientations(maxloc_early(cc,ww));
        cur_init_params = [init_params zeros(1,2)];
        [lfp_gabor_params_early(cc,ww,:),LL_early(cc,ww)] = fit_gabor_params_lfp(cur_init_params,X_resh(cur_used_fixs,:),...
            squeeze(early_avg_amps(cc,:,ww)),[SDIM SDIM],hold_const,LB,UB);
        
        init_params(3) = orientations(maxloc_late(cc,ww));
        cur_init_params = [init_params zeros(1,2)];
        [lfp_gabor_params_late(cc,ww,:),LL_late(cc,ww)] = fit_gabor_params_lfp(cur_init_params,X_resh(cur_used_fixs,:),...
            squeeze(late_avg_amps(cc,:,ww)),[SDIM SDIM],hold_const,LB,UB);
    end
end


%% Compute TBR time-since fix onset
max_tsf = 0.75; nbins = 30;
used_tsf = time_since_fix;
used_tsf(used_tsf > max_tsf) = [];
tax = prctile(used_tsf,linspace(100/nbins,100-100/nbins,nbins));
% tax = 0.01:0.02:0.5;
% tax = logspace(log10(dt),log10(max_tsf),nbins);
Tmat = tbrep(time_since_fix,tax);

%%
xv_frac = 0.;
rperm = randperm(length(cur_used_fixs));
xv_fixs = rperm(1:round(xv_frac*length(cur_used_fixs)));
xv_inds = find(ismember(all_model_fixids,xv_fixs));
% tr_inds = find(~ismember(all_model_fixids,xv_fixs));
tr_inds = 1:length(all_model_fixids);
tr_inds(all_model_blockids(tr_inds) > 3) = [];

% xv_fixs = find(all_stim_filtered(cur_used_fixs)==0);
% xv_inds = find(ismember(all_model_fixids,xv_fixs));
% tr_inds = find(~ismember(all_model_fixids,xv_fixs));
%% COMPUTE SACCADE TRIG AVG LFP AMPS
sac_trg_amp = zeros(24,length(wfreqs),length(tax));
taxb = [0 tax];
for i = 1:length(tax)
    curset = find(time_since_fix(all_used_inds) >= taxb(i) & time_since_fix(all_used_inds) < taxb(i+1));
    sac_trg_amp(:,:,i) = mean(all_ampgrams(:,curset,:),2);
    n_occ(i) = length(curset);
end

%%
% used_freq = 8;
% sac_thresh = 1;
% min_sac_amp = 0;
% micro_sacs = find(all_sac_amps(used_fixs(cur_used_fixs)) < sac_thresh & all_sac_amps(used_fixs(cur_used_fixs)) > min_sac_amp);
% macro_sacs = find(all_sac_amps(used_fixs(cur_used_fixs)) > sac_thresh);
% %
% % all_sacs = randperm(length(cur_used_fixs));
% % micro_sacs = all_sacs(1:length(micro_sacs));
% % macro_sacs = all_sacs(end-length(micro_sacs)+1:end);
%
% % macro_sacs = macro_sacs(randperm(length(macro_sacs)));
% % macro_sacs(length(micro_sacs)+1:end) = [];
% %
% micro_set = find(ismember(all_model_fixids(all_used_inds),micro_sacs));
% macro_set = find(ismember(all_model_fixids(all_used_inds),macro_sacs));
%
%
% micro_sacs = find(all_stim_filtered(used_fixs(cur_used_fixs)) == 1);
% macro_sacs = find(all_stim_filtered(used_fixs(cur_used_fixs)) == 0);
% 
% % all_sacs = randperm(length(cur_used_fixs));
% % micro_sacs = all_sacs(1:length(micro_sacs));
% % macro_sacs = all_sacs(end-length(micro_sacs)+1:end);
% 
% micro_set = find(ismember(all_model_fixids(all_used_inds),micro_sacs));
% macro_set = find(ismember(all_model_fixids(all_used_inds),macro_sacs));

% set1 = find(blockids(used_fixs(cur_used_fixs)) == 1);
% set2 = find(blockids(used_fixs(cur_used_fixs)) == 3);
% micro_set = find(ismember(all_model_fixids(all_used_inds),set1));
% macro_set = find(ismember(all_model_fixids(all_used_inds),set2));


Tmat = round(Tmat);
for cc = 1:24
    % cc = 16;
    cc
    %     ww = 9;
    % for ww = 6:10
    for ww = 2
        %         for ww = 1:length(wfreqs)
        %         disp('Fitting sac hist model')
        %         %first fit sachist model
        %         Xmat = Tmat;
        %
        %         Y = all_ampgrams(cc,:,used_freq)';
        %         XX = [Tmat(all_used_inds,:)];
        %         %     [B_sachist(cc,:),~] = REGRESS(Y,XX);
        %
        %         init_params = randn(length(tax),1);
        %         [B_sachist(cc,:)] = smoothed_regression(Y,XX,init_params,2e3,[1 length(tax)]);
        
        disp('Fitting total stim model')
        X_resh = reshape(X,size(X,1),SDIM^2);
        cur_mask1 = get_pgabor_mask(lfp_gabor_params(cc,ww,1:6),0,[SDIM SDIM]);
        cur_mask2 = get_pgabor_mask(lfp_gabor_params(cc,ww,1:6),pi/2,[SDIM SDIM]);
        cur_gmask = get_pgauss_mask(lfp_gabor_params(cc,ww,1:6),[SDIM SDIM]);
        cur_gmask = cur_gmask(:);
        mask1_out = X_resh*cur_mask1(:);
        mask2_out = X_resh*cur_mask2(:);
        gmask_out = bsxfun(@times,X_resh,cur_gmask');
        loc_mean = mean(gmask_out,2);
        loc_cont = std(gmask_out,[],2);
        
        mask1_out = (mask1_out - loc_mean)./loc_cont;
        mask2_out = (mask2_out - loc_mean)./loc_cont;
        
        loc_mean = zscore(loc_mean);
        loc_cont = zscore(loc_cont);
        energy_out = lfp_gabor_params(cc,ww,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
        total_out = zscore(energy_out);
        
%         
%                  Y = all_ampgrams(cc,micro_set,ww)';
%                 XX = [Tmat(all_used_inds(micro_set),:) bsxfun(@times,Tmat(all_used_inds(micro_set),:),total_out(all_used_inds(micro_set)))];
%                 init_params = randn(2*length(tax),1);
% %                  [B_stimmod_micro(cc,ww,:)] = smoothed_regression(Y,XX,init_params,[2e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)]);
%         %          [B_stimmod_micro(cc,ww,:)] = smoothed_regression(Y,XX,init_params,[0 0],[1 length(tax); length(tax)+1 2*length(tax)]);
%                  [B_stimmod_micro(cc,ww,:),B_stimmod_micro_ci(cc,ww,:,:)] = regress(Y,XX);
%         
%                  Y = all_ampgrams(cc,macro_set,ww)';
%                 XX = [Tmat(all_used_inds(macro_set),:) bsxfun(@times,Tmat(all_used_inds(macro_set),:),total_out(all_used_inds(macro_set)))];
%                 init_params = randn(2*length(tax),1);
% %                  [B_stimmod_macro(cc,ww,:)] = smoothed_regression(Y,XX,init_params,[2e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)]);
%         %          [B_stimmod_macro(cc,ww,:)] = smoothed_regression(Y,XX,init_params,[0 0],[1 length(tax); length(tax)+1 2*length(tax)]);
%                  [B_stimmod_macro(cc,ww,:),B_stimmod_macro_ci(cc,ww,:,:)] = regress(Y,XX);
%         
        %         Y = all_ampgrams(cc,:,ww)';
        %         XX = [Tmat(all_used_inds,:) bsxfun(@times,Tmat(all_used_inds,:),total_out(all_used_inds))];
        %         %     [B_stimmod(cc,:),~] = REGRESS(Y,XX);
        %         init_params = randn(2*length(tax),1);
        % %         [B_stimmod(cc,ww,:)] = smoothed_regression(Y,XX,init_params,[2e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)]);
        %          [B_stimmod_test(cc,:)] = smoothed_regression(Y,XX,init_params,[2e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)]);
        
        for bb = 1:4
            set = find(blockids(used_fixs(cur_used_fixs)) == bb);
            set = find(ismember(all_model_fixids(all_used_inds),set));
        
            Y = all_ampgrams(cc,set,ww)';
            XX = [Tmat(all_used_inds(set),:) bsxfun(@times,Tmat(all_used_inds(set),:),total_out(all_used_inds(set)))];
            init_params = randn(2*length(tax),1);
            [B_stimmod_set(cc,ww,bb,:)] = smoothed_regression(Y,XX,init_params,[2e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)]);
%             [B_stimmod_set(cc,ww,bb,:),B_stimmod_set_ci(cc,ww,bb,:,:)] = regress(Y,XX);
        end
        
%         Y = all_ampgrams(cc,:,ww)';
%         XX = [Tmat(all_used_inds,:) bsxfun(@times,Tmat(all_used_inds,:),total_out(all_used_inds))];
%         init_params = randn(2*length(tax),1);
%         [B_stimmod(cc,ww,:)] = smoothed_regression(Y,XX,init_params,[5e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)]);
%                 [B_stimmod(cc,ww,:),B_stimmod_set_ci(cc,ww,:,:)] = regress(Y,XX);
                
%         
%             disp('Fitting separate stim model')
%     ori_vals = linspace(0,pi,13);
%     ori_vals(end) = [];
%     
%     ori_outs = zeros(length(total_out),length(ori_vals));
%     cur_params = lfp_gabor_params(cc,ww,1:6);
%     for oo = 1:length(ori_vals)
%         cur_params(3) = ori_vals(oo);
%         cur_mask1 = get_pgabor_mask(cur_params,0,[SDIM SDIM]);
%         cur_mask2 = get_pgabor_mask(cur_params,pi/2,[SDIM SDIM]);
%         mask1_out = X_resh*cur_mask1(:);
%         mask2_out = X_resh*cur_mask2(:);
% %         lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
% %         energy_out = gabor_params_fin(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
%         energy_out = sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
% %         ori_outs(:,oo) = lin_out + energy_out;
%         ori_outs(:,oo) = energy_out;
%     end
%     cur_mask1 = get_pgabor_mask(lfp_gabor_params(cc,ww,1:6),0,[SDIM SDIM]);
%     cur_mask2 = get_pgabor_mask(lfp_gabor_params(cc,ww,1:6),pi/2,[SDIM SDIM]);
%     mask1_out = X_resh*cur_mask1(:);
%     mask2_out = X_resh*cur_mask2(:);
% %     lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
%     energy_out = sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
%     total_out = energy_out;
%     
%     contrast_out = sum(ori_outs,2);
%     total_out = total_out./contrast_out;
%     
%     contrast_out = zscore(contrast_out);
%     total_out = zscore(total_out);
%         
%     Xmat = Tmat(all_used_inds,:);
%     Xmat = [Xmat bsxfun(@times,Tmat(all_used_inds,:),total_out(all_used_inds))];
%     Xmat = [Xmat bsxfun(@times,Tmat(all_used_inds,:),contrast_out(all_used_inds))];
%     init_params = randn(3*length(tax),1);
%     [sep_mod] = smoothed_regression(Y,Xmat,init_params,[2e3 5e3 5e3],[1 length(tax); length(tax)+1 2*length(tax); 2*length(tax)+1 3*length(tax)]);
%     sep_mod_sachist(cc,ww,:) = sep_mod(1:length(tax));
%     sep_mod_stim(cc,ww,:) = sep_mod(length(tax)+1:length(tax)*2);
%     sep_mod_cont(cc,ww,:) = sep_mod(length(tax)*2+1:length(tax)*3);
    
    end
    
end

%%
close all
used_ch = 1:24;
for i = 1:length(used_ch)
    subplot(4,6,i)
    pcolor(tax,wfreqs,squeeze(B_stimmod(used_ch(i),:,length(tax)+1:end)));shading flat
    % caxis([-0.25 0.5])
%     caxis([-0.25 2])
%     ylim([5 50])
    xlim([0 0.5])
    set(gca,'yscale','log')
    xlabel('Time (s)','fontsize',14)
    ylabel('Frequency','fontsize',14)
    title(sprintf('Channel %d',used_ch(i)));
end

% figure
% used_ch = 1:2:23;
% for i = 1:length(used_ch)
%     subplot(3,4,i)
%     pcolor(tax,wfreqs,squeeze(B_stimmod(used_ch(i),:,31:end)));shading flat
% % caxis([-0. 0.25])
% %  caxis([-0. 0.45])
% % ylim([5 50])
% xlim([0 0.4])
%    xlabel('Time (s)','fontsize',14)
%     ylabel('Frequency','fontsize',14)
%     title(sprintf('Channel %d',used_ch(i)));
%
% end

% figure
% for i = 1:6
%     subplot(3,2,i)
%     pcolor(tax,wfreqs,squeeze(B_stimmod_set_ci(used_ch(i),:,26:end,1)));shading flat
% caxis([-0. 0.25])
%     xlabel('Time (s)','fontsize',14)
%     ylabel('Frequency','fontsize',14)
%
% end
% %
% fillPage(gcf,'margins',[0 0 0 0],'papersize',[20 13]);
% print('-dpdf','-painters','Fix_trig_avg2');
% close

%%
for i = 1:3
    subplot(3,3,(i-1)*3+1)
    pcolor(tax,1:24,squeeze(B_stimmod_set(:,ww,i,1:25)));shading flat
    caxis([-0.3 0.4])
    xlim([0 0.5])
    xlabel('Time (s)','fontsize',14)
    ylabel('Channel','fontsize',14)
    
    subplot(3,3,(i-1)*3+2)
    pcolor(tax,1:24,squeeze(B_stimmod_set_ci(:,ww,i,1:25)));shading flat
    caxis([-0.3 0.4])
    xlim([0 0.5])
    xlabel('Time (s)','fontsize',14)
    ylabel('Channel','fontsize',14)
    
    subplot(3,3,(i-1)*3+3)
    pcolor(tax,1:24,squeeze(B_stimmod_set_ci(:,ww,i,1:25)));shading flat
    caxis([-0.3 0.4])
    xlim([0 0.5])
    xlabel('Time (s)','fontsize',14)
    ylabel('Channel','fontsize',14)
    
end
shg
%%
ww = 2;
for i = 1:3
    subplot(3,1,i)
    pcolor(tax,1:24,squeeze(B_stimmod_set(:,ww,i,31:end)));shading flat
    caxis([-0. 0.3])
    xlim([0 0.5])
    xlabel('Time (s)','fontsize',14)
    ylabel('Channel','fontsize',14)
%     
%     subplot(3,3,(i-1)*3+2)
%     pcolor(tax,1:24,squeeze(B_stimmod_set(:,ww,i,31:end)));shading flat
%     caxis([-0. 0.4])
%     xlim([0 0.5])
%     xlabel('Time (s)','fontsize',14)
%     ylabel('Channel','fontsize',14)
%     
%     subplot(3,3,(i-1)*3+3)
%     pcolor(tax,1:24,squeeze(B_stimmod_set(:,ww,i,31:end)));shading flat
%     caxis([-0. 0.4])
%     xlim([0 0.5])
%     xlabel('Time (s)','fontsize',14)
%     ylabel('Channel','fontsize',14)
    
end
shg

%%
figure
pcolor(tax,1:24,B_stimmod(:,31:end));shading flat
xlabel('Time (s)','fontsize',14)
ylabel('Channel Number','fontsize',14)

figure
pcolor(tax,1:24,B_stimmod(:,1:30));shading flat
xlabel('Time (s)','fontsize',14)
ylabel('Channel Number','fontsize',14)

%%
close all
ch = 16;
figure
pcolor(tax,wfreqs,squeeze(B_stimmod(ch,:,31:end)));shading flat
xlabel('Time (s)','fontsize',14)
ylabel('Frequency (Hz)','fontsize',14)
caxis([0 0.35])

figure
pcolor(tax,wfreqs,squeeze(B_stimmod(ch,:,1:30)));shading flat
xlabel('Time (s)','fontsize',14)
ylabel('Frequency (Hz)','fontsize',14)

%%
close all
cfreq = 3;
figure
set(gca,'fontname','arial')
pcolor(tax,1:24,squeeze(B_stimmod(:,cfreq,31:end)));shading flat
xlabel('Time (s)','fontsize',14)
ylabel('Frequency (Hz)','fontsize',14)
caxis([0 0.23])
xlim([0.01 0.4])
colorbar

figure
set(gca,'fontname','arial')
pcolor(tax,1:24,squeeze(B_stimmod(:,cfreq,1:30)));shading flat
xlabel('Time (s)','fontsize',14)
ylabel('Frequency (Hz)','fontsize',14)
xlim([0.01 0.4])
caxis([-0.2 0.3])
colorbar

%%
% cd /Users/James/Data/bruce/2_27_12
% save temp_stimgamma tax wfreqs B_stimmod 
