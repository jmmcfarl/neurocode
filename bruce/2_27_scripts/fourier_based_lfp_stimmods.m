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
scales = logspace(log10(2.5),log10(50),30);
% scales = [ 5.7 ];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

use_lfps = 1:2:24;
%%
all_lfp_pow = [];
within_fix_avgs = [];
fix_nums = [];
all_used_inds = [];
all_ampgrams = [];
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
    lfp_samps = lfp_samps(:,use_lfps);
    
    lfp_samps = filtfilt(b,a,lfp_samps);
    lfp_sampsd = downsample(lfp_samps,dsf);
    %     lfp_sampsd = abs(hilbert(lfp_sampsd));
    %     lfp_sampsd = zscore(lfp_sampsd);
    lfp_timed = downsample(lfp_time,dsf);
    
    cur_all_model = find(all_model_blockids==blockid);
    interp_ampgram = zeros(length(use_lfps),length(cur_all_model),length(wfreqs));
    for cc = 1:length(use_lfps)
        %     cc = 17;
        fprintf('Channel %d of %d\n',cc,length(use_lfps));
        temp = cwt(lfp_sampsd(:,cc),scales,'cmor1-1');
        tampgram = abs(temp);
        %         tampgram = bsxfun(@minus,tampgram,mean(tampgram,2));
        %         tampgram = bsxfun(@rdivide,tampgram,std(tampgram,[],2));
        tampgram = tampgram';
        interp_ampgram(cc,:,:) = interp1(lfp_timed,tampgram,all_model_time_axis(cur_all_model));
    end
    clear temp tampgram
    
    %%
    used_tpoints = find(~isnan(interp_ampgram(1,:,1)));
    all_ampgrams = cat(2,all_ampgrams,interp_ampgram(:,used_tpoints,:));
    all_used_inds = [all_used_inds; cur_all_model(used_tpoints)'];
    
end
% all_ampgrams = sqrt(all_ampgrams);
all_ampgrams = bsxfun(@minus,all_ampgrams,mean(all_ampgrams,2));
all_ampgrams = bsxfun(@rdivide,all_ampgrams,std(all_ampgrams,[],2));

%%
cur_used_fixs = unique(all_model_fixids(all_used_inds));
n_fixs = length(cur_used_fixs);
fix_avg_amps = zeros(length(use_lfps),n_fixs,length(wfreqs));
late_avg_amps = zeros(length(use_lfps),n_fixs,length(wfreqs));
early_avg_amps = zeros(length(use_lfps),n_fixs,length(wfreqs));
for i = 1:n_fixs
    cur_set = find(all_model_fixids(all_used_inds) == cur_used_fixs(i));
    fix_avg_amps(:,i,:) = squeeze(mean(all_ampgrams(:,cur_set,:),2));
    cur_set2 = cur_set(time_since_fix(all_used_inds(cur_set)) >= 0.15 & time_since_fix(all_used_inds(cur_set)) < 0.3);
    late_avg_amps(:,i,:) = squeeze(mean(all_ampgrams(:,cur_set2,:),2));
    cur_set2 = cur_set(time_since_fix(all_used_inds(cur_set)) < 0.15);
    early_avg_amps(:,i,:) = squeeze(mean(all_ampgrams(:,cur_set2,:),2));
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
en_out = bsxfun(@minus,en_out,mean(en_out,2));
en_out = bsxfun(@rdivide,en_out,std(en_out,[],2));

rsq = zeros(length(use_lfps),length(wfreqs),12);
rsq_late = zeros(length(use_lfps),length(wfreqs),12);
rsq_early = zeros(length(use_lfps),length(wfreqs),12);
for cc = 1:length(use_lfps)
    for ww = 1:length(wfreqs)
        %     for ww = 9
        for i = 1:12
            temp = [en_out(i,:)' ones(length(cur_used_fixs),1)];
            [B,BINT,R,RINT,STATS] = regress(squeeze(fix_avg_amps(cc,:,ww))',temp);
            rsq(cc,ww,i) = STATS(1);
            %             [B,BINT,R,RINT,STATS] = regress(squeeze(late_avg_amps(cc,:,ww))',temp);
            %             rsq_late(cc,ww,i) = STATS(1);
            %             [B,BINT,R,RINT,STATS] = regress(squeeze(early_avg_amps(cc,:,ww))',temp);
            %             rsq_early(cc,ww,i) = STATS(1);
        end
    end
end

[max_rsq,maxloc] = max(rsq,[],3);
% [max_rsq_late,maxloc_late] = max(rsq_late,[],3);
% [max_rsq_early,maxloc_early] = max(rsq_early,[],3);

%%
tot_nfixs = length(cur_used_fixs);
xv_frac = 0.;
xv_fixs = randperm(tot_nfixs);
xv_fixs(round(tot_nfixs*xv_frac)+1:end) = [];
xv_fixs = cur_used_fixs(xv_fixs);
tr_fixs = setdiff(cur_used_fixs,xv_fixs);

xv_inds = find(ismember(cur_used_fixs,xv_fixs));
tr_inds = find(ismember(cur_used_fixs,tr_fixs));
%%

used_freq = 8;
hold_const = [0 0 0 0 0 1 0 1 1 0];
LB = [-10 -10 0 6 2 0.5 -Inf -Inf -Inf -Inf];
UB = [10 10 pi SDIM(1)/2 SDIM(1)/6 2 Inf Inf Inf Inf];
lfp_gabor_params = zeros(length(use_lfps),length(wfreqs),10);
LL_gabor = zeros(length(use_lfps),length(wfreqs));
for cc = 1:length(use_lfps)
    fprintf('CC %d of %d\n',cc,length(use_lfps));
    for ww = 1:length(wfreqs)
        %     for ww = 9
        fprintf('freq %d of %d\n',ww,length(wfreqs));
        init_params(3) = orientations(maxloc(cc,ww));
        cur_init_params = [init_params zeros(1,4)];
        [lfp_gabor_params(cc,ww,:),LL_gabor(cc,ww)] = fit_gabor_params_lfp(cur_init_params,X_resh(tr_fixs,:),...
            squeeze(fix_avg_amps(cc,tr_inds,ww)),[SDIM SDIM],hold_const,LB,UB);
        
        if xv_frac > 0
        cur_mask1 = get_pgabor_mask(lfp_gabor_params(cc,ww,1:6),0,[SDIM SDIM]);
        cur_mask2 = get_pgabor_mask(lfp_gabor_params(cc,ww,1:6),pi/2,[SDIM SDIM]);
        
        mask1_out = X_resh(xv_fixs,:)*cur_mask1(:);
        mask2_out = X_resh(xv_fixs,:)*cur_mask2(:);
        
        energy_out = lfp_gabor_params(cc,ww,7)*sqrt(mask1_out.^2+mask2_out.^2);
        g = energy_out + lfp_gabor_params(cc,ww,8);
        xvLL_gabor(cc,ww) = sum((g - squeeze(fix_avg_amps(cc,xv_inds,ww))').^2);
        end
        %         init_params(3) = orientations(maxloc_late(cc,ww));
        %         cur_init_params = [init_params zeros(1,4)];
        %         [lfp_gabor_params_late(cc,ww,:),LL_fin_late(cc,ww)] = fit_gabor_params_lfp(cur_init_params,X_resh(cur_used_fixs,:),...
        %             squeeze(late_avg_amps(cc,:,ww)),[SDIM SDIM],hold_const,LB,UB);
        %
        %         init_params(3) = orientations(maxloc_early(cc,ww));
        %         cur_init_params = [init_params zeros(1,4)];
        %         [lfp_gabor_params_early(cc,ww,:),LL_fin_early(cc,ww)] = fit_gabor_params_lfp(cur_init_params,X_resh(cur_used_fixs,:),...
        %             squeeze(early_avg_amps(cc,:,ww)),[SDIM SDIM],hold_const,LB,UB);
    end
end

%%
clear *_pred_out beta*
spatial_scale = 6;
xax = -floor(SDIM/2):floor(SDIM/2);
yax = -floor(SDIM/2):floor(SDIM/2);
xax(SDIM+1:end) = [];
yax(SDIM+1:end) = [];
[XX,YY] = meshgrid(xax,yax);
for cc = 1:length(use_lfps)
    for ww = 1:length(wfreqs)
        cur_params = squeeze(lfp_gabor_params(cc,ww,1:6));
        cur_params(5) = spatial_scale;
        % cur_dist = sqrt((XX-cur_params(1)).^2 + (YY-cur_params(2)).^2);
        uset = find(abs(XX-round(cur_params(1))) <= 2*cur_params(5) & abs(YY-round(cur_params(2))) <= 2*cur_params(5));
        
        gauss_mask = get_pgauss_mask(cur_params,[SDIM SDIM]);
        gauss_mask = gauss_mask(uset);
        X_out = bsxfun(@times,X_resh(:,uset),gauss_mask');
        n_dim = length(uset);
        per_dim = sqrt(n_dim);
        X_out = reshape(X_out,[size(X_out,1) per_dim per_dim]);
        fft_out = nan(size(X_out));
        for i = 1:size(X_out,1)
            fft_out(i,:,:) = abs(fftshift(fft2(squeeze(X_out(i,:,:)))));
        end
        pred_mat = reshape(fft_out,size(X_out,1),per_dim^2);
        pred_mat = zscore(pred_mat);
        
        sparse_lambda = 500;
        smooth_lambda = 2500;
        init_beta = zeros(n_dim,1);
        [beta(cc,ww,:),fft_LL(cc,ww)] = smooth_regress_2d(squeeze(fix_avg_amps(cc,tr_inds,ww))',pred_mat(tr_fixs,:),init_beta,smooth_lambda,sparse_lambda);
        fft_pred_out(cc,ww,:) = pred_mat*squeeze(beta(cc,ww,:));
        
        [beta_early(cc,ww,:),fft_LL_early(cc,ww)] = smooth_regress_2d(squeeze(early_avg_amps(cc,tr_inds,ww))',pred_mat(tr_fixs,:),init_beta,smooth_lambda,sparse_lambda);
        fft_pred_out_early(cc,ww,:) = pred_mat*squeeze(beta_early(cc,ww,:));
        
         [beta_late(cc,ww,:),fft_LL_late(cc,ww)] = smooth_regress_2d(squeeze(late_avg_amps(cc,tr_inds,ww))',pred_mat(tr_fixs,:),init_beta,smooth_lambda,sparse_lambda);
        fft_pred_out_late(cc,ww,:) = pred_mat*squeeze(beta_late(cc,ww,:));

        if xv_frac > 0
            fft_pred_xv = pred_mat(xv_fixs,:)*squeeze(beta(cc,ww,:));
            fft_xvLL(cc,ww) = sum((fft_pred_xv - squeeze(fix_avg_amps(cc,xv_inds,ww))').^2);
        end
    end
end
%% Compute TBR time-since fix onset
max_tsf = 0.75; %nbins = 30;
used_tsf = time_since_fix;
used_tsf(used_tsf > max_tsf) = [];
% tax = prctile(used_tsf,linspace(100/nbins,100-100/nbins,nbins));
tax = 0.01:0.0125:0.5;
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
sac_trg_amp = zeros(length(use_lfps),length(wfreqs),length(tax));
taxb = [0 tax];
for i = 1:length(tax)
    curset = find(time_since_fix(all_used_inds) >= taxb(i) & time_since_fix(all_used_inds) < taxb(i+1));
    sac_trg_amp(:,:,i) = mean(all_ampgrams(:,curset,:),2);
    n_occ(i) = length(curset);
end

%%
clear B_stimmod*

Tmat = round(Tmat);
for cc = 1:length(use_lfps)
    for ww = 1:length(wfreqs)
        
        total_out = zscore(squeeze(fft_pred_out(cc,ww,all_model_fixids(all_used_inds)))');
        total_out_early = zscore(squeeze(fft_pred_out_early(cc,ww,all_model_fixids(all_used_inds)))');
        total_out_late = zscore(squeeze(fft_pred_out_late(cc,ww,all_model_fixids(all_used_inds)))');
        
        Y = all_ampgrams(cc,:,ww)';
        XX = [Tmat(all_used_inds,:) bsxfun(@times,Tmat(all_used_inds,:),total_out')];
        init_params = randn(2*length(tax),1);%
        [B_stimmod(cc,ww,:)] = smoothed_regression(Y,XX,init_params,[1e3 2e3],[1 length(tax); length(tax)+1 2*length(tax)],0);
        
        XX = [Tmat(all_used_inds,:) bsxfun(@times,Tmat(all_used_inds,:),total_out_early')];
        [B_stimmod_early(cc,ww,:)] = smoothed_regression(Y,XX,init_params,[1e3 2e3],[1 length(tax); length(tax)+1 2*length(tax)],0);
        
        XX = [Tmat(all_used_inds,:) bsxfun(@times,Tmat(all_used_inds,:),total_out_late')];
        [B_stimmod_late(cc,ww,:)] = smoothed_regression(Y,XX,init_params,[1e3 2e3],[1 length(tax); length(tax)+1 2*length(tax)],0);
    end
    
end

%%
cd ~/James_scripts/bruce/modelfits/
save fourier_mod_fits_v3 beta* per_dim tax B_stimmod* wfreqs scales lfp_gabor_params LL* *xvLL* use_lfps
%%
cur_set = 1:2:length(use_lfps);
for ww = 1:length(wfreqs);
    fprintf('Freq: %.3f\n',wfreqs(ww));
    for cc = 1:length(cur_set)
        subplot(3,length(cur_set),cc)
        imagesc(reshape(beta(cur_set(cc),ww,:),per_dim,per_dim))
        cax = caxis();
        mc = max(abs(cax));
        caxis([-mc mc]);
        subplot(3,length(cur_set),cc+length(cur_set))
        imagesc(reshape(beta_early(cur_set(cc),ww,:),per_dim,per_dim))
        caxis([-mc mc]);
        subplot(3,length(cur_set),cc+2*length(cur_set))
        imagesc(reshape(beta_late(cur_set(cc),ww,:),per_dim,per_dim))
        caxis([-mc mc]);
    end
    
    pause
    clf
end

%%

%%
for cc = 1:24
    fprintf('Ch %d of %d\n',cc,24);
    for ww = 1:length(wfreqs);
        subplot(5,6,ww)
        imagesc(reshape(beta(cc,ww,:),per_dim,per_dim));
        title(sprintf('%.3f',wfreqs(ww)));
    end
    
    pause
    clf
end

%%
for cc = 1:length(use_lfps)
   subplot(3,1,1)
   pcolor(tax,wfreqs,squeeze(B_stimmod(cc,:,length(tax)+1:end)));shading flat
   colorbar
   set(gca,'yscale','log')
    subplot(3,1,2)
   pcolor(tax,wfreqs,squeeze(B_stimmod_early(cc,:,length(tax)+1:end)));shading flat
   colorbar
      set(gca,'yscale','log')
subplot(3,1,3)
   pcolor(tax,wfreqs,squeeze(B_stimmod_late(cc,:,length(tax)+1:end)));shading flat
   colorbar
      set(gca,'yscale','log')

   cc
  pause
    clf
end