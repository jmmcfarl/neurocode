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

cd ~/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

cd ~/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

SDIM = length(xpatch_inds);
cd ~/Data/bruce/2_27_12/stimrecon
load fixation_image_patches_corrected_raw

NT = length(used_fixs);
X_resh = reshape(X,NT,SDIM^2);
dt = .005;

cd ~/James_scripts/bruce/modelfits
load pref_oris
cd ~/James_scripts/bruce/modelfits/
% load fixation_gabor_models
load gabor_tempmodfits_en_lin

%%
Fs = 1000;
dsf = 2;
Fsd = Fs/dsf;
niqf = Fs/2;
% [b,a] = butter(4,[45 60]/niqf);
[b,a] = butter(2,120/niqf,'low');
% scales = logspace(log10(2.5),log10(35),20);
scales = logspace(log10(2),log10(60),30);
% scales = [5.15 5.7 6.3];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
[b_alpha,a_alpha] = butter(2,[5 12]/niqf);
alpha_ch = 1;
[b_gamma,a_gamma] = butter(2,[25 60]/niqf);
%%
all_lfp_pow = [];
within_fix_avgs = [];
fix_nums = [];
all_used_inds = [];
all_ampgrams = [];
% all_alpha_freq = [];
% all_gamma_freq = [];
all_alpha_phase = [];
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    
    %%
    cd ~/Data/bruce/2_27_12
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
%     interp_gammafreq = zeros(24,length(cur_all_model));
%     interp_alphafreq = zeros(24,length(cur_all_model));
    for cc = 1:24
        %     cc = 17;
        fprintf('Channel %d of %d\n',cc,24);
        temp = cwt(lfp_sampsd(:,cc),scales,'cmor1-1');
        tampgram = abs(temp);
        %         tampgram = bsxfun(@minus,tampgram,mean(tampgram,2));
        %         tampgram = bsxfun(@rdivide,tampgram,std(tampgram,[],2));
        tampgram = tampgram';
        interp_ampgram(cc,:,:) = interp1(lfp_timed,tampgram,all_model_time_axis(cur_all_model));
                
%         gamma_freq = dnf_hilbert_instfreq(lfp_gamma(:,cc),Fsd,5);
%         alpha_freq = dnf_hilbert_instfreq(lfp_alpha(:,cc),Fsd,5);
%         interp_gammafreq(cc,:) = interp1(lfp_timed,gamma_freq.f_smooth,all_model_time_axis(cur_all_model));
%         interp_alphafreq(cc,:) = interp1(lfp_timed,alpha_freq.f_smooth,all_model_time_axis(cur_all_model));
    end
    clear temp tampgram
    
    alpha_phase = angle(hilbert(lfp_alpha(:,alpha_ch)));
    alpha_phase = unwrap_phase_monotonic(alpha_phase);
    interp_alpha_phase = interp1(lfp_timed,alpha_phase,all_model_time_axis(cur_all_model));
    
    %%
    used_tpoints = find(~isnan(interp_ampgram(1,:,1)));
    all_ampgrams = cat(2,all_ampgrams,interp_ampgram(:,used_tpoints,:));
    all_alpha_phase = [all_alpha_phase; interp_alpha_phase(:)];
%     all_gamma_freq = cat(2,all_gamma_freq,interp_gammafreq(:,used_tpoints,:));
%     all_alpha_freq = cat(2,all_alpha_freq,interp_alphafreq(:,used_tpoints,:));
   all_used_inds = [all_used_inds; cur_all_model(used_tpoints)'];
    
end
% all_ampgrams = sqrt(all_ampgrams);
all_ampgrams = bsxfun(@minus,all_ampgrams,mean(all_ampgrams,2));
all_ampgrams = bsxfun(@rdivide,all_ampgrams,std(all_ampgrams,[],2));

%%
cur_used_fixs = unique(all_model_fixids(all_used_inds));
n_fixs = length(cur_used_fixs);
fix_avg_amps = zeros(24,n_fixs,length(wfreqs));
for i = 1:n_fixs
    cur_set = find(all_model_fixids(all_used_inds) == cur_used_fixs(i));
    fix_avg_amps(:,i,:) = squeeze(mean(all_ampgrams(:,cur_set,:),2));
end
fix_avg_amps = bsxfun(@minus,fix_avg_amps,mean(fix_avg_amps,2));
fix_avg_amps = bsxfun(@rdivide,fix_avg_amps,std(fix_avg_amps,[],2));

%%
phase_since_trial = nan(size(all_alpha_phase));
for i = 1:n_fixs
    cur_set = find(all_model_fixids(all_used_inds) == cur_used_fixs(i));
   phase_since_trial(cur_set) = all_alpha_phase(cur_set) - all_alpha_phase(cur_set(1));
end

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
en_out = bsxfun(@minus,en_out,mean(en_out));
en_out = bsxfun(@rdivide,en_out,std(en_out,[],2));

rsq = zeros(24,length(wfreqs),12);
R_vals = zeros(24,length(wfreqs),12,3);

for cc = 1:24
    fprintf('Channel %d of %d\n',cc,24);
    for ww = 1:length(wfreqs)
        %     for ww = 9
        for i = 1:12
            temp = [en_out(i,:)' ones(length(cur_used_fixs),1)];
            [B,BINT,R,RINT,STATS] = regress(squeeze(fix_avg_amps(cc,:,ww))',temp);
            rsq(cc,ww,i) = STATS(1);           
        end
    end
end

%%
[max_rsq,maxloc] = max(rsq,[],3);
cur_oris = linspace(0,pi-pi/8,8);
en_out = nan(8,n_fixs);
X_used = X_resh(cur_used_fixs,:);

hold_const = [0 0 0 0 0 1 0 0];
LB = [-5 -5 0 6 2 0.5 -Inf -Inf];
UB = [10 10 pi SDIM(1)/2 SDIM(1)/6 2 Inf Inf Inf];
lfp_gabor_params = zeros(24,length(wfreqs),8);
LL_fin = zeros(24,length(wfreqs));
for cc = 1:24
    fprintf('CC %d of %d\n',cc,24);
    for ww = 1:length(wfreqs)
        fprintf('freq %d of %d\n',ww,length(wfreqs));
        init_params(3) = orientations(maxloc(cc,ww));
        cur_init_params = [init_params zeros(1,2)];
        [lfp_gabor_params(cc,ww,:),LL_fin(cc,ww)] = fit_gabor_params_lfp(cur_init_params,X_resh(cur_used_fixs,:),...
            squeeze(fix_avg_amps(cc,:,ww)),[SDIM SDIM],hold_const,LB,UB);
    end
end

%%
spatial_scale = 6;
xax = -floor(SDIM/2):floor(SDIM/2);
yax = -floor(SDIM/2):floor(SDIM/2);
xax(SDIM+1:end) = [];
yax(SDIM+1:end) = [];
[XX,YY] = meshgrid(xax,yax);
for cc = 1:24
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
        [beta(cc,ww,:),fft_LL(cc,ww)] = smooth_regress_2d(squeeze(fix_avg_amps(cc,:,ww))',pred_mat(:,:),init_beta,smooth_lambda,sparse_lambda);
        fft_pred_out(cc,ww,:) = pred_mat*squeeze(beta(cc,ww,:));       
    end
end

%% Compute TBR time-since fix onset
max_phase = 2*pi*5;
nbins = 60;
uset = find(phase_since_trial <= max_phase);
tax = linspace(0,max_phase,nbins);
Tmat = tbrep(phase_since_trial,tax);


%%
xv_frac = 0.;
rperm = randperm(length(cur_used_fixs));
xv_fixs = rperm(1:round(xv_frac*length(cur_used_fixs)));
xv_inds = find(ismember(all_model_fixids,xv_fixs));
% tr_inds = find(~ismember(all_model_fixids,xv_fixs));
tr_inds = 1:length(all_model_fixids);
tr_inds = tr_inds(ismember(tr_inds,uset));
% tr_inds(all_model_blockids(tr_inds) > 3) = [];


%%
for cc = 1:24
    for ww = 1:length(wfreqs)
        
        disp('Fitting total stim model')
        X_resh = reshape(X,size(X,1),SDIM^2);
        cur_mask1 = get_pgabor_mask(lfp_gabor_params(cc,ww,1:6),0,[SDIM SDIM]);
        cur_mask2 = get_pgabor_mask(lfp_gabor_params(cc,ww,1:6),pi/2,[SDIM SDIM]);
        mask1_out = X_resh*cur_mask1(:);
        mask2_out = X_resh*cur_mask2(:);
                
        energy_out = lfp_gabor_params(cc,ww,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
        total_out = zscore(energy_out);
        
        
        Y = all_ampgrams(cc,tr_inds,ww)';
        XX = [Tmat(tr_inds,:) bsxfun(@times,Tmat(tr_inds,:),total_out(tr_inds))];
        init_params = randn(2*length(tax),1);
        [B_stimmod_gabor(cc,ww,:)] = smoothed_regression(Y,XX,init_params,[2e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)],0);

    
        total_out = squeeze(zscore(fft_pred_out(cc,ww,all_model_fixids)));
                
        Y = all_ampgrams(cc,tr_inds,ww)';
        XX = [Tmat(tr_inds,:) bsxfun(@times,Tmat(tr_inds,:),total_out(tr_inds))];
        init_params = randn(2*length(tax),1);
        [B_stimmod_fft(cc,ww,:)] = smoothed_regression(Y,XX,init_params,[2e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)],0);

    end
end

%%
save lfp_alphaphase_models_ch1 wfreqs tax B_stimmod* lfp_gabor_params beta 

%%
stim_vals = linspace(-2.5,2.5,10);

close all
uw = find(wfreqs >= 30 & wfreqs <= 60);
for cc = 1:24
    
%     subplot(2,2,1)
%     pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_gabor(cc,:,length(tax)+1:end)));shading flat
%     set(gca,'yscale','log')
% %     colorbar
%     caxis([0 0.3])
%     subplot(2,2,2)
%     pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_gabor(cc,:,1:length(tax))));shading flat
%     set(gca,'yscale','log')
    subplot(3,2,1)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft(cc,:,length(tax)+1:end)));shading flat
    set(gca,'yscale','log')
%     colorbar
    caxis([0 0.3])
    caxis([-0 0.275])
    subplot(3,2,2)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft(cc,:,1:length(tax))));shading flat
    set(gca,'yscale','log')
   caxis([-0.3 0.3])
   
   avg_gam_stim_dep = squeeze(mean(B_stimmod_fft(cc,uw,length(tax)+1:end),2));
   avg_gam_stim_ind = squeeze(mean(B_stimmod_fft(cc,uw,1:length(tax)),2));
   pred_outs = bsxfun(@times,avg_gam_stim_dep,stim_vals);
   pred_outs = bsxfun(@plus,pred_outs,avg_gam_stim_ind);
    subplot(3,2,3)
    plot(tax/(2*pi),avg_gam_stim_dep);shading flat
%     set(gca,'yscale','log')
    caxis([-0 0.275])
    subplot(3,2,4)
    plot(tax/(2*pi),avg_gam_stim_ind) ;shading flat
%     set(gca,'yscale','log')
    subplot(3,2,5)
    plot(tax/(2*pi),pred_outs)
   cc
    pause
    clf
end
    
%%
% for cc = 1:24
%     fprintf('Ch %d of %d\n',cc,24);
%     for ww = 1:length(wfreqs);
%         subplot(5,6,ww)
%         imagesc(reshape(beta(cc,ww,:),per_dim,per_dim));
%         title(sprintf('%.3f',wfreqs(ww)));
%     end
%     
%     pause
%     clf
% end

