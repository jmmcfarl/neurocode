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
[b_alpha,a_alpha] = butter(2,[5 12]/niqf);
alpha_ch = 2;
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
cd ~/James_scripts/bruce/modelfits/
load ./fourier_mod_fits_v3 

% spatial_scale = 6;
% xax = -floor(SDIM/2):floor(SDIM/2);
% yax = -floor(SDIM/2):floor(SDIM/2);
% xax(SDIM+1:end) = [];
% yax(SDIM+1:end) = [];
% [XX,YY] = meshgrid(xax,yax);
% for cc = 1:length(use_lfps)
%     for ww = 1:length(wfreqs)
%         cur_params = squeeze(lfp_gabor_params(cc,ww,1:6));
%         cur_params(5) = spatial_scale;
%         % cur_dist = sqrt((XX-cur_params(1)).^2 + (YY-cur_params(2)).^2);
%         uset = find(abs(XX-round(cur_params(1))) <= 2*cur_params(5) & abs(YY-round(cur_params(2))) <= 2*cur_params(5));
%         
%         gauss_mask = get_pgauss_mask(cur_params,[SDIM SDIM]);
%         gauss_mask = gauss_mask(uset);
%         X_out = bsxfun(@times,X_resh(:,uset),gauss_mask');
%         n_dim = length(uset);
%         per_dim = sqrt(n_dim);
%         X_out = reshape(X_out,[size(X_out,1) per_dim per_dim]);
%         fft_out = nan(size(X_out));
%         for i = 1:size(X_out,1)
%             fft_out(i,:,:) = abs(fftshift(fft2(squeeze(X_out(i,:,:)))));
%         end
%         pred_mat = reshape(fft_out,size(X_out,1),per_dim^2);
%         pred_mat = zscore(pred_mat);
%         
% %         sparse_lambda = 500;
% %         smooth_lambda = 2500;
% %         init_beta = zeros(n_dim,1);
% %         [beta(cc,ww,:),fft_LL(cc,ww)] = smooth_regress_2d(squeeze(fix_avg_amps(cc,tr_inds,ww))',pred_mat(tr_fixs,:),init_beta,smooth_lambda,sparse_lambda);
%         fft_pred_out(cc,ww,:) = pred_mat*squeeze(beta(cc,ww,:));       
%     end
% end



%% Compute TBR time-since fix onset
% max_tsf = 0.75; nbins = 30;
% used_tsf = time_since_fix;
% used_tsf(used_tsf > max_tsf) = [];
% tax = prctile(used_tsf,linspace(100/nbins,100-100/nbins,nbins));
% % tax = 0.01:0.02:0.5;
% % tax = logspace(log10(dt),log10(max_tsf),nbins);
% Tmat = tbrep(time_since_fix,tax);

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

% xv_fixs = find(all_stim_filtered(cur_used_fixs)==0);
% xv_inds = find(ismember(all_model_fixids,xv_fixs));
% tr_inds = find(~ismember(all_model_fixids,xv_fixs));
%% COMPUTE SACCADE TRIG AVG LFP AMPS
% sac_trg_amp = zeros(24,length(wfreqs),length(tax));
% taxb = [0 tax];
% for i = 1:length(tax)
%     curset = find(time_since_fix(all_used_inds) >= taxb(i) & time_since_fix(all_used_inds) < taxb(i+1));
%     sac_trg_amp(:,:,i) = mean(all_ampgrams(:,curset,:),2);
%     n_occ(i) = length(curset);
% end

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
