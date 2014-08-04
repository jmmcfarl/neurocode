
%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt2_compiled_windata_d1p5.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

fix_win_dur = 0.25;
early_win_dur = 0.15;

Pix2Deg = 0.018837;
NT = length(full_fix_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

Fs = 3e4;
dsf = 30;Fsd = Fs/dsf;
scales = logspace(log10(3.2),log10(60),20);
scales = [scales 70 80 100 120 150];
% scales = [5:10 12 14 16 20 25 30 35 40 45 50 55 60];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
use_lfps = [1:8:96];
[b_gamma,a_gamma] = butter(2,[35 80]/(Fsd/2));

%%
n_fixs = length(full_fix_wends);
fix_expt_num = nan(n_fixs,1);
for i = 1:n_fixs
    cur_inds = full_fix_starts(i):full_fix_wends(i);
    fix_expt_num(i) = unique(full_expt_vec(cur_inds));
end
fix_durs = fix_end_times-fix_start_times;
%%
Expt_nu = [13 14 15 16 25 28 29];
n_expts = length(Expt_nu);
fixavg_lfp_amps = nan(length(fix_start_inds),length(wfreqs),length(use_lfps));
earlyavg_lfp_amps = nan(length(fix_start_inds),length(wfreqs),length(use_lfps));
lateavg_lfp_amps = nan(length(fix_start_inds),length(wfreqs),length(use_lfps));
fixavg_gam_amps = nan(length(fix_start_inds),length(use_lfps));
earlyavg_gam_amps = nan(length(fix_start_inds),length(use_lfps));
lateavg_gam_amps = nan(length(fix_start_inds),length(use_lfps));
fixavg_gam_freqs = nan(length(fix_start_inds),length(use_lfps));
earlyavg_gam_freqs = nan(length(fix_start_inds),length(use_lfps));
lateavg_gam_freqs = nan(length(fix_start_inds),length(use_lfps));

full_ampgrams = [];
full_gam_amp = [];
full_gam_freq = [];
for ee = 1:n_expts
    fprintf('Expt %d of %d\n',ee,n_expts);
    clear ampgram all_V* Vgamma_amp Vgamma_freq
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        ampgram(:,:,ll) = abs(cwt(V,scales,'cmor1-1'))';
        Vgamma = filtfilt(b_gamma,a_gamma,V);
        Vgamma_h = hilbert(Vgamma);
        Vgamma_amp(:,ll) = abs(Vgamma_h);
        Vgamma_phase = angle(Vgamma_h);
        Vgamma_phase = unwrap_phase_monotonic(Vgamma_phase);
        cVgamma_freq = smooth(Vgamma_phase,30,'lowess')';
        Vgamma_freq(:,ll) = [nan diff(cVgamma_freq)]/(2*pi)*Fsd;
    end
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    fix_set = find(fix_expt_num == ee);
    for n = 1:length(fix_set)
        if fix_durs(fix_set(n)) > 0.25
            cur_inds = find(t_ax >= fix_start_times(fix_set(n)) & t_ax < fix_start_times(fix_set(n))+0.25);
            fixavg_lfp_amps(fix_set(n),:,:) = mean(ampgram(cur_inds,:,:));
            fixavg_gam_amps(fix_set(n),:) = mean(Vgamma_amp(cur_inds,:));
            fixavg_gam_freqs(fix_set(n),:) = mean(Vgamma_freq(cur_inds,:));
        end
        cur_early_inds = find(t_ax >= fix_start_times(fix_set(n)) & t_ax < (fix_start_times(fix_set(n)) + 0.15));
        earlyavg_lfp_amps(fix_set(n),:,:) = mean(ampgram(cur_early_inds,:,:));
        earlyavg_gam_amps(fix_set(n),:) = mean(Vgamma_amp(cur_early_inds,:));
        earlyavg_gam_freqs(fix_set(n),:) = mean(Vgamma_freq(cur_early_inds,:));
        cur_late_inds = find(t_ax >= (0.15 + fix_start_times(fix_set(n))) & t_ax < (fix_end_times(fix_set(n))));
        if ~isempty(cur_late_inds)
            lateavg_lfp_amps(fix_set(n),:,:) = mean(ampgram(cur_late_inds,:,:),1);
            lateavg_gam_amps(fix_set(n),:) = mean(Vgamma_amp(cur_late_inds,:));
            lateavg_gam_freqs(fix_set(n),:) = mean(Vgamma_freq(cur_late_inds,:));
        end
    end
    fix_inds = find(full_expt_vec == ee);
    ufix_set = fix_inds(~isnan(full_fix_ids(fix_inds)));
    cur_interp_amps = interp1(t_ax,ampgram,full_t(ufix_set));
    cur_interp_gam_amp = interp1(t_ax,Vgamma_amp,full_t(ufix_set));
    cur_interp_gam_freq = interp1(t_ax,Vgamma_freq,full_t(ufix_set));
    full_ampgrams = [full_ampgrams; cur_interp_amps];
    full_gam_amp = [full_gam_amp; cur_interp_gam_amp];
    full_gam_freq = [full_gam_freq; cur_interp_gam_freq];
end

fixavg_lfp_amps = nanzscore(fixavg_lfp_amps);
earlyavg_lfp_amps = nanzscore(earlyavg_lfp_amps);
lateavg_lfp_amps = nanzscore(lateavg_lfp_amps);
fixavg_gam_amps = nanzscore(fixavg_gam_amps);
earlyavg_gam_amps = nanzscore(earlyavg_gam_amps);
lateavg_gam_amps = nanzscore(lateavg_gam_amps);
full_ampgrams = nanzscore(full_ampgrams);
fixavg_gam_freqs = nanzscore(fixavg_gam_freqs);
earlyavg_gam_freqs = nanzscore(earlyavg_gam_freqs);
lateavg_gam_freqs = nanzscore(lateavg_gam_freqs);
full_ampgrams = nanzscore(full_ampgrams);
%% RECONSTRUCT MAP STIMULUS
new_RF_patch = [-0.3 1.2; -1.2 0.3]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);

xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
% fullX_cropped = fullX(:,new_crop);

sdim = length(xpatch_inds);
[XXc,YYc] = meshgrid(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped));
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

load ./fixation_based_corrections_34
resh_X = reshape(resh_all_stims',[sdim sdim n_fixs]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:n_fixs
    d2 = dist_shift2d(resh_X(:,:,ii), -x_cor{end}(ii), 2,0);
    d2 = dist_shift2d(d2,-y_cor{end}(ii),1,0);
    resh_X_sh(:,:,ii) = d2;
end
fullX_sh = reshape(resh_X_sh,sdim^2,n_fixs)';
fullX_sh_cropped = fullX_sh(:,new_crop);
% fullX = reshape(resh_X,sdim^2,n_fixs)';
% fullX_cropped = fullX(:,new_crop);
clear fullX_sh
%% INTERPOLATE RF CENTERS OF LFP TUNING BASED ON UNIT GABOR FITS
% load ./gabor_initfits_d1p25_varrate
load ./gabor_tracking_varmeans
gabor_params_used = gabor_params_f{2};
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

%fit smoothed retinotopic surface
orientations = linspace(0,pi-pi/12,12);
for i = 1:96
    tempx(Y_pos(i),X_pos(i)) = gabor_params_used(i,1);
    tempy(Y_pos(i),X_pos(i)) = gabor_params_used(i,2);
    tempinds(Y_pos(i),X_pos(i)) = i;
end
weights = ones(10,10);
weights(tempx==0) = 0;
xpos_interp = smoothn(tempx,weights,'robust');
ypos_interp = smoothn(tempy,weights,'robust');
used_inds = find(weights == 1);
tempinds = tempinds(used_inds);
interp_x(tempinds) = xpos_interp(used_inds);
interp_y(tempinds) = ypos_interp(used_inds);
gabor_params_used(:,1) = interp_x;
gabor_params_used(:,2) = interp_y;

%% FIT FFT MODEL PARAMS
load ./expt2_lfp_amp_models_full_v2 lfp_gabor_params
all_lfp_params = nan(96,length(wfreqs),9);
all_lfp_params(1:2:end,:,:) = lfp_gabor_params;

clear fft_beta*
width = 0.25;
spatial_scale = 0.2;
sparse_lambda = 200;
smooth_lambda = 20000;
gam_freq = 13;
for cc = 1:length(use_lfps)
    fprintf('Ch %d of %d\n',cc,length(use_lfps));
%     for ww = 1:length(wfreqs)
%         cur_params = squeeze(lfp_gabor_params(cc,ww,1:6));
%         cur_params(5) = spatial_scale;
%         uset = (abs(XXc-cur_params(1)) <= 2*cur_params(5) & ...
%             abs(YYc-cur_params(2)) <= 2*cur_params(5));
%         xdims = max(uset);
%         ydims = max(uset,[],2);
%         if sum(xdims) > sum(ydims)
%             last_col = find(xdims == 1,1,'first');
%             uset(:,last_col) = 0;
%         end
%         if sum(ydims) > sum(xdims)
%             last_col = find(ydims == 1,1,'first');
%             uset(last_col,:) = 0;
%         end
%         uset = find(uset==1);
%         
%         gauss_mask = get_pgauss_mask_v2(XXc,YYc,cur_params);
%         gauss_mask = gauss_mask(uset);
%         X_out = bsxfun(@times,fullX_sh_cropped(:,uset),gauss_mask');
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
%         init_beta = zeros(n_dim,1);
%         yobs = squeeze(fixavg_lfp_amps(:,ww,cc));
%         use_pts = find(~isnan(yobs));
%         [fft_beta{cc,ww},best_LL(cc,ww)] = smooth_regress_2d(yobs(use_pts),pred_mat(use_pts,:),init_beta,smooth_lambda,sparse_lambda);
%                 
%         init_beta = zeros(n_dim,1);
%         yobs = squeeze(earlyavg_lfp_amps(:,ww,cc));
%         use_pts = find(~isnan(yobs));
%         [fft_beta_early{cc,ww},best_LL_early(cc,ww)] = smooth_regress_2d(yobs(use_pts),pred_mat(use_pts,:),init_beta,smooth_lambda,sparse_lambda);
% 
%             init_beta = zeros(n_dim,1);
%         yobs = squeeze(lateavg_lfp_amps(:,ww,cc));
%         use_pts = find(~isnan(yobs));
%         [fft_beta_late{cc,ww},best_LL_late(cc,ww)] = smooth_regress_2d(yobs(use_pts),pred_mat(use_pts,:),init_beta,smooth_lambda,sparse_lambda);
% end
    
    cur_params = squeeze(all_lfp_params(use_lfps(cc),gam_freq,1:6));
    cur_params(5) = spatial_scale;
    uset = (abs(XXc-cur_params(1)) <= 2*width & ...
        abs(YYc-cur_params(2)) <= 2*width);
    xdims = max(uset);
    ydims = max(uset,[],2);
    if sum(xdims) > sum(ydims)
        last_col = find(xdims == 1,1,'first');
        uset(:,last_col) = 0;
    end
    if sum(ydims) > sum(xdims)
        last_col = find(ydims == 1,1,'first');
        uset(last_col,:) = 0;
    end
    uset = find(uset==1);
    
    gauss_mask = get_pgauss_mask_v2(XXc,YYc,cur_params);
    gauss_mask = gauss_mask(uset);
    X_out = bsxfun(@times,fullX_sh_cropped(:,uset),gauss_mask');
    n_dim = length(uset);
    per_dim = sqrt(n_dim);
    X_out = reshape(X_out,[size(X_out,1) per_dim per_dim]);
    fft_out = nan(size(X_out));
    for i = 1:size(X_out,1)
        fft_out(i,:,:) = abs(fftshift(fft2(squeeze(X_out(i,:,:)))));
    end
    pred_mat = reshape(fft_out,size(X_out,1),per_dim^2);
    pred_mat = zscore(pred_mat);
    
    init_beta = zeros(n_dim,1);
    yobs = squeeze(fixavg_gam_amps(:,cc));
    use_pts = find(~isnan(yobs));
    [fft_beta_amp{cc},best_LL_amp(cc)] = smooth_regress_2d(yobs(use_pts),pred_mat(use_pts,:),init_beta,smooth_lambda,sparse_lambda);
    
    yobs = squeeze(fixavg_gam_freqs(:,cc));
    use_pts = find(~isnan(yobs));
    [fft_beta_freq{cc},best_LL_freq(cc)] = smooth_regress_2d(yobs(use_pts),pred_mat(use_pts,:),init_beta,smooth_lambda,sparse_lambda);
    
    yobs = squeeze(lateavg_gam_amps(:,cc));
    use_pts = find(~isnan(yobs));
    [fft_beta_lamp{cc},best_LL_lamp(cc)] = smooth_regress_2d(yobs(use_pts),pred_mat(use_pts,:),init_beta,smooth_lambda,sparse_lambda);
    yobs = squeeze(earlyavg_gam_amps(:,cc));
    use_pts = find(~isnan(yobs));
    [fft_beta_eamp{cc},best_LL_eamp(cc)] = smooth_regress_2d(yobs(use_pts),pred_mat(use_pts,:),init_beta,smooth_lambda,sparse_lambda);
    
    yobs = squeeze(fixavg_gam_freqs(:,cc));
    use_pts = find(~isnan(yobs));
    [fft_beta_freq{cc},best_LL_freq(cc)] = smooth_regress_2d(yobs(use_pts),pred_mat(use_pts,:),init_beta,smooth_lambda,sparse_lambda);
    yobs = squeeze(lateavg_gam_freqs(:,cc));
    use_pts = find(~isnan(yobs));
    [fft_beta_lfreq{cc},best_LL_lfreq(cc)] = smooth_regress_2d(yobs(use_pts),pred_mat(use_pts,:),init_beta,smooth_lambda,sparse_lambda);
    yobs = squeeze(earlyavg_gam_freqs(:,cc));
    use_pts = find(~isnan(yobs));
    [fft_beta_efreq{cc},best_LL_efreq(cc)] = smooth_regress_2d(yobs(use_pts),pred_mat(use_pts,:),init_beta,smooth_lambda,sparse_lambda);

end

%%
%%

for cc = 1:length(use_lfps)
    fprintf('Ch %d of %d\n',cc,length(use_lfps));
%     for ww = 1:length(wfreqs);
%         figure(1)
%         subplot(5,5,ww)
%         cur_pdim = sqrt(length(fft_beta_early{cc,ww}));
%         imagesc(reshape(fft_beta_early{cc,ww},cur_pdim,cur_pdim))
%         cax = caxis();
%         cm = max(abs(cax));
%         caxis([-cm cm])
%         title(sprintf('F: %.1f LL: %.1f',wfreqs(ww),best_LL(cc,ww)));
%         figure(2)
%         subplot(5,5,ww)
%         cur_pdim = sqrt(length(fft_beta_late{cc,ww}));
%         imagesc(reshape(fft_beta_late{cc,ww},cur_pdim,cur_pdim))
%         cax = caxis();
%         cm = max(abs(cax));
%         caxis([-cm cm])
%     end
%         figure(3);clf;
%         plot(wfreqs,best_LL(cc,:),'k')
%         hold on
%         plot(wfreqs,best_LL_early(cc,:))
%         plot(wfreqs,best_LL_late(cc,:),'r')
%         
%         pause
%         
        cur_pdim = sqrt(length(fft_beta_amp{cc}));
    figure(2)
    subplot(3,1,1)
        imagesc(reshape(fft_beta_amp{cc},cur_pdim,cur_pdim))
        cax = caxis();
        cm = max(abs(cax));
        caxis([-cm cm])
    subplot(3,1,2)
        imagesc(reshape(fft_beta_lamp{cc},cur_pdim,cur_pdim))
        cax = caxis();
        cm = max(abs(cax));
        caxis([-cm cm])
    subplot(3,1,3)
        imagesc(reshape(fft_beta_eamp{cc},cur_pdim,cur_pdim))
        cax = caxis();
        cm = max(abs(cax));
        caxis([-cm cm])
    
    figure(3)
    subplot(3,1,1)
        imagesc(reshape(fft_beta_freq{cc},cur_pdim,cur_pdim))
        cax = caxis();
        cm = max(abs(cax));
        caxis([-cm cm])
    subplot(3,1,2)
        imagesc(reshape(fft_beta_lfreq{cc},cur_pdim,cur_pdim))
        cax = caxis();
        cm = max(abs(cax));
        caxis([-cm cm])
    subplot(3,1,3)
        imagesc(reshape(fft_beta_efreq{cc},cur_pdim,cur_pdim))
        cax = caxis();
        cm = max(abs(cax));
        caxis([-cm cm])
    pause
%     figure(2);clf
%     figure(3);clf
% close all
end

