%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt3_fixbased_data.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

fix_win_dur = 0.15;

Pix2Deg = 0.018837;
NT = length(full_stim_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

Fs = 3e4;
dsf = 30;Fsd = Fs/dsf;
[b_alpha,a_alpha] = butter(2,[4 12]/(Fsd/2));
scales = logspace(log10(3.2),log10(60),20);
scales = [scales 70 80 100 120 150];
% scales = logspace(log10(2.5),log10(40),25);
% scales = [5:10 12 14 16 20 25 30 35 40 45 50 55 60];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
use_lfps = [1:2:96];

load ./expt2_lfp_amp_models_full_v2
%%
trial_start_times = full_t(trial_start_inds);
trial_stop_times = full_t(trial_stop_inds);
trial_durs = trial_stop_times - trial_start_times;
min_trial_dur = 0.25;
used_trials = find(trial_durs >= min_trial_dur);
trial_start_inds = trial_start_inds(used_trials);
trial_stop_inds = trial_stop_inds(used_trials);
trial_imnum = trial_imnum(used_trials);
trial_stimnum = trial_stimnum(used_trials);
trial_start_times = trial_start_times(used_trials);
trial_stop_times = trial_stop_times(used_trials);
trial_stop_winds = trial_start_inds + round(min_trial_dur/dt);
trial_stop_wtimes = full_t(trial_stop_winds);

n_trials = length(trial_start_inds);
fix_expt_num = nan(n_trials,1);
for i = 1:n_trials
    cur_inds = trial_start_inds(i):trial_stop_winds(i);
    fix_expt_num(i) = unique(full_expt_vec(cur_inds));
end

%%
% Expt_nu = [13 14 15 16 25 28 29];
Expt_nu = [7:12]; %expt 3 34
n_expts = length(Expt_nu);
fixavg_lfp_amps = nan(length(trial_start_inds),length(wfreqs),length(use_lfps));

full_ampgrams = [];
full_alpha_phase = [];
for ee = 1:n_expts
    fprintf('Expt %d of %d\n',ee,n_expts);
    clear ampgram all_V* alpha_phase
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        V_alpha = filtfilt(b_alpha,a_alpha,V);
        V_alpha_phase = angle(hilbert(V_alpha));
        V_alpha_phase = unwrap_phase_monotonic(V_alpha_phase);
        alpha_phase(:,ll) = V_alpha_phase;
        ampgram(:,:,ll) = abs(cwt(V,scales,'cmor1-1'))';
    end
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    fix_inds = find(full_expt_vec == Expt_nu(ee));
    cur_interp_amps = interp1(t_ax,ampgram,full_t(fix_inds));
    cur_interp_phase = interp1(t_ax,alpha_phase,full_t(fix_inds));
    full_ampgrams = [full_ampgrams; cur_interp_amps];
    full_alpha_phase = [full_alpha_phase; cur_interp_phase];
end

full_ampgrams = nanzscore(full_ampgrams);

%%
phase_since_trial = nan(size(full_alpha_phase));
for i = 1:n_trials
   cur_set = trial_start_inds(i):trial_stop_inds(i);
   phase_since_trial(cur_set,:) = bsxfun(@minus,full_alpha_phase(cur_set,:),full_alpha_phase(cur_set(1),:));
end
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

fullX_cropped = resh_all_stims(trial_stimnum,new_crop);


%%
clear fft_pred_out
spatial_scale = 0.2;
sparse_lambda = 500;
smooth_lambda = 10000;

for cc = 1:length(use_lfps)
    fprintf('Ch %d of %d\n',cc,length(use_lfps));
    for ww = 1:length(wfreqs)
        cur_params = squeeze(lfp_gabor_params(cc,ww,1:6));
        cur_params(5) = spatial_scale;
        uset = (abs(XXc-cur_params(1)) <= 2*cur_params(5) & ...
            abs(YYc-cur_params(2)) <= 2*cur_params(5));
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
        X_out = bsxfun(@times,fullX_cropped(:,uset),gauss_mask');
        n_dim = length(uset);
        per_dim = sqrt(n_dim);
        X_out = reshape(X_out,[size(X_out,1) per_dim per_dim]);
        fft_out = nan(size(X_out));
        for i = 1:size(X_out,1)
            fft_out(i,:,:) = abs(fftshift(fft2(squeeze(X_out(i,:,:)))));
        end
        pred_mat = reshape(fft_out,size(X_out,1),per_dim^2);
        pred_mat = zscore(pred_mat);
        
%         init_beta = zeros(n_dim,1);
%         [fft_beta{cc,ww},best_LL(cc,ww)] = smooth_regress_2d(squeeze(fixavg_lfp_amps(:,ww,cc)),pred_mat,init_beta,smooth_lambda,sparse_lambda);
        fft_pred_out(cc,ww,:) = pred_mat*squeeze(fft_beta{cc,ww});
    end
end

%%
NT = length(full_expt_vec);
max_phase = 2*pi*5;
nbins = 50;
tax = prctile(phase_since_trial(:,1),linspace(100/nbins,100-100/nbins,nbins));


%%
tot_n_stims = length(unique(full_stim_ids));
fft_out_all = nan(size(fft_pred_out,1),size(fft_pred_out,2),tot_n_stims);
fft_out_all(:,:,trial_stimnum) = fft_pred_out;
clear B_stimmod*
uinds = find(ismember(full_stim_ids,trial_stimnum));
for cc = 1:length(use_lfps)
    fprintf('Fitting Cell %d of %d\n',cc,length(use_lfps));
    
    Tmat = tbrep(phase_since_trial(:,cc),tax);

    for ww = 1:length(wfreqs)
        fprintf('Frequency %d of %d\n',ww,length(wfreqs));
        

        total_out = zscore(squeeze(fft_out_all(cc,ww,full_stim_ids(uinds)))');
        stim_dep_X = bsxfun(@times,Tmat(uinds,:),total_out');
        ov_Xmat = [stim_dep_X Tmat(uinds,:)];
        yobs = full_ampgrams(uinds,ww,cc);
        init_params = randn(2*nbins,1);
        [B_stimmod_fft(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[20e3 2e3],[1 nbins; nbins+1 2*nbins],0);
        
    end
end

%%
cd ~/Data/bruce/7_15_12/G034/
save expt3_lfp_amp_models_alphaphase B_stimmod* wfreqs tent_centers dt LL rsq* fft_beta best_LL

%%
for i = 1:length(use_lfps)
    subplot(2,1,1)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft(i,:,1:nbins)));shading flat
    set(gca,'yscale','log')
    caxis([-0.02 0.2])
    
    cax = caxis();
    colorbar
    subplot(2,1,2)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft(i,:,nbins+1:end)));shading flat
    set(gca,'yscale','log')
%     caxis([-0.02 0.25])
    colorbar
    
    pause
    clf
end