%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt3_fixbased_data.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
use_win = [-7 7;-7 7];

min_fix_dur = 0.15;
use_lfps = [1:4:96];

Fs = 3e4;
dsf = 60;Fsd = Fs/dsf;
[b,a] = butter(2,[1 100]/(Fsd/2));
[b2,a2] = butter(2,[4 12]/(Fsd/2));
[b3,a3] = butter(2,[30 70]/(Fsd/2));

forwardlag = round(Fsd*0.6);
backlag = round(Fsd*0.2);
lags = -backlag:forwardlag;

%%
% Expt_nu = [3 8 15 19 26 29];
% Expt_nu = [13 14 15 16 25 28 29];
Expt_nu = [7:12]; %expt 3 34
n_allunits = 96;
pooled_Vbb = [];
pooled_Valpha = [];
pooled_Vgamma = [];
for ee = 1:length(Expt_nu)
    
    %%
    clear ampgram all_V*
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        V_alpha = filtfilt(b2,a2,V);
        V_gamma = filtfilt(b3,a3,V);
        V_bb = filtfilt(b,a,V);
        all_Vbb(:,ll) = V_bb;
        all_Valpha(:,ll) = V_alpha;
        all_Vgamma(:,ll) = abs(hilbert(V_gamma));
    end
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    fix_inds = find(full_expt_vec == Expt_nu(ee));
    cur_interp_Vbb = interp1(t_ax,all_Vbb,full_t(fix_inds));
    cur_interp_Valpha = interp1(t_ax,all_Valpha,full_t(fix_inds));
    cur_interp_Vgamma = interp1(t_ax,all_Vgamma,full_t(fix_inds));
    
    pooled_Vbb = [pooled_Vbb; cur_interp_Vbb];
    pooled_Valpha = [pooled_Valpha; cur_interp_Valpha];
    pooled_Vgamma = [pooled_Vgamma; cur_interp_Vgamma];
    
end
pooled_Vbb = zscore(pooled_Vbb);
pooled_Valpha = zscore(pooled_Valpha);
pooled_Vgamma = zscore(pooled_Vgamma);

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
load ./expt2_lfp_amp_models_full_v2 lfp_gabor_params fft_beta wfreqs
lfp_gabor_params = lfp_gabor_params(1:2:end,:,:);
fft_beta = fft_beta(1:2:end,:);
gam_freq = 14;

clear fft_pred_out
spatial_scale = 0.2;
for cc = 1:length(use_lfps)
    fprintf('Ch %d of %d\n',cc,length(use_lfps));
    cur_params = squeeze(lfp_gabor_params(cc,gam_freq,1:6));
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
    
    fft_pred_out(cc,:) = pred_mat*squeeze(fft_beta{cc,gam_freq});
end

%%
n_trials = length(trial_start_inds);
trig_avg_Vbb = nan(n_trials,length(use_lfps),length(lags));
trig_avg_Valpha = nan(n_trials,length(use_lfps),length(lags));
trig_avg_Vgamma = nan(n_trials,length(use_lfps),length(lags));


