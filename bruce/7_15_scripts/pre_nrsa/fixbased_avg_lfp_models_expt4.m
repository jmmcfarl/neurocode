%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G035/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt4_fixbased_data.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

fix_win_dur = 0.15;

Pix2Deg = 0.018837;
NT = length(full_stim_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

Fs = 3e4;
dsf = 30;Fsd = Fs/dsf;
% scales = logspace(log10(2.5),log10(40),25);
scales = logspace(log10(3.2),log10(60),20);
scales = [scales 70 80 100 120 150];
% scales = [5:10 12 14 16 20 25 30 35 40 45 50 55 60];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
use_lfps = [1:2:96];

load ~/Data/bruce/7_15_12/G034/expt2_lfp_amp_models_full_v2.mat

%%
n_trials = length(trial_start_inds);
fix_expt_num = nan(n_trials,1);
for i = 1:n_trials
    cur_inds = trial_start_inds(i):trial_stop_winds(i);
    fix_expt_num(i) = unique(full_expt_vec(cur_inds));
end
trial_start_times = full_t(trial_start_inds);
trial_stop_wtimes = full_t(trial_stop_winds);
trial_stop_times = full_t(trial_stop_inds);

%%
% Expt_nu = [13 14 15 16 25 28 29];
% Expt_nu = [7:12]; %expt 3 34
Expt_nu = [1:5 8]; %expt 4 35
n_expts = length(Expt_nu);
fixavg_lfp_amps = nan(length(trial_start_inds),length(wfreqs),length(use_lfps));
earlyavg_lfp_amps = nan(length(trial_start_inds),length(wfreqs),length(use_lfps));
lateavg_lfp_amps = nan(length(trial_start_inds),length(wfreqs),length(use_lfps));

full_ampgrams = [];
for ee = 1:n_expts
    fprintf('Expt %d of %d\n',ee,n_expts);
    clear ampgram all_V*
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        ampgram(:,:,ll) = abs(cwt(V,scales,'cmor1-1'))';
    end
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    fix_set = find(fix_expt_num == ee);
    for n = 1:length(fix_set)
        cur_inds = find(t_ax >= trial_start_times(fix_set(n)) & t_ax < trial_stop_times(fix_set(n)));
        fixavg_lfp_amps(fix_set(n),:,:) = mean(ampgram(cur_inds,:,:));
%         cur_early_inds = find(t_ax >= trial_start_times(fix_set(n)) & t_ax < (trial_start_times(fix_set(n)) + 0.15));
%         earlyavg_lfp_amps(fix_set(n),:,:) = mean(ampgram(cur_early_inds,:,:));
%         cur_late_inds = find(t_ax >= (0.15 + trial_start_times(fix_set(n))) & t_ax < (trial_stop_times(fix_set(n))));
%         if ~isempty(cur_late_inds)
%             lateavg_lfp_amps(fix_set(n),:,:) = mean(ampgram(cur_late_inds,:,:),1);
%         end
    end
    fix_inds = find(full_expt_vec == Expt_nu(ee));
    cur_interp_amps = interp1(t_ax,ampgram,full_t(fix_inds));
    full_ampgrams = [full_ampgrams; cur_interp_amps];
end

fixavg_lfp_amps = nanzscore(fixavg_lfp_amps);
% earlyavg_lfp_amps = nanzscore(earlyavg_lfp_amps);
% lateavg_lfp_amps = nanzscore(lateavg_lfp_amps);
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

fullX_cropped = resh_all_stims(trial_stimnum,new_crop);

%%
clear fft_pred_out
spatial_scale = 0.2;
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
        
%         sparse_lambda = 500;
%         smooth_lambda = 5000;
%         init_beta = zeros(n_dim,1);
%         [beta{cc,ww},best_LL(cc,ww)] = smooth_regress_2d(squeeze(fixavg_lfp_amps(:,ww,cc)),pred_mat,init_beta,smooth_lambda,sparse_lambda);
        fft_pred_out(cc,ww,:) = pred_mat*squeeze(fft_beta{cc,ww});
        
    end
end

%%
flen = 80;
flen_t = flen*dt;
NT = length(full_expt_vec);

tent_centers = [0:dt:0.15];
cur_sp = dt;
while max(tent_centers) < flen_t
    tent_centers = [tent_centers (tent_centers(end)+cur_sp)];
    cur_sp = cur_sp + dt/6;
end
tent_centers = round(tent_centers/dt);
if tent_centers(end) >= flen
    tent_centers(end) = [];
end
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

tbmat = [zeros(ntents,tblen-1) tbmat];

sac_inds = zeros(size(full_t));
sac_inds(trial_start_inds) = 1;

sac_Tmat = zeros(length(full_t),ntents);
for i = 1:ntents
    sac_Tmat(:,i) = conv(sac_inds,tbmat(i,:),'same');
end

%%
tot_n_stims = length(unique(full_stim_ids));
fft_out_all = nan(size(fft_pred_out,1),size(fft_pred_out,2),tot_n_stims);
fft_out_all(:,:,trial_stimnum) = fft_pred_out;
clear B_stimmod*
uinds = find(ismember(full_stim_ids,trial_stimnum));
for cc = 1:length(use_lfps)
    fprintf('Fitting Cell %d of %d\n',cc,length(use_lfps));
    for ww = 1:length(wfreqs)
        fprintf('Frequency %d of %d\n',ww,length(wfreqs));
        
        gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,lfp_gabor_params(cc,ww,1:6),0);
        gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,lfp_gabor_params(cc,ww,1:6),pi/2);
        all_gabor_out1 = fullX_cropped*gabor_emp1(:);
        all_gabor_out2 = fullX_cropped*gabor_emp2(:);
        spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);
        cur_spatial_mod_out = nan(tot_n_stims,1);
        cur_spatial_mod_out(trial_stimnum) = spatial_mod_out;
        cur_spatial_mod_out = zscore(cur_spatial_mod_out(full_stim_ids(uinds)));
 
        stim_dep_X = bsxfun(@times,sac_Tmat(uinds,:),cur_spatial_mod_out);
        ov_Xmat = [stim_dep_X sac_Tmat(uinds,:)];
        yobs = full_ampgrams(uinds,ww,cc);
        init_params = randn(2*ntents,1);
        [B_stimmod_gabor(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[50e3 2e3],[1 ntents; ntents+1 2*ntents],0);

        total_out = zscore(squeeze(fft_out_all(cc,ww,full_stim_ids(uinds)))');
        stim_dep_X = bsxfun(@times,sac_Tmat(uinds,:),total_out');
        ov_Xmat = [stim_dep_X sac_Tmat(uinds,:)];
        yobs = full_ampgrams(uinds,ww,cc);
        init_params = randn(2*ntents,1);
        [B_stimmod_fft(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[50e3 2e3],[1 ntents; ntents+1 2*ntents],0);
        
        %         cur_spatial_mod_out = zscore(spatial_mod_out(full_fix_ids(uinds(nat_inds))));
        %         stim_dep_X = bsxfun(@times,sac_Tmat(uinds(nat_inds),:),cur_spatial_mod_out);
        %         ov_Xmat = [stim_dep_X sac_Tmat(uinds(nat_inds),:)];
        %         yobs = full_ampgrams(nat_inds,ww,cc);
        %         init_params = randn(2*ntents,1);
        %         [B_stimmod_nat(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[20e3 2e3],[1 ntents; ntents+1 2*ntents],0);
        %
        %         cur_spatial_mod_out = zscore(spatial_mod_out(full_fix_ids(uinds(noise_inds))));
        %         stim_dep_X = bsxfun(@times,sac_Tmat(uinds(noise_inds),:),cur_spatial_mod_out);
        %         ov_Xmat = [stim_dep_X sac_Tmat(uinds(noise_inds),:)];
        %         yobs = full_ampgrams(noise_inds,ww,cc);
        %         init_params = randn(2*ntents,1);
        %         [B_stimmod_noise(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[20e3 2e3],[1 ntents; ntents+1 2*ntents],0);
    end
end

%%
cd ~/Data/bruce/7_15_12/G035/
save expt4_lfp_amp_models_v4 B_stimmod* wfreqs tent_centers dt LL 

%%
% cur_freq = 20;
% figure
% cnt = 1;
% for i = 1:10
%     for j = 1:10
%         cur = find(X_pos == i & Y_pos == j);
%         if ~isempty(cur)
%             subplot(10,10,cnt)
%             hold on
%             plot(tent_centers,squeeze(B_stimmod(cur,cur_freq,1:ntents)),'b');
% %             hold on
% %             plot(tent_centers,squeeze(B_stimmod(cur,cur_freq,ntents+1:end)),'r');
%                      
%             axis tight;
%             xl = xlim();
%             line(xl,[0 0],'color','k')
%         end
%         cnt = cnt + 1;
%     end
% end

%%
close all
cd ~/Data/bruce/7_15_12/G035/
load ./expt4_lfp_amp_models_v4 
e4_B_stimod = B_stimmod_fft;
e4_tent_centers = tent_centers*dt;
e4_ntents = length(e4_tent_centers);
cd ~/Data/bruce/7_15_12/G034/
load ./expt3_lfp_amp_models_v4
e3_B_stimod = B_stimmod_fft;
e3_tent_centers = tent_centers*dt;
e3_ntents = length(e3_tent_centers);
load ./expt2_lfp_amp_models_full_v2.mat
ntents = length(tent_centers);
B_stimmod = B_stimmod_gabor;
for i = 1:length(use_lfps)
    subplot(2,3,1)
    pcolor(tent_centers,wfreqs,squeeze(B_stimmod_fft(i,:,1:ntents)));shading flat
    set(gca,'yscale','log')
        caxis([-0.02 0.075])

    cax = caxis();
    colorbar
    subplot(2,3,2)
    pcolor(e3_tent_centers,wfreqs,squeeze(e3_B_stimod(i,:,1:e3_ntents)));shading flat
    set(gca,'yscale','log')
        caxis([-0.02 0.2])
    colorbar
    subplot(2,3,3)
    pcolor(e4_tent_centers,wfreqs,squeeze(e4_B_stimod(i,:,1:e4_ntents)));shading flat
    set(gca,'yscale','log')
%         caxis([-0.02 0.15])
    colorbar

        subplot(2,3,4)
    pcolor(tent_centers,wfreqs,squeeze(B_stimmod_fft(i,:,ntents+1:end)));shading flat
    set(gca,'yscale','log')
%     caxis([-0.2 0.4])
    colorbar
    
    subplot(2,3,5)
    pcolor(e3_tent_centers,wfreqs,squeeze(e3_B_stimod(i,:,e3_ntents+1:end)));shading flat
    set(gca,'yscale','log')
    colorbar
    subplot(2,3,6)
    pcolor(e4_tent_centers,wfreqs,squeeze(e4_B_stimod(i,:,e4_ntents+1:end)));shading flat
    set(gca,'yscale','log')
    colorbar
pause
clf
end
