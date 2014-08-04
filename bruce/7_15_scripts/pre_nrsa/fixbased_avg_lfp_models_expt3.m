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
dsf = 100;Fsd = Fs/dsf;
scales = logspace(log10(10),log10(60),20);
scales = [scales 70 80 90 100 110 120 150 180 210 240];
% scales = logspace(log10(2.5),log10(40),25);
% scales = [5:10 12 14 16 20 25 30 35 40 45 50 55 60];
scales = scales*30/dsf;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
use_lfps = [1:2:96];

dir_fit = 1;
% load ./expt2_lfp_amp_models_full_v2
%%
n_trials = length(trial_start_inds);
trial_stop_sacinds = trial_stop_inds;
for i = 1:n_trials
    cur_inds = trial_start_inds(i):trial_stop_inds(i);
    cur_sac_stops = find(full_insac(cur_inds) == 1,1,'first');
    if ~isempty(cur_sac_stops)
        cur_sac_stops = cur_inds(cur_sac_stops);
        trial_stop_sacinds(i) = cur_sac_stops;
    end
end
trial_stop_inds = trial_stop_sacinds;
trial_start_times = full_t(trial_start_inds);
trial_stop_times = full_t(trial_stop_inds);
trial_durs = trial_stop_times - trial_start_times;

min_trial_dur = 0.4;
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
full_intrial = zeros(size(full_t));
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
    
    fix_set = find(fix_expt_num == Expt_nu(ee));
    for n = 1:length(fix_set)
        cur_inds = find(t_ax >= trial_start_times(fix_set(n)) & t_ax < trial_stop_wtimes(fix_set(n)));
        fixavg_lfp_amps(fix_set(n),:,:) = mean(ampgram(cur_inds,:,:));
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

%% INTERPOLATE RF CENTERS OF LFP TUNING BASED ON UNIT GABOR FITS
if dir_fit == 1
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
    
    %% estimate preferred orientations
    rsq = zeros(length(use_lfps),length(wfreqs),12);
    for cc = 1:length(use_lfps)
        fprintf('Ch %d of %d\n',cc,length(use_lfps));
        init_params(1:5) = gabor_params_used(use_lfps(cc),1:5);
        init_params(6) = 1;
        for i = 1:12
            init_params(3) = orientations(i);
            cur_mask1 = get_pgabor_mask_v2(XXc,YYc,init_params,0);
            cur_mask2 = get_pgabor_mask_v2(XXc,YYc,init_params,pi/2);
            
            mask1_out = fullX_cropped*cur_mask1(:);
            mask2_out = fullX_cropped*cur_mask2(:);
            
            en_out(i,:) = sqrt(mask1_out.^2+mask2_out.^2);
        end
        en_out = bsxfun(@minus,en_out,mean(en_out,2));
        en_out = bsxfun(@rdivide,en_out,std(en_out,[],2));
        
        for ww = 1:length(wfreqs)
            %     for ww = 9
            for i = 1:12
                temp = [en_out(i,:)' ones(n_trials,1)];
                [B,BINT,R,RINT,STATS] = regress(squeeze(fixavg_lfp_amps(:,ww,cc)),temp);
                rsq(cc,ww,i) = STATS(1);
            end
        end
    end
    
    %% FIT LFP MODEL PARAMS
    [max_rsq,maxloc] = max(rsq,[],3);
    
    hold_const = [1 1 0 0 0 1 0 0];
    LB = [-0.1 -0.8 0 0.125 0.025 0.2 0 -Inf];
    UB = [0.8 0.1 pi 0.5 0.4 6 Inf Inf];
    lfp_gabor_params = zeros(length(use_lfps),length(wfreqs),9);
    LL = zeros(length(use_lfps),length(wfreqs));
    clear init_params
    for cc = 1:length(use_lfps)
        fprintf('CC %d of %d\n',cc,length(use_lfps));
        for ww = 1:length(wfreqs)
            fprintf('freq %d of %d\n',ww,length(wfreqs));
            init_params(1:5) = gabor_params_used(use_lfps(cc),1:5);
            init_params(6) = 1;
            init_params(3) = orientations(maxloc(cc,ww));
            cur_init_params = [init_params zeros(1,3)];
            [lfp_gabor_params(cc,ww,:),LL(cc,ww)] = fit_gabor_params_lfp_v2(XXc,YYc,cur_init_params,fullX_cropped...
                ,squeeze(fixavg_lfp_amps(:,ww,cc)),hold_const,LB,UB);
            LL(cc,ww) = LL(cc,ww)/n_trials;
        end
    end
end
%%
clear fft_pred_out
spatial_scale = 0.2;
sparse_lambda = 350;
smooth_lambda = 3500;

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
        
        if dir_fit == 1
            init_beta = zeros(n_dim,1);
            [fft_beta{cc,ww},best_LL(cc,ww)] = smooth_regress_2d(squeeze(fixavg_lfp_amps(:,ww,cc)),pred_mat,init_beta,smooth_lambda,sparse_lambda);
            best_LL(cc,ww) = best_LL(cc,ww)/n_trials;
        end
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
        
%         gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,lfp_gabor_params(cc,ww,1:6),0);
%         gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,lfp_gabor_params(cc,ww,1:6),pi/2);
%         all_gabor_out1 = fullX_cropped*gabor_emp1(:);
%         all_gabor_out2 = fullX_cropped*gabor_emp2(:);
%         spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);
%         cur_spatial_mod_out = nan(tot_n_stims,1);
%         cur_spatial_mod_out(trial_stimnum) = spatial_mod_out;
%         cur_spatial_mod_out = zscore(cur_spatial_mod_out(full_stim_ids(uinds)));
%         
%         stim_dep_X = bsxfun(@times,sac_Tmat(uinds,:),cur_spatial_mod_out);
%         ov_Xmat = [stim_dep_X sac_Tmat(uinds,:)];
%         yobs = full_ampgrams(uinds,ww,cc);
%         init_params = randn(2*ntents,1);
%         [B_stimmod_gabor(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[50e3 2e3],[1 ntents; ntents+1 2*ntents],0);
        
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
cd ~/Data/bruce/7_15_12/G034/
save expt3_lfp_amp_models_dirfit_v3 B_stimmod* wfreqs tent_centers dt lfp_gabor_params LL rsq* orientations fft_beta best_LL

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
for i = 1:length(use_lfps)
    pcolor(tent_centers*dt,wfreqs,squeeze(B_stimmod_fft(i,:,1:ntents)));shading flat
    set(gca,'yscale','log')
    caxis([-0.02 0.25])
        
    pause
    clf
end
%%
close all
load ./expt3_lfp_amp_models_dirfit.mat
e3_B_stimod = B_stimmod_fft;
e3_tent_centers = tent_centers*dt;
e3_ntents = length(e3_tent_centers);
load ./expt2_lfp_amp_models_full_v2.mat
ntents = length(tent_centers);
B_stimmod = B_stimmod_gabor;
for i = 1:length(use_lfps)
    subplot(2,2,1)
    pcolor(tent_centers,wfreqs,squeeze(B_stimmod_fft(i,:,1:ntents)));shading flat
    set(gca,'yscale','log')
    caxis([-0.02 0.1])
    
    cax = caxis();
    colorbar
    subplot(2,2,2)
    pcolor(e3_tent_centers,wfreqs,squeeze(e3_B_stimod(i,:,1:e3_ntents)));shading flat
    set(gca,'yscale','log')
    caxis([-0.02 0.25])
    colorbar
    
    subplot(2,2,3)
    pcolor(tent_centers,wfreqs,squeeze(B_stimmod_fft(i,:,ntents+1:end)));shading flat
    set(gca,'yscale','log')
    caxis([-0.2 0.5])
    colorbar
    
    subplot(2,2,4)
    pcolor(e3_tent_centers,wfreqs,squeeze(e3_B_stimod(i,:,e3_ntents+1:end)));shading flat
    set(gca,'yscale','log')
    caxis([-0.2 0.5])
    colorbar
    pause
    clf
end

%%
load ./expt3_lfp_amp_models_dirfit.mat
dirfit_beta = fft_beta;
dirfit_LLfft = best_LL;
dirfit_LL = LL;
load ./expt2_lfp_amp_models_full_v2.mat
for cc = 1:length(use_lfps)
    fprintf('Ch %d of %d\n',cc,length(use_lfps));
    for ww = 1:length(wfreqs);
        figure(1)
        subplot(5,5,ww)
        cur_pdim = sqrt(length(dirfit_beta{cc,ww}));
        imagesc(reshape(dirfit_beta{cc,ww},cur_pdim,cur_pdim))
        cax = caxis();
        cm = max(abs(cax));
        caxis([-cm cm])
        title(sprintf('F: %.1f LL: %.1f',wfreqs(ww),best_LL(cc,ww)));
        figure(2)
        subplot(5,5,ww)
        cur_pdim = sqrt(length(fft_beta{cc,ww}));
        imagesc(reshape(fft_beta{cc,ww},cur_pdim,cur_pdim))
        cax = caxis();
        cm = max(abs(cax));
        caxis([-cm cm])
        title(sprintf('F: %.1f LL: %.1f',wfreqs(ww),best_LL(cc,ww)));
    end
    %     figure(2)
    %     hold on
    %     plot(wfreqs,best_LL(cc,:),'g')
    %     plot(wfreqs,LL(cc,:),'k')
    pause
    figure(1)
    clf
    figure(2)
    clf
end

%%
% cur_freq = 13;
% close all
% load ./expt3_lfp_amp_models_v3
% e3_B_stimod = B_stimmod_gabor;
% e3_tent_centers = tent_centers*dt;
% e3_ntents = length(e3_tent_centers);
% load ./expt2_lfp_amp_models_full.mat
% ntents = length(tent_centers);
% B_stimmod = B_stimmod_gabor;
% for i = 1:length(use_lfps)
%     subplot(2,1,1)
%     plot(tent_centers,squeeze(B_stimmod(i,cur_freq,1:ntents)))
%     hold on
%     plot(e3_tent_centers,squeeze(e3_B_stimod(i,cur_freq,1:e3_ntents)),'r')
%     xlim([0 0.5])
%     subplot(2,1,2)
%     plot(tent_centers,squeeze(B_stimmod(i,cur_freq,ntents+1:end)))
%     hold on
%     plot(e3_tent_centers,squeeze(e3_B_stimod(i,cur_freq,e3_ntents+1:end)),'r')
%     xlim([0 0.5])
% pause
% clf
% end