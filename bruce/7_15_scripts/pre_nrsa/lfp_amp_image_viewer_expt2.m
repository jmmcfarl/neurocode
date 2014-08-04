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
[b_alpha,a_alpha] = butter(2,[5 15]/(Fsd/2));
[b_gamma,a_gamma] = butter(2,[35 70]/(Fsd/2));
use_lfps = [1:96];

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
trial_stop_winds = trial_start_inds + round(0.15/dt);
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
fixavg_peakalpha = nan(length(trial_start_inds),length(use_lfps));
fixavg_peakloc = nan(length(trial_start_inds),length(use_lfps));
fixavg_avgalpha = nan(length(trial_start_inds),length(use_lfps));
fixavg_avggamma = nan(length(trial_start_inds),length(use_lfps));
full_alpha = [];
full_gamma = [];
full_alphaphase = [];
full_V = [];
for ee = 1
    fprintf('Expt %d of %d\n',ee,n_expts);
    clear all_alpha all_V all_alphaphase all_gamma
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        all_alpha(:,ll) = filtfilt(b_alpha,a_alpha,V);
        all_gamma(:,ll) = filtfilt(b_gamma,a_gamma,V);
        all_V(:,ll) = V;
        all_alphaphase(:,ll) = angle(hilbert(all_alpha(:,ll)));
    end
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    all_alpha_env = abs(hilbert(all_alpha));
    all_gamma_env = abs(hilbert(all_gamma));
    
    fix_set = find(fix_expt_num == Expt_nu(ee));
    for n = 1:length(fix_set)
        cur_inds = find(t_ax >= trial_start_times(fix_set(n)) & t_ax < trial_stop_wtimes(fix_set(n)));
        [fixavg_peakalpha(fix_set(n),:),fixavg_peakloc(fix_set(n),:)] = min(all_alpha(cur_inds,:));
        cur_inds = find(t_ax >= trial_start_times(fix_set(n)) & t_ax < trial_stop_times(fix_set(n)));
        fixavg_avgalpha(fix_set(n),:) = mean(all_alpha_env(cur_inds,:));
        cur_inds = find(t_ax >= trial_start_times(fix_set(n)) & t_ax < trial_stop_times(fix_set(n)));
        fixavg_avggamma(fix_set(n),:) = mean(all_gamma_env(cur_inds,:));
    end
    fix_inds = find(full_expt_vec == Expt_nu(ee));
    cur_interp_alpha = interp1(t_ax,all_alpha,full_t(fix_inds));
    cur_interp_alphaphase = interp1(t_ax,all_alphaphase,full_t(fix_inds));
    cur_interp_V = interp1(t_ax,all_V,full_t(fix_inds));
    full_alpha = [full_alpha; cur_interp_alpha];
    full_alphaphase = [full_alphaphase; cur_interp_alphaphase];
    full_V = [full_V; cur_interp_V];
end

fixavg_peakalpha = nanzscore(fixavg_peakalpha);
fixavg_peakloc = fixavg_peakloc/Fsd;
fixavg_avgalpha = nanzscore(fixavg_avgalpha);
fixavg_avggamma = nanzscore(fixavg_avggamma);
% earlyavg_lfp_amps = nanzscore(earlyavg_lfp_amps);
% lateavg_lfp_amps = nanzscore(lateavg_lfp_amps);
full_alpha = nanzscore(full_alpha);

%%
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

array_peakalpha = nan(size(fixavg_avgalpha,1),10,10);
% array_peakalpha = nan(size(fixavg_avgalpha,1),10,10);
for i = 1:96
    %     array_peakalpha(:,Y_pos(i),X_pos(i)) = fixavg_peakalpha(:,i);
    array_peakalpha(:,Y_pos(i),X_pos(i)) = fixavg_peakloc(:,i);
    %     array_peakalpha(:,Y_pos(i),X_pos(i)) = fixavg_avggamma(:,i);
end

xi = linspace(min(interp_x),max(interp_x),50);
yi = linspace(min(interp_y),max(interp_y),50);
[Xi,Yi] = meshgrid(xi,yi);

%%
new_RF_patch = [0.1 0.7; -0.7 -0.1]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);

xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
% fullX_cropped = fullX(:,new_crop);

sdim = length(xpatch_inds);
sdimc = length(xpatch_inds_cropped);
[XXc,YYc] = meshgrid(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped));
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

fullX_cropped = resh_all_stims(trial_stimnum,new_crop);
%%
uset = 1:96;
uset([9 10 11 16 35 92]) = [];
close all
rect_pos = [0.1 -0.7 0.6 0.6]
xx = xax(xpatch_inds);
yy = yax(ypatch_inds);
xxc = xax(xpatch_inds_cropped); 
yyc = yax(ypatch_inds_cropped);
for i = 1:n_trials
    cur_image = reshape(resh_all_stims(trial_stimnum(i),:),sdim,sdim);
    cur_image_c = reshape(fullX_cropped(i,:),sdimc,sdimc);
   
%     Vq = griddata(interp_x,interp_y,fixavg_peakalpha(i,:),Xi,Yi);
%     Vq = griddata(interp_x,interp_y,fixavg_avggamma(i,:),Xi,Yi);
%     Vq = griddata(interp_x,interp_y,fixavg_avgalpha(i,:),Xi,Yi);
    Vq = griddata(interp_x(uset),interp_y(uset),fixavg_peakloc(i,uset),Xi,Yi);

%     subplot(2,1,1)
%     imagesc(squeeze(array_peakalpha(i,:,:))); colorbar; %caxis([-1.5 1.5]); 
%     subplot(2,1,2)
%     imagesc(xi,yi,Vq); %caxis([-1 1]); colorbar; 
    pcolor(xi,yi,Vq); shading flat; colorbar; 
    
%     subplot(2,2,2)
%     imagesc(xx,yy,cur_image);colormap(jet)
%     rectangle('Position',rect_pos,'edgecolor','r')
%     subplot(2,2,4)
%     imagesc(xxc,yyc,cur_image_c);colormap(jet)
    pause
    clf
    
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
        
            %     for ww = 9
            for i = 1:12
                temp = [en_out(i,:)' ones(n_trials,1)];
                [B,BINT,R,RINT,STATS] = regress(squeeze(fixavg_peakalpha(:,cc)),temp);
                rsq(cc,i) = STATS(1);
            end
        end
    
    %% FIT LFP MODEL PARAMS
    [max_rsq,maxloc] = max(rsq,[],3);
    
    hold_const = [1 1 0 0 0 1 0 0];
    LB = [-0.1 -0.8 0 0.125 0.025 0.2 0 -Inf];
    UB = [0.8 0.1 pi 0.5 0.4 6 Inf Inf];
    lfp_gabor_params = zeros(length(use_lfps),9);
   LL = zeros(length(use_lfps),1);
    clear init_params
    for cc = 1:length(use_lfps)
        fprintf('CC %d of %d\n',cc,length(use_lfps));
            init_params(1:5) = gabor_params_used(use_lfps(cc),1:5);
            init_params(6) = 1;
            init_params(3) = orientations(maxloc(cc));
            cur_init_params = [init_params zeros(1,3)];
            [lfp_gabor_params(cc,:),LL(cc)] = fit_gabor_params_lfp_v2(XXc,YYc,cur_init_params,fullX_cropped...
                ,squeeze(fixavg_peakalpha(:,cc)),hold_const,LB,UB);
            LL(cc) = LL(cc)/n_trials;
    end
%%
clear fft_pred_out
spatial_scale = 0.2;
sparse_lambda = 100;
smooth_lambda = 2000;

for cc = 1:length(use_lfps)
    fprintf('Ch %d of %d\n',cc,length(use_lfps));
        cur_params = squeeze(lfp_gabor_params(cc,1:6));
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
        
            init_beta = zeros(n_dim,1);
            [fft_beta{cc},best_LL(cc)] = smooth_regress_2d(squeeze(fixavg_peakalpha(:,cc)),pred_mat,init_beta,smooth_lambda,sparse_lambda);
            best_LL(cc) = best_LL(cc)/n_trials;
        fft_pred_out(cc,:) = pred_mat*squeeze(fft_beta{cc});

            [fft_beta2{cc},best_LL2(cc)] = smooth_regress_2d(squeeze(fixavg_avgalpha(:,cc)),pred_mat,init_beta,smooth_lambda,sparse_lambda);
            best_LL2(cc) = best_LL2(cc)/n_trials;
        fft_pred_out(cc,:) = pred_mat*squeeze(fft_beta2{cc});

end


%%
%%
close all
for cc = 1:length(use_lfps)
    fprintf('Ch %d of %d\n',cc,length(use_lfps));
        cur_pdim = sqrt(length(fft_beta2{cc}));
        imagesc(reshape(fft_beta2{cc},cur_pdim,cur_pdim))
        cax = caxis();
        cm = max(abs(cax));
        caxis([-cm cm])
    pause
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