%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/
% load ./eye_calibration_data
% load ./G029Expts.mat
% cd G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

% load ./Expt1_compiled_data_fixedlag.mat
% load ./Expt1_compiled_data_fixedlag.mat
cd ~/Data/bruce/7_15_12/G029/
% load ./Expt1_compiled_data_fixedlag_driftcor
% load ./Expt1_compiled_data_fixedlag_driftcor_leftonly
% load ./Expt1_compiled_data_fixedlag_driftcor_leftonly
load ./Expt1_newcompiled_data_fixedlag_leftonly_d1p5.mat
fullX = fullX/std(fullX(:));

% fullX = bsxfun(@minus,fullX,mean(fullX));
% fullX = bsxfun(@rdivide,fullX,std(fullX));

Pix2Deg = 0.018837;
%%
new_RF_patch = [-0.1 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);
new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
fullX = fullX(:,new_crop);

xpatch_inds = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

%% PARSE DATA INTO FIXATIONS
[NT,klen] = size(fullX);
min_fix_dur = 0.1;
trial_flips = 1 + find(diff(full_trial_vec) ~= 0);
sac_inds = find(full_insac==1);
all_break_inds = unique([trial_flips; sac_inds]);
fix_start_inds = [1; all_break_inds+1];
fix_stop_inds = [all_break_inds-1; NT];
fix_durs = (fix_stop_inds-fix_start_inds)*dt;
used_fixs = find(fix_durs > min_fix_dur);
fprintf('Using %d of %d fixations\n',length(used_fixs),length(fix_durs));

used_ims = [];
for i = 1:length(used_fixs)
    used_ims = [used_ims fix_start_inds(used_fixs(i)):fix_stop_inds(used_fixs(i))];
end

%% REFIT GABOR PARAMS
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
cd ~/Data/bruce/7_15_12/
cd G029/
load ./grating_mu_data
gr_oris = unique_oris;
[~,sf_inds] = max(avg_mu_sf_profile,[],2);
pref_sfs = unique_sfs(sf_inds);
load ./oned_fixation_fits_v3.mat

nfold = 5;
nparts = 60;
NT = length(used_ims);
[tr_samp,xv_samp] = create_xv_set(NT,nfold,nparts);

tr_samp = used_ims(tr_samp);
xv_samp = used_ims(xv_samp);

%%
for t = 1:96;
    fprintf('SU %d of %d\n',t,96);
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/pref_sfs(t); %lambda
    init_params(5) = 0.4*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    init_params(7) = 0; %weight of quad term
    init_params(8) = 0; %const offset
    LB = [-0.1 -0.8 0 0.12 0.025 0.5 0 -Inf];
    UB = [0.8 0.1 pi 0.4 0.25 2 Inf Inf];

%     hold_const = [0 0 1 1 1 1 0 0];
%     [gabor_params_f(t,:),LL(t)] = fit_gabor_energy_mod(XX,YY,init_params,fullX(tr_samp,:),full_binned_spks(tr_samp,t),hold_const,LB,UB);
%     xv_LL_f(t) = get_gabor_energy_LL(gabor_params_f(t,:),fullX(xv_samp,:),full_binned_spks(xv_samp,t),XX,YY);
% 
%     hold_const = [0 0 1 1 0 1 0 0];
%     [gabor_params_f2(t,:),LL2(t)] = fit_gabor_energy_mod(XX,YY,init_params,fullX(tr_samp,:),full_binned_spks(tr_samp,t),hold_const,LB,UB);
%     xv_LL_f2(t) = get_gabor_energy_LL(gabor_params_f2(t,:),fullX(xv_samp,:),full_binned_spks(xv_samp,t),XX,YY);

    hold_const = [0 0 1 0 0 0 0 0];
    [gabor_params_f3(t,:),LL3(t)] = fit_gabor_energy_mod(XX,YY,init_params,fullX(tr_samp,:),full_binned_spks(tr_samp,t),hold_const,LB,UB);
    xv_LL_f3(t) = get_gabor_energy_LL(gabor_params_f3(t,:),fullX(xv_samp,:),full_binned_spks(xv_samp,t),XX,YY);
    
    hold_const = [1 1 1 1 1 1 0 0];
    [gabor_params_n(t,:),LL_n(t)] = fit_gabor_energy_mod(XX,YY,init_params,fullX(tr_samp,:),full_binned_spks(tr_samp,t),hold_const,LB,UB);
    xv_LL_n(t) = get_gabor_energy_LL(gabor_params_n(t,:),fullX(xv_samp,:),full_binned_spks(xv_samp,t),XX,YY);
    
end

%%

cd ~/Data/bruce/7_15_12/
cd G029/
load ./grating_mu_data
gr_oris = unique_oris;
[~,sf_inds] = max(avg_mu_sf_profile,[],2);
pref_sfs = unique_sfs(sf_inds);
load ./oned_fixation_fits_v3.mat
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));


close all
for t = 1:96
    t
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    
    cur_params = init_params;
    cur_params(1:2) = gabor_params_f(t,1:2);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,cur_params,0);

    cur_params = init_params;
    cur_params([1 2 5]) = gabor_params_f2(t,[1 2 5]);
    gabor_emp3 = get_pgabor_mask_v2(XX,YY,cur_params,0);
    
    cur_params = init_params;
    cur_params([1 2 4 5 6]) = gabor_params_f3(t,[1 2 4 5 6]);
    gabor_emp4 = get_pgabor_mask_v2(XX,YY,cur_params,0);
    
% %     figure
%     subplot(2,2,1)
%     imagesc(gabor_emp1); set(gca,'ydir','normal'); colormap(gray);
%     title(sprintf('xvLL: %.3f',xv_LL_n(t)));
%     subplot(2,2,2)
%     imagesc(gabor_emp2); set(gca,'ydir','normal'); colormap(gray);
%     title(sprintf('xvLL: %.3f',xv_LL_f(t)));
%     subplot(2,2,3)
%     imagesc(gabor_emp3); set(gca,'ydir','normal'); colormap(gray);
%     title(sprintf('xvLL: %.3f',xv_LL_f2(t)));
%     subplot(2,2,4)
%     imagesc(gabor_emp4); set(gca,'ydir','normal'); colormap(gray);
%     title(sprintf('xvLL: %.3f',xv_LL_f3(t)));

subplot(2,2,1)
imagesc(gabor_emp1); set(gca,'ydir','normal'); colormap(gray);
title(sprintf('xvLL: %.3f',xv_LL_n(t)));

subplot(2,2,2)
imagesc(gabor_emp4); set(gca,'ydir','normal'); colormap(gray);
title(sprintf('xvLL: %.3f',xv_LL_f3(t)));

    subplot(2,2,3)
    plot(unique_oris,avg_mu_ori_profile(t,:))
    axis tight
    subplot(2,2,4)
    plot(parallel_profile(t,:));
hold on
plot(orth_profile(t,:),'r')
    axis tight
    
    if ismember(t,single_units)
        title('SU')
    end
    
    pause
    clf 
end

