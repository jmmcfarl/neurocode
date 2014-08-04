%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/
load ./eye_calibration_data
% load ./G029Expts.mat

% cd ../G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

% load ./Expt1_compiled_data_fixedlag.mat
% load ./Expt1_compiled_data_fixedlag.mat

cd ~/Data/bruce/7_15_12/G029/
% load ./Expt1_compiled_data_fixedlag_leftonly_d1_rfshift.mat
load ./Expt1_compiled_data_fixedlag_leftonly_d1p5_flipx_revx.mat
fullX_wr = fullX;

load ./Expt1_compiled_data_fixedlag_leftonly_d1p5.mat

% load ./Expt1_compiled_data_fixedlag_driftcor
% load ./Expt1_compiled_data_fixedlag_driftcor_leftonly_xshiftn1.mat
% fullX_wr = fullX;
% load ./Expt1_compiled_data_fixedlag_driftcor_leftonly.mat

% cd ~/Data/bruce/7_15_12/G034/
% load ./Expt1_compiled_data_fixedlag_driftcor_leftonly_rptframes.mat

Pix2Deg = 0.018837;
%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

cd ~/Data/bruce/7_15_12/
cd G029/

% cd G034/

load ./grating_mu_data
gr_oris = unique_oris;
[~,sf_inds] = max(avg_mu_sf_profile,[],2);
pref_sfs = unique_sfs(sf_inds);
% load ./oned_fixation_fits_v2.mat
load ./oned_fixation_fits_v3.mat

% mean_x = mean_x_old;
% mean_y = mean_y_old;

clear gabor_emp*
for t = 1:96
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    %     init_params(4) = 1/pref_sfs(t);
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
end

for t = 1:96
    init_params(1) = mean_x(t); %x0
%     init_params(1) = mean_x(ceil(rand*96)); %y0
    init_params(2) = mean_y(t); %y0
%     init_params(2) = mean_y(ceil(rand*96)); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
%     init_params(3) = rand*2*pi; %theta
    init_params(4) = 1/6; %lambda
%         init_params(4) = 1/pref_sfs(t);
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
    
    gabor_emp1_filt_wr(t,:) = gabor_emp1(:);
    gabor_emp2_filt_wr(t,:) = gabor_emp2(:);
end

[NT,klen] = size(fullX);

gabor_emp_X1 = fullX*gabor_emp1_filt';
gabor_emp_X2 = fullX*gabor_emp2_filt';
gabor_emp_X = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);

gabor_emp_X1 = fullX_wr*gabor_emp1_filt';
gabor_emp_X2 = fullX_wr*gabor_emp2_filt';
gabor_emp_X_wr = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
% gabor_emp_X1 = fullX*gabor_emp1_filt_wr';
% gabor_emp_X2 = fullX*gabor_emp2_filt_wr';
% gabor_emp_X_wr = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);

%% PARSE DATA INTO FIXATIONS
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

%% DIVIDE INTO TRAINING AND XV SETS
xv_frac = 0.2;
NT = length(used_ims);
xv_NT = round(xv_frac*NT);
% xv_samp = randperm(NT);
xv_samp = 1:NT;
xv_samp(xv_NT+1:end) = [];
tr_samp = setdiff(1:NT,xv_samp);

tr_samp = used_ims(tr_samp);
xv_samp = used_ims(xv_samp);

%% FIT AND VALIDATE GabEn MODELS

for t = 1:96;
    fprintf('SU %d of %d\n',t,96);
    cur_binned_spks = full_binned_spks(:,t);
    
    [beta(t,:),dev(t),stats(t)] = glmfit(squeeze(gabor_emp_X(tr_samp,t)),cur_binned_spks(tr_samp),'poisson');
    p_vals(t,:) = stats.p;
    pred_r = squeeze(gabor_emp_X(tr_samp,t))*beta(t,2) + beta(t,1);
    pred_r = exp(pred_r);
    LL(t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
    pred_r = squeeze(gabor_emp_X(xv_samp,t))*beta(t,2) + beta(t,1);
    pred_r = exp(pred_r);
    xvLL(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
    
    avg_rate(t) = mean(cur_binned_spks(tr_samp));
    null_LL(t) = -sum(bsxfun(@minus,cur_binned_spks(tr_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(tr_samp));
    null_xvLL(t) = -sum(bsxfun(@minus,cur_binned_spks(xv_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(xv_samp));
    
    [beta_wr(t,:),dev_wr(t),stats_wr(t)] = glmfit(squeeze(gabor_emp_X_wr(tr_samp,t)),cur_binned_spks(tr_samp),'poisson');
    p_vals_wr(t,:) = stats.p;
    pred_r = squeeze(gabor_emp_X_wr(tr_samp,t))*beta_wr(t,2) + beta_wr(t,1);
    pred_r = exp(pred_r);
    LL_wr(t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
    pred_r = squeeze(gabor_emp_X_wr(xv_samp,t))*beta_wr(t,2) + beta_wr(t,1);
    pred_r = exp(pred_r);
    xvLL_wr(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));

end

%%
expt_nums = unique(full_expt_vec);
max_len = 3915;
for t = 1:96;
    fprintf('SU %d of %d\n',t,96);
    cur_binned_spks = full_binned_spks(:,t);
    
    [beta(t,:),dev(t),stats(t)] = glmfit(squeeze(gabor_emp_X(tr_samp,t)),cur_binned_spks(tr_samp),'poisson');
    p_vals(t,:) = stats.p;
    pred_r = squeeze(gabor_emp_X(tr_samp,t))*beta(t,2) + beta(t,1);
    pred_r = exp(pred_r);
    LL(t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
    pred_r = squeeze(gabor_emp_X(xv_samp,t))*beta(t,2) + beta(t,1);
    pred_r = exp(pred_r);
    xvLL(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
    
    avg_rate(t) = mean(cur_binned_spks(tr_samp));
    null_LL(t) = -sum(bsxfun(@minus,cur_binned_spks(tr_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(tr_samp));
    null_xvLL(t) = -sum(bsxfun(@minus,cur_binned_spks(xv_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(xv_samp));
    
    for ee = 1:length(expt_nums)
        
        cur_set = tr_samp(full_expt_vec(tr_samp) == expt_nums(ee));
        cur_set(max_len+1:end) = [];
        [beta_e(t,ee,:),dev_e(t,ee),stats(t)] = glmfit(squeeze(gabor_emp_X(cur_set,t)),cur_binned_spks(cur_set),'poisson');
        cur_len(ee) = length(cur_set);
    end
    
end

