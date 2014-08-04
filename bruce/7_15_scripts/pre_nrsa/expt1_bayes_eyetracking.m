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
% load ./Expt1_newcompiled_data_fixedlag_leftonly_d1p5.mat
% load ./Expt1_newcompiled_data_fixedlag_withdrift_leftonly_d1p5
load ./Expt1_newcompiled_data_fixedlag_d1p5
fullX = fullX/std(fullX(:));

% fullX = bsxfun(@minus,fullX,mean(fullX));
% fullX = bsxfun(@rdivide,fullX,std(fullX));

Pix2Deg = 0.018837;
[NT,klen] = size(fullX);

%%
new_RF_patch = [-0.1 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);
new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
fullX_cropped = fullX(:,new_crop);

xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

%% PARSE DATA INTO FIXATIONS

n_fixs = length(fix_start_inds);

% avg_x_pos = 0.5*full_eyepos(:,1) + 0.5*full_eyepos(:,3);
% avg_y_pos = 0.5*full_eyepos(:,2) + 0.5*full_eyepos(:,4);

%use only left eye position
avg_x_pos = full_eyepos(:,1);
avg_y_pos = full_eyepos(:,2);
% avg_x_pos = full_eyepos(:,3);
% avg_y_pos = full_eyepos(:,4);

diff_used_inds = diff(used_inds);
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_fixs = length(rel_fix_start_inds);
x_drifts = nan(size(avg_x_pos));
y_drifts = nan(size(avg_y_pos));
medx_pos = nan(n_fixs,1);
medy_pos = nan(n_fixs,1);
for i = 1:n_fixs
    cur_set = rel_fix_start_inds(i):rel_fix_stop_inds(i);
    medx_pos(i) = median(avg_x_pos(cur_set));
    medy_pos(i) = median(avg_y_pos(cur_set));
    x_drifts(cur_set) = avg_x_pos(cur_set) - medx_pos(i);
    y_drifts(cur_set) = avg_y_pos(cur_set) - medy_pos(i);
end
fix_deltaxs = [0; diff(medx_pos)];
fix_deltays = [0; diff(medy_pos)];

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
sdim = length(xpatch_inds);
[XXc,YYc] = meshgrid(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped));
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

cd ~/Data/bruce/7_15_12/
cd G029/
load ./grating_mu_data
gr_oris = unique_oris;
[~,sf_inds] = max(avg_mu_sf_profile,[],2);
pref_sfs = unique_sfs(sf_inds);
load ./oned_fixation_fits_v3.mat

LB = [-0.1 -0.8 0 0.12 0.025 0.5 0 -Inf];
UB = [0.8 0.1 pi 0.4 0.25 2 Inf Inf];
hold_const = [0 0 1 0 0 0 0 0];
all_const = [1 1 1 1 1 1 0 0];

clear gabor_emp*
for t = 1:96
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.4*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    init_params(7) = 0; %weight of quad term
    init_params(8) = 0; %const offset
    
    [gabor_params_f(t,:),LL(t)] = fit_gabor_energy_mod(XXc,YYc,init_params,fullX_cropped,full_binned_spks(:,t),hold_const,LB,UB);
    [gabor_params_n(t,:),LLn(t)] = fit_gabor_energy_mod(XXc,YYc,init_params,fullX_cropped,full_binned_spks(:,t),all_const,LB,UB);
    
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params_f(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params_f(t,1:6),pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
end

gabor_emp_X1 = fullX*gabor_emp1_filt';
gabor_emp_X2 = fullX*gabor_emp2_filt';
gabor_emp_X = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);


%% DIVIDE INTO TRAINING AND XV SETS
% nfold = 5;
% nparts = 50;
% [tr_samp,xv_samp] = create_xv_set(NT,nfold,nparts);
tr_samp = 1:NT;
xv_samp = [];

%% FIT AND VALIDATE GabEn MODELS

for t = 1:96;
    fprintf('SU %d of %d\n',t,96);
    cur_binned_spks = full_binned_spks(:,t);
    
    [beta(t,:),dev(t),stats(t)] = glmfit(squeeze(gabor_emp_X(tr_samp,t)),cur_binned_spks(tr_samp),'poisson');
    p_vals(t,:) = stats.p;
    pred_r = squeeze(gabor_emp_X(tr_samp,t))*beta(t,2) + beta(t,1);
    pred_r = exp(pred_r);
    LL(t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
%     pred_r = squeeze(gabor_emp_X(xv_samp,t))*beta(t,2) + beta(t,1);
%     pred_r = exp(pred_r);
%     xvLL(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
    
    avg_rate(t) = mean(cur_binned_spks(tr_samp));
    null_LL(t) = -sum(bsxfun(@minus,cur_binned_spks(tr_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(tr_samp));
%     null_xvLL(t) = -sum(bsxfun(@minus,cur_binned_spks(xv_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(xv_samp));
    
end
pvals = [stats(:).p];
sig_fits = find(pvals(2,:) < .01);

%% COMPUTE BETA FOR EACH MULT FACTOR

gain_mults = [0.75 1 1.5 2];
for t = 1:96;
    fprintf('SU %d of %d\n',t,96);
    spikebins = convert_to_spikebins(full_binned_spks(:,t));
        
    [fitp,grad] = GLMsolve_jmm( gabor_emp_X(:,t), spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
    beta(t,:) = fitp.k;
    LL(t) = fitp.LL;
    for i = 1:length(gain_mults)
        [fitp,grad] = GLMsolve_jmm( gabor_emp_X(:,t), spikebins, [gain_mults(i)*beta(t,1); 0], 1, [], [], [], [], [], [1], 0);
        beta_m(t,i,:) = fitp.k;
        LL_m(t,i) = fitp.LL;
    end
end

%%
max_shift = 16;
dshift = 1;
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;
[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
SH = [Xsh(:) Ysh(:)];
n_shifts = size(SH,1);

fix_LLs = zeros(n_fixs,96,length(gain_mults),n_shifts);

for fix = 1:n_fixs
    
    fprintf('Fixation %d of %d\n',fix,n_fixs); pause(0.1)
    
    tot_LLs = zeros(96,length(gain_mults),n_shifts);
    cur_ims = rel_fix_start_inds(fix):rel_fix_stop_inds(fix);
    
    %     base_LLs = zeros(n_shifts,1);
    
    %     l_like_table = zeros(length(gain_mults),n_shifts);
    for ii = 1:length(cur_ims);
        %         fprintf('Image %d of %d\n',ii,length(cur_ims));
        cur_im = reshape(fullX(cur_ims(ii),:),sdim,sdim);
        Robs = full_binned_spks(cur_ims(ii),:)';
        
        shifted_ims = nan(length(x_shifts)*length(y_shifts),sdim^2);
        cnt = 1;
        for xx = 1:length(x_shifts)
            for yy = 1:length(y_shifts)
                d2 = dist_shift2d(cur_im, x_shifts(xx), 2,0);
                d2 = dist_shift2d(d2,y_shifts(yy),1,0);
                shifted_ims(cnt,:) = d2(:);
                cnt = cnt + 1;
            end
        end
        
        gabor_outs1 = gabor_emp1_filt*shifted_ims';
        gabor_outs2 = gabor_emp2_filt*shifted_ims';
        gabor_outs = gabor_outs1.^2 + gabor_outs2.^2;
        gabor_outs = sqrt(gabor_outs);
        gabor_outs = reshape(gabor_outs,[96 1 n_shifts]);
%         gabor_outs = repmat(reshape(gabor_outs,[96 1 n_shifts]),[1 length(gain_mults) 1]);
        pred_Rs = bsxfun(@times,gabor_outs,beta_m(:,:,1));
        pred_Rs = bsxfun(@plus,pred_Rs,beta_m(:,:,2));
        pred_Rs = log(1+exp(pred_Rs));
%         LLs = bsxfun(@times,log(pred_Rs),reshape(Robs,[96 1 1]));
        LLs = bsxfun(@times,log(pred_Rs),Robs);
        LLs = LLs - pred_Rs;
        %constant with respect to model params (theta and eps)
%         LLs = bsxfun(@minus,LLs,reshape(log(factorial(Robs)),[96 1 1])); 

        
        fix_LLs(fix,:,:,:) = fix_LLs(fix,:,:,:) + reshape(LLs,[1 size(LLs)]);
    end
end



%%
% sig_fits = 1:96;
NSIG = length(sig_fits);
xv_frac = 0.2;
tr_set = randperm(NSIG);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NSIG,tr_set);
tr_set = sig_fits(tr_set);
xv_set = sig_fits(xv_set);
bad_set = find(~ismember(1:96,sig_fits));
%%
null_shift = find(SH(:,1) == 0 & SH(:,2) == 0);

clear lpeps lptheta
n_iter = 2;
lpeps(1,:,:) = nan(n_fixs,n_shifts);

% lptheta(1,:,:) = log(ones(96,length(gain_mults))/length(gain_mults)); %uniform prior
mval = 1e-50;
lptheta(1,:,:) = mval*ones(96,length(gain_mults));
lptheta(1,:,2) = 1;
lptheta(1,:,:) = log(lptheta(1,:,:));

lptheta(1,:,:) = bsxfun(@minus,lptheta(1,:,:),logsumexp(lptheta(1,:,:),3));

eps_prior_sigma = 0.25;
leps_prior = zeros(n_shifts,1);
leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior));

fix_deltax_eps = round(fix_deltaxs*Fsd/dshift);
fix_deltay_eps = round(fix_deltays*Fsd/dshift);
fix_delta_d = sqrt(fix_deltaxs.^2 + fix_deltays.^2);
max_delta_d = 0.5;

deps_sigma = 0.15*Fsd;
NT = n_fixs;
for it = 2:n_iter + 1
    
    fprintf('Iteration %d of %d...',it-1,n_iter);
    
    %% COMPUTE POSTERIOR ON EPS GIVEN P(theta)
    temp = reshape(lptheta(it-1,:,:),[1 96 length(gain_mults)]);
    lposterior_eps = bsxfun(@plus,fix_LLs,temp); %multiply L(theta,eps) by p(theta)
    lposterior_eps = squeeze(logsumexp(lposterior_eps,3)); %integrate out theta
    
    lposterior_eps = squeeze(sum(lposterior_eps(:,tr_set,:),2)); %now sum across cells
    lposterior_eps = bsxfun(@plus,lposterior_eps,leps_prior'); %prior on eps
    lposterior_eps = bsxfun(@minus,lposterior_eps,logsumexp(lposterior_eps,2));
    lpeps(it,:,:) = lposterior_eps;
    
    
%     lB = squeeze(sum(lposterior_eps(:,tr_set,:),2)); %now sum across cells
%     if any(isnan(lB))
%         error('NANS DETECTED in LL!')
%     end
%     
%     fprintf(' Forward messages ')
%     lalpha=zeros(NT,n_shifts);
%     lscale=zeros(NT,1); %initialize rescaling parameters
%     %compute rescaled forward messages
%     lalpha(1,:) = leps_prior' + lB(1,:);
%     lscale(1)=logsumexp(lalpha(1,:));
%     lalpha(1,:) = lalpha(1,:) - lscale(1);
%     for t=2:NT
%         if fix_delta_d(t) < max_delta_d
%             cdist = pdist2(SH,bsxfun(@plus,SH,[fix_deltax_eps(t) fix_deltay_eps(t)]));
%             lA = -cdist.^2/(2*deps_sigma^2);
%             lA = bsxfun(@minus,lA,logsumexp(lA,2));
%         else
%             lA = repmat(leps_prior',n_shifts,1);
%         end
%         
%         %         lalpha(t,:) = log((lalpha(t-1,:))*A) + lB(t,:);
%         lalpha(t,:) = logmulexp(lalpha(t-1,:),lA) + lB(t,:);
%         lscale(t) = logsumexp(lalpha(t,:));
%         lalpha(t,:)= lalpha(t,:) - lscale(t);
%     end
%     
%     fprintf(' Backwards messages ')
%     %compute rescaled backward messages
%     lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
%     for t=NT-1:-1:1
%         if fix_delta_d(t) < max_delta_d
%             cdist = pdist2(SH,bsxfun(@plus,SH,[fix_deltax_eps(t) fix_deltay_eps(t)]));
%             lA = -cdist.^2/(2*deps_sigma^2);
%             lA = bsxfun(@minus,lA,logsumexp(lA,2));
%         else
%             lA = repmat(leps_prior',n_shifts,1);
%         end
%         %         beta(t,:)=(beta(t+1,:).*B(t+1,:))*(A')/scale(t);
%         lf1 = lbeta(t+1,:) + lB(t+1,:);
%         lbeta(t,:) = logmulexp(lf1,lA') - lscale(t);
%         %         lbeta(t,:) = logsumexp(bsxfun(@plus,lf1,lf2)) - lscale(t);
%         %         lbeta(t,:)= (lbeta(t+1,:) + lB(t+1,:)) *(A') -lscale(t);
%     end
%     %             beta(t,:)=(beta(t+1,:).*B(t+1,:))*(A')/scale(t);
%     
%     %compute posteriors over hidden states
%     lgamma= lalpha + lbeta;
%     lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
%     lpeps(it,:,:) = lgamma;
    
    %% NOW COMPUTE P(THETA) GIVEN P(EPS)
    fprintf(' Computing p(theta)\n');
    temp = reshape(lpeps(it,:,:),[NT 1 1 n_shifts]);
    lposterior_theta = bsxfun(@plus,fix_LLs,temp); %multiply like by current estimates of p(eps)
    lposterior_theta = logsumexp(lposterior_theta,4); %sum over epsilon
    %     lposterior_theta = bsxfun(@minus,lposterior_theta,logsumexp(lposterior_theta,3)); %normalize to posteriors on theta
    lposterior_theta = squeeze(sum(lposterior_theta)); %sum log posteriors over fixations
%     llike_theta(it,:,:) = lposterior_theta;
    lposterior_theta = bsxfun(@minus,lposterior_theta,logsumexp(lposterior_theta,2)); %normalize p(thetas)
    lptheta(it,:,:) = lposterior_theta;
    
    %%
    
end

% temp_lpeps = ones(length(used_fixs),n_shifts);
% temp_lpeps = bsxfun(@rdivide,temp_lpeps,sum(temp_lpeps,2));
% temp_lpeps = log(temp_lpeps);
% temp = reshape(temp_lpeps,[length(used_fixs) 1 1 n_shifts]);
% lposterior_theta = bsxfun(@plus,fix_LLs,temp); %multiply like by current estimates of p(eps)
% lposterior_theta = logsumexp(lposterior_theta,4); %sum over epsilon
% lposterior_theta = squeeze(sum(lposterior_theta)); %sum log posteriors over fixations
% null_llike = lposterior_theta;
% 
% temp_lpeps = zeros(length(used_fixs),n_shifts);
% temp_lpeps(:,null_shift) = 1;
% temp_lpeps = log(temp_lpeps);
% temp = reshape(temp_lpeps,[length(used_fixs) 1 1 n_shifts]);
% lposterior_theta = bsxfun(@plus,fix_LLs,temp); %multiply like by current estimates of p(eps)
% lposterior_theta = logsumexp(lposterior_theta,4); %sum over epsilon
% lposterior_theta = squeeze(sum(lposterior_theta)); %sum log posteriors over fixations
% orig_llike = lposterior_theta;

% FIND MAP SEQUENCE OF ERRORS
clear eps
for i = 1:n_iter
    [max_post,max_loc] = max(lpeps(i+1,:,:),[],3);
    eps(i,:,:) = SH(max_loc,:);
end

%% create new estimate of stimulus
clear beta beta_sh LL_sh

% REFIT GE MODELS
gabor_emp_X1 = fullX*gabor_emp1_filt';
gabor_emp_X2 = fullX*gabor_emp2_filt';
gabor_emp_X = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
for ii = 1:n_iter
    fullX_sh = zeros(size(fullX));
    for fix = 1:n_fixs
        fprintf('Fixation %d of %d\n',fix,n_fixs);
        cur_ims = rel_fix_start_inds(fix):rel_fix_stop_inds(fix);
        for jj = 1:length(cur_ims);
            cur_im = reshape(fullX(cur_ims(jj),:),sdim,sdim);
            d2 = dist_shift2d(cur_im, eps(ii,fix,1), 2,0);
            d2 = dist_shift2d(d2,eps(ii,fix,2),1,0);
            fullX_sh(cur_ims(jj),:) = d2(:);
        end
    end
    
    gabor_emp_X1 = fullX_sh*gabor_emp1_filt';
    gabor_emp_X2 = fullX_sh*gabor_emp2_filt';
    gabor_emp_X_sh = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
    for t = 1:96;
        fprintf('SU %d of %d\n',t,96);
        spikebins = convert_to_spikebins(full_binned_spks(:,t)); 
        avg_rate(t) = mean(full_binned_spks(:,t));
        n_spks = sum(full_binned_spks(:,t));
        null_LL(t) = -sum(bsxfun(@minus,full_binned_spks(:,t)*log(avg_rate(t)),avg_rate(t)))/n_spks;
                
        [fitp,grad] = GLMsolve_jmm(gabor_emp_X_sh(:,t), spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
        beta_sh(ii,t,:) = fitp.k;
        LL_sh(ii,t) = fitp.LL;        
        
        [fitp,grad] = GLMsolve_jmm(gabor_emp_X(:,t), spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
        beta(t,:) = fitp.k;
        LL(t) = fitp.LL;               
    end
end
LL_imp = bsxfun(@minus,LL_sh,LL);

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
fullX_sh = zeros(size(fullX));
use_it = 2;
for fix = 1:n_fixs
    fprintf('Fixation %d of %d\n',fix,n_fixs);
    cur_ims = rel_fix_start_inds(fix):rel_fix_stop_inds(fix);
    for jj = 1:length(cur_ims);
        cur_im = reshape(fullX(cur_ims(jj),:),sdim,sdim);
        d2 = dist_shift2d(cur_im, eps(use_it,fix,1), 2,0);
        d2 = dist_shift2d(d2,eps(use_it,fix,2),1,0);
        fullX_sh(cur_ims(jj),:) = d2(:);
    end
end
fullX_cropped = fullX_sh(:,new_crop);

LB = [-0.1 -0.8 0 0.12 0.025 0.16 0 -Inf];
UB = [0.8 0.1 pi 0.4 0.25 6.25 Inf Inf];
hold_const = [0 0 1 0 0 0 0 0];
all_const = [1 1 1 1 1 1 0 0];

clear gabor_emp*
for t = 1:96
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    
    init_params = gabor_params_f(t,:);
    [gabor_params_r(t,:),LLr(t)] = fit_gabor_energy_mod(XXc,YYc,init_params,fullX_cropped,full_binned_spks(:,t),hold_const,LB,UB);
    
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params_r(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params_r(t,1:6),pi/2);
    
    gabor_emp1_filt_r(t,:) = gabor_emp1(:);
    gabor_emp2_filt_r(t,:) = gabor_emp2(:);
end


%%
fullX_cropped_new = fullX_sh(:,new_crop);
fullX_cropped_old = fullX(:,new_crop);
LB = [-0.1 -0.8 0 0.12 0.025 0.5 0 -Inf -Inf -Inf];
UB = [0.8 0.1 pi 0.4 0.25 2 Inf Inf Inf Inf];
hold_const = [0 0 1 0 0 0 0 0 0 0];
for t = 1:96
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    
    cur_gabor_0 = get_pgabor_mask_v2(XXc,YYc,gabor_params_r(t,1:6),0);
    cur_gabor_1 = get_pgabor_mask_v2(XXc,YYc,gabor_params_r(t,1:6),pi/2);
    gabor_out_0 = fullX_cropped_new*cur_gabor_0(:);
    gabor_out_1 = fullX_cropped_new*cur_gabor_1(:);
    gabor_energy = sqrt(gabor_out_0.^2 + gabor_out_1.^2);
    
    [beta_new(t,:),dev_new(t),stats_new(t)] = glmfit([gabor_out_0 gabor_out_1 gabor_energy],full_binned_spks(:,t),'Poisson'); 

    cur_gabor_0 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f(t,1:6),0);
    cur_gabor_1 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f(t,1:6),pi/2);
    gabor_out_0 = fullX_cropped_old*cur_gabor_0(:);
    gabor_out_1 = fullX_cropped_old*cur_gabor_1(:);
    gabor_energy = sqrt(gabor_out_0.^2 + gabor_out_1.^2);
    
    [beta_old(t,:),dev_old(t),stats_old(t)] = glmfit([gabor_out_0 gabor_out_1 gabor_energy],full_binned_spks(:,t),'Poisson'); 
end
pvals_new = [stats_new(:).p];
pvals_old = [stats_old(:).p];
%%
close all
for t = 1:96;
    gabor_emp1_r = get_pgabor_mask_v2(XXc,YYc,gabor_params_r(t,1:6),0);
    gabor_emp2_r = get_pgabor_mask_v2(XXc,YYc,gabor_params_r(t,1:6),pi/2);
    gabor_emp_X1 = fullX_cropped*gabor_emp1_r(:);
    gabor_emp_X2 = fullX_cropped*gabor_emp2_r(:);
    gabor_emp_X_r = zscore(sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2));
    beta_r = glmfit(gabor_emp_X_r,full_binned_spks(:,t));
    
    gabor_emp1_f = get_pgabor_mask_v2(XXc,YYc,gabor_params_f(t,1:6),0);
    gabor_emp2_f = get_pgabor_mask_v2(XXc,YYc,gabor_params_f(t,1:6),pi/2);
    gabor_emp_X1 = fullX_cropped*gabor_emp1_f(:);
    gabor_emp_X2 = fullX_cropped*gabor_emp2_f(:);
    gabor_emp_X_f = zscore(sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2));
    beta_f = glmfit(gabor_emp_X_f,full_binned_spks(:,t));

    gabor_emp1_n = get_pgabor_mask_v2(XXc,YYc,gabor_params_n(t,1:6),0);
    gabor_emp2_n = get_pgabor_mask_v2(XXc,YYc,gabor_params_n(t,1:6),pi/2);
    gabor_emp_X1 = fullX_cropped*gabor_emp1_n(:);
    gabor_emp_X2 = fullX_cropped*gabor_emp2_n(:);
    gabor_emp_X_n = zscore(sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2));
    beta_n = glmfit(gabor_emp_X_n,full_binned_spks(:,t));
    
    subplot(3,1,1)
    imagesc(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped),gabor_emp1_n);set(gca,'ydir','normal');
    if ~ismember(t,xv_set) colormap(gray); else colormap(jet); end
    title(sprintf('Weight %.3f',beta_n(2)));
    subplot(3,1,2)
    imagesc(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped),gabor_emp1_f);set(gca,'ydir','normal');
    if ~ismember(t,xv_set) colormap(gray); else colormap(jet); end
    title(sprintf('Weight %.3f',beta_f(2)));
    subplot(3,1,3)
    imagesc(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped),gabor_emp1_r);set(gca,'ydir','normal');
    if ~ismember(t,xv_set) colormap(gray); else colormap(jet); end
    title(sprintf('Weight %.3f',beta_r(2)));
    pause
    clf
end

%% REFIT GABOR PARAMS
xv_frac = 0.2;
NT = length(used_ims);
xv_NT = round(xv_frac*NT);
% xv_samp = randperm(NT);
xv_samp = 1:NT;
xv_samp(xv_NT+1:end) = [];
tr_samp = setdiff(1:NT,xv_samp);

tr_samp = used_ims(tr_samp);
xv_samp = used_ims(xv_samp);

for t = single_units;
    fprintf('SU %d of %d\n',t,96);
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    init_params(7) = 0; %weight of quad term
    init_params(8) = 0; %weight of lin term 0 phase
    init_params(9) = 0; %weight of lin term pi/2 phase
    init_params(10) = 0; %const offset
    hold_const = [0 0 1 1 1 1 0 1 1 0];
    LB = [0.1 -0.7 0 0.125 0.02 0.5 -Inf -Inf -Inf -Inf];
    UB = [0.7 -0.1 pi 0.25 0.25 2 Inf Inf Inf Inf];
    
    [gabor_params_f(t,:),LL_o(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX(tr_samp,:),full_binned_spks(tr_samp,t),hold_const,LB,UB);
    gabor_bank_f(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_f(t,1:6),0);
    xv_LL_f(t) = get_pgabor_LL_v2(gabor_params_f(t,:),fullX(xv_samp,:),full_binned_spks(xv_samp,t),XX,YY);
    
    hold_const = [1 1 1 1 1 1 0 1 1 0];
    [gabor_params_n(t,:),LL_n(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX(tr_samp,:),full_binned_spks(tr_samp,t),hold_const,LB,UB);
    gabor_bank_n(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_n(t,1:6),0);
    xv_LL_n(t) = get_pgabor_LL_v2(gabor_params_n(t,:),fullX(xv_samp,:),full_binned_spks(xv_samp,t),XX,YY);
    
end

%% FIT RF POSITION MAPs
nsh1 = sqrt(n_shifts);
cur_lpeps = squeeze(lpeps(end,:,:));
lposmap = zeros(n_shifts,96);
for cid = 1:96;
    fprintf('Cell %d of %d\n',cid,96);
    cur_lptheta = squeeze(lptheta(end,cid,:));
    %     cur_lptheta = zeros(1,length(gain_mults));
    %     cur_lptheta(4) = 1;
    %     cur_lptheta = log(cur_lptheta)';
    for t = 1:NT
        %         fprintf('Fixation %d of %d\n',t,NT);
        temp = bsxfun(@plus,squeeze(fix_LLs(t,cid,:,:)),cur_lptheta);
        cur_LL = reshape(logsumexp(temp,1),nsh1,nsh1);
        cur_lprior = reshape(cur_lpeps(t,:),nsh1,nsh1);
        max_L = max(cur_LL(:));
        %         max_p = max(cur_lprior(:));
        cur_L = exp(cur_LL-max_L);
        cur_prior = exp(cur_lprior);
        
        %         cur_prior(cur_prior < max(cur_prior(:))) = 1e-100; cur_prior = bsxfun(@rdivide,cur_prior,sum(cur_prior(:)));
        
        %         cur_L = exp(bsxfun(@minus,cur_LL,max_L));
        %         cur_prior = exp(bsxfun(@minus,cur_lprior,max_p));
        cur_conv = conv2(cur_L,fliplr(flipud(cur_prior)),'same');
        cur_conv_temp = conv2(ones(nsh1,nsh1),fliplr(flipud(cur_prior)),'same');
        cur_conv = cur_conv./cur_conv_temp;
        cur_conv(isnan(cur_conv)) = 0;
        %         lcur_conv = bsxfun(@plus,log(cur_conv),max_L);
        lcur_conv = log(cur_conv);
        lposmap(:,cid) = lposmap(:,cid) + lcur_conv(:);
    end
end
%
prior_sigma = 0.3*Fsd;
poserror_prior = exp(-sum(SH.^2,2)./(2*prior_sigma^2));
poserror_prior = poserror_prior/sum(poserror_prior);

lposmap_wprior = bsxfun(@plus,lposmap,log(poserror_prior));

% lposmap = bsxfun(@minus,lposmap,logsumexp(lposmap,1));
lposmap_wprior = bsxfun(@minus,lposmap_wprior,logsumexp(lposmap_wprior,1));
[a,b] = max(lposmap_wprior);
% [a,b] = max(lposmap);
rf_x = SH(b,1);
rf_y = SH(b,2);
dd = sqrt(rf_x.^2+rf_y.^2)/Fsd;

%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

%% REFIT GABOR PARAMS to see if RF correction improves fits
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
    
    init_params(1) = init_params(1) - rf_x(t)/Fsd;
    init_params(2) = init_params(2) - rf_y(t)/Fsd;
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
    gabor_emp1_filtc(t,:) = gabor_emp1(:);
    gabor_emp2_filtc(t,:) = gabor_emp2(:);
end

% [NT,klen] = size(fullX);

gabor_emp_X1 = fullX_sh*gabor_emp1_filt';
gabor_emp_X2 = fullX_sh*gabor_emp2_filt';
gabor_emp_X = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);

gabor_emp_X1 = fullX_sh*gabor_emp1_filtc';
gabor_emp_X2 = fullX_sh*gabor_emp2_filtc';
gabor_emp_Xc = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);

for t = 1:96;
    fprintf('SU %d of %d\n',t,96);
    cur_binned_spks = full_binned_spks(:,t);
    
    [beta(t,:),dev(t),stats(t)] = glmfit(squeeze(gabor_emp_X(used_ims,t)),cur_binned_spks(used_ims),'poisson');
    [betac(t,:),devc(t),statsc(t)] = glmfit(squeeze(gabor_emp_Xc(used_ims,t)),cur_binned_spks(used_ims),'poisson');
    
end

%% FIT RF POSITION MAPs for cross-val check
nsh1 = sqrt(n_shifts);
cur_lpeps = squeeze(lpeps(end,:,:));
lposmaptr = zeros(n_shifts,96);
lposmapxv = zeros(n_shifts,96);
xv_set = randperm(NT);
xv_set(round(NT/2):end) = [];
tr_set = setdiff(1:NT,xv_set);

for cid = 1:96;
    fprintf('Cell %d of %d\n',cid,96);
    cur_lptheta = squeeze(lptheta(end,cid,:));
    for t = 1:length(tr_set)
        temp = bsxfun(@plus,squeeze(fix_LLs(tr_set(t),cid,:,:)),cur_lptheta);
        cur_LL = reshape(logsumexp(temp,1),nsh1,nsh1);
        cur_lprior = reshape(cur_lpeps(tr_set(t),:),nsh1,nsh1);
        max_L = max(cur_LL(:));
        cur_L = exp(cur_LL-max_L);
        cur_prior = exp(cur_lprior);
        
        cur_conv = conv2(cur_L,fliplr(flipud(cur_prior)),'same');
        cur_conv_temp = conv2(ones(nsh1,nsh1),fliplr(flipud(cur_prior)),'same');
        cur_conv = cur_conv./cur_conv_temp;
        cur_conv(isnan(cur_conv)) = 0;
        lcur_conv = log(cur_conv);
        lposmaptr(:,cid) = lposmaptr(:,cid) + lcur_conv(:);
    end
end
for cid = 1:96;
    fprintf('Cell %d of %d\n',cid,96);
    cur_lptheta = squeeze(lptheta(end,cid,:));
    for t = 1:length(xv_set)
        temp = bsxfun(@plus,squeeze(fix_LLs(xv_set(t),cid,:,:)),cur_lptheta);
        cur_LL = reshape(logsumexp(temp,1),nsh1,nsh1);
        cur_lprior = reshape(cur_lpeps(xv_set(t),:),nsh1,nsh1);
        max_L = max(cur_LL(:));
        cur_L = exp(cur_LL-max_L);
        cur_prior = exp(cur_lprior);
        
        cur_conv = conv2(cur_L,fliplr(flipud(cur_prior)),'same');
        cur_conv_temp = conv2(ones(nsh1,nsh1),fliplr(flipud(cur_prior)),'same');
        cur_conv = cur_conv./cur_conv_temp;
        cur_conv(isnan(cur_conv)) = 0;
        lcur_conv = log(cur_conv);
        lposmapxv(:,cid) = lposmapxv(:,cid) + lcur_conv(:);
    end
end

prior_sigma = 0.05*Fsd;
poserror_prior = exp(-sum(SH.^2,2)./(2*prior_sigma^2));
poserror_prior = poserror_prior/sum(poserror_prior);

lposmap_wpriortr = bsxfun(@plus,lposmaptr,log(poserror_prior));
lposmap_wpriorxv = bsxfun(@plus,lposmapxv,log(poserror_prior));

[a,b] = max(lposmaptr);
% [a,b] = max(lposmap_wpriortr);
rf_xtr = SH(b,1);
rf_ytr = SH(b,2);
ddtr = sqrt(rf_xtr.^2+rf_ytr.^2)/Fsd;

[a,b] = max(lposmapxv);
% [a,b] = max(lposmap_wpriorxv);
rf_xxv= SH(b,1);
rf_yxv = SH(b,2);
ddxv = sqrt(rf_xxv.^2+rf_yxv.^2)/Fsd;

%% CONSTRUCT ADJUSTED LL ARRAY USING NEW RF POSITIONS
fix_LLs_adj = fix_LLs;
adj_set = find(dd ~= 0);
NT = size(fix_LLs,1);
for c = 1:length(adj_set)
    fprintf('Adjusting cell %d of %d\n',c,length(adj_set));
    
    cur_x_adj = -rf_x(adj_set(c));
    cur_y_adj = -rf_y(adj_set(c));
    
    cur_array = reshape(fix_LLs(:,adj_set(c),:,:),[NT length(gain_mults) nsh1 nsh1]);
    nullvals = squeeze(mean(fix_LLs(:,adj_set(c),:,:),4));
    mapper = mod((0:(nsh1-1))-cur_x_adj,nsh1)+1;
    cur_array_sh = cur_array(:,:,:,mapper);
    if cur_x_adj > 0
        cur_array_sh(:,:,:,1:cur_x_adj) = repmat(nullvals,[1 1 nsh1 cur_x_adj]);
    elseif cur_x_adj < 0
        cur_array_sh(:,:,:,(nsh1-cur_x_adj+1:end)) = repmat(nullvals,[1 1 nsh1 cur_x_adj]);
    end
    
    mapper = mod((0:(nsh1-1))-cur_y_adj,nsh1)+1;
    cur_array_sh = cur_array_sh(:,:,mapper,:);
    if cur_y_adj > 0
        cur_array_sh(:,:,1:cur_y_adj,:) = repmat(nullvals,[1 1 cur_y_adj nsh1]);
    elseif cur_y_adj < 0
        cur_array_sh(:,:,(nsh1-cur_y_adj+1:end),:) = repmat(nullvals,[1 1 cur_y_adj nsh1]);
    end
    
    fix_LLs_adj(:,adj_set(c),:,:) = reshape(cur_array_sh,[NT 1 length(gain_mults) n_shifts]);
end

%%
n_iter = 2;
lpepsr(1,:,:) = nan(length(used_fixs),n_shifts);

% lptheta(1,:,:) = log(ones(96,length(gain_mults))/length(gain_mults)); %uniform prior
mval = 1e-50;
lpthetar(1,:,:) = mval*ones(96,length(gain_mults));
lpthetar(1,:,3) = 1;
lpthetar(1,:,:) = log(lpthetar(1,:,:));
lpthetar(1,:,:) = bsxfun(@minus,lpthetar(1,:,:),logsumexp(lpthetar(1,:,:),3));

deps_sigma = 0.2*Fsd;
NT = length(used_fixs);
for it = 2:n_iter + 1
    
    fprintf('Iteration %d of %d...',it-1,n_iter);
        
    %%
    temp = reshape(lpthetar(it-1,:,:),[1 96 length(gain_mults)]);
    lposterior_eps = bsxfun(@plus,fix_LLs_adj,temp); %multiply L(theta,eps) by p(theta)
%     lposterior_eps(isnan(lposterior_eps)) = -Inf; %strange problem with infs converted to nans
    lposterior_eps_s = squeeze(logsumexp(lposterior_eps,3)); %integrate out theta
    
%     lposterior_eps = squeeze(sum(lposterior_eps_s(:,tr_set,:),2)); %now sum across cells
%     lposterior_eps = bsxfun(@plus,lposterior_eps,leps_prior'); %prior on eps
%     lposterior_eps = bsxfun(@minus,lposterior_eps,logsumexp(lposterior_eps,2));
%     lpeps(it,:,:) = lposterior_eps;
    
    
    lB = squeeze(sum(lposterior_eps_s(:,tr_set,:),2)); %now sum across cells
    if any(isnan(lB))
        error('NANS DETECTED in LL!')
    end
    
    fprintf(' Forward messages ')
    lalpha=zeros(NT,n_shifts);
    lscale=zeros(NT,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior' + lB(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:NT
        if fix_delta_d(t) < max_delta_d
            cdist = pdist2(SH,bsxfun(@plus,SH,[fix_deltax_eps(t) fix_deltay_eps(t)]));
            lA = -cdist.^2/(2*deps_sigma^2);
            lA = bsxfun(@minus,lA,logsumexp(lA,2));
        else
            lA = repmat(leps_prior',n_shifts,1);
        end
        
        %         lalpha(t,:) = log((lalpha(t-1,:))*A) + lB(t,:);
        lalpha(t,:) = logmulexp(lalpha(t-1,:),lA) + lB(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:)= lalpha(t,:) - lscale(t);
    end
    
    fprintf(' Backwards messages ')
    %compute rescaled backward messages
    lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
    for t=NT-1:-1:1
        if fix_delta_d(t) < max_delta_d
            cdist = pdist2(SH,bsxfun(@plus,SH,[fix_deltax_eps(t) fix_deltay_eps(t)]));
            lA = -cdist.^2/(2*deps_sigma^2);
            lA = bsxfun(@minus,lA,logsumexp(lA,2));
        else
            lA = repmat(leps_prior',n_shifts,1);
        end
        %         beta(t,:)=(beta(t+1,:).*B(t+1,:))*(A')/scale(t);
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,lA') - lscale(t);
        %         lbeta(t,:) = logsumexp(bsxfun(@plus,lf1,lf2)) - lscale(t);
        %         lbeta(t,:)= (lbeta(t+1,:) + lB(t+1,:)) *(A') -lscale(t);
    end
    %             beta(t,:)=(beta(t+1,:).*B(t+1,:))*(A')/scale(t);
    
    %compute posteriors over hidden states
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    lpepsr(it,:,:) = lgamma;
    
    %% NOW COMPUTE P(THETA) GIVEN P(EPS)
    fprintf(' Computing p(theta)\n');
    temp = reshape(lpepsr(it,:,:),[length(used_fixs) 1 1 n_shifts]);
    lposterior_theta = bsxfun(@plus,fix_LLs_adj,temp); %multiply like by current estimates of p(eps)
    lposterior_theta = logsumexp(lposterior_theta,4); %sum over epsilon
    lposterior_theta = squeeze(sum(lposterior_theta)); %sum log posteriors over fixations
    lposterior_theta = bsxfun(@minus,lposterior_theta,logsumexp(lposterior_theta,2)); %normalize p(thetas)
    lpthetar(it,:,:) = lposterior_theta;
    
end

% FIND MAP SEQUENCE OF ERRORS
for i = 1:n_iter
    [max_post,max_loc] = max(lpepsr(i+1,:,:),[],3);
    epsr(i,:,:) = SH(max_loc,:);
end

%%
for ii = 1:n_iter
    fullX_sh = zeros(size(fullX));
    for fix = 1:length(used_fixs)
        fprintf('Fixation %d of %d\n',fix,length(used_fixs));
        cur_ims = fix_start_inds(used_fixs(fix)):fix_stop_inds(used_fixs(fix));
        for jj = 1:length(cur_ims);
            cur_im = reshape(fullX(cur_ims(jj),:),sdim,sdim);
            d2 = dist_shift2d(cur_im, epsr(ii,fix,1), 2,0);
            d2 = dist_shift2d(d2,epsr(ii,fix,2),1,0);
            fullX_sh(cur_ims(jj),:) = d2(:);
        end
    end
    
    gabor_emp_X1 = fullX_sh*gabor_emp1_filt';
    gabor_emp_X2 = fullX_sh*gabor_emp2_filt';
    gabor_emp_X_sh = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
    for t = 1:96;
        fprintf('SU %d of %d\n',t,96);
        spikebins = convert_to_spikebins(full_binned_spks(used_ims,t)); 
                
        [fitp,grad] = GLMsolve_jmm(gabor_emp_X_sh(used_ims,t), spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
        beta_sh(ii,t,:) = fitp.k;
        LL_shr(ii,t) = fitp.LL;        
    end
end
LL_impr = bsxfun(@minus,LL_shr,LL);
