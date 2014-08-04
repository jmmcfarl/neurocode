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

cd ~/Data/bruce/7_15_12/G029/
load ./Expt1_newcompiled_data_fixedlag_d1p25_full.mat
fullX = fullX/std(fullX(:));

Pix2Deg = 0.018837;
[NT,klen] = size(fullX);

%% crop stimulus for the purpose of faster gabor function fitting
new_RF_patch = [-0.11 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);

xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

% if length(ypatch_inds_cropped) > length(xpatch_inds_cropped)
%     ypatch_inds_cropped(end) = [];
% elseif length(ypatch_inds_cropped) < length(xpatch_inds_cropped)
%     xpatch_inds_cropped(end) = [];
% end

new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
fullX_cropped = fullX(:,new_crop);

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

% PARSE DATA INTO FIXATIONS
diff_used_inds = [1; diff(used_inds)];
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_fixs = length(rel_fix_start_inds);

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA

gab_priors(1).type = 'gauss';
gab_priors(1).theta(2) = 0.25; %prior std

gab_priors(2).type = 'gauss';
gab_priors(2).theta(2) = 0.25;

gab_priors(4).type = 'gam';
gab_priors(4).theta(1) = 8; %shape
gab_priors(4).theta(2) = 0.03; %scale

gab_priors(5).type = 'gam';
gab_priors(5).theta(1) = 8; %shape
gab_priors(5).theta(2) = 0.011; %scale

gab_priors(6).type = 'gam';
gab_priors(6).theta(1) = 2; %shape
gab_priors(6).theta(2) = 2; %scale

LB = [-0.1 -0.8 0 0.1 0.02 0.2 0 -Inf];
UB = [0.8 0.1 pi 0.4 0.25 6 Inf Inf];
hold_const = [0 0 1 0 0 0 0 0];
all_const = [1 1 1 1 1 1 0 0];

load ./expt1_eyecor_d1p25_nosac_v2.mat
resh_X = reshape(fullX',[sdim sdim NT]);
resh_X_sh  = zeros(size(resh_X));
for ii = 1:NT
    d2 = dist_shift2d(resh_X(:,:,ii), -x_cor{end}(ii), 2,0);
    d2 = dist_shift2d(d2,-y_cor{end}(ii),1,0);
    resh_X_sh(:,:,ii) = d2;
end
fullX_sh_old = reshape(resh_X_sh,sdim^2,NT)';
fullX_cropped = fullX_sh_old(:,new_crop);

for t = 1:96
    gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f{end}(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f{end}(t,1:6),pi/2);
    gabor_emp1_filt_c(t,:) = gabor_emp1(:);
    gabor_emp2_filt_c(t,:) = gabor_emp2(:);

    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params_f{end}(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params_f{end}(t,1:6),pi/2);
    gabor_emp1_filt_f(t,:) = gabor_emp1(:);
    gabor_emp2_filt_f(t,:) = gabor_emp2(:);
end
gabor_out1_c = fullX_cropped*gabor_emp1_filt_c';
gabor_out2_c = fullX_cropped*gabor_emp2_filt_c';
gabor_en_out_c = sqrt(gabor_out1_c.^2 + gabor_out2_c.^2);

%%
lam = 100;
llist = [lam 1 2 3];
for t = 1:96;
    spikebins = convert_to_spikebins(full_binned_spks(:,t));
    X = [gabor_out1_c(:,t) gabor_out2_c(:,t) gabor_en_out_c(:,t)];
    [fitp,grad] = GLMsolve_jmm(X, spikebins, [0; 0], 1, [], [], [], [], llist, [], 0);
    gabor_weights(t,:) = fitp.k(1:3);
    gabor_offset(t) = fitp.k(end);
    LL_full(t) = fitp.LL;
    
    [beta(t,:),dev(t),stats(t)] = glmfit(X,full_binned_spks(:,t),'poisson');
end

%% ESTIMATE LL for each shift in each stimulus frame
max_shift = 25;
dshift = 1;
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;
[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
SH = [Xsh(:) Ysh(:)];
n_shifts = size(SH,1);

frame_LLs = zeros(NT,n_shifts);
Robs = full_binned_spks(:,tr_set);

n_tr_cells = length(tr_set);
gabor_filt_bank1 = reshape(gabor_emp1_filt_f(tr_set,:)',[sdim sdim n_tr_cells]);
gabor_filt_bank2 = reshape(gabor_emp2_filt_f(tr_set,:)',[sdim sdim n_tr_cells]);
shifted_gabor_bank1 = nan(sdim^2,n_tr_cells);
shifted_gabor_bank2 = nan(sdim^2,n_tr_cells);

shift_cnt = 1;
for xx = 1:length(x_shifts)
    for yy = 1:length(y_shifts)
        fprintf('Shift %d of %d\n',shift_cnt,n_shifts);
        d2 = dist_shift3d(gabor_filt_bank1,x_shifts(xx),2);
        d2 = dist_shift3d(d2,y_shifts(yy),1);
        shifted_gabor_bank1 = reshape(d2,sdim^2,n_tr_cells);
        d2 = dist_shift3d(gabor_filt_bank2,x_shifts(xx),2);
        d2 = dist_shift3d(d2,y_shifts(yy),1);
        shifted_gabor_bank2 = reshape(d2,sdim^2,n_tr_cells);
        
        gabor_outs1 = fullX*shifted_gabor_bank1;
        gabor_outs2 = fullX*shifted_gabor_bank2;
        energy_out = sqrt(gabor_outs1.^2 + gabor_outs2.^2);
        
        gfun_e = bsxfun(@times,energy_out,gabor_weights(tr_set,3)');
        gfun_s = bsxfun(@times,gabor_outs1,gabor_weights(tr_set,1)') + ...
            bsxfun(@times,gabor_outs2,gabor_weights(tr_set,2)');
        gfun = bsxfun(@plus,gfun_e+gfun_s,gabor_offset(tr_set));
        too_large = gfun > 50;
        pred_rate = log(1+exp(gfun));
        pred_rate(too_large) = gfun(too_large);
        pred_rate(pred_rate < 1e-20) = 1e-20;
        
        LLs = Robs.*log(pred_rate) - pred_rate;
        frame_LLs(:,shift_cnt) = sum(LLs,2);
        shift_cnt = shift_cnt + 1;
    end
end


%% HMM for inferring sequence of stimulus translations
chunk_dur = 4;

%overall prior on shifts
eps_prior_sigma = 0.3; %0.2
leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
% cur_cent = [ov_eps(cur_fix,1) ov_eps(cur_fix,2)];
% cdist = sum((bsxfun(@minus,SH,cur_cent)/Fsd).^2,2);
% leps_prior = -cdist/(2*eps_prior_sigma^2);
% leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior));

%state transition matrix (includes a 'constant' prior)
deps_sigma = 0.04; %0.06
cdist = squareform(pdist(SH/Fsd));
lA = -cdist.^2/(2*deps_sigma^2);
lA = bsxfun(@plus,lA,leps_prior'); %factor in constant prior
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

ov_lgamma = nan(NT,n_shifts);

for cf = 1:n_fixs
    fprintf('Fixation %d of %d\n',cf,n_fixs);
    cur_im_nums = rel_fix_start_inds(cf):rel_fix_stop_inds(cf);
    
    n_chunks = ceil(length(cur_im_nums)/chunk_dur);
    
    chunk_starts = (0:n_chunks-1)*chunk_dur + 1;
    chunk_stops = chunk_starts + chunk_dur;
    chunk_stops(chunk_stops > length(cur_im_nums)) = length(cur_im_nums);
    
    chunk_assignments = nan(length(cur_im_nums),1);
    lB = nan(n_chunks,n_shifts);
    for i = 1:n_chunks
        chunk_assignments(chunk_starts(i):chunk_stops(i)) = i;
        lB(i,:) = sum(frame_LLs(cur_im_nums(chunk_starts(i):chunk_stops(i)),:));
    end
    
    lalpha=zeros(n_chunks,n_shifts);
    lbeta = zeros(n_chunks,n_shifts);
    lscale=zeros(n_chunks,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior' + lB(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:n_chunks
        lalpha(t,:) = logmulexp(lalpha(t-1,:),lA) + lB(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:)= lalpha(t,:) - lscale(t);
    end
    
    %compute rescaled backward messages
    lbeta(n_chunks,:)=log(ones(1,n_shifts)) - lscale(n_chunks);
    for t=n_chunks-1:-1:1
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,lA') - lscale(t);
    end
    
    %compute posteriors over hidden states
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    ov_lgamma(cur_im_nums,:) = lgamma(chunk_assignments,:);
    
end

%% RECONSTRUCT MAP STIMULUS
[max_post,max_loc] = max(ov_lgamma,[],2);
x_cor_new = SH(max_loc,1);
y_cor_new = SH(max_loc,2);

%% RECONSTRUCT NEW STIMULUS MATRIX
resh_X = reshape(fullX',[sdim sdim NT]);
resh_X_sh_new = zeros(size(resh_X));
for ii = 1:NT
    %     if mod(ii,100)==0 fprintf('%d of %d\n',ii,NT); end
    d2 = dist_shift2d(resh_X(:,:,ii), -x_cor_new(ii), 2,0);
    d2 = dist_shift2d(d2,-y_cor_new(ii),1,0);
    resh_X_sh_new(:,:,ii) = d2;

end

%%
fullX_sh_new = reshape(resh_X_sh_new,sdim^2,NT)';
fullX_cropped_new = fullX_sh_new(:,new_crop);
gabor_out1_new = fullX_cropped_new*gabor_emp1_filt_c';
gabor_out2_new = fullX_cropped_new*gabor_emp2_filt_c';
gabor_en_out_new = sqrt(gabor_out1_new.^2 + gabor_out2_new.^2);

%%
for t = 1:96;
    spikebins = convert_to_spikebins(full_binned_spks(:,t));
    
    X = [gabor_out1_new(:,t) gabor_out2_new(:,t) gabor_en_out_new(:,t)];
    [fitp,grad] = GLMsolve_jmm(X, spikebins, [0; 0], 1, [], [], [], [], llist, [], 0);
    gabor_weights_new(t,:) = fitp.k(1:3);
    gabor_offset_new(t) = fitp.k(end);
    LL_full_new(t) = fitp.LL;
    
end

%%
nfold = 5;
nparts = 60;
[tr_samp,xv_samp] = create_xv_set(NT,nfold,nparts);

tr_samp = tr_samp{1};
xv_samp = xv_samp{1};
%%
gabor_out1 = fullX_cropped*gabor_emp1_filt_c';
gabor_out2 = fullX_cropped*gabor_emp2_filt_c';
gabor_en_out = sqrt(gabor_out1.^2 + gabor_out2.^2);

for t = 1:96;
    cur_binned_spikes_tr = full_binned_spks(tr_samp,t);
    cur_binned_spikes_xv = full_binned_spks(xv_samp,t);
    
    X = [gabor_out1(:,t) gabor_out2(:,t) gabor_en_out(:,t)];
    [beta_sh2(t,:),dev_sh2(t),stats2(t)] = glmfit(X(tr_samp,:),cur_binned_spikes_tr,'poisson');
    pred_r =X(xv_samp,1)*beta_sh2(t,2) + X(xv_samp,2)*beta_sh2(t,3) + X(xv_samp,3)*beta_sh2(t,4) + beta_sh2(t,1);
    pred_r = exp(pred_r);
    LL_sh2(t) = -sum(cur_binned_spikes_xv.*log(pred_r) - pred_r)/sum(cur_binned_spikes_xv);
    
    X = [gabor_en_out(:,t)];
    [beta_sh3(t,:),dev_sh3(t),stats3(t)] = glmfit(X(tr_samp,:),cur_binned_spikes_tr,'poisson');
    pred_r =X(xv_samp,1)*beta_sh3(t,2)+beta_sh3(t,1);
    pred_r = exp(pred_r);
    LL_sh3(t) = -sum(cur_binned_spikes_xv.*log(pred_r) - pred_r)/sum(cur_binned_spikes_xv);
end

%%
save expt1_eyecor_d1p25_nosac_v4 gabor_* xv_set tr_set x_cor_new y_cor_new ov_lgamma
