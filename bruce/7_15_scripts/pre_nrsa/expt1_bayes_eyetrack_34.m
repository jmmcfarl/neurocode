%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/
load ./eye_calibration_data
% load ./G029Expts.mat
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt1_newcompiled_data_fixeddelay_d1p25_34_v2.mat
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

new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
fullX_cropped = fullX(:,new_crop);

sdim = length(xpatch_inds);
[XXc,YYc] = meshgrid(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped));
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

cd ~/Data/bruce/7_15_12/G034/
load ./grating_mu_data
gr_oris = unique_oris;
[~,sf_inds] = max(avg_mu_sf_profile,[],2);
pref_sfs = unique_sfs(sf_inds);
cd ~/Data/bruce/7_15_12/G029/
load ./oned_fixation_fits_v3.mat

dt = dt*2;

% SET UP XV CELL SET OF CELLS
NSIG = 96;
xv_frac = 0.2;
tr_set = randperm(NSIG);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NSIG,tr_set);
n_tr_cells = length(tr_set);

% PARSE DATA INTO FIXATIONS
diff_used_inds = [1; diff(used_inds)];
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_fixs = length(rel_fix_start_inds);

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA

[~,~,full_expt_inds] = unique(full_expt_vec); 

gab_priors(1).type = 'gauss';
gab_priors(1).theta(2) = 0.3; %prior std

gab_priors(2).type = 'gauss';
gab_priors(2).theta(2) = 0.3;

gab_priors(4).type = 'gam';
gab_priors(4).theta(1) = 8; %shape
gab_priors(4).theta(2) = 0.03; %scale

gab_priors(5).type = 'gam';
gab_priors(5).theta(1) = 8; %shape
gab_priors(5).theta(2) = 0.011; %scale

gab_priors(6).type = 'gam';
gab_priors(6).theta(1) = 2; %shape
gab_priors(6).theta(2) = 2; %scale

LB = [-0.1 -0.8 0 0.125 0.025 0.2 0];
UB = [0.8 0.1 pi 0.4 0.25 6 Inf];
hold_const = [0 0 1 0 0 0 0];
all_const = [1 1 1 1 1 1 0];

for t = 1:96
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    
    gab_priors(1).theta(1) = mean_x(t);
    gab_priors(2).theta(1) = mean_y(t);

    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.4*init_params(4); %sigma
    init_params(6) = 2; %eccentricity of gaussian ellipse
    init_params(7) = 0; %weight of quad term
%     init_params(8) = 0; %const offset
    
    [gabor_params_f{1}(t,:),LL(1,t)] = fit_gabor_energy_mod_varmeanrate(XXc,YYc,init_params,fullX_cropped,...
        full_binned_spks(:,t),full_expt_inds,hold_const,LB,UB,gab_priors);
    
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params_f{1}(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params_f{1}(t,1:6),pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
    
end

%% SAVE THESE INITIAL GABOR MODELS
cd ~/Data/bruce/7_15_12/G034/
save gabor_initfits_d1p25_varrate gabor_params_* LL gabor_emp*_filt gab_priors hold_const all_const

%%
% [un_expts,~,full_expt_inds] = unique(full_expt_vec);

load ./all_eyedata_expt1_34
orig_eyedt = all_t(2)-all_t(1);
interp_eyespeeds = interp1(all_t,all_eyespeed,full_t);
avg_interp_eyespeed = interp_eyespeeds;
interp_eyepos = interp1(all_t,all_eyepos,full_t);

max_shift = 22;
dshift = 1;
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;
[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
SH = [Xsh(:) Ysh(:)];
n_shifts = size(SH,1);

chunk_dur = 1;

%overall prior on shifts
eps_prior_sigma = 0.2; %0.2
leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior',n_shifts,1);

%state transition matrix (includes a 'constant' prior)
min_deps_sigma = 0.005*chunk_dur; 
max_deps_sigma = 0.4;
deps_offset = 5;
deps_slope = 2;
sac_deps_sigmas = (max_deps_sigma-min_deps_sigma)./(1+exp(-(avg_interp_eyespeed-deps_offset)*deps_slope))+min_deps_sigma;

cdist = squareform(pdist(SH/Fsd));
% lA = -cdist.^2/(2*deps_sigma^2);
% lA = bsxfun(@plus,lA,leps_prior'); %factor in constant prior
% lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

ov_lgamma = nan(NT,n_shifts);

%%
n_iter = 3;
for it = 1:n_iter
    
    %% ESTIMATE LL for each shift in each stimulus frame    
    frame_LLs = zeros(NT,n_shifts);
    Robs = full_binned_spks(:,tr_set);
    
    gabor_filt_bank1 = reshape(gabor_emp1_filt(tr_set,:)',[sdim sdim n_tr_cells]);
    gabor_filt_bank2 = reshape(gabor_emp2_filt(tr_set,:)',[sdim sdim n_tr_cells]);
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
            
            gfun = bsxfun(@times,energy_out,gabor_params_f{it}(tr_set,7)');
            gfun = gfun + gabor_params_f{it}(tr_set,full_expt_inds+7)';
             
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
    fprintf('Processing Chunked Data...\n');
    chunk_assignments = [];
    lB = [];
    chunk_labels = [];
    chunk_fix_nums = [];
    chunk_eyespeed = [];
    chunk_cnt = 0;
    for cf = 1:n_fixs
        cur_im_nums = rel_fix_start_inds(cf):rel_fix_stop_inds(cf);        
        n_chunks = ceil(length(cur_im_nums)/chunk_dur);
        chunk_starts = (0:n_chunks-1)*chunk_dur + 1;
        chunk_stops = chunk_starts + chunk_dur-1;
        chunk_stops(chunk_stops > length(cur_im_nums)) = length(cur_im_nums);
        
        temp_chunk_assignments = nan(length(cur_im_nums),1);
        temp_lB = nan(n_chunks,n_shifts);
        temp_chunk_eyespeeds = nan(n_chunks,1);
        for i = 1:n_chunks
            temp_chunk_assignments(chunk_starts(i):chunk_stops(i)) = i;
            temp_lB(i,:) = sum(frame_LLs(cur_im_nums(chunk_starts(i):chunk_stops(i)),:),1);
            temp_chunk_eyespeeds(i) = mean(avg_interp_eyespeed(cur_im_nums(chunk_starts(i):chunk_stops(i))));
        end
        chunk_assignments = [chunk_assignments; temp_chunk_assignments+chunk_cnt];
        chunk_fix_nums = [chunk_fix_nums; cf*ones(n_chunks,1)];
        chunk_eyespeed = [chunk_eyespeed; temp_chunk_eyespeeds];
        lB = [lB; temp_lB];
        chunk_labels = [chunk_labels; 1; zeros(n_chunks-1,1)];
        chunk_cnt = chunk_cnt + n_chunks;
    end
    
    tot_n_chunks = size(lB,1);

    lalpha=zeros(tot_n_chunks,n_shifts);
    lbeta = zeros(tot_n_chunks,n_shifts);
    lscale=zeros(tot_n_chunks,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior' + lB(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:tot_n_chunks
        fprintf('%d of %d\n',t,tot_n_chunks);
        if chunk_labels(t)==0
            cur_lA = -cdist.^2/(2*sac_deps_sigmas(t-1)^2);
%             cur_lA = bsxfun(@plus,cur_lA,leps_prior'); %factor in constant prior
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2)); %normalize
        elseif chunk_labels(t)==1
            cur_lA = lA_tflip;
        end
        lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:)= lalpha(t,:) - lscale(t);
    end
    
    %compute rescaled backward messages
    lbeta(tot_n_chunks,:)=log(ones(1,n_shifts)) - lscale(tot_n_chunks);
    for t=tot_n_chunks-1:-1:1
        fprintf('%d\n',t);
        if chunk_labels(t+1)==0
            cur_lA = -cdist.^2/(2*sac_deps_sigmas(t)^2);
%             cur_lA = bsxfun(@plus,cur_lA,leps_prior'); %factor in constant prior
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2)); %normalize
        elseif chunk_labels(t+1)==1
            cur_lA = lA_tflip;
         end
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
    end

    %compute posteriors over hidden states
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    ov_lgamma = lgamma(chunk_assignments,:);
    
    %% RECONSTRUCT MAP STIMULUS
    [max_post,max_loc] = max(ov_lgamma,[],2);
    x_cor{it} = SH(max_loc,1);
    y_cor{it} = SH(max_loc,2);
    
    %% RECONSTRUCT NEW STIMULUS MATRIX
    resh_X = reshape(fullX',[sdim sdim NT]);
    resh_X_sh = zeros(size(resh_X));
    for ii = 1:NT
        %     if mod(ii,100)==0 fprintf('%d of %d\n',ii,NT); end
        d2 = dist_shift2d(resh_X(:,:,ii), -x_cor{it}(ii), 2,0);
        d2 = dist_shift2d(d2,-y_cor{it}(ii),1,0);
        resh_X_sh(:,:,ii) = d2;
    end
    
    %% REFIT GEM PARAMETERS
    fullX_sh = reshape(resh_X_sh,sdim^2,NT)';
    fullX_cropped = fullX_sh(:,new_crop);
        
    for t = 1:96
        fprintf('Fitting GEM: Cell %d of %d\n',t,96);
        gab_priors(1).theta(1) = gabor_params_f{it}(t,1);
        gab_priors(2).theta(1) = gabor_params_f{it}(t,2);
        
        init_params = gabor_params_f{it}(t,:);
        [gabor_params_f{it+1}(t,:),LL(it+1,t)] = fit_gabor_energy_mod_varmeanrate(XXc,YYc,init_params,fullX_cropped,...
            full_binned_spks(:,t),full_expt_inds,hold_const,LB,UB,gab_priors);
        
        gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params_f{it+1}(t,1:6),0);
        gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params_f{it+1}(t,1:6),pi/2);
        
        gabor_emp1_filt(t,:) = gabor_emp1(:);
        gabor_emp2_filt(t,:) = gabor_emp2(:);
    end
end

%%
cd ~/Data/bruce/7_15_12/G034/
save gabor_tracking_varmeans_v2 gabor_* LL xv_set tr_set x_cor y_cor 

