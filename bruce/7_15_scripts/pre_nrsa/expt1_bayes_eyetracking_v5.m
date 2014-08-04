%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

cd ~/Data/bruce/7_15_12/G029/
load ./Expt1_newcompiled_data_fixedlag_d1p5_new.mat
fullX = fullX/std(fullX(:));

load ./eye_calibration_data

Pix2Deg = 0.018837;
[NT,klen] = size(fullX);

%% crop stimulus for the purpose of faster gabor function fitting
new_RF_patch = [-0.11 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);
new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
fullX_cropped = fullX(:,new_crop);

xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

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

% SET UP XV CELL SET
NSIG = 96;
xv_frac = 0.2;
tr_set = randperm(NSIG);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NSIG,tr_set);
n_tr_cells = length(tr_set);

% PARSE DATA INTO FIXATIONS

diff_used_inds = diff(used_inds);
rel_fix_start_inds = [1; 1+find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)); NT];
n_fixs = length(rel_fix_start_inds);

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
load ./expt1_eyecor_d1p25_nosac_v4 gabor*
gabor_params = gabor_params_f{end};
clear gabor_emp*
for t = 1:96
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
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
        
        gfun = bsxfun(@times,energy_out,gabor_params(tr_set,7)');
        gfun = bsxfun(@plus,gfun,gabor_params(tr_set,8)');
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
chunk_dur = 2;

%overall prior on shifts
eps_prior_sigma = 0.15; %0.2
leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize

eps_loose_prior_sigma = 0.4; %0.2
leps_loose_prior = -sum((SH/Fsd).^2,2)/(2*eps_loose_prior_sigma^2);
leps_loose_prior = bsxfun(@minus,leps_loose_prior,logsumexp(leps_loose_prior)); %normalize

lA_tflip = repmat(leps_prior',n_shifts,1);

%state transition matrix (includes a 'constant' prior)
deps_sigma = 0.01*chunk_dur; 
cdist = squareform(pdist(SH/Fsd));
lA = -cdist.^2/(2*deps_sigma^2);
lA = bsxfun(@plus,lA,leps_loose_prior'); %factor in constant prior
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

sac_deps_sigma = 0.05;

ov_lgamma = nan(NT,n_shifts);

fprintf('Processing Chunked Data...\n');
chunk_assignments = [];
lB = [];
chunk_labels = [];
chunk_fix_nums = [];
chunk_cnt = 0;
for cf = 1:n_fixs
    cur_im_nums = rel_fix_start_inds(cf):rel_fix_stop_inds(cf);
    
    n_chunks = ceil(length(cur_im_nums)/chunk_dur);
    
    chunk_starts = (0:n_chunks-1)*chunk_dur + 1;
    chunk_stops = chunk_starts + chunk_dur;
    chunk_stops(chunk_stops > length(cur_im_nums)) = length(cur_im_nums);
    
    temp_chunk_assignments = nan(length(cur_im_nums),1);
    temp_lB = nan(n_chunks,n_shifts);
    for i = 1:n_chunks
        temp_chunk_assignments(chunk_starts(i):chunk_stops(i)) = i;
        temp_lB(i,:) = sum(frame_LLs(cur_im_nums(chunk_starts(i):chunk_stops(i)),:));
    end
    chunk_assignments = [chunk_assignments; temp_chunk_assignments+chunk_cnt];
    chunk_fix_nums = [chunk_fix_nums; cf*ones(n_chunks,1)];
    lB = [lB; temp_lB];
    if from_sac(cf) == 0
        chunk_labels = [chunk_labels; 1; zeros(n_chunks-1,1)];
    else
        chunk_labels = [chunk_labels; 2; zeros(n_chunks-1,1)];
    end
    chunk_cnt = chunk_cnt + n_chunks;
end

tot_n_chunks = size(lB,1);
used_sac_dx = sac_amp_left(:,1);
used_sac_dy = sac_amp_left(:,2);

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
        cur_lA = lA;
    elseif chunk_labels(t)==1
        cur_lA = lA_tflip;
    else
        cdist = pdist2(SH/Fsd,bsxfun(@minus,SH/Fsd,[used_sac_dx(chunk_fix_nums(t)) used_sac_dy(chunk_fix_nums(t))]));
%         cur_lA = -cdist.^2/(2*sac_deps_sigmas(chunk_fix_nums(t))^2);
        cur_lA = -cdist.^2/(2*sac_deps_sigma^2);
        cur_lA = bsxfun(@plus,cur_lA,leps_loose_prior'); %factor in constant prior
        cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2)); %normalize
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
        cur_lA = lA;
    elseif chunk_labels(t+1)==1
        cur_lA = lA_tflip;
    else
        cdist = pdist2(SH/Fsd,bsxfun(@minus,SH/Fsd,[used_sac_dx(chunk_fix_nums(t+1)) used_sac_dy(chunk_fix_nums(t+1))]));
%         cur_lA = -cdist.^2/(2*sac_deps_sigmas(chunk_fix_nums(t+1))^2^2);
        cur_lA = -cdist.^2/(2*sac_deps_sigma^2);
        cur_lA = bsxfun(@plus,cur_lA,leps_loose_prior'); %factor in constant prior
        cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2)); %normalize
    end
    lf1 = lbeta(t+1,:) + lB(t+1,:);
    lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
end

%compute posteriors over hidden states
lgamma= lalpha + lbeta;
lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));

ov_lgamma_left = lgamma(chunk_assignments,:);

%% NOW WITH RIGHT EYE
used_sac_dx = sac_amp_right(:,1);
used_sac_dy = sac_amp_right(:,2);

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
        cur_lA = lA;
    elseif chunk_labels(t)==1
        cur_lA = lA_tflip;
    else
        cdist = pdist2(SH/Fsd,bsxfun(@minus,SH/Fsd,[used_sac_dx(chunk_fix_nums(t)) used_sac_dy(chunk_fix_nums(t))]));
%         cur_lA = -cdist.^2/(2*sac_deps_sigmas(chunk_fix_nums(t))^2);
        cur_lA = -cdist.^2/(2*sac_deps_sigma^2);
        cur_lA = bsxfun(@plus,cur_lA,leps_loose_prior'); %factor in constant prior
        cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2)); %normalize
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
        cur_lA = lA;
    elseif chunk_labels(t+1)==1
        cur_lA = lA_tflip;
    else
        cdist = pdist2(SH/Fsd,bsxfun(@minus,SH/Fsd,[used_sac_dx(chunk_fix_nums(t+1)) used_sac_dy(chunk_fix_nums(t+1))]));
%         cur_lA = -cdist.^2/(2*sac_deps_sigmas(chunk_fix_nums(t+1))^2^2);
        cur_lA = -cdist.^2/(2*sac_deps_sigma^2);
        cur_lA = bsxfun(@plus,cur_lA,leps_loose_prior'); %factor in constant prior
        cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2)); %normalize
    end
    lf1 = lbeta(t+1,:) + lB(t+1,:);
    lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
end

%compute posteriors over hidden states
lgamma= lalpha + lbeta;
lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));

ov_lgamma_right = lgamma(chunk_assignments,:);

%% RECONSTRUCT MAP STIMULUS
[max_post_left,max_loc] = max(ov_lgamma_left,[],2);
x_cor_left = SH(max_loc,1);
y_cor_left = SH(max_loc,2);

[max_post_right,max_loc] = max(ov_lgamma_right,[],2);
x_cor_right = SH(max_loc,1);
y_cor_right = SH(max_loc,2);

%%
load ./all_eyedata_expt1_v2
interp_eyepos = interp1(all_t,all_eyepos,full_t);

