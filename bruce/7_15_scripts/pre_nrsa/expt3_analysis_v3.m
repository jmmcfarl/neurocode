%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

cd ~/Data/bruce/7_15_12/G029/
load ./Expt3_newcompiled_data_d1p25_htres.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

load ./all_eyedata_expt3
Pix2Deg = 0.018837;
NT = length(full_stim_ids);

%% crop stimulus for the purpose of faster gabor function fitting
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

%%
load ./expt1_eyecor_d1p25_nosac_v2.mat gabor*
gabor_params = gabor_params_f{end};
clear gabor_*filt
for t = 1:96
    
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);    
end

idier = [full_expt_vec full_trial_vec];
unique_ids = unique(idier,'rows');
n_trials = size(unique_ids,1);
xv_frac =0;
xv_num = round(xv_frac*n_trials);
xv_set = randperm(n_trials);
xv_set(xv_num+1:end) = [];
tr_set = setdiff(1:n_trials,xv_set);
xv_inds = [];
for ii = 1:xv_num
   cur_inds = find(full_expt_vec == unique_ids(xv_set(ii),1) & full_trial_vec == unique_ids(xv_set(ii),2));
   xv_inds = [xv_inds; cur_inds];
end
tr_inds = setdiff(1:NT,xv_inds);

%%
% PARSE DATA INTO FIXATIONS
diff_used_inds = [1; diff(used_inds)];
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_trials = length(rel_fix_start_inds);

flen = 88;

im_flip_inds = [1 (1+find(diff(full_stim_ids)~=0)) NT];
im_flip_log = zeros(size(full_t));
im_flip_log(im_flip_inds) = 1;
time_since_imflip = zeros(size(full_t));
for i = 1:length(im_flip_inds)-1
    cur_set = im_flip_inds(i):(im_flip_inds(i+1)-1);
    time_since_imflip(cur_set) = (0:dt:(dt*length(cur_set)-dt));
end
im_flip_inds(end) = [];

im_tax = [0:dt:0.1];
cur_sp = dt;
while max(im_tax) < max(time_since_imflip)
    im_tax = [im_tax (im_tax(end)+cur_sp)];
    cur_sp = cur_sp + dt/3;
end
if im_tax(end) >= max(time_since_imflip)
    im_tax(end) = [];
end
nbins = length(im_tax);
imflip_Tmat = tbrep(time_since_imflip,im_tax);

sac_start_inds = [(1+find(diff(full_insac)>0))];
sac_stop_inds = (1+find(diff(full_insac) < 0));
time_since_sac = Inf(size(full_t));
for i = 1:length(sac_start_inds)-1
    next_trial = rel_fix_start_inds(find(rel_fix_start_inds > sac_stop_inds(i),1));
    next_sac = sac_start_inds(find(sac_start_inds > sac_stop_inds(i),1));
    cur_set = sac_stop_inds(i):min(next_sac,next_trial);
    time_since_sac(cur_set) = (0:dt:(dt*length(cur_set)-dt));
end
sac_Tmat = tbrep(time_since_sac,im_tax);

flen = nbins;

% ov_im_X = [];
% ov_sac_X = [];
% for i = 1:n_trials
%     cur_inds = rel_fix_start_inds(i):rel_fix_stop_inds(i);
%     cur_im_flips = im_flip_log(cur_inds);
%     cur_im_X = makeStimRows(cur_im_flips(:),flen);
%     ov_im_X = [ov_im_X; cur_im_X];
% 
%     cur_in_sac = full_insac(cur_inds);
%     cur_sac_X = makeStimRows(cur_in_sac(:),flen);
%     ov_sac_X = [ov_sac_X; cur_sac_X];
% end

[un_expts,~,full_expt_inds] = unique(full_expt_vec);
n_un_expts = length(un_expts);
linX = zeros(NT,n_un_expts-1);
for i = 1:n_un_expts-1
    linX(full_expt_inds==i,i) = 1;
end

load ./all_eyedata_expt3.mat
orig_eyedt = all_t(2)-all_t(1);
interp_eyespeeds = interp1(all_t,all_eyespeed,full_t);
interp_eyepos = interp1(all_t,all_eyepos,full_t);
avg_interp_eyespeed = dt/orig_eyedt*mean(interp_eyespeeds,2);

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
all_gabor_out1 = resh_all_stims*gabor_emp1_filt';
all_gabor_out2 = resh_all_stims*gabor_emp2_filt';
spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);

l2 = 100;
silent = 1;

clear stim_dep_kern stim_ind_kern sac_kern sac_stimkern 
for t = 1:96
    fprintf('Fitting Cell %d of 96\n',t);
    cur_spatial_mod_out = spatial_mod_out(full_stim_ids,t);
    stim_dep_X = bsxfun(@times,imflip_Tmat,cur_spatial_mod_out);
%     stim_dep_X = bsxfun(@times,ov_im_X,cur_spatial_mod_out);
    stim_dep_sacX = bsxfun(@times,sac_Tmat,cur_spatial_mod_out);
    
    ov_Xmat = [stim_dep_X imflip_Tmat stim_dep_sacX sac_Tmat linX];
%     ov_Xmat = [stim_dep_X ov_im_X linX];
%     ov_Xmat = [stim_dep_X ov_im_X ov_sac_X linX];
% %     ov_Xmat = [stim_dep_X ov_im_X stim_dep_sacX ov_sac_X linX];
    
    Robs = full_binned_spks(tr_inds,t);
    Robs_xv = full_binned_spks(xv_inds,t);

    poss_spkcnts = unique(Robs(Robs > 0));
    spksN = [];
    for i = 1:length(poss_spkcnts)
        cur_set = find(Robs == poss_spkcnts(i));
        spksN = [spksN; repmat(cur_set,poss_spkcnts(i),1)];
    end
    spksN = sort(spksN);
    
    klen = size(ov_Xmat,2);
    K0 = zeros(klen,1);
%     lamrange2 = [l2 1 flen;l2 (flen+1) 2*flen];
%     lamrange2 = [l2 1 flen;l2 (flen+1) 2*flen;l2 (2*flen+1) 3*flen];
%     lamrange2 = [l2 1 flen;l2 (flen+1) 2*flen;l2 (2*flen+1) 3*flen];
    lamrange2 = [l2 1 flen;l2 (flen+1) 2*flen;l2 (2*flen+1) 3*flen;l2 (3*flen+1) 4*flen];
    llist = [];
    [fitp,grad] = GLMsolve_jmm(ov_Xmat(tr_inds,:), spksN, K0, silent, [], lamrange2, [], [],llist, [], 0);
    stim_dep_kern(t,:) = fitp.k(1:flen);
    stim_ind_kern(t,:) = fitp.k((flen+1):2*flen);
%     sac_kern(t,:) = fitp.k((2*flen+1):3*flen);
    sac_stimkern(t,:) = fitp.k((2*flen+1):3*flen);
    sac_kern(t,:) = fitp.k((3*flen+1):4*flen);
%     block_kern(t,:) = fitp.k((2*flen+1):end-1);
%     block_kern(t,:) = fitp.k((3*flen+1):end-1);
    block_kern(t,:) = fitp.k((4*flen+1):end-1);
    ov_const(t) = fitp.k(end);
    %     LL(t) = fitp.LL;
    gfun = ov_Xmat(tr_inds,:)*fitp.k(1:end-1) + fitp.k(end);
    too_large = find(gfun > 50);
    predrate = log(1+exp(gfun));
    predrate(too_large) = gfun(too_large);
    LL(t) = -sum(Robs.*log(predrate)-predrate)/sum(Robs);
    
    if xv_frac > 0
        gfun = ov_Xmat(xv_inds,:)*fitp.k(1:end-1) + fitp.k(end);
        too_large = find(gfun > 50);
        xv_predrate = log(1+exp(gfun));
        xv_predrate(too_large) = gfun(too_large);
        xv_LL(t) = -sum(Robs_xv.*log(xv_predrate)-xv_predrate)/sum(Robs_xv);
    end
end

%%
for t = single_units
    subplot(2,1,1)
    plot(im_tax,stim_dep_kern(t,:))
    hold on
    plot(im_tax,sac_stimkern(t,:),'r')
    xlim([0 0.51])
    subplot(2,1,2)
    plot(im_tax,stim_ind_kern(t,:))
    hold on
    plot(im_tax,sac_kern(t,:),'r')
    xlim([0 0.51])
    t
    pause
    clf
end

%%
% avg_eyespeed = mean(all_eyespeed,2);
% int_eyespeed = interp1(all_t,avg_eyespeed,full_t);
% 
% for t = 1:96
%     obs_rate = full_binned_spks(:,t);
%     cur_spatial_mod_out = spatial_mod_out(full_stim_ids,t);
%     stim_dep_X = bsxfun(@times,ov_im_X,cur_spatial_mod_out);
%     ov_Xmat = [stim_dep_X ov_im_X linX];
%     
%     fitK = [stim_dep_kern(t,:)  stim_ind_kern(t,:)  block_kern(t,:)];
%     pred_rate = ov_Xmat*fitK'+ov_const(t);
%     pred_rate = log(1+exp(pred_rate));
%     
%     figure
%     plot(full_t,smooth(obs_rate,5))
%     hold on
%     plot(full_t(im_flip_inds),ones(size(im_flip_inds)),'ko','linewidth',3)
%     plot(full_t(rel_fix_start_inds),ones(size(rel_fix_start_inds)),'go','linewidth',3)
%     plot(full_t,10*int_eyespeed,'c')
%     
%     plot(full_t,smooth(mean(full_binned_spks,2),5)*4,'k')
%     
%     plot(full_t,pred_rate,'r')
%     pause
%     close all
% end

%% SET UP XV CELL SET
NSIG = 96;
xv_frac = 0.2;
tr_set = randperm(NSIG);
tr_set = sort(tr_set(1:round(length(tr_set)*(1-xv_frac))));
xv_set = setdiff(1:NSIG,tr_set);
n_tr_cells = length(tr_set);

%%
max_shift = 22;
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

stim_dep_im_X = zeros(NT,length(tr_set));
stim_dep_sac_X = zeros(NT,length(tr_set));

shift_cnt = 1;
for xx = 1:length(x_shifts)
    for yy = 1:length(y_shifts)
        
        fprintf('Shift %d of %d\n',shift_cnt,n_shifts);
        d2 = dist_shift3d(gabor_filt_bank1,x_shifts(xx),2);
        d2 = dist_shift3d(d2,y_shifts(yy),2);
        shifted_gabor_bank1 = reshape(d2,sdim^2,n_tr_cells);
        d2 = dist_shift3d(gabor_filt_bank2,x_shifts(xx),1);
        d2 = dist_shift3d(d2,y_shifts(yy),2);
        shifted_gabor_bank2 = reshape(d2,sdim^2,n_tr_cells);
        
        all_gabor_out1 = resh_all_stims*shifted_gabor_bank1;
        all_gabor_out2 = resh_all_stims*shifted_gabor_bank2;
        spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);
        
        cur_spatial_mod_out = spatial_mod_out(full_stim_ids,:);
        for t = 1:length(tr_set)
            stim_dep_im_X(:,t) = bsxfun(@times,imflip_Tmat,cur_spatial_mod_out(:,t))*stim_dep_kern(tr_set(t),:)';
            stim_dep_sac_X(:,t) = bsxfun(@times,sac_Tmat,cur_spatial_mod_out(:,t))*sac_stimkern(tr_set(t),:)';
        end
        
        gfun = stim_dep_im_X + stim_dep_sac_X + imflip_Tmat*stim_ind_kern(tr_set,:)' + ...
            sac_Tmat*sac_kern(tr_set,:)' + linX*block_kern(tr_set,:)';
        gfun = bsxfun(@plus,gfun,ov_const(tr_set));
        
        too_large = gfun > 50;
        pred_rate = log(1+exp(gfun));
        pred_rate(too_large) = gfun(too_large);
        pred_rate(pred_rate < 1e-20) = 1e-20;
        
        LLs = Robs.*log(pred_rate) - pred_rate;
        frame_LLs(:,shift_cnt) = sum(LLs,2);
        shift_cnt = shift_cnt + 1;
    end
end

%%
save temp_data frame_LLs

%% HMM for inferring sequence of stimulus translations
chunk_dur = 88;

%overall prior on shifts
eps_prior_sigma = 0.15; %0.2
leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior',n_shifts,1);

deps_sigma = 0.05; %0.06
cdist = squareform(pdist(SH/Fsd));
lA = -cdist.^2/(2*deps_sigma^2);
lA = bsxfun(@plus,lA,leps_prior'); %factor in constant prior
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

fprintf('Processing Chunked Data...\n');
chunk_assignments = [];
lB = [];
chunk_labels = [];
chunk_fix_nums = [];
chunk_eyespeed = [];
chunk_cnt = 0;
chunk_inds = [];

min_seg_dur = 10;

all_im_starts = 1+find(diff(full_stim_ids) ~= 0);
all_im_stops = all_im_starts-1;
% all_seg_starts = sort(unique([all_im_starts'; rel_fix_start_inds]));
% all_seg_starts = sort([all_im_starts'; rel_fix_start_inds; sac_stop_inds']);
% all_seg_stops = sort([all_im_stops'; rel_fix_stop_inds; sac_start_inds']);
all_seg_starts = sort([1; all_im_starts'; sac_stop_inds']);
all_seg_stops = sort([all_im_stops'; sac_start_inds'; NT]);
% all_seg_starts = rel_fix_start_inds;
% all_seg_stops = rel_fix_stop_inds;

seg_durs = all_seg_stops-all_seg_starts;
bad_segs = find(seg_durs < min_seg_dur);
all_seg_starts(bad_segs) = [];
all_seg_stops(bad_segs) = [];

n_fixs = length(all_seg_stops);
for cf = 1:n_fixs
    cur_im_nums = all_seg_starts(cf):all_seg_stops(cf);
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
    chunk_inds = [chunk_inds; cur_im_nums(:)];
    chunk_cnt = chunk_cnt + n_chunks;
end

tot_n_chunks = size(lB,1);

lPost = bsxfun(@plus,lB,leps_prior');
lPost = bsxfun(@minus,lPost,logsumexp(lPost,2));
ov_lgamma = lPost(chunk_assignments,:);

% lalpha=zeros(tot_n_chunks,n_shifts);
% lbeta = zeros(tot_n_chunks,n_shifts);
% lscale=zeros(tot_n_chunks,1); %initialize rescaling parameters
% %compute rescaled forward messages
% lalpha(1,:) = leps_prior' + lB(1,:);
% lscale(1)=logsumexp(lalpha(1,:));
% lalpha(1,:) = lalpha(1,:) - lscale(1);
% for t=2:tot_n_chunks
%     fprintf('%d of %d\n',t,tot_n_chunks);
% %     if chunk_labels(t) == 1
%         cur_lA = lA_tflip;
% %     else
% %         cur_lA = lA;
% %     end
%     lalpha(t,:) = logmulexp(lalpha(t-1,:),lA) + lB(t,:);
%     lscale(t) = logsumexp(lalpha(t,:));
%     lalpha(t,:)= lalpha(t,:) - lscale(t);
% end
% 
% %compute rescaled backward messages
% lbeta(tot_n_chunks,:)=log(ones(1,n_shifts)) - lscale(tot_n_chunks);
% for t=tot_n_chunks-1:-1:1
%     fprintf('%d\n',t);
% %     if chunk_labels(t+1)==1
%         cur_lA = lA_tflip;
% %     else
% %         cur_lA = lA;
% %     end
%     lf1 = lbeta(t+1,:) + lB(t+1,:);
%     lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
% end
% 
% %compute posteriors over hidden states
% lgamma= lalpha + lbeta;
% lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
% 
% ov_lgamma = lgamma(chunk_assignments,:);

%% RECONSTRUCT MAP STIMULUS
[max_post,max_loc] = max(ov_lgamma,[],2);
x_cor = SH(max_loc,1);
y_cor = SH(max_loc,2);

full_xcor = zeros(size(full_t));
full_ycor = zeros(size(full_t));
full_xcor(chunk_inds) = x_cor;
full_ycor(chunk_inds) = y_cor;

n_stims = size(resh_all_stims,1);
resh_X = reshape(resh_all_stims',[sdim sdim n_stims]);
resh_X_sh = zeros(sdim,sdim,NT);
n_chunk_inds = length(chunk_inds);
chunk_log = zeros(size(full_t));
chunk_log(chunk_inds) = chunk_inds;
for ii = 1:NT
    if chunk_log(ii) ~= 0
        d2 = dist_shift2d(resh_X(:,:,full_stim_ids(ii)), -full_xcor(ii), 2,0);
        d2 = dist_shift2d(d2,-full_ycor(ii),1,0);
        resh_X_sh(:,:,ii) = d2;
    else
        resh_X_sh(:,:,ii) = resh_X(:,:,full_stim_ids(ii));
    end
end
fullX_sh = reshape(resh_X_sh,sdim^2,NT)';
 
%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
l2 = 100;
silent = 1;

all_gabor_out1 = fullX_sh*gabor_emp1_filt';
all_gabor_out2 = fullX_sh*gabor_emp2_filt';
spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);

for t = 1:96
    fprintf('Fitting Cell %d of 96\n',t);
    cur_spatial_mod_out = spatial_mod_out(:,t);
    stim_dep_X = bsxfun(@times,imflip_Tmat,cur_spatial_mod_out);
%     stim_dep_X = bsxfun(@times,ov_im_X,cur_spatial_mod_out);
    stim_dep_sacX = bsxfun(@times,sac_Tmat,cur_spatial_mod_out);
    
    ov_Xmat = [stim_dep_X imflip_Tmat stim_dep_sacX sac_Tmat linX];
    
    Robs = full_binned_spks(:,t);
    poss_spkcnts = unique(Robs(Robs > 0));
    spksN = [];
    for i = 1:length(poss_spkcnts)
        cur_set = find(full_binned_spks(:,t) == poss_spkcnts(i));
        spksN = [spksN; repmat(cur_set,poss_spkcnts(i),1)];
    end
    spksN = sort(spksN);
    
    klen = size(ov_Xmat,2);
    K0 = zeros(klen,1);
%     lamrange2 = [l2 1 flen;l2 (flen+1) 2*flen];
    lamrange2 = [l2 1 flen;l2 (flen+1) 2*flen;l2 (2*flen+1) 3*flen;l2 (3*flen+1) 4*flen];
    llist = [];
    [fitp,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],...
        llist, [], 0);
    stim_dep_kern2(t,:) = fitp.k(1:flen);
    stim_ind_kern2(t,:) = fitp.k((flen+1):2*flen);
    sac_stimkern2(t,:) = fitp.k((2*flen+1):3*flen);
    sac_kern2(t,:) = fitp.k((3*flen+1):4*flen);
    block_kern2(t,:) = fitp.k((4*flen+1):end-1);
    ov_const2(t) = fitp.k(end);
%     LL2(t) = fitp.LL;
    gfun = ov_Xmat(tr_inds,:)*fitp.k(1:end-1) + fitp.k(end);
    too_large = find(gfun > 50);
    predrate = log(1+exp(gfun));
    predrate(too_large) = gfun(too_large);
    LL2(t) = -sum(Robs.*log(predrate)-predrate)/sum(Robs);    
    
end
