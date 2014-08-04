%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

cd ~/Data/bruce/7_15_12/G029/
load ./Expt2_compiled_data_d1p5.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

Pix2Deg = 0.018837;
NT = length(full_fix_ids);
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

infix_ids = find(~isnan(full_fix_ids));
n_fixs = length(unique(full_fix_ids(infix_ids)));
xv_frac =0.;
xv_num = round(xv_frac*n_fixs);
xv_set = randperm(n_fixs);
xv_set(xv_num+1:end) = [];
tr_set = setdiff(1:n_fixs,xv_set);
xv_inds = find(ismember(full_fix_ids,xv_set));
tr_inds = setdiff(find(~isnan(full_fix_ids)),xv_inds);

%%
flen = 80;
flen_t = flen*dtd;

tent_centers = [0:dtd:0.1];
cur_sp = dtd;
while max(tent_centers) < flen_t
    tent_centers = [tent_centers (tent_centers(end)+cur_sp)];
    cur_sp = cur_sp + dtd/4;
end
tent_centers = round(tent_centers/dtd)*dtd;
if tent_centers(end) >= flen_t
    tent_centers(end) = [];
end
tbmat = construct_tent_bases(tent_centers,dtd);
[ntents,tblen] = size(tbmat);

tbmat = [zeros(ntents,tblen-1) tbmat];

sac_inds = zeros(size(full_t));
sac_inds(full_sac_inds) = 1;

sac_Tmat = zeros(length(full_t),ntents);
for i = 1:ntents
    sac_Tmat(:,i) = conv(sac_inds,tbmat(i,:),'same');
end

[un_expts,~,full_expt_inds] = unique(full_expt_vec);
n_un_expts = length(un_expts);
linX = zeros(NT,n_un_expts-1);
for i = 1:n_un_expts-1
    linX(full_expt_inds==i,i) = 1;
end

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
all_gabor_out1 = resh_all_stims*gabor_emp1_filt';
all_gabor_out2 = resh_all_stims*gabor_emp2_filt';
spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);

l2_ind = 100;
l2_dep = 500;
silent = 1;

clear stim_dep_kern stim_ind_kern
for t = 1:96
    fprintf('Fitting Cell %d of 96\n',t);
    cur_spatial_mod_out = spatial_mod_out(full_fix_ids(tr_inds),t);
    stim_dep_X = bsxfun(@times,sac_Tmat(tr_inds,:),cur_spatial_mod_out);    
    ov_Xmat = [stim_dep_X sac_Tmat(tr_inds,:) linX(tr_inds,:)];
    
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
    lamrange2 = [l2_dep 1 ntents;l2_ind (ntents+1) 2*ntents];
%     lamrange2 = [l2 1 flen;l2 (flen+1) 2*flen;l2 (2*flen+1) 3*flen];
%     lamrange2 = [l2 1 flen;l2 (flen+1) 2*flen;l2 (2*flen+1) 3*flen];
%     lamrange2 = [l2 1 flen;l2 (flen+1) 2*flen;l2 (2*flen+1) 3*flen;l2 (3*flen+1) 4*flen];
    llist = [];
    [fitp,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], 0);
    stim_dep_kern(t,:) = fitp.k(1:ntents);
    stim_ind_kern(t,:) = fitp.k((ntents+1):2*ntents);
    block_kern(t,:) = fitp.k((2*ntents+1):end-1);
    ov_const(t) = fitp.k(end);

    gfun = ov_Xmat*fitp.k(1:end-1) + fitp.k(end);
    too_large = find(gfun > 50);
    predrate = log(1+exp(gfun));
    predrate(too_large) = gfun(too_large);
    LL(t) = -sum(Robs.*log(predrate)-predrate)/sum(Robs);
    
    if xv_frac > 0
        cur_spatial_mod_out = spatial_mod_out(full_fix_ids(xv_inds),t);
        stim_dep_X = bsxfun(@times,sac_Tmat(xv_inds,:),cur_spatial_mod_out);
        ov_Xmat = [stim_dep_X sac_Tmat(xv_inds,:) linX(xv_inds,:)];
        
        gfun = ov_Xmat*fitp.k(1:end-1) + fitp.k(end);
        too_large = find(gfun > 50);
        xv_predrate = log(1+exp(gfun));
        xv_predrate(too_large) = gfun(too_large);
        xv_LL(t) = -sum(Robs_xv.*log(xv_predrate)-xv_predrate)/sum(Robs_xv);
    end
end

%%
for t = 1:96
    ax = plotyy(tent_centers,stim_dep_kern(t,:),tent_centers,stim_ind_kern(t,:));
    set(ax(1),'xlim',[0 0.5])
    set(ax(2),'xlim',[0 0.5])
    t
    pause
    clf
end

%% SET UP XV CELL SET
NSIG = 96;
xv_frac = 0.2;
tr_set = randperm(NSIG);
tr_set = sort(tr_set(1:round(length(tr_set)*(1-xv_frac))));
xv_set = setdiff(1:NSIG,tr_set);
n_tr_cells = length(tr_set);

%%
NT = length(tr_inds);
max_shift = 35;
dshift = 2;
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;
[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
SH = [Xsh(:) Ysh(:)];
n_shifts = size(SH,1);

frame_LLs = zeros(NT,n_shifts);
Robs = full_binned_spks(tr_inds,tr_set);

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
        
        cur_spatial_mod_out = spatial_mod_out(full_fix_ids(tr_inds),:);
        for t = 1:length(tr_set)
            stim_dep_im_X(:,t) = bsxfun(@times,sac_Tmat(tr_inds,:),cur_spatial_mod_out(:,t))*stim_dep_kern(tr_set(t),:)';
        end
        
        gfun = stim_dep_im_X + sac_Tmat(tr_inds,:)*stim_ind_kern(tr_set,:)' + linX(tr_inds,:)*block_kern(tr_set,:)';
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

%overall prior on shifts
eps_prior_sigma = 0.3; %0.2
leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize

used_fixnums = full_fix_ids(tr_inds);

lB = nan(n_fixs,n_shifts);
for cf = 1:n_fixs
    cur_inds = find(used_fixnums==cf);
    lB(cf,:) = sum(frame_LLs(cur_inds,:),1);
end

lPost = bsxfun(@plus,lB,leps_prior');
lPost = bsxfun(@minus,lPost,logsumexp(lPost,2));


%% RECONSTRUCT MAP STIMULUS
[max_post,max_loc] = max(lPost,[],2);
x_cor = SH(max_loc,1);
y_cor = SH(max_loc,2);

full_xcor = x_cor(used_fixnums);
full_ycor = y_cor(used_fixnums);

resh_X = reshape(resh_all_stims',[sdim sdim n_fixs]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:n_fixs
    d2 = dist_shift2d(resh_X(:,:,ii), -full_xcor(ii), 2,0);
    d2 = dist_shift2d(d2,-full_ycor(ii),1,0);
    resh_X_sh(:,:,ii) = d2;
end
fullX_sh = reshape(resh_X_sh,sdim^2,n_fixs)';
 
%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA

all_gabor_out1 = fullX_sh*gabor_emp1_filt';
all_gabor_out2 = fullX_sh*gabor_emp2_filt';
spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);

for t = 1:96
    fprintf('Fitting Cell %d of 96\n',t);
    cur_spatial_mod_out = spatial_mod_out(full_fix_ids(tr_inds),t);
    stim_dep_X = bsxfun(@times,sac_Tmat(tr_inds,:),cur_spatial_mod_out);    
    ov_Xmat = [stim_dep_X sac_Tmat(tr_inds,:) linX(tr_inds,:)];
    
    Robs = full_binned_spks(tr_inds,t);

    poss_spkcnts = unique(Robs(Robs > 0));
    spksN = [];
    for i = 1:length(poss_spkcnts)
        cur_set = find(Robs == poss_spkcnts(i));
        spksN = [spksN; repmat(cur_set,poss_spkcnts(i),1)];
    end
    spksN = sort(spksN);
    
    klen = size(ov_Xmat,2);
    K0 = zeros(klen,1);
    lamrange2 = [l2_dep 1 ntents;l2_ind (ntents+1) 2*ntents];
%     lamrange2 = [l2 1 flen;l2 (flen+1) 2*flen;l2 (2*flen+1) 3*flen];
%     lamrange2 = [l2 1 flen;l2 (flen+1) 2*flen;l2 (2*flen+1) 3*flen];
%     lamrange2 = [l2 1 flen;l2 (flen+1) 2*flen;l2 (2*flen+1) 3*flen;l2 (3*flen+1) 4*flen];
    llist = [];
    [fitp,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], 0);
    stim_dep_kern2(t,:) = fitp.k(1:ntents);
    stim_ind_kern2(t,:) = fitp.k((ntents+1):2*ntents);
    block_kern2(t,:) = fitp.k((2*ntents+1):end-1);
    ov_const2(t) = fitp.k(end);

    gfun = ov_Xmat*fitp.k(1:end-1) + fitp.k(end);
    too_large = find(gfun > 50);
    predrate = log(1+exp(gfun));
    predrate(too_large) = gfun(too_large);
    LL2(t) = -sum(Robs.*log(predrate)-predrate)/sum(Robs);
    
end

