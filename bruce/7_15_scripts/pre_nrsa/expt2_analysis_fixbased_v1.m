%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt2_compiled_windata_d1p5.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

fix_win_dur = 0.15;

Pix2Deg = 0.018837;
NT = length(full_fix_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

%%
n_fixs = length(full_fix_wends);
fix_binned_spks = nan(n_fixs,96);
fix_expt_num = nan(n_fixs,1);
for i = 1:n_fixs
    cur_inds = full_fix_starts(i):full_fix_wends(i);
    fix_binned_spks(i,:) = sum(full_binned_spks(cur_inds,:));
    
    fix_expt_num(i) = unique(full_expt_vec(cur_inds));
end
[un_expts,~,full_expt_inds] = unique(fix_expt_num);
n_un_expts = length(un_expts);
linX = zeros(n_fixs,n_un_expts-1);
for i = 1:n_un_expts-1
    linX(fix_expt_num==i,i) = 1;
end

%%
% load ./expt1_eyecor_d1p25_nosac_v2.mat gabor*
load ./gabor_tracking_varmeans.mat gabor*
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

% %% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
% all_gabor_out1 = resh_all_stims*gabor_emp1_filt';
% all_gabor_out2 = resh_all_stims*gabor_emp2_filt';
% spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);
% 
% silent = 1;
% 
% for t = 1:96
%     fprintf('Fitting Cell %d of 96\n',t);
%     cur_spatial_mod_out = spatial_mod_out(:,t);
%     ov_Xmat = [cur_spatial_mod_out linX];
%     
%     Robs = fix_binned_spks(:,t);
%     
%     poss_spkcnts = unique(Robs(Robs > 0));
%     spksN = [];
%     for i = 1:length(poss_spkcnts)
%         cur_set = find(Robs == poss_spkcnts(i));
%         spksN = [spksN; repmat(cur_set,poss_spkcnts(i),1)];
%     end
%     spksN = sort(spksN);
%     
%     klen = size(ov_Xmat,2);
%     K0 = zeros(klen,1);
%     lamrange2 = [];
%     llist = [];
%     [fitp,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], 0);
%     stim_dep_kern(t) = fitp.k(1);
%     block_kern(t,:) = fitp.k(2:end-1);
%     ov_const(t) = fitp.k(end);
%     
%     gfun = ov_Xmat*fitp.k(1:end-1) + fitp.k(end);
%     too_large = find(gfun > 50);
%     predrate = log(1+exp(gfun));
%     predrate(too_large) = gfun(too_large);
%     LL(t) = -sum(Robs.*log(predrate)-predrate)/sum(Robs);
%     
%     if xv_frac > 0
%         cur_spatial_mod_out = spatial_mod_out(full_fix_ids(xv_inds),t);
%         stim_dep_X = bsxfun(@times,sac_Tmat(xv_inds,:),cur_spatial_mod_out);
%         ov_Xmat = [stim_dep_X sac_Tmat(xv_inds,:) linX(xv_inds,:)];
%         
%         gfun = ov_Xmat*fitp.k(1:end-1) + fitp.k(end);
%         too_large = find(gfun > 50);
%         xv_predrate = log(1+exp(gfun));
%         xv_predrate(too_large) = gfun(too_large);
%         xv_LL(t) = -sum(Robs_xv.*log(xv_predrate)-xv_predrate)/sum(Robs_xv);
%     end
% end


% %% SET UP XV CELL SET
% NSIG = 96;
% xv_frac = 0.2;
% tr_set = randperm(NSIG);
% tr_set = sort(tr_set(1:round(length(tr_set)*(1-xv_frac))));
% xv_set = setdiff(1:NSIG,tr_set);
% n_tr_cells = length(tr_set);
% 
% %%
% all_stim_dep_kern{1} = stim_dep_kern;
% all_block_kern{1} = block_kern;
% all_ov_const{1} = ov_const;
% all_LL{1} = LL;
% 
% n_iter = 5;
% for it = 1:n_iter
%     NT = length(tr_inds);
%     max_shift = 35;
%     dshift = 2;
%     x_shifts = -max_shift:dshift:max_shift;
%     y_shifts = -max_shift:dshift:max_shift;
%     [Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
%     SH = [Xsh(:) Ysh(:)];
%     n_shifts = size(SH,1);
%     
%     fix_LLs = zeros(n_fixs,n_shifts);
%     Robs = fix_binned_spks(:,tr_set);
%     
%     gabor_filt_bank1 = reshape(gabor_emp1_filt(tr_set,:)',[sdim sdim n_tr_cells]);
%     gabor_filt_bank2 = reshape(gabor_emp2_filt(tr_set,:)',[sdim sdim n_tr_cells]);
%     shifted_gabor_bank1 = nan(sdim^2,n_tr_cells);
%     shifted_gabor_bank2 = nan(sdim^2,n_tr_cells);
%     
%     stim_dep_X = zeros(n_fixs,length(tr_set));
%     
%     shift_cnt = 1;
%     for xx = 1:length(x_shifts)
%         for yy = 1:length(y_shifts)
%             
%             fprintf('Shift %d of %d\n',shift_cnt,n_shifts);
%             d2 = dist_shift3d(gabor_filt_bank1,x_shifts(xx),2);
%             d2 = dist_shift3d(d2,y_shifts(yy),2);
%             shifted_gabor_bank1 = reshape(d2,sdim^2,n_tr_cells);
%             d2 = dist_shift3d(gabor_filt_bank2,x_shifts(xx),1);
%             d2 = dist_shift3d(d2,y_shifts(yy),2);
%             shifted_gabor_bank2 = reshape(d2,sdim^2,n_tr_cells);
%             
%             all_gabor_out1 = resh_all_stims*shifted_gabor_bank1;
%             all_gabor_out2 = resh_all_stims*shifted_gabor_bank2;
%             spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);
%             
%             stim_dep_X = bsxfun(@times,spatial_mod_out,all_stim_dep_kern{it}(tr_set));
%             
%             gfun = stim_dep_X + linX*all_block_kern{it}(tr_set,:)';
%             gfun = bsxfun(@plus,gfun,all_ov_const{it}(tr_set));
%             
%             too_large = gfun > 50;
%             pred_rate = log(1+exp(gfun));
%             pred_rate(too_large) = gfun(too_large);
%             pred_rate(pred_rate < 1e-20) = 1e-20;
%             
%             LLs = Robs.*log(pred_rate) - pred_rate;
%             fix_LLs(:,shift_cnt) = sum(LLs,2);
%             shift_cnt = shift_cnt + 1;
%         end
%     end
%     
%     %% HMM for inferring sequence of stimulus translations
%     
%     %overall prior on shifts
%     eps_prior_sigma = 0.5; %0.2
%     leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
%     leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
%     
%     % leps_prior = ones(n_shifts,1)/n_shifts;
%     
%     lPost = bsxfun(@plus,fix_LLs,leps_prior');
%     lPost = bsxfun(@minus,lPost,logsumexp(lPost,2));
%     
%     % for i = 500:n_fixs
%     %    imagesc(reshape(lPost(i,:),length(x_shifts),length(y_shifts)))
%     %    pause
%     %    clf
%     % end
%     %% RECONSTRUCT MAP STIMULUS
%     [max_post,max_loc] = max(lPost,[],2);
%     x_cor{it} = SH(max_loc,1);
%     y_cor{it} = SH(max_loc,2);
%     
%     resh_X = reshape(resh_all_stims',[sdim sdim n_fixs]);
%     resh_X_sh = zeros(size(resh_X));
%     for ii = 1:n_fixs
%         d2 = dist_shift2d(resh_X(:,:,ii), -x_cor{it}(ii), 2,0);
%         d2 = dist_shift2d(d2,-y_cor{it}(ii),1,0);
%         resh_X_sh(:,:,ii) = d2;
%     end
%     fullX_sh = reshape(resh_X_sh,sdim^2,n_fixs)';
%     
%     %% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
%     
%     all_gabor_out1 = fullX_sh*gabor_emp1_filt';
%     all_gabor_out2 = fullX_sh*gabor_emp2_filt';
%     spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);
%     
%     for t = 1:96
%         fprintf('Fitting Cell %d of 96\n',t);
%         cur_spatial_mod_out = spatial_mod_out(:,t);
%         ov_Xmat = [cur_spatial_mod_out linX];
%         
%         cur_Robs = fix_binned_spks(:,t);
%         
%         poss_spkcnts = unique(cur_Robs(cur_Robs > 0));
%         spksN = [];
%         for i = 1:length(poss_spkcnts)
%             cur_set = find(cur_Robs == poss_spkcnts(i));
%             spksN = [spksN; repmat(cur_set,poss_spkcnts(i),1)];
%         end
%         spksN = sort(spksN);
%         
%         klen = size(ov_Xmat,2);
%         K0 = zeros(klen,1);
%         lamrange2 = [];
%         llist = [];
%         [fitp,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], 0);
%         all_stim_dep_kern{it+1}(t) = fitp.k(1);
%         all_block_kern{it+1}(t,:) = fitp.k(2:end-1);
%         all_ov_const{it+1}(t) = fitp.k(end);
%         
%         gfun = ov_Xmat*fitp.k(1:end-1) + fitp.k(end);
%         too_large = find(gfun > 50);
%         predrate = log(1+exp(gfun));
%         predrate(too_large) = gfun(too_large);
%         all_LL{it+1}(t) = -sum(cur_Robs.*log(predrate)-predrate)/sum(cur_Robs);
%         
%     end
%     
% end
%%

%%
% save fixation_based_corrections_34 LL* x_cor y_cor all_* 
% % save fixation_based_corrections LL* x_cor y_cor stim_dep* block_* ov_*

%%
load ./fixation_based_corrections_34
it = 5;
resh_X = reshape(resh_all_stims',[sdim sdim n_fixs]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:n_fixs
    d2 = dist_shift2d(resh_X(:,:,ii), -x_cor{it}(ii), 2,0);
    d2 = dist_shift2d(d2,-y_cor{it}(ii),1,0);
    resh_X_sh(:,:,ii) = d2;
end
fullX_sh = reshape(resh_X_sh,sdim^2,n_fixs)';

%%
infix_ids = find(~isnan(full_fix_ids));
n_fixs = length(unique(full_fix_ids(infix_ids)));
xv_frac =0.;
xv_num = round(xv_frac*n_fixs);
xv_f = randperm(n_fixs);
xv_f(xv_num+1:end) = [];
tr_f = setdiff(1:n_fixs,xv_f);
xv_inds = find(ismember(full_fix_ids,xv_f));
tr_inds = setdiff(find(~isnan(full_fix_ids)),xv_inds);

flen = 80;
flen_t = flen*dtd;
NT = length(full_expt_vec);

tent_centers = [0:dtd:0.15];
cur_sp = dtd;
while max(tent_centers) < flen_t
    tent_centers = [tent_centers (tent_centers(end)+cur_sp)];
    cur_sp = cur_sp + dtd/6;
end

% tent_centers = [0:dtd:0.1];
% cur_sp = dtd;
% while max(tent_centers) < flen_t
%     tent_centers = [tent_centers (tent_centers(end)+cur_sp)];
%     cur_sp = cur_sp + dtd/4;
% end
tent_centers = round(tent_centers/dtd)*dtd;
if tent_centers(end) >= flen_t
    tent_centers(end) = [];
end
tbmat = construct_tent_bases(tent_centers,dtd);
[ntents,tblen] = size(tbmat);

tbmat = [zeros(ntents,tblen-1) tbmat];

sac_inds = zeros(size(full_t));
sac_inds(full_sac_inds(:,2)) = 1;
sac_sinds = zeros(size(full_t));
sac_sinds(full_sac_inds(:,1)) = 1;

sac_Tmat = zeros(length(full_t),ntents);
sac_sTmat = zeros(length(full_t),ntents);
for i = 1:ntents
    sac_Tmat(:,i) = conv(sac_inds,tbmat(i,:),'same');
    sac_sTmat(:,i) = conv(sac_sinds,tbmat(i,:),'same');
end

[un_expts,~,full_expt_inds] = unique(full_expt_vec);
n_un_expts = length(un_expts);
linX = zeros(NT,n_un_expts-1);
for i = 1:n_un_expts-1
    linX(full_expt_inds==i,i) = 1;
end

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
all_gabor_out1 = fullX_sh*gabor_emp1_filt';
all_gabor_out2 = fullX_sh*gabor_emp2_filt';
spatial_mod_out = zscore(sqrt(all_gabor_out1.^2 + all_gabor_out2.^2));

% all_gabor_out1 = resh_all_stims*gabor_emp1_filt';
% all_gabor_out2 = resh_all_stims*gabor_emp2_filt';
% ospatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);

l2_ind = 100;
l2_dep = 2000;
silent = 1;
clear fstim_* ostim_*
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
    llist = [];
    [fitp,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], 0);
    fstim_dep_kern(t,:) = fitp.k(1:ntents);
    fstim_ind_kern(t,:) = fitp.k((ntents+1):2*ntents);
    fblock_kern(t,:) = fitp.k((2*ntents+1):end-1);
    fov_const(t) = fitp.k(end);
    
    gfun = ov_Xmat*fitp.k(1:end-1) + fitp.k(end);
    too_large = find(gfun > 50);
    predrate = log(1+exp(gfun));
    predrate(too_large) = gfun(too_large);
    LLf(t) = -sum(Robs.*log(predrate)-predrate)/sum(Robs);
    
    cur_lambda = 50;
    ov_Xmat = [sac_Tmat(tr_inds,:) linX(tr_inds,:)];
    lamrange2 = [cur_lambda 1 ntents];
    llist = [];
    [fitp,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], 0);
    trig_avg_kern = fitp.k(1:ntents);
    offset = fitp.k(end) + mean(fitp.k(ntents+1:end-1));
    trig_avg_rate(t,:) = log(1+exp(trig_avg_kern + offset));

    
%     cur_spatial_mod_out = ospatial_mod_out(full_fix_ids(tr_inds),t);
%     stim_dep_X = bsxfun(@times,sac_Tmat(tr_inds,:),cur_spatial_mod_out);
%     ov_Xmat = [stim_dep_X sac_Tmat(tr_inds,:) linX(tr_inds,:)];
%     
%     [fitp,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], 0);
%     ostim_dep_kern(t,:) = fitp.k(1:ntents);
%     ostim_ind_kern(t,:) = fitp.k((ntents+1):2*ntents);
%     oblock_kern(t,:) = fitp.k((2*ntents+1):end-1);
%     oov_const(t) = fitp.k(end);
%     
%     gfun = ov_Xmat*fitp.k(1:end-1) + fitp.k(end);
%     too_large = find(gfun > 50);
%     predrate = log(1+exp(gfun));
%     predrate(too_large) = gfun(too_large);
%     LLo(t) = -sum(Robs.*log(predrate)-predrate)/sum(Robs);
    
end

%%
% save expt2_unit_tempmods LLo LLf fstim* fblock* fov_* ostim* oblock* oov* tent_centers ntents dtd trig_avg_rate
save expt2_unit_tempmods_v2 LLf fstim* fblock* fov_* tent_centers ntents dtd trig_avg_rate
%%
for t = 1:96
    subplot(3,1,1)
%     plot(tent_centers,ostim_dep_kern(t,:));
    hold on
    plot(tent_centers,fstim_dep_kern(t,:),'r')
    xlim([0 0.5])
    subplot(3,1,2)
%     plot(tent_centers,ostim_ind_kern(t,:));
    hold on
    plot(tent_centers,fstim_ind_kern(t,:),'r')
    xlim([0 0.5])
    subplot(3,1,3)
    plot(tent_centers,trig_avg_rate(t,:),'k')
    xlim([0 0.5])
    t
    pause
    clf
end

%%
im_members = nan(n_fixs,1);
im_members(fix_imnum < 686) = 1;
im_members(fix_imnum >= 686 & fix_imnum <= 913) = 2;
im_members(fix_imnum > 913 & fix_imnum <= 1141) = 3;
im_members(fix_imnum > 1141) = 4;
full_im_members = nan(size(full_expt_vec));
full_im_members(tr_inds) = im_members(full_fix_ids(tr_inds));
full_im_members_xv = im_members(full_fix_ids(xv_inds));

back_win = round(0.05/dtd);
for i = 1:n_fixs
    cur_set = (full_fix_starts(i)-back_win):full_fix_starts(i);
    full_im_members(cur_set) = 0;
end
uset = find(~isnan(full_im_members));
full_fix_ids(full_im_members==0) = 1;

%%
% zspatial_mod_out = zscore(spatial_mod_out);
t = 4;

cur_spatial_mod_out = spatial_mod_out(full_fix_ids(uset),t);
cur_spatial_mod_out(full_im_members(uset) == 0) = 0;

% net_spatial_mod_out = mean(spatial_mod_out(full_fix_ids(tr_inds),:),2);
stim_dep_X = bsxfun(@times,sac_Tmat(uset,:),cur_spatial_mod_out);
% stim_dep_X2 = bsxfun(@times,sac_sTmat(uset,:),cur_spatial_mod_out);
% clear stim_dep_netX
% for i = 1:4
%     stim_dep_netX(i,:,:) = bsxfun(@times,sac_Tmat(tr_inds,:),cur_spatial_mod_out);
%     stim_dep_netX(i,full_im_members~=i,:) = 0;
% end
% stim_dep_netX = reshape(stim_dep_netX,length(tr_inds),ntents*4);

% ov_linX = [sac_sTmat(tr_inds,:) linX(tr_inds,:)];
ov_linX = [sac_Tmat(uset,:) sac_sTmat(uset,:) linX(uset,:)];
% ov_linX = [sac_Tmat(uset,:) linX(uset,:)];
% ov_linX = [sac_Tmat(tr_inds,:)];

stim_dep_X = [stim_dep_X];

stim_params.spatial_dims = 0;
stim_params.sdim = 1;
stim_params.flen = ntents;

Robs = full_binned_spks(uset,t);
Robs_xv = full_binned_spks(xv_inds,t);

poss_spkcnts = unique(Robs(Robs > 0));
spksN = [];
for i = 1:length(poss_spkcnts)
    cur_set = find(Robs == poss_spkcnts(i));
    spksN = [spksN; repmat(cur_set,poss_spkcnts(i),1)];
end
spksN = sort(spksN);

nmods = 1;
init_signs = [1];
kern_types{1} = 'lin';
init_kerns = zeros(ntents*stim_params.sdim,nmods);
defmod.lambda_d2T = 500;
defmod.lambda_L1x = 0;
nlk = 2;
% init_linK = zeros(ntents+n_un_expts-1,1);
init_linK = zeros(nlk*ntents+n_un_expts-1,1);
% init_linK = zeros(ntents,1);
gnm = createGNM_v2(init_kerns,init_signs,kern_types,init_linK,defmod,stim_params);
gnm.lambda_d2T_lin = 100;
gnm.lambda_L1_lin = 5;
gnm.nlk = nlk;
gnm = fitGNM_filters_v2(gnm,stim_dep_X,ov_linX,spksN,'none',[],[],[]);

% nmods = 2;
% init_signs = [1 -1];
% kern_types{1} = 'threshlin';
% kern_types{2} = 'threshlin';
% init_kerns = zeros(ntents*2,nmods);
% defmod.lambda_d2T = 500;
% init_linK = zeros(ntents+n_un_expts-1,1);
% % init_linK = zeros(ntents,1);
% gnm2 = createGNM_v2(init_kerns,init_signs,kern_types,init_linK,defmod,stim_params);
% gnm2.lambda_d2T_lin = 500;
% gnm2 = fitGNM_filters_v2(gnm2,stim_dep_X,ov_linX,spksN,'none',[],[],[]);

cur_spatial_mod_out = spatial_mod_out(full_fix_ids(xv_inds),t);
stim_dep_X = bsxfun(@times,sac_Tmat(xv_inds,:),cur_spatial_mod_out);
% stim_dep_X2 = bsxfun(@times,sac_sTmat(uset,:),cur_spatial_mod_out);
% clear stim_dep_netX
% for i = 1:4
%     stim_dep_netX(i,:,:) = bsxfun(@times,sac_Tmat(tr_inds,:),cur_spatial_mod_out);
%     stim_dep_netX(i,full_im_members~=i,:) = 0;
% end
% stim_dep_netX = reshape(stim_dep_netX,length(tr_inds),ntents*4);

% ov_linX = [sac_sTmat(tr_inds,:) linX(tr_inds,:)];
ov_linX = [sac_Tmat(xv_inds,:) sac_sTmat(xv_inds,:) linX(xv_inds,:)];
% ov_linX = [sac_Tmat(xv_inds,:) linX(xv_inds,:)];
% ov_linX = [sac_Tmat(tr_inds,:)];

stim_dep_X = [stim_dep_X];

poss_spkcnts = unique(Robs_xv(Robs_xv > 0));
spksN = [];
for i = 1:length(poss_spkcnts)
    cur_set = find(Robs_xv == poss_spkcnts(i));
    spksN = [spksN; repmat(cur_set,poss_spkcnts(i),1)];
end
spksN = sort(spksN);

xvLL(t) = getLL_GNM_v2(gnm,stim_dep_X,ov_linX,spksN,'none');
% xvLL2(t) = getLL_GNM_v2(gnm2,stim_dep_X,ov_linX,spksN,'none');

