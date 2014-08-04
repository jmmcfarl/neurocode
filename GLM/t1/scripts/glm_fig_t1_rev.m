%% Load Data
clear all;
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/Timm/tim1/functions/')
addpath(genpath('~/James_Scripts'));

cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

cd ~/Data/blanche/matlabdata/
% load ./spks7576-all.mat;

load ./sub_Xmat2
stype='all';

tcell = 13;
tsbs      = 1+floor(aselspks{tcell}/dt);

sdim = 20;
fsdim = sdim^2;
flen = 8;
stim_params.spatial_dims = 2;
stim_params.sdim = sdim;
stim_params.flen = flen;
dt = 19.9920031987*1e-3; %true sweep time.

X = X/std(X(:));

%% create XV data
nfold = 5;
[stimlen,k_len] = size(X);
expt_len = 6000;
num_expts = stimlen/expt_len;
expt_inds = nan(stimlen,1);
xvset_inds = nan(stimlen,1);
for i = 1:num_expts
    cur_inds = (i-1)*expt_len + (1:expt_len);
    expt_inds(cur_inds) = i;
    xvset_inds(cur_inds) = floor(1:nfold/expt_len:(1+nfold-nfold/expt_len));
end
for i = 1:num_expts
    xv_sets(i,:) = randperm(nfold);
end
for i = 1:nfold
    xv_inds{i} = [];
    for j = 1:num_expts
        cur_set = find(xvset_inds == xv_sets(j,i) & expt_inds == j);
        xv_inds{i} = [xv_inds{i} cur_set'];
    end
    tr_inds{i} = setdiff(1:stimlen,xv_inds{i});
end


%%
for xv = 1:nfold;
    xv
    compids   = 1:k_len;
    
    tsbs      = 1+floor(aselspks{tcell}/dt);
    tr_spbs = find(ismember(tr_inds{xv},tsbs));
    spikebins = tr_spbs;
    xv_spbs = find(ismember(xv_inds{xv},tsbs));
    
    X_tr = X(tr_inds{xv},:);
    X_xv = X(xv_inds{xv},:);
    
    stimlen = length(tr_inds{xv});
    stimlenxv = length(xv_inds{xv});
    Robs = zeros(1,stimlen);
    ftable = tabulate(spikebins);
    Robs(ftable(:,1)) = ftable(:,2);
    Robsxv = zeros(1,stimlenxv);
    ftable = tabulate(xv_spbs);
    Robsxv(ftable(:,1)) = ftable(:,2);
    
    %FIT NULL MODEL
    avg_rate = length(spikebins)/stimlen;
    xvpred_rate = ones(1,stimlenxv)*avg_rate;
    trpred_rate = ones(1,stimlen)*avg_rate;
    null_LL(xv) = -sum(Robs.*log(trpred_rate) - trpred_rate)/sum(Robs)
    null_xvLL(xv) = -sum(Robsxv.*log(xvpred_rate) - xvpred_rate)/sum(Robsxv)
    
    %%
    compids = 500;
    [COEFF, SCORE, LATENT] = princomp(X_tr);
    whitemat = diag(1./sqrt(LATENT)); %whitening (rescaling) transformation from training set
    wstim = SCORE(:,1:compids)*whitemat(1:compids,1:compids); %whitened training data
    % STAC WITH WHITENING
    nneg = 0;
    npos = 4;
    spike_cond_stim = wstim(spikebins,:);
    sta = mean(spike_cond_stim) - mean(wstim);
    sta = sta/norm(sta);
    proj_mat = sta'/(sta*sta')*sta;
    %     stim_proj = wstim - wstim*proj_mat; %project out STA
    stim_proj = wstim; %don't project out STA
    stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
    [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
    stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
    %
    stcs = [sta' stcs];
    pix_conv_mat = COEFF(:,1:compids)';
    stcs = (stcs'*pix_conv_mat)';
    
    %%
%     kern_types{1} = 'lin';
%     sta_mod = createGNM(stcs(:,1),1,kern_types,[],stim_params);
%     g_mat = X_tr*stcs(:,1);
%     sta_mod = fitGNM_weights(sta_mod,g_mat,spikebins,1);
%     [~, ~, ~, ~, g] = getLL_GNM(sta_mod,X_tr,spikebins,'none');
%     sta_mod = fitGNM_spkNL(sta_mod,g,spikebins,0);
%     sta_mod_xvLL(xv) = getLL_GNM(sta_mod,X(xv_inds{xv},:),xv_spbs,'none');
%     
%     %      stcs = stcs(:,[1 2 4]);
%     g_mat = X_tr*stcs;
%     for i = 2:size(stcs,2)
%         kern_types{i} = 'quad';
%     end
%     stc_mod = createGNM(stcs,ones(size(stcs,2)),kern_types,[],stim_params);
%     stc_mod = fitGNM_weights(stc_mod,g_mat,spikebins,1);
%     [~, ~, ~, ~, g] = getLL_GNM(stc_mod,X_tr,spikebins,'none');
%     stc_mod = fitGNM_spkNL(stc_mod,g,spikebins,0);
%     stc_mod_xvLL(xv) = getLL_GNM(stc_mod,X(xv_inds{xv},:),xv_spbs,'none');
    
    
    %%
%     %initialize model
%     defmod.lambda_dT = 0;
%     defmod.lambda_d2X = 6000;
%     defmod.lambda_L1x = 12;
%     
%     cur_ndims = 3;
%     init_signs = ones(cur_ndims,1);
%     init_kerns = 0.01*randn(k_len,cur_ndims);
%     kern_types{1} = 'lin';
%     for i = 2:cur_ndims
%         kern_types{i} = 'quad';
%     end
%     quad_mod(xv) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
%     %     quad_mod(xv).mods(1).lambda_L1x = 75;
%     %     quad_mod(xv).mods(1).lambda_d2X = 5000;
%     quad_mod(xv) = fitGNM_filters(quad_mod(xv),X_tr,spikebins,'none',[],1e-4,1e-6,0);
%     quad_xvLL(xv) = getLL_GNM(quad_mod(xv),X(xv_inds{xv},:),xv_spbs,'none');
%     
%     quadr_mod(xv,1) = quad_mod(xv);
%     for k = 2
%         [~, ~, ~, ~, g] = getLL_GNM(quadr_mod(xv,k-1),X_tr,spikebins,'none');
%         quadr_mod(xv,k) = fitGNM_spkNL(quadr_mod(xv,k-1),g,spikebins,0);
%         quadr_mod(xv,k) = fitGNM_filters(quadr_mod(xv,k),X_tr,spikebins,'none',[],1e-4,1e-6,0);
%         quadr_xvLL(xv,k) = getLL_GNM(quadr_mod(xv,k),X(xv_inds{xv},:),xv_spbs,'none');
%     end
%     
%     %initialize model
%     defmod.lambda_dT = 0;
%     defmod.lambda_d2X = 6000;
%     defmod.lambda_L1x = 0;
%     cur_ndims = 3;
%     init_signs = ones(cur_ndims,1);
%     init_kerns = 0.01*randn(k_len,cur_ndims);
%     kern_types{1} = 'lin';
%     for i = 2:cur_ndims
%         kern_types{i} = 'quad';
%     end
%     quad_mod_ls(xv) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
%     quad_mod_ls(xv) = fitGNM_filters(quad_mod_ls(xv),X_tr,spikebins,'none',[],1e-4,1e-6,0);
%     quad_ls_xvLL(xv) = getLL_GNM(quad_mod_ls(xv),X(xv_inds{xv},:),xv_spbs,'none');
    
    %%
    quad_basis = get_k_mat(quad_mod_ls(xv));
    quad_out = X_tr*quad_basis;
    quad_xv_out = X(xv_inds{xv},:)*quad_basis;
    quad_stim_params.spatial_dims = 1;
    quad_stim_params.sdim = size(quad_basis,2);
    quad_stim_params.flen = 1;
    quad_klen = size(quad_basis,2);
    cur_ndims = 3;
    init_signs = ones(cur_ndims,1);
    clear kenr_types
    for i = 1:cur_ndims
        kern_types{i} = 'threshlin';
    end
    gnm_init_attempts = 10;
    test_xv = nan(gnm_init_attempts,1);
    clear test_gnm
    for i = 1:gnm_init_attempts
        init_kerns = 0.01*randn(quad_klen,cur_ndims);
        test_gnm(i) = createGNM(init_kerns,init_signs,kern_types,[],quad_stim_params);
        test_gnm(i) = fitGNM_filters(test_gnm(i),quad_out,spikebins,'none',[],1e-4,1e-6);
        test_LL(xv,i) = getLL_GNM(test_gnm(i),quad_out,spikebins,'none');
        test_xvLL(xv,i) = getLL_GNM(test_gnm(i),quad_xv_out,xv_spbs,'none');
    end
    [~,best_loc] = min(test_LL(xv,:));
    basis_mod = test_gnm(best_loc);
    [~, ~, ~, ~, g] = getLL_GNM(basis_mod,quad_out,spikebins,'none');
    basis_mod = fitGNM_spkNL(basis_mod,g,spikebins,0);
    
    %%
    %initialize model
    defmod.lambda_dT = 0;
    defmod.lambda_d2X = 6000; %3000
    defmod.lambda_L1x = 50; %40
    init_kerns = quad_basis*get_k_mat(basis_mod);
    gnm(xv) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
    gnm(xv).spk_alpha = basis_mod.spk_alpha;
    gnm(xv).spk_beta = basis_mod.spk_beta;
    gnm(xv).spk_theta = basis_mod.spk_theta;
    gnm(xv)  = fitGNM_filters(gnm(xv),X_tr,spikebins,'none',[],1e-4,1e-6,0);
    gnm_xvLL(xv) = getLL_GNM(gnm(xv),X(xv_inds{xv},:),xv_spbs,'none');
    
    gnmr(xv) = adjust_all_reg(gnm(xv),'nlx',linspace(-3,3,25));
    gnmr(xv) = setGNM_NLBFs(gnmr(xv),X_tr);
    gnmr(xv) = adjust_all_reg(gnmr(xv),'nltype','uncon');
    gnmr(xv) = adjust_all_reg(gnmr(xv),'nlmon',1);
    gnmr(xv) = adjust_all_reg(gnmr(xv),'lnl2',400);
    gnmr1(xv) = fitGNM_internal_NLs(gnmr(xv),X_tr,spikebins,0,0);
    gnmr1_xvLL(xv) = getLL_GNM(gnmr1(xv),X(xv_inds{xv},:),xv_spbs,'none');
    
    gnmr2(xv) = fitGNM_filters(gnmr1(xv),X_tr,spikebins,'none',[],1e-4,1e-6,0);
    gnmr2_xvLL(xv) = getLL_GNM(gnmr2(xv),X(xv_inds{xv},:),xv_spbs,'none');
    
%     [~, ~, ~, ~, g] = getLL_GNM(gnmr2(xv),X_tr,spikebins,'none');
%     gnmr3(xv) = fitGNM_spkNL(gnmr2(xv),g,spikebins,0);
%     gnmr3_xvLL(xv) = getLL_GNM(gnmr3(xv),X(xv_inds{xv},:),xv_spbs,'none');
%     
    
    
end
%%
% cd ~/James_scripts/GLM/t1/
% save gnm_fig_data_rev_v2 gnm* quad* null* test* stc* sta* xv_sets

%%
sta_rel_xvLL = (null_xvLL - sta_mod_xvLL)/log(2);
stc_rel_xvLL = (null_xvLL - stc_mod_xvLL)/log(2);
% glm_rel_xvLL = (null_xvLL - glm_xvLL)/log(2);
quad_rel_xvLL = (null_xvLL - quadr_xvLL(:,2)')/log(2);
gnmr_rel_xvLL = (null_xvLL - gnmr2_xvLL)/log(2);
% gnmr_rel_xvLL = (null_xvLL - gnm_xvLL)/log(2);

Y = [sta_rel_xvLL(:); stc_rel_xvLL(:); quad_rel_xvLL(:); gnmr_rel_xvLL(:)];
% Y = [sta_rel_xvLL(:); glm_rel_xvLL(:); stc_rel_xvLL(:); quad_rel_xvLL(:); gnmr_rel_xvLL(:)];
G = [ones(nfold,1); 2*ones(nfold,1); 3*ones(nfold,1); 4*ones(nfold,1);];
boxplot(Y,G)

hold on
% Y = [sta_rel_xvLL(:) glm_rel_xvLL(:) stc_rel_xvLL(:) quad_rel_xvLL(:) gnmr_rel_xvLL(:)];
Y = [sta_rel_xvLL(:) stc_rel_xvLL(:) quad_rel_xvLL(:) gnmr_rel_xvLL(:)];
for i = 1:nfold
    plot(1:4,Y(i,:),'ko-')
end
set(gca,'fontsize',16,'fontname','arial')

%% NOW FOR OVERALL MODEL FITS
[stimlen,k_len] = size(X);
ov_tr_inds = 1:stimlen;
tsbs      = 1+floor(aselspks{tcell}/dt);
tr_spbs = find(ismember(ov_tr_inds,tsbs));
spikebins = tr_spbs;
X_tr = X(ov_tr_inds,:);

Robs = zeros(1,stimlen);
ftable = tabulate(spikebins);
Robs(ftable(:,1)) = ftable(:,2);

%FIT NULL MODEL
avg_rate = length(spikebins)/stimlen;
trpred_rate = ones(1,stimlen)*avg_rate;
ov_null_LL = -sum(Robs.*log(trpred_rate) - trpred_rate)/sum(Robs)

%initialize model
clear defmod
defmod.lambda_dT = 0;
defmod.lambda_d2X = 6000;
defmod.lambda_L1x = 12;

cur_ndims = 3;
init_signs = ones(cur_ndims,1);
init_kerns = 0.01*randn(k_len,cur_ndims);
kern_types{1} = 'lin';
for i = 2:cur_ndims
    kern_types{i} = 'quad';
end
quad_mod_ov = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
% quad_mod_ov.mods(1).lambda_L1x = 75;
% quad_mod_ov.mods(1).lambda_d2X = 5000;
quad_mod_ov = fitGNM_filters(quad_mod_ov,X_tr,spikebins,'none',[],1e-4,1e-6);
[~,~,~,~,g] = getLL_GNM(quad_mod_ov,X_tr,spikebins,'none');
quadr_mod_ov = fitGNM_spkNL(quad_mod_ov,g,spikebins,0);
[~, ~, ~, ~, g] = getLL_GNM(quadr_mod_ov,X_tr,spikebins,'none');
quadr_mod_ov = fitGNM_spkNL(quadr_mod_ov,g,spikebins,0);
quadr_mod_ov = fitGNM_filters(quadr_mod_ov,X_tr,spikebins,'none',[],1e-4,1e-6,0);

% quad_k = get_k_mat(quadr_mod_ov);
% quad_k = bsxfun(@rdivide,quad_k,sqrt(sum(quad_k.^2)));
% quad_ov_norm = set_k_mat(quadr_mod_ov,quad_k);
% g_mat = X_tr*quad_k;
% quad_ov_norm = fitGNM_weights(quad_ov_norm,g_mat,spikebins,1);

%% FULL GNM FIT
quad_basis = get_k_mat(quadr_mod_ov);
quad_out = X_tr*quad_basis;
quad_stim_params.spatial_dims = 1;
quad_stim_params.sdim = size(quad_basis,2);
quad_stim_params.flen = 1;
quad_klen = size(quad_basis,2);
cur_ndims = 4;
init_signs = ones(cur_ndims,1);
clear kenr_types
for i = 1:cur_ndims
    kern_types{i} = 'threshlin';
end
gnm_init_attempts = 50;
clear test_LL
for i = 1:gnm_init_attempts
    init_kerns = 0.01*randn(quad_klen,cur_ndims);
    test_gnm(i) = createGNM(init_kerns,init_signs,kern_types,[],quad_stim_params);
    test_gnm(i) = fitGNM_filters(test_gnm(i),quad_out,spikebins,'none',[],1e-4,1e-6);
    test_LL(i) = getLL_GNM(test_gnm(i),quad_out,spikebins,'none');
end
[~,best_loc] = min(test_LL);
basis_mod = test_gnm(best_loc);
[~, ~, ~, ~, g] = getLL_GNM(basis_mod,quad_out,spikebins,'none');
basis_mod = fitGNM_spkNL(basis_mod,g,spikebins,0);

%initialize model
defmod.lambda_dT = 0;
defmod.lambda_d2X = 6000; %3000
defmod.lambda_L1x = 50; %40
init_kerns = quad_basis*get_k_mat(basis_mod);
gnm_ov = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
gnm_ov.spk_alpha = basis_mod.spk_alpha;
gnm_ov.spk_beta = basis_mod.spk_beta;
gnm_ov.spk_theta = basis_mod.spk_theta;
gnm_ov  = fitGNM_filters(gnm_ov,X_tr,spikebins,'none',[],1e-4,1e-6,0);

gnm_ov = adjust_all_reg(gnm_ov,'nlx',linspace(-3,3,25));
gnm_ov = setGNM_NLBFs(gnm_ov,X_tr);
gnm_ov = adjust_all_reg(gnm_ov,'nltype','uncon');
gnm_ov = adjust_all_reg(gnm_ov,'nlmon',1);
gnm_ov = adjust_all_reg(gnm_ov,'lnl2',400);
gnm_ov = fitGNM_internal_NLs(gnm_ov,X_tr,spikebins,0,0);

gnm_ov2 = fitGNM_filters(gnm_ov,X_tr,spikebins,'none',[],1e-4,1e-6,0);
[~, ~, ~, ~, g] = getLL_GNM(gnm_ov2,X_tr,spikebins,'none');
gnm_ov2 = fitGNM_spkNL(gnm_ov2,g,spikebins,0);

% gnm_k = get_k_mat(gnm_ov2);
% gnm_k = bsxfun(@rdivide,gnm_k,sqrt(sum(gnm_k.^2)));
% gnm_ov_norm = set_k_mat(gnm_ov2,gnm_k);
% g_mat = X_tr*gnm_k;
% gnm_ov_norm = fitGNM_weights(gnm_ov_norm,g_mat,spikebins,1);

%%
cd ~/James_scripts/GLM/t1/
save gnm_fig_data_ov_mods *_ov*


%%
% used_mod = quad_mod_ov
% used_mod = gnm_ov2
% used_mod = gnmr2(1);
used_mod = quadr_mod(1,2);

pix_mat = get_k_mat(used_mod);
n_filts = length(used_mod.mods);

smooth_mat = epanechikov(3)'*epanechikov(3);
for i = 1:n_filts
    cur_filt_kmat = reshape(pix_mat(:,i),flen,sdim^2);
    cur_filt_smkmat = zeros(size(cur_filt_kmat));
    for j = 1:flen
        temp_k = reshape(cur_filt_kmat(j,:),sdim,sdim);
        cur_sm_kmat = conv2(temp_k,smooth_mat,'same');
        cur_filt_smkmat(j,:) = cur_sm_kmat(:);
    end
    pix_mat(:,i) = cur_filt_smkmat(:);
end

clear gabor_fits gabor_fitvals nly_mat
for i = 1:n_filts
    nly_mat(i,:) = used_mod.mods(i).nly;
    [gabor_fits(i,:),gabor_fitvals(i,:,:)] = get_gaborfits_allslices(pix_mat(:,i),flen,sdim);
end

pix_mat = get_k_mat(used_mod);

nlx = used_mod.mods(1).nlx;
[~,best_slice_ids] = get2dMaxSlices(pix_mat,flen,sdim,1:(sdim^2));

best_thetas = nan(n_filts,1);
x_cents = nan(n_filts,1);
y_cents = nan(n_filts,1);
for i = 1:n_filts
    best_thetas(i) = gabor_fits(i,best_slice_ids(i)).theta;
    x_cents(i) = gabor_fits(i,best_slice_ids(i)).x0;
    y_cents(i) = gabor_fits(i,best_slice_ids(i)).y0;
end
best_thetas = mod(best_thetas,pi);
avg_best_theta = circ_mean(best_thetas);
avg_best_x = mean(x_cents);
avg_best_y = mean(y_cents);

% f1 = figure;
figure
for ii = 1:n_filts
    
    %compute space-time projections
    gabor_filtmat = squeeze(gabor_fitvals(ii,:,:));
    cur_filtmat = reshape(pix_mat(:,ii),flen,fsdim);
    
    %     filt_projprof = project_2drf(cur_filtmat(:),gabor_fits(ii,best_slice_ids(ii)).theta,sdim);
    filt_projprof = project_2drf(cur_filtmat(:),avg_best_theta,sdim);
    %     filt_projprof = project_2drf(cur_filtmat(:),best_thetas(ii),sdim);
    
    subplot(n_filts,4,(ii-1)*4+1)
    imagesc(reshape(cur_filtmat(best_slice_ids(ii),:),sdim,sdim));colormap(jet);
    zmax = max(abs(cur_filtmat(best_slice_ids(ii),:))); caxis([-zmax zmax]);
    title('Best time-slice')
    set(gca,'ydir','normal')
    hold on
    plot(avg_best_x,avg_best_y,'ro')
    xax = xlim();
    yax = ylim();
    cur_linex = linspace(xax(1),xax(end),50);
    cur_liney = tan(avg_best_theta)*(cur_linex - avg_best_x)+avg_best_y;
    plot(cur_linex,cur_liney,'color','w')
    xlim(xax);ylim(yax);
    
    subplot(n_filts,4,(ii-1)*4+2)
    imagesc(reshape(gabor_filtmat(best_slice_ids(ii),:),sdim,sdim));colormap(jet);
    hold on
    set(gca,'ydir','normal')
    xax = xlim();
    yax = ylim();
    cur_linex = linspace(xax(1),xax(end),50);
    cur_liney = tan(avg_best_theta)*(cur_linex - avg_best_x)+avg_best_y;
    plot(cur_linex,cur_liney,'color','w')
    xm = max(abs(gabor_filtmat(:)));
    xlim(xax);ylim(yax);
    caxis([-xm xm]);
    title('Space-time projection')
    
    subplot(n_filts,4,(ii-1)*4+3)
    imagesc(1:sdim,-(0:5)*0.02,filt_projprof);
    xm = max(abs(filt_projprof(:)));
    caxis([-xm xm]/1.5);
    title('Space-time projection')
    
    subplot(n_filts,4,(ii-1)*4+4)
    plot(used_mod.mods(ii).fox,used_mod.mods(ii).foy,'k')
    hold on
    plot(nlx,nly_mat(ii,:),'b')
    axis tight; %xlim([-3.1 3.1])
%     xlim(nlx([1 end]))
xlim([-1 1])
    title('Module non-linearity')
    
end
colormap(gray)

%%
% used_mod = quad_mod(xv);
% % used_mod = stc_mod;
% pix_mat = get_k_mat(used_mod);
% n_filts = length(used_mod.mods);
%
% smooth_mat = epanechikov(3)'*epanechikov(3);
% for i = 1:n_filts
%     cur_filt_kmat = reshape(pix_mat(:,i),flen,sdim^2);
%     cur_filt_smkmat = zeros(size(cur_filt_kmat));
%     for j = 1:flen
%         temp_k = reshape(cur_filt_kmat(j,:),sdim,sdim);
%         cur_sm_kmat = conv2(temp_k,smooth_mat,'same');
%         cur_filt_smkmat(j,:) = cur_sm_kmat(:);
%     end
%     pix_mat(:,i) = cur_filt_smkmat(:);
% end
%
% clear gabor_fits gabor_fitvals nly_mat
% for i = 1:n_filts
%     nly_mat(i,:) = used_mod.mods(i).nly;
%     [gabor_fits(i,:),gabor_fitvals(i,:,:)] = get_gaborfits_allslices(pix_mat(:,i),flen,sdim);
% end
%
% nlx = used_mod.mods(1).nlx;
% [~,best_slice_ids] = get2dMaxSlices(pix_mat,flen,sdim,1:(sdim^2));
%
% best_thetas = nan(n_filts,1);
% for i = 1:n_filts
%     best_thetas(i) = gabor_fits(i,best_slice_ids(i)).theta;
% end
% avg_best_theta = circ_mean(best_thetas);
%
% f1 = figure;
% for ii = 1:n_filts
%
%     %compute space-time projections
%     gabor_filtmat = squeeze(gabor_fitvals(ii,:,:));
%     cur_filtmat = reshape(pix_mat(:,ii),flen,fsdim);
%
%     %     filt_projprof = project_2drf(cur_filtmat(:),gabor_fits(ii,best_slice_ids(ii)).theta,sdim);
%     filt_projprof = project_2drf(cur_filtmat(:),avg_best_theta,sdim);
%     subplot(n_filts,3,(ii-1)*3+1)
%     imagesc(reshape(cur_filtmat(best_slice_ids(ii),:),sdim,sdim));colormap(jet);
%     zmax = max(abs(cur_filtmat(best_slice_ids(ii),:))); caxis([-zmax zmax]);
%     title('Best time-slice')
%
%     subplot(n_filts,3,(ii-1)*3+2)
%     imagesc(1:sdim,-(0:5)*0.02,filt_projprof);
%     xm = max(abs(filt_projprof(:)));
%     %     caxis([-xm xm]);
%     title('Space-time projection')
%     subplot(n_filts,3,(ii-1)*3+3)
%     plot(used_mod.mods(ii).fox,used_mod.mods(ii).foy,'k')
%     hold on
%     axis tight; xlim([-2 2])
%     plot(nlx,nly_mat(ii,:),'b')
%     hold on
%     title('Module non-linearity')
%
% end
% colormap(gray)
