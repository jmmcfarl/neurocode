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
used_expts = [1:9];
% used_expts = [10:15];
n_use_expts = length(used_expts);
for i = 1:nfold
    xv_inds{i} = [];
    for j = 1:n_use_expts
        cur_set = find(xvset_inds == xv_sets(j,i) & expt_inds == used_expts(j));
        xv_inds{i} = [xv_inds{i} cur_set'];
    end
    tr_inds{i} = setdiff(find(ismember(expt_inds,used_expts)),xv_inds{i});
end

%%
for xv = 1:nfold;
    xv
    compids   = 1:k_len;
    
    tsbs      = 1+floor(aselspks{tcell}/dt);
    tr_spbs = find(ismember(tr_inds{xv},tsbs));
    spikebins = tr_spbs;
%     xv_spbs = find(ismember(xv_inds{xv},tsbs));
    
    X_tr = X(tr_inds{xv},:);
%     X_xv = X(xv_inds{xv},:);
    
    stimlen = length(tr_inds{xv});
    stimlenxv = length(xv_inds{xv});
    Robs = zeros(1,stimlen);
    ftable = tabulate(spikebins);
    Robs(ftable(:,1)) = ftable(:,2);
%     Robsxv = zeros(1,stimlenxv);
%     ftable = tabulate(xv_spbs);
%     Robsxv(ftable(:,1)) = ftable(:,2);
    
    %FIT NULL MODEL
    avg_rate = length(spikebins)/stimlen;
    xvpred_rate = ones(1,stimlenxv)*avg_rate;
    trpred_rate = ones(1,stimlen)*avg_rate;
    null_LL(xv) = -sum(Robs.*log(trpred_rate) - trpred_rate)/sum(Robs)
%     null_xvLL(xv) = -sum(Robsxv.*log(xvpred_rate) - xvpred_rate)/sum(Robsxv)
    
%     %% TRY FITTING SEQUENCE OF STC MODELS WTIH INCREASING DIMENSIONALITY
%     %initialize model
%     defmod.lambda_dT = 0;
%     defmod.lambda_d2X = 5000;
%     defmod.lambda_L1x = 75;
%     
%     cur_ndims = 1;
%     init_signs = ones(cur_ndims,1);
%     init_kerns = 0.01*randn(k_len,cur_ndims);
%     kern_types{1} = 'lin';
%     glm(xv) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
%     glm(xv) = fitGNM_filters(glm(xv),X_tr,spikebins,'none',[],1e-4,1e-6);
%     glm_xvLL(xv) = getLL_GNM(glm(xv),X(xv_inds{xv},:),xv_spbs,'none');
%     [~,~,~,~,g] = getLL_GNM(glm(xv),X_tr,spikebins,'none');
%     glmr(xv) = fitGNM_spkNL(glm(xv),g,spikebins,0);
%     glmr_xvLL(xv) = getLL_GNM(glmr(xv),X(xv_inds{xv},:),xv_spbs,'none');
    
    %% TRY FITTING SEQUENCE OF STC MODELS WTIH INCREASING DIMENSIONALITY
    %initialize model
    defmod.lambda_dT = 0;
    defmod.lambda_d2X = 1500;
    defmod.lambda_L1x = 5;
    
    cur_ndims = 3;
    init_signs = ones(cur_ndims,1);
    init_kerns = 0.01*randn(k_len,cur_ndims);
    kern_types{1} = 'lin';
    for i = 2:cur_ndims
        kern_types{i} = 'quad';
    end
    quad_mod(xv) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
    quad_mod(xv).mods(1).lambda_L1x = 75;
    quad_mod(xv).mods(1).lambda_d2X = 5000;
    quad_mod(xv) = fitGNM_filters(quad_mod(xv),X_tr,spikebins,'none',[],1e-4,1e-6);
    quad_xvLL(xv) = getLL_GNM(quad_mod(xv),X(xv_inds{xv},:),xv_spbs,'none');
    [~,~,~,~,g,int_g] = getLL_GNM(quad_mod(xv),X_tr,spikebins,'none');
    quadr_mod(xv) = fitGNM_spkNL_reg(quad_mod(xv),X_tr,spikebins);
    quadr_mod(xv) = fitGNM_filters(quadr_mod(xv),X_tr,spikebins,'none',[],1e-4,1e-6);
    quadr_xvLL(xv) = getLL_GNM(quadr_mod(xv),X(xv_inds{xv},:),xv_spbs,'none');
   
        
    %%
    %initialize model
    defmod.lambda_dT = 0;
    defmod.lambda_d2X = 3000;
    defmod.lambda_L1x = 30;
    
    cur_ndims = 4;
    init_signs = ones(cur_ndims,1);
    init_kerns = 0.01*randn(k_len,cur_ndims);
    clear kenr_types
    for i = 1:cur_ndims
        kern_types{i} = 'threshlin';
    end
    gnm(xv) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
    gnm(xv) = fitGNM_filters(gnm(xv),X_tr,spikebins,'none',[],1e-4,1e-6);
%     gnm_LL(xv) = getLL_GNM(gnm(xv),X_tr,spikebins,'none');
%     gnm_xvLL(xv) = getLL_GNM(gnm(xv),X_xv,xv_spbs,'none');

    [~,~,~,~,g] = getLL_GNM(gnm(xv),X_tr,spikebins,'none');
    temp = fitGNM_spkNL(gnm(xv),g,spikebins,1);
%     [~,~,~,prate,g] = getLL_GNM(gnm(xv),X_xv,xv_spbs,'none');
%     [~,~,~,prate2,g] = getLL_GNM(temp,X_xv,xv_spbs,'none');
    
%     [~,~,~,~,g] = getLL_GNM(gnm(xv),X_xv,xv_spbs,'none');
%     temp = fitGNM_spkNL(gnm(xv),g,xv_spbs,1);

    
    gnmr(xv) = fitGNM_spkNL_reg(gnm(xv),X_tr,spikebins);
    gnmr_LL(xv) = getLL_GNM(gnmr(xv),X(tr_inds{xv},:),spikebins,'none');
%     gnmr_xvLL(xv) = getLL_GNM(gnmr(xv),X(xv_inds{xv},:),xv_spbs,'none')
    
    gnmr2(xv) = adjust_all_reg(gnmr(xv),'nltype','uncon');
    gnmr2(xv) = adjust_all_reg(gnmr2(xv),'nlmon',1);
    gnmr2(xv) = adjust_all_reg(gnmr2(xv),'lnl2',200);
    gnmr2(xv) = fitGNM_internal_NLs(gnmr2(xv),X_tr,spikebins,0,0);
%     gnmr2_xvLL(xv) = getLL_GNM(gnmr2(xv),X(xv_inds{xv},:),xv_spbs,'none')
    
%     gnmr3(xv) = fitGNM_filters(gnmr2(xv),X_tr,spikebins,'none',[],1e-4,1e-6);
%     gnmr3_xvLL(xv) = getLL_GNM(gnmr3(xv),X(xv_inds{xv},:),xv_spbs,'none')
    
end


%%
% sta_rel_xvLL = (null_xvLL - sta_mod_xvLL)/log(2);
% stc_rel_xvLL = (null_xvLL - stc_mod_xvLL)/log(2);
glm_rel_xvLL = (null_xvLL - glm_xvLL)/log(2);
quad_rel_xvLL = (null_xvLL - quad_xvLL)/log(2);
gnmr_rel_xvLL = (null_xvLL - gnmr2_xvLL)/log(2);

% Y = [sta_rel_xvLL(:); stc_rel_xvLL(:); quad_rel_xvLL(:); gnmr_rel_xvLL(:)];
Y = [glm_rel_xvLL(:); quad_rel_xvLL(:); gnmr_rel_xvLL(:)];
G = [ones(nfold,1); 2*ones(nfold,1); 3*ones(nfold,1)];
boxplot(Y,G)

hold on
Y = [glm_rel_xvLL(:) quad_rel_xvLL(:) gnmr_rel_xvLL(:)];
for i = 1:nfold
    plot(1:3,Y(i,:),'ko-')
end
set(gca,'fontsize',16,'fontname','arial')


%%
cd ~/James_scripts/GLM/t1/
save gnm_fig_data_v5 gnm* glm* quad* null*

%%
% xv = 4;
used_mod = gnmr2(xv);
% used_mod = quad_mod(xv);
% 
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
avg_best_theta = circ_mean(best_thetas);
avg_best_x = mean(x_cents);
avg_best_y = mean(y_cents);

f1 = figure;
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
    axis tight; xlim([-3.1 3.1])
    title('Module non-linearity')
    
end
colormap(gray)

% %%
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
