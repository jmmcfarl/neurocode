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


% poss_reg_params = [0.1 0.5 2 5 10 25 50 100]/2;
poss_reg_params = [0 0.05 0.25 1 2 5 10 20];
cent_reg_val = 4;
defmod.lambda_dT = 0;
defmod.lambda_d2X = 6000;
defmod.lambda_L1x = 50;
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
    cur_ndims = 3;
    init_signs = ones(cur_ndims,1);
    init_kerns = 0.001*randn(k_len,cur_ndims);
    kern_types{1} = 'lin';
    for i = 2:cur_ndims
        kern_types{i} = 'quad';
    end
    for pp = 1:length(poss_reg_params)
        fprintf('Reg params %d of %d\n',pp,length(poss_reg_params));
        if pp == 1 %start from scratch if first reg param
            quad_mod(xv,pp) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
        else %otherwise initialize filters based on previous reg params
            quad_mod(xv,pp) = quad_mod(xv,pp-1);
        end
        quad_mod(xv,pp) = adjust_all_reg(quad_mod(xv,pp),'lambda_d2X',defmod.lambda_d2X*poss_reg_params(pp));
        quad_mod(xv,pp) = adjust_all_reg(quad_mod(xv,pp),'lambda_L1x',0);
        %         quad_mod(xv,pp) = adjust_all_reg(quad_mod(xv,pp),'lambda_L1x',defmod.lambda_L1x*poss_reg_params(pp));
        
        quad_mod(xv,pp) = fitGNM_filters(quad_mod(xv,pp),X_tr,spikebins,'none',[],1e-4,1e-6,1);
        quad_xvLL(xv,pp) = getLL_GNM(quad_mod(xv,pp),X(xv_inds{xv},:),xv_spbs,'none');
        
        [~, ~, ~, ~, g] = getLL_GNM(quad_mod(xv,pp),X_tr,spikebins,'none');
        quadr_mod(xv,pp) = fitGNM_spkNL(quad_mod(xv,pp),g,spikebins,0);
        quadr_mod(xv,pp) = fitGNM_filters(quadr_mod(xv,pp),X_tr,spikebins,'none',[],1e-4,1e-6);
        quadr_xvLL(xv,pp) = getLL_GNM(quadr_mod(xv,pp),X(xv_inds{xv},:),xv_spbs,'none');
    end
    
    %%
    quad_basis = get_k_mat(quadr_mod(xv,cent_reg_val));
    quad_out = X_tr*quad_basis;
    quad_xv_out = X(xv_inds{xv},:)*quad_basis;
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
    gnm_init_attempts = 10;
    test_xv = nan(gnm_init_attempts,1);
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
    init_kerns = quad_basis*get_k_mat(basis_mod);
    for pp = 1:length(poss_reg_params)
        fprintf('Reg params %d of %d\n',pp,length(poss_reg_params));
        
        gnm(xv,pp) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
        
        gnm(xv,pp) = adjust_all_reg(gnm(xv,pp),'lambda_d2X',defmod.lambda_d2X*poss_reg_params(pp));
        gnm(xv,pp) = adjust_all_reg(gnm(xv,pp),'lambda_L1x',0);
        %         gnm(xv,pp) = adjust_all_reg(gnm(xv,pp),'lambda_L1x',defmod.lambda_L1x*poss_reg_params(pp));
        
        gnm(xv,pp) = fitGNM_filters(gnm(xv,pp),X_tr,spikebins,'none',[],1e-4,1e-6);
        gnm(xv,pp).spk_alpha = basis_mod.spk_alpha;
        gnm(xv,pp).spk_beta = basis_mod.spk_beta;
        gnm(xv,pp).spk_theta = basis_mod.spk_theta;
        gnm(xv,pp)  = fitGNM_filters(gnm(xv,pp),X_tr,spikebins,'none',[],1e-4,1e-6,1);
        gnm_xvLL(xv,pp) = getLL_GNM(gnm(xv,pp),X(xv_inds{xv},:),xv_spbs,'none');
        
        gnmr(xv,pp) = adjust_all_reg(gnm(xv,pp),'nlx',linspace(-3,3,25));
        gnmr(xv,pp) = setGNM_NLBFs(gnmr(xv,pp),X_tr);
        gnmr(xv,pp) = adjust_all_reg(gnmr(xv,pp),'nltype','uncon');
        gnmr(xv,pp) = adjust_all_reg(gnmr(xv,pp),'nlmon',1);
        gnmr(xv,pp) = adjust_all_reg(gnmr(xv,pp),'lnl2',400);
        gnmr1(xv,pp) = fitGNM_internal_NLs(gnmr(xv,pp),X_tr,spikebins,0,0);
        gnmr1_xvLL(xv,pp) = getLL_GNM(gnmr1(xv,pp),X(xv_inds{xv},:),xv_spbs,'none');
        
        gnmr2(xv,pp) = fitGNM_filters(gnmr1(xv,pp),X_tr,spikebins,'none',[],1e-4,1e-6,1);
        gnmr2_xvLL(xv,pp) = getLL_GNM(gnmr2(xv,pp),X(xv_inds{xv},:),xv_spbs,'none');
        
    end
    
cd ~/James_scripts/GLM/t1/
save gnm_fig_data_rev_v1_regtest_smooth2 gnm* quad* null* 
end
%%
cd ~/James_scripts/GLM/t1/
save gnm_fig_data_rev_v1_regtest_smooth2 gnm* quad* null* 

%%
use_xv = 1:5;
gnm_xv_imp = -bsxfun(@minus,gnmr2_xvLL(use_xv,:),null_xvLL(use_xv)')/log(2);
quad_xv_imp = -bsxfun(@minus,quadr_xvLL(use_xv,:),null_xvLL(use_xv)')/log(2);

reg_eps = 0.01;

figure
% errorbar(poss_reg_params*50+1,mean(gnm_xv_imp),std(gnm_xv_imp)/sqrt(length(use_xv)));
% hold on
% errorbar(poss_reg_params*50+1,mean(quad_xv_imp),std(quad_xv_imp)/sqrt(length(use_xv)),'r');
errorbar(poss_reg_params + reg_eps,mean(gnm_xv_imp),std(gnm_xv_imp)/sqrt(length(use_xv)));
hold on
errorbar(poss_reg_params + reg_eps,mean(quad_xv_imp),std(quad_xv_imp)/sqrt(length(use_xv)),'r');
set(gca,'xscale','log')
yl = ylim();
ylim([0 yl(2)])
% xlim([0.5 2000])
xlim([0.005 40])
set(gca,'fontname','arial','fontsize',16)
xlabel('Smoothness penalty','fontsize',18,'fontname','arial')
ylabel('Log-likelihood improvement (bits/spk)','fontsize',18,'fontname','arial')

%%
% used_mod = gnmr2(2,end-4);
% used_mod = quadr_mod(4,end-5);
used_mod = gnmr2(2,8); %4,x
% used_mod = gnmr2(4,end);

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
use_thetas = find(~isnan(best_thetas));
avg_best_theta = circ_mean(best_thetas(use_thetas));
avg_best_x = mean(x_cents);
avg_best_y = mean(y_cents);

% f1 = figure;
figure
for ii = 1:n_filts
    
    %compute space-time projections
    gabor_filtmat = squeeze(gabor_fitvals(ii,:,:));
    cur_filtmat = reshape(pix_mat(:,ii),flen,fsdim);
    
    if max(abs(cur_filtmat(:))) > 0
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
    axis tight; %xlim([-3.1 3.1])
    xlim(nlx([1 end]))
    title('Module non-linearity')
    end
end
colormap(gray)

