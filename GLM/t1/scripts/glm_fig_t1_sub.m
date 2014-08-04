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

load ./sub_Xmat
stype='all';

sdim = 18;
fsdim = 18^2;
pids =1:fsdim;
flen = 7;
%% create XV data
[stimlen,k_len] = size(X);
nparts = 100;
nfold = 10;
partlen = floor(stimlen/nparts);
nxvparts = nparts/nfold;

%boundaries of parts
pbounds  = [(0:nparts-1)*partlen+1;(1:nparts)*partlen]';

%create Nfold different sets of XV and TR data
for i = 1:nfold
    xv_inds{i} = [];
    tr_inds{i} = [];
    
    cur_perm = randperm(nparts);
    cur_xv_parts = sort(cur_perm(1:nxvparts));
    cur_tr_parts = setdiff(1:nparts,cur_xv_parts);
    
    xv_new_inds = nan(stimlen,1);
    for j = 1:length(cur_xv_parts)
        cur_start = pbounds(cur_xv_parts(j),1);
        cur_stop = pbounds(cur_xv_parts(j),2);
        
        cur_inds = (cur_start:cur_stop) - cur_start + length(xv_inds{i}) + 1;
        xv_new_inds(cur_start:cur_stop) = cur_inds;
        
        xv_inds{i} = [xv_inds{i} cur_start:cur_stop];
    end
    
    tr_new_inds = nan(stimlen,1);
    for j = 1:length(cur_tr_parts)
        cur_start = pbounds(cur_tr_parts(j),1);
        cur_stop = pbounds(cur_tr_parts(j),2);
        
        cur_inds = (cur_start:cur_stop) - cur_start + length(tr_inds{i}) + 1;
        tr_new_inds(cur_start:cur_stop) = cur_inds;
        
        tr_inds{i} = [tr_inds{i} cur_start:cur_stop];
    end
end

xv = 1;
compids   = 1:k_len;

t = 13;
tcell = t;
tsbs      = 1+floor(aselspks{tcell}/dt);
tr_spbs = find(ismember(tr_inds{xv},tsbs));
spikebins = tr_spbs;
xv_spbs = find(ismember(xv_inds{xv},tsbs));

X_tr = X(tr_inds{xv},:);
%%
compids = 500;
[COEFF, SCORE, LATENT] = princomp(X_tr);
whitemat = diag(1./sqrt(LATENT)); %whitening (rescaling) transformation from training set
wstim = SCORE(:,1:compids)*whitemat(1:compids,1:compids); %whitened training data
%%
nneg =5;
npos =5;
spike_cond_stim = wstim(spikebins,:);
sta = mean(spike_cond_stim) - mean(wstim);
sta = sta/norm(sta);
proj_mat = sta'/(sta*sta')*sta;
stim_proj = wstim - wstim*proj_mat;
stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

stcs = [sta' stcs];
pix_conv_mat = COEFF(:,1:compids)';
stcs = (stcs'*pix_conv_mat)';
figure
plotfilterbank(stcs,sdim,1:fsdim);

%%
% nneg =5;
% npos =5;
% spike_cond_stim = X(spikebins,:);
% sta = mean(spike_cond_stim) - mean(X);
% sta = sta/norm(sta);
% proj_mat = sta'/(sta*sta')*sta;
% stim_proj = X - X*proj_mat;
% stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
% [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
% stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
% 
% stcs = [sta' stcs];
% figure
% plotfilterbank(stcs,sdim,1:fsdim);

%% TRY FITTING SEQUENCE OF STC MODELS WTIH INCREASING DIMENSIONALITY
%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 100;
defmod.lambda_L1x = 0;
defmod.lambda_dT = 100;

defmod.fsdim = fsdim;
defmod.pids = 1:fsdim;
defmod.SDIM = sdim;
basis = 'pix';
cur_ndims = 4;

init_signs = ones(cur_ndims,1);
init_kerns = randn(k_len,cur_ndims);

for xv = 1
    X_tr = X(tr_inds{xv},:); %whitened training data
    %spike bins during training data
    tr_spbs = find(ismember(tr_inds{xv},tsbs));
    spikebins = tr_spbs;
    
    [glm_stcb] = createGLM_quad(init_kerns,init_signs,defmod,basis);
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,X_tr);
    g_mat = X_tr*get_k_mat(glm_stcb);
    temp = fitWeights_full(glm_stcb,g_mat,spikebins,1);
    glm_stcb.const = temp.const;
    
    % glm_stcb.spk_nl = 'logexp';
        glm_stcb.image_type = '2d';
    quad_glm(xv) = fitGLM_lexp(glm_stcb,X_tr,spikebins,'tots',50,1e-4,1e-6);
end


%%
quad_glm2(xv) = adjust_all_reg(quad_glm(xv),'lambda_L1x',5);
quad_glm2(xv) = fitGLM_lexp(quad_glm2(xv),X_tr,spikebins,'tots',50,1e-4,1e-6);

%% TRY FITTING SEQUENCE OF STC MODELS WTIH INCREASING DIMENSIONALITY
%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 100;
defmod.lambda_L1x = 1;
defmod.lambda_dT = 100;

defmod.fsdim = fsdim;
defmod.pids = 1:fsdim;
defmod.SDIM = sdim;
basis = 'pix';
cur_ndims = 4;

init_signs = ones(cur_ndims,1);
init_kerns = randn(k_len,cur_ndims);
init_betas = 2*ones(cur_ndims,1);
init_thetas = zeros(cur_ndims,1);
for xv = 1
    X_tr = X(tr_inds{xv},:); %whitened training data
    %spike bins during training data
    tr_spbs = find(ismember(tr_inds{xv},tsbs));
    spikebins = tr_spbs;
    
    [glm_stcb] = createGLM_lexp(init_kerns,init_signs,init_betas,init_thetas,defmod,basis);
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,X_tr);
    g_mat = X_tr*get_k_mat(glm_stcb);
    temp = fitWeights_full(glm_stcb,g_mat,spikebins,1);
    glm_stcb.const = temp.const;
    
    % glm_stcb.spk_nl = 'logexp';
        glm_stcb.image_type = '2d';
    lexp_glm(xv) = fitGLM_lexp(glm_stcb,X_tr,spikebins,'tots',50,1e-4,1e-6);
end

%%
lexp_glm2.mods(2) = [];
% lexp_glm2 = adjust_all_reg(lexp_glm(xv),'lambda_dX',100);
% lexp_glm2 = adjust_all_reg(lexp_glm2,'lambda_dT',100);
% lexp_glm2 = adjust_all_reg(lexp_glm2,'lambda_L1x',2);
[lexp_glm2,norm_vals] = normalizeRFs_full(lexp_glm2,X_tr);
g_mat = X_tr*get_k_mat(lexp_glm2);
lexp_glm2 = fitWeights_lexp(lexp_glm2,g_mat,spikebins,0);
for i = 1:length(lexp_glm2.mods)
    lexp_glm2.mods(i).lambda_dX = defmod.lambda_dX*lexp_glm2.mods(i).w^2;
    lexp_glm2.mods(i).lambda_dT = defmod.lambda_dT*lexp_glm2.mods(i).w^2;
    glm_quadf.mods(i).lambda_L1x = defmod.lambda_L1x*lexp_glm2.mods(i).w;
end
lexp_glm2 = fitGLM_lexp(lexp_glm2,X_tr,spikebins,'tots',50,1e-4,1e-6);
%%
lexp_glm2(xv) = adjust_all_reg(lexp_glm(xv),'lambda_dX',5);
quad_glm2(xv) = fitGLM_lexp(quad_glm2(xv),X_tr,spikebins,'tots',50,1e-4,1e-6);



%%
lexp_xvLL = getLLGLM_lexp(lexp_glm2,X(xv_inds{xv},:),xv_spbs,'tots')
quad_xvLL = getLLGLM_lexp(quad_glm,X(xv_inds{xv},:),xv_spbs,'tots')
%%
% used_mod = quad_glm(xv);
used_mod = lexp_glm2;
pix_mat = get_pix_mat(used_mod);
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

for i = 1:n_filts
    nly_mat(i,:) = used_mod.mods(i).nly;
    [gabor_fits(i,:),gabor_fitvals(i,:,:)] = get_gaborfits_allslices(pix_mat(:,i),flen,sdim);
end

nlx = used_mod.mods(1).nlx;
[~,best_slice_ids] = get2dMaxSlices(pix_mat,flen,sdim,pids);

f1 = figure;
for ii = 1:n_filts
    
    %compute space-time projections
    gabor_filtmat = squeeze(gabor_fitvals(ii,:,:));
    cur_filtmat = reshape(pix_mat(:,ii),flen,fsdim);
    
    filt_projprof = project_2drf(cur_filtmat(:),gabor_fits(ii,best_slice_ids(ii)).theta,sdim);
    subplot(n_filts,3,(ii-1)*3+1)
    imagesc(reshape(cur_filtmat(best_slice_ids(ii),:),sdim,sdim));colormap(jet);
    zmax = max(abs(cur_filtmat(best_slice_ids(ii),:))); caxis([-zmax zmax]);
title('Best time-slice')

    subplot(n_filts,3,(ii-1)*3+2)
    imagesc(1:sdim,-(0:5)*0.02,filt_projprof);
    xm = max(abs(filt_projprof(:)));
%     caxis([-xm xm]);
    title('Space-time projection')
    subplot(n_filts,3,(ii-1)*3+3)
    plot(nlx,nly_mat(ii,:),'k')
    axis tight; xlim([-3.1 3.1])
    title('Module non-linearity')
    
end
colormap(gray)
% 
% set(f1,'PaperUnits','centimeters');
% set(f1, 'PaperSize', [30 7*length(rotated_mod{t}.mods)]);
% set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
% 

