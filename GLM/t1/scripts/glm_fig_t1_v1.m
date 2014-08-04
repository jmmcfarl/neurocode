%% Load Data

clear all;
addpath('~/James_scripts/GLM')
addpath('~/James_scripts/GLM/t1')
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/Timm/tim1/functions/')
addpath(genpath('~/James_Scripts'));

cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

% load spks7576-ns.mat; load PCScores-ns.mat; stype='ns'
% load spks7576-pn.mat; load PCScores-pn.mat; stype='pn'
% load spks7576-ps.mat; load PCScores-ps.mat; stype='ps'

%load z-scored PC data
cd ~/Data/blanche/matlabdata/
load spks7576-all.mat;

% cd ~/Data/blanche/rec_76/matlabdata/
% load PCScores_f9_z-all.mat;
% stype='all';
load PCScores_z-all.mat;
stype='all';

%throw away extra stimulus dimensions that we definitely won't be using
scorevars(1001:end) = [];
coefs(:,1001:end) = [];
scores(:,1001:end) = [];

pids =1:1024;

%% create XV data
stimlen = size(scores,1);
nparts = 20;
partlen = floor(stimlen/nparts);
nfold = 10;
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
ncomps = 300;
compids   = 1:ncomps;
scorevars = var(scores(tr_inds{xv},:)); %variance of PCA dimensions within the training set
whitemat = diag(1./sqrt(scorevars(compids))); %whitening (rescaling) transformation from training set
WX        = scores(tr_inds{xv},compids); %whitened training data
pix_conv_mat = coefs(:,compids)';
kern_conv_mat = coefs(:,compids);
X_avg = mean(WX); %average whitened training stim
NT_x = size(WX,1); %number of training samples

%% Cross-validation
WX_xv  = scores(xv_inds{xv},compids);

%% Compute whitened data, and STA/STCs
cd '/Users/James/James_scripts/GLM/t1/allcells_fits_mlstc_v4/'
% for t = 1:20

t = 13;
tcell = t;
tsbs      = 1+floor(aselspks{tcell}/dt);

%spike bins during training data
tr_spbs = find(ismember(tr_inds{xv},tsbs));
spikebins = tr_spbs(tr_spbs>flen & tr_spbs<(size(WX,1)+flen-1))-flen+1 ;

%spike bins for XV set
xv_spbs = find(ismember(xv_inds{xv},tsbs));
xv_spikebins = xv_spbs(xv_spbs>flen & xv_spbs<(size(WX_xv,1)+flen-1))-flen+1 ;

% compute STE and STA
WS        = WX(spikebins,:);
rsta      = mean(WS) - X_avg;
rsta = rsta*whitemat; %correct for unequal variance bias in STA

%project STA out of STE
proj_mat = rsta'*inv(rsta*rsta')*rsta;
WXp = WX - WX*proj_mat;
WSp = WXp(spikebins,:);

%Compute STC
stvcv = cov(WSp); %covaraince of STE
%[evecs,evals] = eig(stvcv-eye(ncomps)); %assume exact whitening...
[evecs,evals] = eig(stvcv-cov(WXp));
evs{t} = diag(evals);

npos=4; nneg=4; %retain a handful of positive and negative dimensions to check
stcs  = evecs(:,[1:nneg,length(evs{t})-npos+1:end]); stcs  = stcs(:,end:-1:1);
rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
STCbvs = [stcs(:,1:npos) rstcs(:,1:nneg)];
STCbvs = (STCbvs'*whitemat)'; %correct for unequal variance bias in STC dims
STCbvs = [rsta' STCbvs];

kimages = STCbvs'*pix_conv_mat;
f1 = figure('Name',stype); plotfilterbank(kimages',32,pids)
set(f1,'PaperUnits','centimeters');
set(f1, 'PaperSize', [30 70]);
set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
%     fname = sprintf('Cell%d_Init_STC',t);
%     print(f1,'-dpng',fname);close all


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
defmod.lambda_dX = 50;
defmod.lambda_L1x = 2;
defmod.lambda_dT = 5;
% defmod.lambda_dX = 15000;
% defmod.lambda_L1x = 25;
% defmod.lambda_dT = 1500;
% defmod.lambda_dX = 4000;
% defmod.lambda_L1x = 20;
% defmod.lambda_dT = 400;
defmod.lambda_dX = 6000;
defmod.lambda_L1x = 25;
defmod.lambda_dT = 600;

defmod.fsdim = fsdim;
defmod.pids = 1:fsdim;
defmod.SDIM = sdim;
basis = 'white';

% clear stc_glm init_nls nltypes
%first fit sta model
cur_ndims = 4;
cur_basis = STCbvs(:,1:cur_ndims); %just use STA
STCcf_0 = eye(cur_ndims);

init_signs = ones(cur_ndims,1);
% [glm_stcb] = createGLM_tlin(cur_basis*STCcf_0,init_signs,defmod,basis,pix_conv_mat,kern_conv_mat)

init_betas = 2*ones(cur_ndims,1);
init_thetas = zeros(cur_ndims,1);
[glm_stcb] = createGLM_lexp(cur_basis*STCcf_0,init_signs,init_betas,init_thetas,defmod,basis,pix_conv_mat,kern_conv_mat)
[glm_stcb,norm_vals1] = normalizeRFs_full(glm_stcb,WX);
g_mat = WX*get_k_mat(glm_stcb);
temp = fitWeights_full(glm_stcb,g_mat,spikebins,1);
glm_stcb.const = temp.const;

% glm_stcb.spk_nl = 'logexp';
glm_stcb.image_type = '2d';
tlin_glm2 = fitGLM_lexp(glm_stcb,WX,spikebins,'tots',70,1e-4,1e-6);
tlin_xvLL = getLLGLM_FULL2d(tlin_glm,WX_xv,xv_spikebins,'tots') %determine XVLL
% lexp_glm = fitGLM_lexp(glm_stcb,WX,spikebins,'tots',70,1e-4,1e-6);
% lexp_xvLL = getLLGLM_FULL2d(lexp_glm,WX_xv,xv_spikebins,'tots') %determine XVLL

%%
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 5000;
defmod.lambda_L1x = 25;
defmod.lambda_dT = 500;
% defmod.lambda_dX = 1000;
% defmod.lambda_L1x = 20;
% defmod.lambda_dT = 100;
% 
defmod.fsdim = fsdim;
defmod.pids = 1:fsdim;
defmod.SDIM = sdim;
basis = 'white';

% clear stc_glm init_nls nltypes
%first fit sta model
cur_ndims = 3;
cur_basis = STCbvs(:,1:cur_ndims); %just use STA
STCcf_0 = eye(cur_ndims);

init_signs = ones(cur_ndims,1);
[glm_stcb] = createGLM_quad(cur_basis*STCcf_0,init_signs,defmod,basis,pix_conv_mat,kern_conv_mat)

[glm_stcb,norm_vals1] = normalizeRFs_full(glm_stcb,WX);
g_mat = WX*get_k_mat(glm_stcb);
temp = fitWeights_full(glm_stcb,g_mat,spikebins,1);
glm_stcb.const = temp.const;

% glm_stcb.spk_nl = 'logexp';
glm_stcb.image_type = '2d';
glm_stcb.mods(1).dX = 5000000;
glm_stcb.mods(1).L1x = 7000;
glm_stcb.mods(1).dT = 500000;

quad_glm = fitGLM_lexp(glm_stcb,WX,spikebins,'tots',70,1e-4,1e-6);
quad_xvLL = getLLGLM_FULL2d(quad_glm,WX_xv,xv_spikebins,'tots') %determine XVLL

%%
cd /Users/James/Documents/GNM_paper/fig_info
save t1_gnm_fig_v3 quad* tlin* *_parts 

%% Try fitting overcomplete threshlin  model
used_mod = stc_glm5;
sdim = sqrt(SDIM);
basis_vecs = get_k_mat(used_mod);
% FF = mod(1:(flen*fsdim),flen)+1;
% basis_pixs = get_pix_mat(used_mod);
% basis_pixs(FF ~= 1,:) = 0;
% basis_pixs1 = basis_pixs; basis_pixs1(basis_pixs1 < 0) = 0;
% basis_pixs2 = basis_pixs; basis_pixs2(basis_pixs2 > 0) = 0;
% basis_pixs = [basis_pixs1 basis_pixs2];
% 
% basis_vecs = basis_pixs;

n_bvs = size(basis_vecs,2);


%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 0; %350
defmod.lambda_L1x = 0; %40
defmod.lambda_dT = 0;
defmod.pids = 1:fsdim;
defmod.SDIM = SDIM;
defmod.fsdim = fsdim;
defmod.image_type = '2d';
white_props.basis = 'white';
white_props.pix_conv_mat = pix_conv_mat;
white_props.kern_conv_mat = kern_conv_mat;
% white_props.basis = 'pix';
% white_props.pix_conv_mat = [];
% white_props.kern_conv_mat = [];

nmods = [3:6];
LL_vals = [];
LP_vals = [];
rotbv_mod = [];
xv_vals = [];
for n = 1:length(nmods)
    
    %define NL initializations: "lin, threshlin, pquad, nquad"
    clear init_nls nltypes
    for i = 1:nmods(n); init_nls{i} = 'threshlin'; end;
    %define NL types: "uncon, lin, threshlin, quad"
    for i = 1:nmods(n); nltypes{i} = 'threshlin'; end;
    
    %points uniformly distributed on the set of Nmod p-spheres (p=n_bvs)
    STCcf_0 = randn(n_bvs,nmods(n));
    for i = 1:nmods(n); STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    
    glm_stcb = create2dGLM0_stcb(basis_vecs,STCcf_0,defmod,ones(nmods(n),1),ones(nmods(n),1),nltypes,init_nls,'test',white_props); %initialize
    glm_stcb.lambdaW = 0;
    glm_stcb.spk_nl = 'logexp';
    rotbv_mod = [rotbv_mod; fitNLHI_stcb2d_nonlpsc(glm_stcb,WX,spikebins,'tots',8,2)];
    LL_vals = [LL_vals rotbv_mod(end).LL];
    LP_vals = [LP_vals rotbv_mod(end).LP];
    xv_vals = [xv_vals getLLGLM_FULL2d(rotbv_mod(end),WX_xv,xv_spikebins,'none')];
    
end

%%
cnmods = 6;
for i = 1:cnmods; init_nls{i} = 'threshlin'; end;
%define NL types: "uncon, lin, threshlin, quad"
for i = 1:cnmods; nltypes{i} = 'threshlin'; end;

%points uniformly distributed on the set of Nmod p-spheres (p=n_bvs)
STCcf_0 = randn(n_bvs,cnmods);
for i = 1:cnmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
glm_stcb = create2dGLM0_stcb(basis_vecs,STCcf_0,defmod,ones(cnmods,1),ones(cnmods,1),nltypes,init_nls,'test',white_props); %initialize
glm_stcb.lambdaW = 0;
glm_stcb.spk_nl = 'logexp';

npos = 4; nneg = 2;
for i = npos+1:cnmods
    glm_stcb.mods(i).w = -1;
end
new_mod = fitNLHI_stcb2d_nonlpsc(glm_stcb,WX,spikebins,'tots',8,2);
new_mod_xvLL = getLLGLM_FULL2d(new_mod,WX_xv,xv_spikebins,'none');

%%
used_mod = tlin_glm(xv);
pix_mat = get_pix_mat(used_mod);
n_filts = length(used_mod.mods);

smooth_mat = epanechikov(5)'*epanechikov(5);
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
    [gabor_fits(i,:),gabor_fitvals(i,:,:)] = get_gaborfits_allslices(pix_mat(:,i),flen,32);
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
