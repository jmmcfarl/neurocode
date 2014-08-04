%% Load Data

clear all;
addpath('~/James_scripts/GLM')
addpath('~/James_scripts/GLM/t1')
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/Timm/tim1/functions/')

cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

% load spks7576-ns.mat; load PCScores-ns.mat; stype='ns'
% load spks7576-pn.mat; load PCScores-pn.mat; stype='pn'
% load spks7576-ps.mat; load PCScores-ps.mat; stype='ps'

%load z-scored PC data
cd ~/Data/blanche/matlabdata/
load spks7576-all.mat;
load PCScores_z-all.mat;
stype='all';

stcmods = 4;
ncomps = 400;
scorevars(1001:end) = [];
coefs(:,1001:end) = [];
scores(:,1001:end) = [];

%% Compute whitened data, and STA/STCs
cd '/Users/James/James_scripts/GLM/t1/allcell_fits/'
% for t = 1:27
t = 12;
tcell = t;
pids =1:1024;
tsbs      = 1+floor(aselspks{tcell}/dt);

compids   = 1:ncomps;

% scalmat = diag(1./sqrt(scorevars(compids)));
scalmat = diag(ones(1,ncomps));
WX        = scores(:,compids)*scalmat;
% WXw = WX*diag(1./sqrt(scorevars(compids)));

spikebins = tsbs(tsbs>flen & tsbs<(size(WX,1)+flen-1))-flen+1 ;
WS        = WX(spikebins,:);
rsta      = mean(WS) - mean(WX);

stvcv = cov(WS);  utvcv = cov(WX);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);

npos=10; nneg=10;
stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
nposdims = 6; nnegdims =6;
posdims = 1:nposdims; negdims = 1:nnegdims;
STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace

% pix_conv_mat = diag(sqrt(scorevars(compids)))*coefs(:,compids)';
% kern_conv_mat = coefs(:,compids)*diag(1./sqrt(scorevars(compids)));
pix_conv_mat = coefs(:,compids)';
kern_conv_mat = coefs(:,compids);

kimages = [rsta',STCbvs]'*pix_conv_mat;

% STCbvs = [rsta',stcs(:,posdims) rstcs(:,negdims)]';
% STCbvsw = STCbvs*diag(1./sqrt(scorevars(compids)));
% pix_conv_mat = diag(sqrt(scorevars(compids)))*coefs(:,compids)';
% kern_conv_mat = coefs(:,compids)*diag(1./sqrt(scorevars(compids)));
% kimages = STCbvsw'*pix_conv_mat;

f1 = figure('Name',stype); plotfilterbank(kimages',32,pids)
figure
plot(diag(evals),'.')

% %add in STA
% STCbvs = [rsta' STCbvs];
% STCbvsw = STCbvsw';
% [klen,Nstcbvs] = size(STCbvsw);
% flen = 6;
% stimlen = size(WX,1);
% 
%% REFINE STC (ML STC)
used_stc_dims = [1:6];
STCbasis = STCbvs(:,used_stc_dims);
[klen,Nstcbvs] = size(STCbasis);
%initialize on STC dims
STCcf_0 = eye(Nstcbvs);
nmods = Nstcbvs;

%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;

%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 500; %350
defmod.lambda_L1x = 10; %40
defmod.lambda_dT = 10;
defmod.pids = pids;

%define NL initializations: "lin, threshlin, pquad, nquad"
clear init_nls nltypes
init_nls{1} = 'lin';
for i = 2:6; init_nls{i} = 'pquad'; end;
% for i = 8:nmods; init_nls{i} = 'nquad'; end;
nltypes{1} = 'lin';
for i = 2:nmods; nltypes{i} = 'quad'; end;

basis = 'pix';
defmod_t = defmod;
defmod_t.fsdim = ncomps;
defmod_t.pids = 1:ncomps;
defmod_t.lambda_dX = 0;
defmod_t.lambda_dT = 0;
defmod_t.lambda_L1x = 0;
glm_stcb = createGLM2d_fullbf(STCbasis,STCcf_0,[],[],defmod_t,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
glm_stcb.image_type = '2d';
[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
stc_glm = fitstc_fullbf(glm_stcb,WX,spikebins,'tots');
ker_filts = get_k_mat(stc_glm);

ker_filts(:,5:6) = [];
ker_filts = [ker_filts -ker_filts];
ker_pix = ker_filts'*pix_conv_mat;
nmods = size(ker_filts,2);
for i = 1:nmods; init_nls{i} = 'threshlin'; end;
for i = 1:nmods; nltypes{i} = 'threshlin'; end;

basis = 'white';
glm_stcb = createGLM2d_fullbf(STCbasis,ones(size(STCbasis,2),nmods),pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
glm_stcb.image_type = '2d';
for i = 1:nmods
    temp_k = ker_filts(:,i);
    temp_pix = ker_pix(i,:);
    glm_stcb.mods(i).k = temp_k(:);
    glm_stcb.mods(i).pix = temp_pix(:);
end
[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
% pix_stc_glm = fitstc_fullbf(glm_stcb,WX,spikebins,'tots',2);
% [nstc_glm,norm_vals] = normalizeRFs_full(stc_glm,WX);

full_glm = fitNLHI2d_fullbf(glm_stcb,WX,spikebins,'tots',2);

%% REFINE STC (THRESHLIN MODEL)
% used_stc_dims = [2:5 7];
% STCbasis = STCbvs(:,used_stc_dims);
% STCbasis = [STCbasis -STCbasis(:,2:end)];
% 
% [klen,Nstcbvs] = size(STCbasis);
% %initialize on STC dims
% STCcf_0 = eye(Nstcbvs);
% nmods = Nstcbvs;
% 
% defmod.h(1:end-1) = []; %eliminate PSC
% defmod.lnl = 0;
% defmod.lh = 0;
% defmod.lnl2 = 0;
% defmod.lh2 = 0;
% defmod.nlcon = 0;
% defmod.nlmon = 0;
% defmod.locLambda = 0;
% defmod.lambda_dX = 500; %350
% defmod.lambda_L1x = 10; %40
% defmod.lambda_dT = 0;
% defmod.pids = pids;
% 
% % clear init_nls nltypes
% % init_nls{1} = 'lin'; for i = 2:nmods; init_nls{i} = 'threshlin'; end;
% % nltypes{1} = 'lin'; for i = 2:nmods; nltypes{i} = 'threshlin'; end;
% clear init_nls nltypes
% for i = 1:nmods; init_nls{i} = 'threshlin'; end;
% for i = 1:nmods; nltypes{i} = 'threshlin'; end;
% 
% 
% % basis = 'pix';
% basis = 'white';
% glm_stcb = createGLM2d_fullbf(STCbasis,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
% [glm_stcb,sub_medians] = median_sub_filters(glm_stcb);
% 
% [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
% % glm_stcb.lambdaW = 10; %sparseness on model weights
% glm_stcb.image_type = '2d';
% full_glm = fitNLHI2d_fullbf(glm_stcb,WX,spikebins,'tots');


%%
%get rid of bad filter
full_glm_3filt = full_glm; full_glm_3filt.mods([4 5 8 9]) = [];
nmods = length(full_glm_3filt.mods);
for i = 1:nmods; full_glm_3filt.mods(i).nltype='uncon'; full_glm_3filt.mods(i).lnl2 = 100; end;

full_glm2 = fitNLHI2d_fullbf(full_glm_3filt,WX,spikebins,'tots');
f2 = plot2d_mod(full_glm2,0);

%%
new_mod_2 = full_glm2.mods(2); 
new_mod_2.k = -new_mod_2.k; new_mod_2.pix = -new_mod_2.pix;
new_mod_2.nly = fliplr(new_mod_2.nly);
new_mod_3 = full_glm2.mods(3); 
new_mod_3.k = -new_mod_3.k; new_mod_3.pix = -new_mod_3.pix;
new_mod_3.nly = fliplr(new_mod_3.nly);
neg_vals = find(full_glm2.mods(2).nlx < 0);
pos_vals = find(full_glm2.mods(2).nlx > 0);
full_glm3 = full_glm2;
full_glm3.mods(2).nly(neg_vals) = 0; 
full_glm3.mods(3).nly(neg_vals) = 0; 
full_glm3.mods = [full_glm3.mods new_mod_2 new_mod_3];
full_glm3.mods(4).nly(neg_vals) = 0;
full_glm3.mods(5).nly(neg_vals) = 0;

%%
for i = 1:length(full_glm3.mods); full_glm3.mods(i).nlmon = 1; end;
full_glm4 = fitNLHI2d_fullbf(full_glm3,WX,spikebins,'tots');
f2 = plot2d_mod(full_glm4,0);

%%
used_mod = pix_stc_glm;
used_mod.mods([5 6]) = [];
basis_vecs = get_k_mat(used_mod);
n_bvs = size(basis_vecs,2);
nmods = 4;
mod_signs = ones(nmods,1);
dim_signs = ones(n_bvs,1);
unused_stcs = (nmods+1):n_bvs;

flen = 6;
%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 100;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 1;
defmod.locLambda = 500;
defmod.lambda_dX = 200; %350
defmod.lambda_L1x = 0; %40
defmod.lambda_dT = 10;
defmod.pids = 1:fsdim;
defmod.SDIM = SDIM;
defmod.fsdim = fsdim;

%define NL initializations: "lin, threshlin, pquad, nquad"
clear init_nls nltypes
for i = 1:nmods; init_nls{i} = 'threshlin'; end;
%define NL types: "uncon, lin, threshlin, quad"
for i = 1:nmods; nltypes{i} = 'uncon'; end;

Nstcbf = size(basis_vecs,2);
klen = size(pix_conv_mat,2);
flen = klen/fsdim;
kern_output = WX*basis_vecs;

clear init_vals all_filtproj all_initproj cur_LL cur_LP fin_vals dist_trav *_lp rotbv_mod
max_reps = 300;
min_reps = 10;
min_LM_fract = 0.95;
eps = 0.002;
cur_reps = 0;
used_cfs = [];
LL_vals = [];
LP_vals = [];
smallest_dists = [];
is_optimized = [];
rotbv_mod = [];
for r = 1:10
r
    %points uniformly distributed on the set of Nmod p-spheres (p=n_bvs)
    STCcf_0 = randn(n_bvs,nmods);
    for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    
    defmod.image_type = '2d';
    white_props.basis = 'white';
    white_props.pix_conv_mat = pix_conv_mat;
    white_props.kern_conv_mat = kern_conv_mat;
    glm_stcb = create2dGLM0_stcb(basis_vecs,STCcf_0,defmod,mod_signs,dim_signs,nltypes,init_nls,'test',white_props); %initialize
    glm_stcb.lambdaW = 50;
    
    %determine LL and LP at current filter point
    glm_stcb = fitWeights_stcb_nonlpsc(glm_stcb,kern_output*STCcf_0,spikebins,1,1e-3);
    [ll0, ll0p] = getLLGLM2d_STCBF_nonlpsc(glm_stcb,kern_output,spikebins,'none');
    
    rotbv_mod = [rotbv_mod; fitNLHI_stcb2d_nopsc(glm_stcb,WX,spikebins,'none',6,2)];
    LL_vals = [LL_vals rotbv_mod(end).LL];
    LP_vals = [LP_vals rotbv_mod(end).LP];
end
[~,best_mod] = min(LL_vals);

%%
used_mod = rotbv_mod(1);
w_vals = arrayfun(@(x) x.w,used_mod.mods);
used_mod.mods(w_vals==0) = [];
%%
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
    [gabor_fits(i,:),gabor_fitvals(i,:,:)] = get_gaborfits_allslices(pix_mat(:,i),6,32);
end

nlx = used_mod.mods(1).nlx;
[~,best_slice_ids] = get2dMaxSlices(pix_mat,flen,sdim,pids);
f1 = figure;
for ii = 1:n_filts
      
    %compute space-time projections
    gabor_filtmat = squeeze(gabor_fitvals(ii,:,:));
    cur_filtmat = reshape(pix_mat(:,ii),6,fsdim);
    
    
    filt_projprof = project_2drf(cur_filtmat(:),gabor_fits(ii,best_slice_ids(ii)).theta,sdim);    
    subplot(n_filts,3,(ii-1)*3+1)
    imagesc(reshape(cur_filtmat(best_slice_ids(ii),:),sdim,sdim));colormap(jet);
    zmax = max(abs(cur_filtmat(best_slice_ids(ii),:))); caxis([-zmax zmax]);
    title('Best time-slice')
    
%     subplot(n_filts,3,(ii-1)*3+2)
%     imagesc(reshape(gabor_filtmat(best_slice_ids(ii),:),sdim,sdim));colormap(jet);
%     zmax = max(abs(gabor_filtmat(best_slice_ids(ii),:))); caxis([-zmax zmax]);
%     title('Best time-slice (Gabor-fit)')    
        
    subplot(n_filts,3,(ii-1)*3+2)
    imagesc(1:sdim,-(0:5)*0.02,filt_projprof);
    xm = max(abs(filt_projprof(:)));
    caxis([-xm xm]);
    title('Space-time projection')
    
    subplot(n_filts,3,(ii-1)*3+3)
    plot(nlx,nly_mat(ii,:),'k')
    axis tight; xlim([-3.1 3.1])
    title('Module non-linearity')
    
end