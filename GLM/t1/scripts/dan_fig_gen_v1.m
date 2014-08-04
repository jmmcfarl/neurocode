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
ncomps = 500;
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
pix_conv_mat = diag(sqrt(scorevars(compids)))*coefs(:,compids)';
kern_conv_mat = coefs(:,compids)*diag(1./sqrt(scorevars(compids)));

WX        = scores(:,compids)*diag(1./sqrt(scorevars(compids)));
spikebins = tsbs(tsbs>flen & tsbs<(size(WX,1)+flen-1))-flen+1 ;
WS        = WX(spikebins,:);
rsta      = mean(WS) - mean(WX);

stvcv = cov(WS);  utvcv = cov(WX);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);

npos=10; nneg=10;
stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
nposdims = stcmods; nnegdims = 0;
% nposdims = 5; nnegdims = 5;
posdims = 1:nposdims-1; negdims = 1:nnegdims;
% STCbvs = [sta' stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
        kimages = [rsta',STCbvs]'*diag(sqrt(scorevars(compids)))*coefs(:,compids)';
f1 = figure('Name',stype); plotfilterbank(kimages',32,pids)

%add in STA
STCbvs = [rsta' STCbvs];

[klen,Nstcbvs] = size(STCbvs);
flen = 6;
stimlen = size(WX,1);

%% INITIALIZE MODEL
nmods = Nstcbvs;
mod_signs = ones(nmods,1);
dim_signs = ones(Nstcbvs,1);

unused_stcs = (nmods+1):Nstcbvs;
%initialize on STC dims
STCcf_0 = eye(Nstcbvs);
STCcf_0(:,unused_stcs) = [];

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
defmod.lambda_L1x = 30; %40
defmod.lambda_dT = 10;
defmod.pids = pids;

%define NL initializations: "lin, threshlin, pquad, nquad"
clear init_nls nltypes
init_nls{1} = 'lin'; for i = 2:nmods; init_nls{i} = 'pquad'; end;
%define NL types: "uncon, lin, threshlin, quad"
nltypes{1} = 'lin'; for i = 2:nmods; nltypes{i} = 'quad'; end;
% %define NL initializations: "lin, threshlin, pquad, nquad"
% for i = 1:nmods; init_nls{i} = 'threshlin'; end;
% %define NL types: "uncon, lin, threshlin, quad"
% for i = 1:nmods; nltypes{i} = 'threshlin'; end;

% basis = 'pix';
basis = 'white';
glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
[glm_stcb,sub_medians] = median_sub_filters(glm_stcb);

[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
glm_stcb.lambdaW = 10; %sparseness on model weights
glm_stcb.image_type = '2d';

%%
full_glm = fitNLHI2d_fullbf(glm_stcb,WX,spikebins,'tots');
f2 = plot2d_mod(full_glm,0);

%%
%get rid of bad filter
full_glm_3filt = full_glm; full_glm_3filt.mods(4) = [];
for i = 1:nmods-1; full_glm_3filt.mods(i).nltype='uncon'; full_glm_3filt.mods(i).lnl2 = 100; end;

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
used_mod = full_glm4;

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