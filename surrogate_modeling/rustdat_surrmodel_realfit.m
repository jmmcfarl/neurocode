clear all;
%% load data
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/James_scripts/GLM/')

pars = load('~/Data/rust/infos/StandardParametersRust.mat');
defmod = pars.defmod;
defmod.h(1:end-1) = [];%restrict PSC term to delta function
flen = pars.flen;
hilen = length(defmod.h);
foff = flen + pars.hilen;
ooptions = optimset('MaxFunEvals',100000,'MaxIter',10000);

% datdir = '/Users/timm/myrepo/workspace/DataRepository/rust/stcbar/DATA/';
tuid = '44-29'; hassta=0;  npos=12; nneg=12;
datdir = '~/Data/rust/stcbar/Data/';
cd ~/Data/rust
dataname = sprintf('stcbf_data-%s.mat',tuid)
eval(['load ' dataname]);
dt = 0.01;
sdim = 24;
flen = 14;

%% load model
cd ~/James_scripts/surrogate_modeling/
load ./rust4429_finmod_14filt.mat

%% create XV data
[stimlen,sdim] = size(stim);
nparts = 18;
partlen = floor(stimlen/nparts);
nfold = 3;
nxvparts = nparts/nfold;
ntrparts = nparts - nxvparts;
NXV = nxvparts*partlen;
NTR = ntrparts*partlen;

%boundaries of parts
pbounds  = [(0:nparts-1)*partlen+1;(1:nparts)*partlen]';
for i = 1:nfold
    xv_inds{i} = [];
    xv_spkbns{i} = [];
    tr_inds{i} = [];
    tr_spkbns{i} = [];
    
    xv_spks = [];
    xv_new_inds = nan(stimlen,1);
    for j = 1:length(cur_xv_parts)
        cur_start = pbounds(cur_xv_parts(j),1);
        cur_stop = pbounds(cur_xv_parts(j),2);
        xv_spks = [xv_spks; spikebins(spikebins >= cur_start & spikebins < cur_stop)];
        
        cur_inds = (cur_start:cur_stop) - cur_start + length(xv_inds{i}) + 1;
        xv_new_inds(cur_start:cur_stop) = cur_inds;
        
        xv_inds{i} = [xv_inds{i} cur_start:cur_stop];
    end
    xv_spkbns{i} = xv_new_inds(xv_spks);
    
    tr_spks = [];
    tr_new_inds = nan(stimlen,1);
    for j = 1:length(cur_tr_parts)
        cur_start = pbounds(cur_tr_parts(j),1);
        cur_stop = pbounds(cur_tr_parts(j),2);
        tr_spks = [tr_spks; spikebins(spikebins >= cur_start & spikebins < cur_stop)];
        
        cur_inds = (cur_start:cur_stop) - cur_start + length(tr_inds{i}) + 1;
        tr_new_inds(cur_start:cur_stop) = cur_inds;
        
        tr_inds{i} = [tr_inds{i} cur_start:cur_stop];
    end
    tr_spkbns{i} = tr_new_inds(tr_spks);
end

xv = 1; %use only one XV set
cur_tr_stim = stim(tr_inds{xv},:);
cur_tr_stimemb = makeStimRows(cur_tr_stim,flen);
cur_tr_spkbns = tr_spkbns{xv};

%% generate surrogate data
pix_mat = get_pix_mat(fin_glm);

%internal filter outputs of each module
gen_functions = cur_tr_stimemb*pix_mat;

kx    = fin_glm.const; %initialize kx with model constant term
nmods = length(fin_glm.mods);
NT = size(cur_tr_stimemb,1);

for imod = 1:nmods %loop across NL modules
    
    fg = gen_functions(:,imod);
    fg = nlin_proc_stim(fg, fin_glm.mods(imod).nly, fin_glm.mods(imod).nlx );
    kx = kx + fg*fin_glm.mods(imod).w;
    
end
r = log(1+exp(kx)); %apply spiking NL tp get predicted rate

spikes = poissrnd(r);
rbins = (find(spikes>0.5));
nsp = spikes(rbins);
spk_vals = unique(spikes); spk_vals(spk_vals==0) = [];
surr_tr_spikebins = [];
for i = 1:length(spk_vals)
    cur_set = find(spikes == spk_vals(i));
    surr_tr_spikebins = [surr_tr_spikebins; repmat(cur_set(:),spk_vals(i),1)];
end
fprintf('Nspks: %d\n',length(surr_tr_spikebins));

%% STC analysis
%for surr data
nfilts = 10;
npos = nfilts; nneg = nfilts;
spike_cond_stim = cur_tr_stimemb(surr_tr_spikebins,:);
sta      = mean(spike_cond_stim) - mean(cur_tr_stimemb);
sta = sta/norm(sta);
stvcv = cov(spike_cond_stim);  utvcv = cov(cur_tr_stimemb);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
figure
subplot(3,nfilts,1)
imagesc(reshape(sta,flen,sdim));
for i = 1:nfilts
    subplot(3,nfilts,nfilts+i)
    imagesc(reshape(stcs(:,i),flen,sdim));
end
for i = 1:nfilts
    subplot(3,nfilts,2*nfilts+i)
    imagesc(reshape(stcs(:,nfilts+i),flen,sdim));
    colormap(gray)
end
%
%for real data
spike_cond_stim = cur_tr_stimemb(cur_tr_spkbns,:);
sta      = mean(spike_cond_stim) - mean(cur_tr_stimemb);
sta = sta/norm(sta);
stvcv = cov(spike_cond_stim);  utvcv = cov(cur_tr_stimemb);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
figure
subplot(3,nfilts,1)
imagesc(reshape(sta,flen,sdim));
for i = 1:nfilts
    subplot(3,nfilts,nfilts+i)
    imagesc(reshape(stcs(:,i),flen,sdim));
end
for i = 1:nfilts
    subplot(3,nfilts,2*nfilts+i)
    imagesc(reshape(stcs(:,nfilts+i),flen,sdim));
    colormap(gray)
end

%%
rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
nposdims = 6; nnegdims = 3;
posdims = 1:nposdims; negdims = 1:nnegdims;
STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace

stimlen = size(cur_tr_stim,1);
Nstcbf = size(STCbvs,2);
kern_output = cur_tr_stimemb*STCbvs;

% % for XV data
% [xvstimlen,sdim] = size(cur_xv_stim);
%
% %precompute the stimulus filtered by each STC kernel
% [klen,cur_Nstcbf] = size(STCbvs);
% flen = klen/sdim;
% xvkern_output = cur_xv_stimemb*STCbvs;

%% NOW FIND BEST OBLIQUE ROTATION WITHIN THE NEW SUBSPACE
basis_vecs = STCbvs;
nSTCbvs = size(basis_vecs,2);
nmods = 10;
mod_signs = ones(nmods,1);
dim_signs = ones(nSTCbvs,1);
unused_stcs = (nmods+1):nSTCbvs;

cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
flen = 14; sdim = 24;
%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 100;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 1;
defmod.locLambda = 200;
defmod.lambda_dX = 20; %350
defmod.lambda_L1x = 10; %40
defmod.lambda_dT = 10;
defmod.pids = 1:sdim;
defmod.SDIM = sdim;
defmod.fsdim = sdim;

%define NL initializations: "lin, threshlin, pquad, nquad"
clear init_nls nltypes
for i = 1:nmods; init_nls{i} = 'threshlin'; end;
%define NL types: "uncon, lin, threshlin, quad"
for i = 1:nmods; nltypes{i} = 'threshlin'; end;

[klen,Nstcbf] = size(basis_vecs);
flen = klen/sdim;
kern_output = cur_tr_stimemb*basis_vecs;

%determine distribution of random interpoint distances
rand_reps = 500;
init_vals = zeros(rand_reps,nSTCbvs,nmods);
for r = 1:rand_reps
    % compute average separation between NN initial points
    STCcf_0 = randn(nSTCbvs,nmods);
    %normalize
    for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    init_vals(r,:,:) = STCcf_0;
end
cur_dists = zeros(1,rand_reps*(rand_reps-1)/2);
for nn = 1:nmods
    cur_dists = cur_dists + pdist(init_vals(:,:,nn),'cosine');
end
rdmean = 0.75*mean(cur_dists);
rdscale = 2*std(cur_dists);
mod_signs = ones(nmods,1);

max_reps = 300;
min_reps = 10;
min_LM_fract = 0.95;
eps = 0.002;
n_its = 25;
cur_reps = 0;
used_cfs = [];
LL_vals = [];
xvLL_vals = [];
LP_vals = [];
smallest_dists = [];
is_optimized = [];
rotbv_mod = [];
for r = 1:n_its
    fprintf('ITERATION: %d\n\n\n\n',r);
    % random initialization
    STCcf_0 = randn(nSTCbvs,nmods);
    %normalize
    for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    
    %initialize model
    glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
    glm_stcb.image_type = '1d';
    glm_stcb.basis = 'pix';
    [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);
    glm_stcb.lambdaW = 0; %sparseness on model weights
    
    %determine LL and LP at current filter point
    glm_stcb = fitWeights_stcb_nonlpsc(glm_stcb,kern_output*STCcf_0,surr_tr_spikebins,1,1e-3);
    [ll0, ll0p] = getLLGLM_STCBF_nonlpsc(glm_stcb,kern_output,cur_tr_spkbns,'none');
%     xvll = getLLGLM_STCBF(glm_stcb,xvkern_output,cur_xv_spkbns,'none');
    
    %store starting point
    used_cfs = cat(3,used_cfs,STCcf_0);
    LL_vals = [LL_vals; ll0]; LP_vals = [LP_vals; ll0p];
%     xvLL_vals = [xvLL_vals; xvll];
    is_optimized = [is_optimized; 0];
    
    
    %compute distance in filter space to nearest better value
    better_vals = find(LP_vals < ll0p);
    n_better_vals = length(better_vals);
    if n_better_vals > 0
        fprintf('%d better values found\n',n_better_vals);
        
        %compute pairwise distances to better points (along the filter
        %n-spheres) as sum of cosine distances
        cur_dists = zeros(n_better_vals,1);
        for n = 1:nmods
            cur_dists = cur_dists + (1-squeeze(sum(repmat(STCcf_0(:,n),[1 1 n_better_vals]).*used_cfs(:,n,better_vals))));
        end
        
        %determine acceptance probaiblity based on smallest distance
        fprintf('Smallest dist %.4f\n',min(cur_dists));
        smallest_dists = [smallest_dists; min(cur_dists)];
        %         acc_prob = 2./(1+exp(-smallest_dists(end)/dscale))-1;
        acc_prob = normcdf(smallest_dists(end),rdmean,rdscale);
    else
        acc_prob = 1;
    end
    
    
    %use this starting point as a seed with probability acc_prob
    fprintf('Will start with probability %.5f\n',acc_prob);
    if rand < acc_prob
        disp('Starting iteration');
        %local optimization
        %             rotbv_mod{nn} = [rotbv_mod{nn}; fitNLHI_stcb_nonlpsc(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'none',6,2)];
        rotbv_mod = [rotbv_mod; fitNLHI_stcb_nopsc(glm_stcb,cur_tr_stimemb,surr_tr_spikebins,'none',6)];
        %         xvLL = getLLGLM_STCBF(rotbv_mod(end),xvkern_output,cur_xv_spkbns,'none');
        
        [spatial_profiles, temporal_profiles, weights, mod_type, space_COM, temp_COM] = ...
            compute_mod_stats(rotbv_mod(end));
        pos_inds = find(weights > 0); neg_inds = find(weights < 0);
        [~,space_com_ord_pos] = sort(space_COM(pos_inds));
        [~,space_com_ord_neg] = sort(space_COM(neg_inds));
        used_ord = [pos_inds(space_com_ord_pos); neg_inds(space_com_ord_neg)];
        rotbv_mod(end).mods = rotbv_mod(end).mods(used_ord);
        
        %store final point
        STCcf_fin = get_STCcf_mat(rotbv_mod(end));
        %normalize filter weights to unit length for distance comparisons
        for i = 1:nmods
            STCcf_fin(:,i) = STCcf_fin(:,i)/norm(STCcf_fin(:,i));
            dist_trav(r,i) = 1-dot(STCcf_fin(:,i),STCcf_0(:,i));
        end
        used_cfs = cat(3,used_cfs,STCcf_fin);
        LL_vals = [LL_vals; rotbv_mod(end).LL]; LP_vals = [LP_vals; rotbv_mod(end).LP];
        %         xvLL_vals = [xvLL_vals; xvLL];
        is_optimized = [is_optimized; 1];
        
        cur_opt_LP_vals = LP_vals(is_optimized==1);
        cur_cluster_assignments = get_eps_ball_clusters(cur_opt_LP_vals,eps);
        w = length(unique(cur_cluster_assignments));
        cur_est_fract = (r-w-1)*(r+w)/(r*(r-1));
        fprintf('Estimating: %d local minima found\n',w);
        if cur_est_fract > min_LM_fract & r > min_reps
            disp('Stopping criteria satisfied.  Halting.');
            break
        end
    end
end

% now refine best rotated model
opt_set = find(is_optimized==1)'
% [~,best_mod] = min(xvLL_vals(is_optimized==1));
[~,best_mod] = min(LL_vals(is_optimized==1));
% xv_min = xvLL_vals(opt_set(best_mod));
ll_min = LL_vals(opt_set(best_mod));
lp_min = LP_vals(opt_set(best_mod));
STCcf_0 = nan(nSTCbvs,nmods);
for i = 1:nmods
    STCcf_0(:,i) = rotbv_mod(best_mod).mods(i).STCcf;
end
basis = 'pix';
flen = 14;
%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 200;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 100;
defmod.locSigma = 3;
defmod.maxLocPen = 250;
defmod.lambda_dX = 100; %350
defmod.lambda_L1x = 0; %40
defmod.lambda_dT = 10;
defmod.pids = 1:sdim;
defmod.SDIM = sdim;
defmod.fsdim = sdim;

for i = 1:nmods; init_nls{i} = 'threshlin'; end;
%define NL types: "uncon, lin, threshlin, quad"
for i = 1:nmods; nltypes{i} = 'uncon'; end;
glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize

%copy over internal NLs from previous fits
for i = 1:nmods
    glm_stcb.mods(i).nly = rotbv_mod(best_mod).mods(i).nly;
    glm_stcb.mods(i).w = rotbv_mod(best_mod).mods(i).w;
end
glm_stcb.const = rotbv_mod(best_mod).const;

[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,cur_tr_stimemb);
glm_stcb.image_type = '1d';
fit_fin_glm = fitNLHI2d_fullbf(glm_stcb,cur_tr_stimemb,surr_tr_spikebins,'tots',2,2);
% fit_fin_xvLL = getLLGLM_FULL2d(fit_fin_glm,cur_xv_stimemb,cur_xv_spkbns,'none');
% true_xvLL = getLLGLM_FULL2d(fin_glm,cur_xv_stimemb,cur_xv_spkbns, 'none')

%% NOW FIT "STC" model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 200;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.locSigma = 0;
defmod.maxLocPen = 0;
defmod.lambda_dX = 0; %350
defmod.lambda_L1x = 0; %40
defmod.lambda_dT = 0;
STCcf_0 = eye(nSTCbvs);
mod_signs(8:end) = -1;
for i = 1:nSTCbvs; init_nls{i} = 'pquad'; end;
%define NL types: "uncon, lin, threshlin, quad"
for i = 1:nSTCbvs; nltypes{i} = 'uncon'; end;
glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize

glm_stcb.image_type = '1d';
glm_stcb.basis = 'pix';
[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,cur_tr_stimemb);
glm_stcb.lambdaW = 0; %sparseness on model weights
[glm_stcb,COMs] = get_filter_coms_1d(glm_stcb);

%determine LL and LP at current filter point
stc_mod = fitNLw_alt_full(glm_stcb,cur_tr_stimemb,surr_tr_spikebins);
[stc_LL, stc_LP] = getLLGLM_FULL2d(stc_mod,cur_tr_stimemb,surr_tr_spikebins,'tots');

% stc_xvLL = getLLGLM_FULL2d(stc_mod,cur_xv_stimemb,cur_xv_spkbns,'none');


%% generate XV surrogate data
NT = NXV;
for rr = 1:50
    rr
    cur_rand_xv_stim = round(rand(NT,sdim)); cur_rand_xv_stim(cur_rand_xv_stim==0)=-1;
    cur_xv_stimemb = makeStimRows(cur_rand_xv_stim,flen);
    
    %internal filter outputs of each module
    xv_gen_functions = cur_xv_stimemb*pix_mat;
    
    kx    = fin_glm.const; %initialize kx with model constant term
    nmods = length(fin_glm.mods);
    
    for imod = 1:nmods %loop across NL modules
        
        fg = xv_gen_functions(:,imod);
        fg = nlin_proc_stim(fg, fin_glm.mods(imod).nly, fin_glm.mods(imod).nlx );
        kx = kx + fg*fin_glm.mods(imod).w;
        
    end
    r = log(1+exp(kx)); %apply spiking NL tp get predicted rate
    
    spikes = poissrnd(r);
    rbins = (find(spikes>0.5));
    nsp = spikes(rbins);
    spk_vals = unique(spikes); spk_vals(spk_vals==0) = [];
    surr_xv_spikebins = [];
    for i = 1:length(spk_vals)
        cur_set = find(spikes == spk_vals(i));
        surr_xv_spikebins = [surr_xv_spikebins; repmat(cur_set(:),spk_vals(i),1)];
    end
    fprintf('Nspks: %d\n',length(surr_xv_spikebins));
    
    fit_fin_xvLL(rr) = getLLGLM_FULL2d(fit_fin_glm,cur_xv_stimemb,surr_xv_spikebins,'none');
    true_xvLL(rr) = getLLGLM_FULL2d(fin_glm,cur_xv_stimemb,surr_xv_spikebins,'none');
    stc_xvLL(rr) = getLLGLM_FULL2d(stc_mod,cur_xv_stimemb,surr_xv_spikebins,'none');   
end

cur_xv_stim = stim(xv_inds{xv},:);
cur_xv_stimemb = makeStimRows(cur_xv_stim,flen);
cur_xv_spkbns = xv_spkbns{xv};
real_fit_fin_xvLL = getLLGLM_FULL2d(fit_fin_glm,cur_xv_stimemb,cur_xv_spkbns,'none');
real_true_xvLL = getLLGLM_FULL2d(fin_glm,cur_xv_stimemb,cur_xv_spkbns,'none');
real_stc_xvLL = getLLGLM_FULL2d(stc_mod,cur_xv_stimemb,cur_xv_spkbns,'none');
