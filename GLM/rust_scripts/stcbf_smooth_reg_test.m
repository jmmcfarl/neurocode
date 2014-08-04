clear all;

addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/James_scripts/GLM')

pars = load('~/Data/rust/infos/StandardParametersRust.mat');
defmod = pars.defmod;
defmod.h(1:end-1) = [];%restrict PSC term to delta function
flen = pars.flen;
hilen = length(defmod.h);
foff = flen + pars.hilen;
ooptions = optimset('MaxFunEvals',100000,'MaxIter',10000);

% datdir = '/Users/timm/myrepo/workspace/DataRepository/rust/stcbar/DATA/';
datdir = '~/Data/rust/stcbar/Data/';
cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles);
fnames = arrayfun(@(x)x.name,cfiles,'UniformOutput',0);
uids   = cellfun(@(x)[x(2:3),'-',x(6:7)],fnames,'UniformOutput',0);

% tuid = '52-18'; hassta=1;  npos=1; nneg=1;
tuid = '44-29'; hassta=0;  npos=12; nneg=12;
% tuid = '43-03'; hassta=0;  npos=6; nneg=6;
% tuid = '43-21'; hassta=0;  npos=6; nneg=6; %was 6,6
% tuid = '33-27';  %% complicated simple
% tuid = '43-03'; %%??
% tuid = '52-18'; %% simple scell
% tuid = '44-29'; %% the complex cell
% tuid = '33-44';  %% ??

%% load STC data
cd ~/Data/rust
dataname = sprintf('stcbf_data-%s.mat',tuid)
eval(['load ' dataname]);
tid =  find(strcmp(tuid,uids));
npos = 6; nneg = 6;

rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
nposdims = npos; nnegdims = nneg;
posdims = 1:nposdims; negdims = 1:nnegdims;
% STCbvs = [sta' stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
negdims = negdims + length(posdims);

%% create XV data
[stimlen,sdim] = size(stim);
nparts = 18;
partlen = floor(stimlen/nparts);
nfold = 3;
nxvparts = nparts/nfold;

%boundaries of parts
pbounds  = [(0:nparts-1)*partlen+1;(1:nparts)*partlen]';

for i = 1:nfold
    xv_inds{i} = [];
    xv_spkbns{i} = [];
    tr_inds{i} = [];
    tr_spkbns{i} = [];
    
    cur_perm = randperm(nparts);
    cur_xv_parts = sort(cur_perm(1:nxvparts));
    cur_tr_parts = setdiff(1:nparts,cur_xv_parts);
    
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

% for xv = 1:nfold
xv = 1;

cur_tr_stim = stim(tr_inds{xv},:);
cur_tr_spkbns = tr_spkbns{xv};

cur_xv_stim = stim(xv_inds{xv},:);
cur_xv_spkbns = xv_spkbns{xv};

%% precompute the stimulus filtered by each STC kernel for training
%% stim
[klen,Nstcbf] = size(STCbvs);
flen = klen/sdim;
stimlen = size(cur_tr_stim,1);
Nstcbf = size(STCbvs,2);
kern_output = zeros(stimlen-flen+1,Nstcbf); %initialize filtered output to be (NT-N_tau+1)xNmods.
for ikern = 1:Nstcbf;
    %for the given NL module, store the output of the stimulus filtered by
    %the internal kernel
    kern_output(:,ikern)= kfilterInput(cur_tr_stim,STCbvs(:,ikern));
end

%% for XV data
[xvstimlen,sdim] = size(cur_xv_stim);

%precompute the stimulus filtered by each STC kernel
[klen,cur_Nstcbf] = size(STCbvs);
flen = klen/sdim;
xvkern_output = zeros(xvstimlen-flen+1,1); %initialize filtered output to be (NT-N_tau+1)xNmods.
for ikern = 1:Nstcbf;
    %for the given NL module, store the output of the stimulus filtered by
    %the internal kernel
    xvkern_output(:,ikern)= kfilterInput(cur_xv_stim,STCbvs(:,ikern));
end


%% Initialize model
defmod.SDIM=sdim;
defmod.locLambda = 50;
defmod.lh = 0;
defmod.lh2 = 0;
defmod.lnl = 0;
defmod.lnl2 = 0;
defmod.hcon = 0;defmod.hmon = 0;defmod.nlcon = 0;defmod.nlmon = 0;
nSTCbvs = size(STCbvs,2);

foff = flen + length(defmod.h);
tr_cspkbs = cur_tr_spkbns(cur_tr_spkbns>foff & cur_tr_spkbns<stimlen)-foff+2;

dim_signs = nan(size(STCbvs,2),1);
dim_signs(posdims) = 1;
dim_signs(negdims) = -1;


%%
%     for jj = 1:length(nmods)
cur_nmods = 2;
%         cur_nmods = nmods(jj);
mod_signs = nan(cur_nmods,1);

posmods = 1:(cur_nmods/2); negmods = (cur_nmods/2+1):cur_nmods;
posmod_inds = 1:length(posmods);
if cur_nmods <= Nstcbf
    negmod_inds = (length(posdims)+1):(length(posdims)+length(negmods));
    unused_stcs = setdiff(1:nSTCbvs,[posmod_inds negmod_inds]);
    mod_signs(posmods) = 1;
    mod_signs(negmods) = -1;
else
    extra_mods = (Nstcbf+1):cur_nmods;
    mod_signs(posdims) = 1;
    mod_signs(negdims) = -1;
    mod_signs(extra_mods(mod(extra_mods,2)==0)) = 1;
    mod_signs(extra_mods(mod(extra_mods,2)==1)) = -1;
    unused_stcs = [];
end

%initialize on STC dims
STCcf_0 = eye(nSTCbvs);
STCcf_0(:,unused_stcs) = [];
if cur_nmods > nSTCbvs
    n_extra = cur_nmods-nSTCbvs;
    STCcf_0 = [STCcf_0 randn(nSTCbvs,n_extra)];
end
%random initialization
%         n_extra = cur_nmods;
%         STCcf_0 = randn(nSTCbvs,n_extra);

%make sure expansive subunits start out in expansive subspace, etc
STCcf_0(negdims,mod_signs==1) = 0;
STCcf_0(posdims,mod_signs==-1) = 0;

%double space
STCcf_0 = [STCcf_0 -STCcf_0];
mod_signs = [mod_signs; mod_signs];

%normalize
for i = 1:2*cur_nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;

%initialize model
glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
[glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);

glm_stcb.lambdaW = 0; %sparseness on model weights

%fit model
stc_posneg_mod{xv} = fitNLHI_stcb_nonlpsc(glm_stcb,cur_tr_stim,tr_cspkbs,'tots');

%%
cur_nmods = 4;
glmod = stc_posneg_mod{xv};
%compute prior covariance matrix on ks
kern_len = length(glmod.mods(1).k);
kern_t = kern_len/sdim;

space_sigma = 3;
time_sigma =0.1;
prior_scale = -5;
[Xspace,Xtime] = meshgrid(1:sdim,1:kern_t);
Xspace_vec = Xspace(:);
Xtime_vec = Xtime(:);
space_distmat = squareform(pdist(Xspace_vec));
time_distmat = squareform(pdist(Xtime_vec));
dist_metric = space_distmat.^2/space_sigma^2 + time_distmat.^2/time_sigma^2;
prior_cov = exp(-1/2*space_distmat.^2/space_sigma^2 - 1/2*time_distmat.^2/time_sigma^2 - prior_scale);
% [V,D] = eig(prior_cov);
% Y = V*cur_kern;
% eigvals = diag(D);
% inv_eigvals = 1./eigvals;
% inv_eigvals(eigvals < 1e-20) = 0;
% Dinv = diag(inv_eigvals);
% glmod.V = V;
% glmod.Dinv = Dinv;
% [U,S,V] = svd(prior_cov);
% temp = 1./diag(S);
% temp(250:end) = 0;
% Sp = diag(temp);
% % prior_precisionmat = V*Sp*U';
prior_precisionmat = inv(prior_cov);
glmod.prior_precisionmat = prior_precisionmat;
% glmod.prior_cov = prior_cov;
% [glmod.R,p] = cholcov(prior_cov);
% glmod.dist_metric = dist_metric;
glmod.dist_pen = 0;

%
tolF = 1e-4;
tolX = 1e-6;
cs = [];

silent = 0;
if silent == 0
    opts = optimset('Algorithm','active-set','GradObj','on','LargeScale','off','Display','iter','MaxIter',20,'TolFun',tolF,'TolX',tolX);
else
    opts = optimset('Algorithm','active-set','GradObj','on','LargeScale','off','Display','off','MaxIter',20,'TolFun',tolF,'TolX',tolX);    
end

lamrange = [];
ustim = cur_tr_stim(flen:end,:);
stimlen = size(cur_tr_stim,1);
%compute binned spike vector
foff = flen + length(glmod.mods(1).h) - 2;
Robs = zeros(1,stimlen - foff);
ftable = tabulate(tr_cspkbs);
Robs(ftable(:,1)) = ftable(:,2);

X    = zeros(stimlen-flen+1,kern_len);
for i = flen:stimlen; X(i-flen+1,1:kern_len) = reshape(cur_tr_stim((i-flen+1):i,:),1,kern_len); end

%%%%%%%% THE FIT %%%%%%%%
%if you don't want exc and sup filters to mix
% [params LL eflag] = fminunc( @(K) STCBF_LLinternal_nonlpsc(K,Robs,ustim,glmod, lamrange), initial_params, opts );

noise = 0.1*randn(kern_len,1);

tempmod = glmod;
tempmod.dist_pen = 1e4;
tempmod.mods(1).k = tempmod.mods(1).k + noise;
targets = 1:4;
initial_params = [];
for m = 1:length(targets)
    initial_params = [initial_params; tempmod.mods(targets(m)).k]; %add STCBcoefs to initial param vector
end
initial_params(end+1) = tempmod.const; %add constant offset term to params

[params LL eflag] = fminunc( @(K) FULLBF_LLinternal_nonlpsc(K,Robs,X,tempmod, lamrange,targets), initial_params, opts );

% Reassign variables
fit1 = tempmod;
fit1.const = params(end);
for n = 1:length(targets)
    cur_kern = params((n-1)*kern_len+(1:kern_len));
    kern_out = X*cur_kern;
    fit1.mods(targets(n)).k = cur_kern/std(kern_out);
end


%%
xvspkbns = cur_xv_spkbns(cur_xv_spkbns>foff & cur_xv_spkbns<xvstimlen)-foff;
stc_posneg_xvLL(xv) = getLLGLM_STCBF(stc_posneg_mod{xv,jj},xvkern_output,xvspkbns,'none');



%     end

% end