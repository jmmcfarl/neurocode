%% Load Data

clear all;
addpath('~/James_scripts/GLM')
addpath('~/James_scripts/GLM/2d')
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/Timm/tim1/functions/')


cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

% load spks7576-ns.mat; load PCScores-ns.mat; stype='ns'
% load spks7576-pn.mat; load PCScores-pn.mat; stype='pn'
% load spks7576-ps.mat; load PCScores-ps.mat; stype='ps'

cd ~/Data/blanche/matlabdata/
load spks7576-all.mat; 
load PCScores_z-all.mat; 
stype='all';


%% Compute whitened data, and STA/STCs

tcell = 26;
pids =1:1024;
tsbs      = 1+floor(aselspks{tcell}/dt);

% figure; ecdf(aselspks{tcell})

ncomps    = 1000; compids   = 1:ncomps;
pix_conv_mat = diag(sqrt(scorevars(compids)))*coefs(:,compids)';
kern_conv_mat = coefs(:,compids)*diag(1./sqrt(scorevars(compids)));

WX        = scores(:,compids)*diag(1./sqrt(scorevars(compids)));
% mWX = mean(WX*pix_conv_mat)*kern_conv_mat;
% WX = WX - repmat(mWX,size(WX,1),1);
spikebins = tsbs(tsbs>flen & tsbs<(size(WX,1)+flen-1))-flen+1 ;
WS        = WX(spikebins,:);
rsta      = mean(WS) - mean(WX);

stvcv = cov(WS);  utvcv = cov(WX);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);

npos=10; nneg=6;
stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
nposdims = 4; nnegdims = 0;
posdims = 1:nposdims; negdims = 1:nnegdims;
% STCbvs = [sta' stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
kimages = [rsta',STCbvs]'*diag(sqrt(scorevars(compids)))*coefs(:,compids)';
f1 = figure('Name',stype); plotfilterbank(kimages',32,pids)

%add in STA
STCbvs = [rsta' STCbvs];

[klen,Nstcbvs] = size(STCbvs);
flen = 6;
stimlen = size(WX,1);
kern_output = WX*STCbvs;

%% INITIALIZE MODEL
addpath('~/James_scripts/GLM/2d')
nmods = Nstcbvs;
mod_signs = ones(nmods,1);
dim_signs = ones(Nstcbvs,1);

unused_stcs = (nmods+1):Nstcbvs;
%initialize on STC dims
STCcf_0 = eye(Nstcbvs);
STCcf_0(:,unused_stcs) = [];
if nmods > Nstcbvs
    n_extra = nmods-Nstcbvs;
    STCcf_0 = [STCcf_0 randn(Nstcbvs,n_extra)];
end
%random initialization
% n_extra = nmods;
% STCcf_0 = randn(Nstcbvs,n_extra);

%make sure expansive subunits start out in expansive subspace, etc
STCcf_0(negdims,mod_signs==1) = 0;
STCcf_0(posdims,mod_signs==-1) = 0;

% %double space
% STCcf_0 = [STCcf_0 -STCcf_0];
% mod_signs = [mod_signs; mod_signs];
% nmods = 2*nmods;

%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;

%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 200;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 50;
defmod.lambda_L1x = 3e-3;
defmod.lambda_dT = 0;
defmod.pids = pids;
nltype = 'uncon';
init_nls{1} = 'l'; for i = 2:nmods; init_nls{i} = 'pq'; end;
% basis = 'pix'; 
basis = 'white';
glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltype,init_nls,basis,'test'); %initialize
[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
% glm_stcb.basis = 'white';
glm_stcb.lambdaW = 0; %sparseness on model weights

figure
plot2d_mod(glm_stcb)
% stc_pix = STCbvs'*pix_conv_mat;
% figure
% plotfilterbank(stc_pix',32,pids);

%% Fit model
% clear WX scores scorevars coefs
% load StimXmat_z-all
% full_glm = fitNLHI2d_fullbf(glm_stcb,X,spikebins,'tots');

%%
full_glm = fitNLHI2d_fullbf(glm_stcb,WX,spikebins,'tots');
figure
plot2d_mod(full_glm);
figure
plot2d_mod(glm_stcb);

%%
skern_len = length(glm_stcb.mods(1).k);
pkern_len = size(pix_conv_mat,2);
kern_t = pkern_len/fsdim;

lamrange = [];
stimlen = size(WX,1);
%compute binned spike vector
Robs = zeros(1,stimlen);
ftable = tabulate(spikebins);
Robs(ftable(:,1)) = ftable(:,2);

tempmod = glm_stcb;
tempmod.kderiv_pen = [1e3 0];
tempmod.k2_pen = 0;
targets = 1:nmods;
initial_params = [];
for m = 1:length(targets)
    cur_pkern = tempmod.mods(targets(m)).k'*tempmod.pix_conv_mat;
    initial_params = [initial_params cur_pkern]; %add STCBcoefs to initial param vector
end
initial_params(end+1) = tempmod.const; %add constant offset term to params
tempmod.kern_conv_mat = kern_conv_mat;

[r,gen_funs] = get_predicted_rate_stcbf(tempmod,WX*STCbvs);
NL_gen_funs = gen_funs;
NL_gen_funs(NL_gen_funs < 0) = 0;

addpath('~/James_scripts/minFunc/')
options = [];
% options.display = 'none';
options.maxFunEvals = 30;
% options.Method = 'newton0';
% options.Method = 'csd';
options.Method = 'lbfgs';
[params LL] = minFunc( @(K) FULLBF2d_LLinternal_nonlpsc(K,Robs,WX,tempmod, NL_gen_funs, lamrange,targets), initial_params', options );

% Reassign variables
fit1 = tempmod;
fit1.const = params(end);
for n = 1:length(targets)
    cur_kern = params((n-1)*pkern_len+(1:pkern_len));
    kern_out = WX * (cur_kern'*kern_conv_mat)';
    norm_kern = cur_kern/std(kern_out);
    fit1.mods(targets(n)).pix_k = norm_kern;
    fit1.mods(targets(n)).k = (norm_kern'*kern_conv_mat);
end
K_MAT_ref = get_k_mat(fit1);
stcbf_pix = STCbvs'*pix_conv_mat;
ref_pix = K_MAT_ref'*pix_conv_mat;
figure
plotfilterbank(stcbf_pix',32,pids);
figure
plotfilterbank(ref_pix',32,pids);



