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
load spks7576-all.mat; load PCScores-all.mat; stype='all'


%%
tcell = 12;
pids =1:1024;
tsbs      = 1+floor(aselspks{tcell}/dt);

% figure; ecdf(aselspks{tcell})

ncomps    = 600; compids   = 1:ncomps;
WX        = scores(:,compids)*diag(1./sqrt(scorevars(compids)));
spikebins = tsbs(tsbs>flen & tsbs<(size(WX,1)+flen-1))-flen+1 ;
WS        = WX(spikebins,:);
rsta      = mean(WS) - mean(WX);

stvcv = cov(WS);  utvcv = cov(WX);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);

%%
npos=50; nneg=6;
stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
nposdims = 10; nnegdims = 0;
posdims = 1:nposdims; negdims = 1:nnegdims;
% STCbvs = [sta' stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
% kimages = [rsta',STCbvs]'*diag(sqrt(scorevars(compids)))*coefs(:,compids)';
f1 = figure('Name',stype); plotfilterbank(kimages',32,pids)

%add in STA
STCbvs = [rsta' STCbvs];
nposdims = 11;

[klen,Nstcbvs] = size(STCbvs);
flen = 6;
stimlen = size(WX,1);
kern_output = WX*STCbvs;

pix_conv_mat = diag(sqrt(scorevars(compids)))*coefs(:,compids)';
kern_conv_mat = coefs(:,compids)*diag(1./sqrt(scorevars(compids)));

%%
addpath('~/James_scripts/GLM/2d')
nmods = 4;
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

%double space
STCcf_0 = [STCcf_0 -STCcf_0];
mod_signs = [mod_signs; mod_signs];
nmods = 2*nmods;

%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;

%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.locLambda = 30;
defmod.pids = pids;
glm_stcb = createGLM2d_stcb(STCbvs,STCcf_0,pix_conv_mat,defmod,mod_signs,dim_signs,'test'); %initialize
[glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);

glm_stcb.lambdaW = 0; %sparseness on model weights

%fit model
stcbf_mod = fitNLHI2d_stcb_nonlpsc(glm_stcb,WX,spikebins,'tots');

stcbf_cfs = get_STCcf_mat(stcbf_mod);
stcbf_pix = (stcbf_mod.STCbasis*stcbf_cfs)'*pix_conv_mat;
stc_pix = STCbvs'*pix_conv_mat;
figure
plotfilterbank(stc_pix',32,pids);
figure
plotfilterbank(stcbf_pix',32,pids);

%%
stimlen = size(WX,1);
Robs = zeros(1,stimlen);
ftable = tabulate(spikebins);
Robs(ftable(:,1)) = ftable(:,2);

[r,gen_funs] = get_predicted_rate_stcbf(stcbf_mod,WX*STCbvs);
NL_gen_funs = gen_funs;
NL_gen_funs(NL_gen_funs < 0) = 0;

C = corr(NL_gen_funs) - eye(nmods);
figure
imagesc(C);colorbar

DKL = nan(nmods,1);
for i = 1:nmods
    ov_gen = gen_funs(:,i);
    gen_spike = ov_gen(spikebins);
    ov_gen_mean = mean(ov_gen);
    ov_gen_var = var(ov_gen);
    gen_spike_mean = mean(gen_spike);
    gen_spike_var = var(gen_spike);
    
    DKL(i) = 1/2*(gen_spike_var/ov_gen_var + (ov_gen_mean - gen_spike_mean).^2/ov_gen_var - log(gen_spike_var/ov_gen_var) - 1);
end

n_rows = 4;
n_cols = ceil(nmods/n_rows);
figure
row_cnt = 1;
col_cnt = 1;
ov_max = max(stcbf_cfs(:));
for i = 1:nmods
   subplot(n_rows,n_cols,(row_cnt-1)*n_cols + col_cnt)
%    polar(linspace(0,360,Nstcbvs),abs(stcbf_cfs(:,i))','.-')
   plot(1:Nstcbvs,abs(stcbf_cfs(:,i)),'o-'), ylim([0 ov_max])
   title(sprintf('Weight %.5f',stcbf_mod.mods(i).w));
   row_cnt = row_cnt + 1;
   if row_cnt > n_rows
       row_cnt = 1;
       col_cnt = col_cnt + 1;
   end
end
  
figure
plot(mean(abs(stcbf_cfs')),'o-')
xlabel('STC BF','fontsize',14)
ylabel('Average weighting','fontsize',14)

w = arrayfun(@(x) x.w,stcbf_mod.mods);
figure
plot(DKL,w,'.')
xlabel('Module Information','fontsize',14)
ylabel('Module weight','fontsize',14)


%%
cur_basevecs = zeros(Nstcbvs,2);
cur_basevecs(1,1) = 1; cur_basevecs(2,2) = 1;
norm_cfs = normvecsl2(stcbf_cfs);
pcofs = norm_cfs'*cur_basevecs;
figure
plot(pcofs(:,1),pcofs(:,2),'o')
xlim([-1 1]), ylim([-1 1])
hold on
circle([0 0],1,100,'k--')
xlabel('Projection onto STC1','fontsize',14)
ylabel('Projection onto STC2','fontsize',14)
%%

cur_nmods = 8;
glmod = stcbf_mod;

skern_len = length(glmod.mods(1).k);
pkern_len = size(pix_conv_mat,2);
kern_t = pkern_len/fsdim;

space_sigma = 3;
time_sigma =0.1;
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
stimlen = size(WX,1);
%compute binned spike vector
Robs = zeros(1,stimlen);
ftable = tabulate(spikebins);
Robs(ftable(:,1)) = ftable(:,2);

tempmod = glmod;
tempmod.kderiv_pen = [0 0];
tempmod.k2_pen = 1e3;
targets = 1:4;
initial_params = [];
for m = 1:length(targets)
    cur_pkern = tempmod.mods(targets(m)).k'*tempmod.pix_conv_mat;
    initial_params = [initial_params cur_pkern]; %add STCBcoefs to initial param vector
end
initial_params(end+1) = tempmod.const; %add constant offset term to params
s
tempmod.kern_conv_mat = kern_conv_mat;

% X = WX*pix_conv_mat;
[r,gen_funs] = get_predicted_rate_stcbf(tempmod,WX*STCbvs);
NL_gen_funs = gen_funs;
NL_gen_funs(NL_gen_funs < 0) = 0;

% [X,grad_X] = FULLBF2d_LLinternal_nonlpsc(initial_params',Robs,WX,tempmod, NL_gen_funs, lamrange,targets);
% [params LL eflag] = fminunc( @(K) FULLBF2d_LLinternal_nonlpsc(K,Robs,WX,tempmod, NL_gen_funs, lamrange,targets), initial_params, opts );


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
K_MAT0 = get_k_mat(glm_stcb);
K_MAT = get_k_mat(stcbf_mod);
K_MAT_ref = get_k_mat(fit1);
init_pix = K_MAT0'*pix_conv_mat;
stcbf_pix = STCbvs'*pix_conv_mat;
pix = K_MAT'*pix_conv_mat;
ref_pix = K_MAT_ref'*pix_conv_mat;
figure
plotfilterbank(stcbf_pix',32,pids);
% figure
% plotfilterbank(init_pix',32,pids);
figure
plotfilterbank(pix',32,pids);
figure
plotfilterbank(ref_pix',32,pids);



