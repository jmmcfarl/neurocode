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
npos = 8; nneg = 8;

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

%%
nmods = 6;
xv = 1; %first xv set

cur_tr_stim = stim(tr_inds{xv},:);
cur_tr_spkbns = tr_spkbns{xv};

foff = flen + 1;
tr_cspkbs = cur_tr_spkbns(cur_tr_spkbns>foff & cur_tr_spkbns<stimlen)-foff+2;

cur_xv_stim = stim(xv_inds{xv},:);
cur_xv_spkbns = xv_spkbns{xv};

%% create X matrix
NT   = size(cur_tr_stim,1); NeK = flen*sdim;
X    = zeros(NT-flen+1,NeK);
for i = flen:NT; X(i-flen+1,1:NeK) = reshape(cur_tr_stim((i-flen+1):i,:),1,NeK); end


%% Initialize model
Nstcbvs = size(STCbvs,2);
mod_signs = ones(nmods,1);
mod_signs(npos+1:end) = -1;
dim_signs = ones(Nstcbvs,1);
dim_signs(npos+1:end) = -1;

unused_stcs = (nmods+1):Nstcbvs;
%initialize on STC dims
% STCcf_0 = eye(Nstcbvs);
% STCcf_0(:,unused_stcs) = [];
% if nmods > Nstcbvs
%     n_extra = nmods-Nstcbvs;
%     STCcf_0 = [STCcf_0 randn(Nstcbvs,n_extra)];
% end
%random initialization
n_extra = nmods;
STCcf_0 = randn(Nstcbvs,n_extra);

%make sure expansive subunits start out in expansive subspace, etc
STCcf_0(negdims,mod_signs==1) = 0;
STCcf_0(posdims,mod_signs==-1) = 0;

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
defmod.SDIM = sdim;
defmod.fsdim = sdim;
nltype = 'uncon';
for i = 1:npos; init_nls{i} = 'pq'; end;
for i = (npos+1):nmods; init_nls{i} = 'nq'; end;
basis = 'pix';
% basis = 'white';

defmod.lambda_L1x = 100;
defmod.lambda_dX = 50;
defmod.lambda_dT = 50;
defmod.locLambda = 0;

glm_stcb = createGLM1d_fullbf(STCbvs,STCcf_0,defmod,nltype,init_nls,basis,'test'); %initialize
[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,X);
% glm_stcb.basis = 'white';
glm_stcb.lambdaW = 0; %sparseness on model weights
glm_stcb.image_type = '1d';
plotfo1d_nopsc(glm_stcb,3)
%%
addpath('~/James_scripts/GLM/2d/')
full_glm = fitNLHI2d_fullbf(glm_stcb,X,tr_cspkbs,'tots');
plotfo1d_nopsc(glm_stcb,3)
plotfo1d_nopsc(full_glm,3)
