clear all;
close all;

addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath(genpath('~/James_Scripts'))

pars = load('~/Data/rust/infos/StandardParametersRust.mat');
defmod = pars.defmod;
defmod.h(1:end-1) = [];%restrict PSC term to delta function
flen = pars.flen;
hilen = length(defmod.h);
foff = flen + pars.hilen;

datdir = '~/Data/rust/stcbar/Data/';
cfiles = dir([datdir,'*stc.mat']); ncells = length(cfiles);
fnames = arrayfun(@(x)x.name,cfiles,'UniformOutput',0);
uids   = cellfun(@(x)[x(2:3),'-',x(6:7)],fnames,'UniformOutput',0);
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

cd /Users/James/James_scripts/stc_sparse_test_figs/stac_allcells
load ./used_stcims.mat


ncells = length(uids);
% for nn = 1:ncells
nn = 1;
cd /Users/James/James_scripts/stc_sparse_test_figs/revised_modfits2/

fprintf('ANALYZING CELL %d OF %d\n\n',nn,ncells);

%% load data
eval(['load ',['~/Data/rust/stcbar/Data/',fnames{nn}]]);

psth = spikes_per_frm(:);
rbins = (find(psth>0.5));
nsp = psth(rbins);
spikebins =[];
uvals = unique(psth);
for i = 1:length(uvals)
    cur_set = find(psth==uvals(i));
    spikebins = [spikebins; repmat(cur_set(:),uvals(i),1)];
end
spikebins = sort(spikebins);

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
[NT,sdim] = size(stim);

%%
xv = 1; %use only one XV set

cur_xv_stim = stim(xv_inds{xv},:);
cur_xv_stimemb = makeStimRows(cur_xv_stim,flen);
cur_xv_spkbns = xv_spkbns{xv};

cur_tr_stim = stim(tr_inds{xv},:);
cur_tr_stimemb = makeStimRows(cur_tr_stim,flen);
cur_tr_spkbns = tr_spkbns{xv};

%compute STA
sta = mean(cur_tr_stimemb(cur_tr_spkbns,:)) - mean(cur_tr_stimemb);

%project out STA
proj_mat = sta'*inv(sta*sta')*sta;
stim_proj = cur_tr_stimemb - cur_tr_stimemb*proj_mat;
stvcv = cov(stim_proj(cur_tr_spkbns,:));
utvcv = cov(stim_proj);
[evecs,evals] = eig(stvcv-utvcv);
evs{nn} = diag(evals);
STCbvs  = evecs;
sta = sta';
nSTCbvs = size(STCbvs,2);

npos = 8; nneg = 8;
stcs_compareset  = evecs(:,[1:nneg,length(evs{nn})-npos+1:end]);
stcs_compareset  = stcs_compareset(:,end:-1:1);
rstcs = fliplr(stcs_compareset); %reversed STC kernels (suppressive first)
stcs_compareset = [stcs_compareset(:,1:npos) rstcs(:,1:nneg)];

f1 = figure;
nfilts = 8;
subplot(3,nfilts,1)
imagesc(reshape(sta,flen,sdim));
for i = 1:nfilts
    subplot(3,nfilts,nfilts+i)
    imagesc(reshape(stcs_compareset(:,i),flen,sdim));
end
for i = 1:nfilts
    subplot(3,nfilts,2*nfilts+i)
    imagesc(reshape(stcs_compareset(:,2*nfilts+1-i),flen,sdim));
end
colormap(gray)

%%
%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 0;
defmod.lambda_L1x = 0; %2
defmod.lambda_dT = 0;
defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
basis = 'pix';

nmods = 5;
init_signs = [-1*ones(1,5)];
% init_signs = [ones(1,4) -1*ones(1,6)];
init_betas = 2*ones(1,nmods);

% init_kerns = randn(size(stcs_compareset,1),nmods);
% for i = 1:nmods
%     init_kerns(:,i) = init_kerns(:,i)/norm(init_kerns(:,i));
% end

cur_basis = stcs_compareset;
STCcf_0 = randn(size(cur_basis,2),nmods);
%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
init_kerns = cur_basis*STCcf_0;

glm_stcb = createGLM_lexp(init_kerns,init_signs,init_betas,defmod,basis);
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
glm_lexp2 = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots');
% glm_lexp2 = normalizeRFs_full(glm_lexp2,cur_tr_stimemb);
% 
% k_mat = get_k_mat(glm_lexp2);
% g_mat = cur_tr_stimemb*k_mat;
% glm_lexp2 = fitWeights_lexp(glm_lexp2,g_mat,cur_tr_spkbns,1);

% glm_lexp4 = fitNLw_alt_lexp(glm_lexp4,cur_tr_stimemb,cur_tr_spkbns,0);
% glm_lexp2 = fitGLM_lexp(glm_lexp6,cur_tr_stimemb,cur_tr_spkbns,'tots');

% glm_lexp6.lambdaW = 1000;
% glm_lexp6 = fitNLw_alt_lexp(glm_lexp6,cur_tr_stimemb,cur_tr_spkbns,0);
