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
nn = 27;

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
evs = diag(evals);
STCbvs  = evecs;
sta = sta';
nSTCbvs = size(STCbvs,2);

npos = 12; nneg = 12;
stcs_compareset  = evecs(:,[1:nneg,length(evs)-npos+1:end]);
stcs_compareset  = stcs_compareset(:,end:-1:1);
rstcs = fliplr(stcs_compareset); %reversed STC kernels (suppressive first)
stcs_compareset = [stcs_compareset(:,1:npos) rstcs(:,1:nneg)];

f1 = figure;
nfilts = 12;
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

%% %initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 20;
defmod.lambda_L1x = 20; %2
defmod.lambda_dT = 20;
defmod.lambda_L2x = 0;
defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
basis = 'pix';

nmods2 = [2:2:10 14:4:24];

npos = 12; nneg = 12;
cur_basis  = evecs(:,[1:nneg,length(evs)-npos+1:end]);
cur_basis  = cur_basis(:,end:-1:1);
cur_basis = [cur_basis(:,1:npos) rstcs(:,1:nneg)];

n_bvs = size(cur_basis,2);

for n = 8:length(nmods)
    fprintf('fitting model with %d subunits\n',nmods(n));
    STCcf_0 = randn(n_bvs,nmods(n));
    %normalize
    for i = 1:nmods(n); STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    init_kerns = cur_basis*STCcf_0;
    init_signs = [1*ones(1,nmods(n)/2) -1*ones(1,nmods(n)/2)];
    basis = 'pix';
    glm_stcb = createGLM_quad(init_kerns,init_signs,defmod,basis,[],[]);
    glm_stcb.image_type = '1d';
    glm_stcb.spk_nl = 'logexp';
    glm_quad(n) = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots');
    xvLL_quad(n) = getLLGLM_lexp(glm_quad(n),cur_xv_stimemb,cur_xv_spkbns,'none');
end

%%
nmods2 = [2:2:10 14:4:42];
defmod.locLambda = 0;
defmod.locSigma = 4;
defmod.maxLocPen = 1e4;
for n = 10:length(nmods2)
    n_bvs = size(cur_basis,2);
    STCcf_0 = randn(n_bvs,nmods2(n));
    %normalize
    for i = 1:nmods2(n); STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    init_kerns = cur_basis*STCcf_0;
    init_signs = [1*ones(1,nmods2(n)/2) -1*ones(1,nmods2(n)/2)];
    %     init_betas = 4*ones(1,nmods);
    basis = 'pix';
    glm_stcb = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
    glm_stcb.image_type = '1d';
    glm_stcb.spk_nl = 'logexp';
    if n > 1
        for ii = 1:nmods2(n-1)
            glm_stcb.mods(ii).k = glm_tlin(n-1).mods(ii).k;
            glm_stcb.mods(ii).pix = glm_tlin(n-1).mods(ii).pix;
        end
    end
    glm_tlin(n) = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots');
    xvLL_lexp(n) = getLLGLM_lexp(glm_tlin(n),cur_xv_stimemb,cur_xv_spkbns,'none');
end

%%
nmods = 30;

n_bvs = size(cur_basis,2);
STCcf_0 = randn(n_bvs,nmods);
%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
init_kerns = cur_basis*STCcf_0;
init_signs = [1*ones(1,nmods/2) -1*ones(1,nmods/2)];
%     init_betas = 4*ones(1,nmods);
basis = 'pix';
glm_stcb = createGLM_rquad(init_kerns,init_signs,defmod,basis,[],[]);
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
if n > 1
    for ii = 1:nmods
        glm_stcb.mods(ii).k = glm_tlin(n-1).mods(ii).k;
        glm_stcb.mods(ii).pix = glm_tlin(n-1).mods(ii).pix;
    end
end
glm_rquad = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots');


%%
n = 8;
for i = 1:length(glm_tlin(n).mods)
    glm_tlin(n).mods(i).locLambda = 500;
    glm_tlin(n).mods(i).locSigma = 3;
    glm_tlin(n).mods(i).maxLocPen = 1e5;
end
glm_tlin_loc = fitGLM_lexp(glm_tlin(n),cur_tr_stimemb,cur_tr_spkbns,'tots');
xvLL_lexp_loc = getLLGLM_lexp(glm_tlin_loc,cur_xv_stimemb,cur_xv_spkbns,'none');

%%
cd /Users/James/James_scripts/stc_sparse_test_figs/revised_modfits2
save rustcell_gnmfigtest_1 *_rot rotxvLL glm_quad glm_tlin xvLL* nn cur_*

%%
cur_basis = get_k_mat(glm_quad(end-1));
n_bvs = size(cur_basis,2);

nmods = [20:4:52];
for n = 7:length(nmods)
    fprintf('%d subunits\n',nmods(n));
    STCcf_0 = randn(n_bvs,nmods(n));
    %normalize
    for i = 1:nmods(n); STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    init_kerns = cur_basis*STCcf_0;
    init_signs = [1*ones(1,nmods(n)/2) -1*ones(1,nmods(n)/2)];
    % init_signs = ones(1,nmods);
%     init_betas = 2*ones(1,nmods(n));
    
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
    
    glm_stcb = createGLM_lexp_rot(cur_basis,STCcf_0,init_signs,ones(size(init_signs)),defmod,basis);
    glm_stcb.image_type = '1d';
    for i = 1:length(glm_stcb.mods)
        glm_stcb.mods(i).nltype = 'threshlin';
        glm_stcb.mods(i).nly = glm_stcb.mods(i).nlx;
        glm_stcb.mods(i).nly(glm_stcb.mods(i).nlx < 0) = 0;
    end
    glm_stcb.spk_nl = 'logexp';
    % glm_stcb.spk_nl = 'exp';pl
    [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,cur_tr_stimemb*cur_basis);
    kern_output = cur_tr_stimemb*cur_basis;
    for i = 1:length(glm_lexp_rot(n).mods)
        glm_lexp_rot(n).mods(i).locLambda = 0;
    end
    glm_lexp_rot2(n) = fitGLM_lexp_rot(glm_lexp_rot(n),kern_output,cur_tr_spkbns,'tots');
    rotxvLL(n) = getLLGLM_lexp_rot(glm_lexp_rot2(n),cur_xv_stimemb*cur_basis,cur_xv_spkbns,'none'); 
    
    loc_lexp(n) = glm_lexp_rot2(n);
    for i = 1:length(loc_lexp(n).mods)
        loc_lexp(n).mods(i).locLambda = 40;
        loc_lexp(n).mods(i).locSigmaX = 4;
        loc_lexp(n).mods(i).locSigmaT = 3;
    end
    loc_lexp(n) = fitGLM_lexp_rot(loc_lexp(n),kern_output,cur_tr_spkbns,'tots');
    loc_xvLL(n) = getLLGLM_lexp_rot(loc_lexp(n),cur_xv_stimemb*cur_basis,cur_xv_spkbns,'none');
    
end

%%
n = 6;
cur_mod = glm_tlin(n);
for i = 1:length(cur_mod.mods)
    cur_mod.mods(i).lambda_dX = 20;
    cur_mod.mods(i).lambda_dT = 20;
    cur_mod.mods(i).lambda_L1x = 20;
    cur_mod.mods(i).locLambda = 20;
end
cur_mod = fitGLM_lexp(cur_mod,cur_tr_stimemb,cur_tr_spkbns,'tots');
cur_mod_xvLL = getLLGLM_lexp(cur_mod,cur_xv_stimemb,cur_xv_spkbns,'none');
%%
nmods = [12:4:72];
for n = 15:length(nmods)
    fprintf('%d of %d mods\n',nmods(n),length(nmods));
    STCcf_0 = randn(n_bvs,nmods(n));
    %normalize
    for i = 1:nmods(n); STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    init_kerns = cur_basis*STCcf_0;
    init_signs = [1*ones(1,nmods(n)/2) -1*ones(1,nmods(n)/2)];
    defmod = pars.defmod;
    defmod.h(1:end-1) = [];%restrict PSC term to delta function
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
    
    glm_stcb = createGLM_lexp_rot(cur_basis,STCcf_0,init_signs,ones(size(init_signs)),defmod,basis);
    glm_stcb.image_type = '1d';
    for i = 1:length(glm_stcb.mods)
        glm_stcb.mods(i).nltype = 'threshlin';
        glm_stcb.mods(i).nly = glm_stcb.mods(i).nlx;
        glm_stcb.mods(i).nly(glm_stcb.mods(i).nlx < 0) = 0;
    end
    glm_stcb.spk_nl = 'logexp';
    % glm_stcb.spk_nl = 'exp';pl
    [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,cur_tr_stimemb*cur_basis);
    kern_output = cur_tr_stimemb*cur_basis;
    cur_lexp_rot_tlin(n) = fitGLM_lexp_rot(glm_stcb,kern_output,cur_tr_spkbns,'tots');
    cur_lexp_rotxvLL_tlin(n) = getLLGLM_lexp_rot(cur_lexp_rot_tlin(n),cur_xv_stimemb*cur_basis,cur_xv_spkbns,'tots');
    
    loc_lexp_tlin(n) = cur_lexp_rot_tlin(n);
    for i = 1:length(loc_lexp_tlin(n).mods)
        loc_lexp_tlin(n).mods(i).locLambda = 40;
        loc_lexp_tlin(n).mods(i).locSigmaX = 4;
        loc_lexp_tlin(n).mods(i).locSigmaT = 3;
    end
    loc_lexp_tlin(n) = fitGLM_lexp_rot(loc_lexp_tlin(n),kern_output,cur_tr_spkbns,'tots');
    loc_xvLL_tlin(n) = getLLGLM_lexp_rot(loc_lexp_tlin(n),cur_xv_stimemb*cur_basis,cur_xv_spkbns,'none');
end

%%
temp = loc_lexp(end-1);
[temp,norm_vals] = normalizeRFs_STCB(temp,kern_output);
s_mat = get_STCcf_mat(temp);
temp = fitWeights_lexp(temp,kern_output*s_mat,cur_tr_spkbns,0);
temp2 = fitNL_lexp(temp,kern_output*s_mat,cur_tr_spkbns,0);
for i = 1:length(temp.mods)
    temp2.mods(i).nltype = 'uncon';
    temp2.mods(i).locLambda = 40*mean(norm_vals);
end
temp3 = fitGLM_lexp_rot(temp2,kern_output,cur_tr_spkbns,'tots');


%%
% nmods = 8;
% cur_basis = get_k_mat
% n_bvs = size(cur_basis,2);
% STCcf_0 = randn(n_bvs,nmods);
% %normalize
% for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
% % init_kerns = cur_basis*STCcf_0;
% init_signs = [1*ones(1,nmods/2) -1*ones(1,nmods/2)];
% init_betas = 3*ones(1,nmods);
%
% glm_stcb = createGLM_lexp_rot(cur_basis,STCcf_0,init_signs,init_betas,defmod,'pix',[],[]);
% glm_stcb.image_type = '1d';
% glm_stcb.spk_nl = 'logexp';
% kern_output = cur_tr_stimemb*cur_basis;
% glm_lexp = fitGLM_lexp_rot(glm_stcb,kern_output,cur_tr_spkbns,'tots',100);
%
% %%
% for i = 1:5
%     glm_lexp2 = glm_lexp;
%     glm_lexp2 = update_ARD_priors(glm_lexp2);
%     pens(i,:) = arrayfun(@(x) x.lambda_L2x,glm_lexp2.mods);
%     glm_lexp2 = fitGLM_lexp(glm_lexp2,cur_tr_stimemb,cur_tr_spkbns,'full');
% end
% %%
% clear weights
% lambda_vals = logspace(log10(1),log10(2e3),20);
% for i = 1:length(lambda_vals)
%     i
%     glm_lexp2(i) = glm_lexp;
%     glm_lexp2(i).lambdaW = lambda_vals(i);
%     glm_lexp2(i) = fitWeights_lexp_2(glm_lexp2(i),cur_tr_stimemb*get_k_mat(glm_lexp2(i)),cur_tr_spkbns);
%     cur_LL(i) = glm_lexp2(i).LL;
%     weights(i,:) = arrayfun(@(x) x.w,glm_lexp2(i).mods);
%     % glm_lexp4 = fitGLM_lexp(glm_lexp3,stim_emb,spikebins,'tots',100);
% end
% %%
% glm_lexp3 = glm_lexp;
% glm_lexp3.lambdaW = 150;
% glm_lexp3 = fitWeights_lexp_2(glm_lexp3,cur_tr_stimemb*get_k_mat(glm_lexp3),cur_tr_spkbns);
% cur_w = arrayfun(@(x) x.w,glm_lexp3.mods);
% glm_lexp3.mods(cur_w==0) = [];
% glm_lexp3 = fitGLM_lexp(glm_lexp3,cur_tr_stimemb,cur_tr_spkbns,'tots');
% cur_w = arrayfun(@(x) x.w,glm_lexp3.mods);
% glm_lexp3.mods(cur_w==0) = [];
%
% %%
% avg_rate = length(cur_tr_spkbns)/length(tr_inds{xv});
%
% null_pred = ones(size(tr_inds{xv}))*avg_rate;
%
% null_LL = -(sum(log(null_pred(cur_tr_spkbns))) - sum(null_pred))/length(cur_tr_spkbns);
