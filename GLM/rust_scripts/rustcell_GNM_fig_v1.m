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
nn = 5; flen = 14;

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
nparts = 20;
partlen = floor(stimlen/nparts);
nfold = 10;
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
% used_stim_dims = [1:sdim-3];
used_stim_dims = [1:sdim];
stim = stim(:,used_stim_dims);
sdim = length(used_stim_dims);

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

n_pos_used = 1;
n_neg_used = 6;
used_stcs = [sta stcs_compareset(:,1:n_pos_used) stcs_compareset(:,nfilts+(1:n_neg_used))];

%% %initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
% defmod.lambda_dX = 50;
% defmod.lambda_L1x = 5; %2
% defmod.lambda_dT = 50;

% for cell 5
defmod.lambda_dX = 0;
defmod.lambda_L1x = 20; %2
defmod.lambda_dT = 0;
defmod.lambda_d2XT = 40;
defmod.lambda_d2X = 0;
defmod.lambda_L2x = 0;
defmod.kscale = 1;
% %for cell 13
% defmod.lambda_dX = 100;
% defmod.lambda_L1x = 10; %2
% defmod.lambda_dT = 100;

% defmod.lambda_dX = 100;
% defmod.lambda_L1x = 5; %2
% defmod.lambda_dT = 100;

% %for cell 34
% defmod.lambda_dX = 300;
% defmod.lambda_L1x = 20; %2
% defmod.lambda_dT = 300;
% 
defmod.lambda_L2x = 0;
defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
basis = 'pix';

nmods = [8];

npos = 12; nneg = 12;
cur_basis  = evecs(:,[1:nneg,length(evs)-npos+1:end]);
cur_basis  = cur_basis(:,end:-1:1);
cur_basis = [cur_basis(:,1:npos) rstcs(:,1:nneg)];

n_bvs = size(cur_basis,2);

for n = length(nmods)
% n = 4;
    fprintf('fitting model with %d subunits\n',nmods(n));
    STCcf_0 = randn(n_bvs,nmods(n));
    %normalize
    for i = 1:nmods(n); STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    init_kerns = cur_basis*STCcf_0;
    % init_signs = [1*ones(1,nmods(n)/2) -1*ones(1,nmods(n)/2)];
%     init_signs = [1*ones(1,nmods(n)/2-1) -1*ones(1,nmods(n)/2+1)];
    init_signs = [1 1 1 -1 -1 -1 -1 -1];
    basis = 'pix';
    glm_stcb = createGLM_quad(init_kerns,init_signs,defmod,basis,[],[]);
    glm_stcb.image_type = '1d';
    glm_stcb.spk_nl = 'logexp';
    glm_quad(n) = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-4,1e-6,[]);
    xvLL_quad(n) = getLLGLM_lexp(glm_quad(n),cur_xv_stimemb,cur_xv_spkbns,'none');
end

% [quad_mod,norm_vals] = normalizeRFs_full(glm_quad(n),cur_tr_stimemb);
% g_mat = cur_tr_stimemb*get_k_mat(quad_mod);
% quad_mod = fitWeights_full(quad_mod,g_mat,cur_tr_spkbns,0);
% for i = 1:length(quad_mod.mods)
%     quad_mod.mods(i).lambda_dX = quad_mod.mods(i).lambda_dX*norm_vals(i)^2;
%     quad_mod.mods(i).lambda_L1x = quad_mod.mods(i).lambda_L1x*norm_vals(i);
%     quad_mod.mods(i).lambda_dT = quad_mod.mods(i).lambda_dT*norm_vals(i)^2;
% end
% quad_mod = fitGLM_lexp(quad_mod,cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-4,1e-6,[]);

%%
quad_old = glm_quad;
qo = glm_quad;
% qo = add_null_filter(qo,'zero',1,'quad');

qo = adjust_all_reg(qo,'lambda_L1x',30);
qo = adjust_all_reg(qo,'lambda_d2XT',150);
qo = adjust_all_reg(qo,'lambda_dX',300);
qo = adjust_all_reg(qo,'lambda_dT',300);
qo = fitGLM_lexp(qo,cur_tr_stimemb,cur_tr_spkbns,'tots',50,1e-4,1e-6,[]);
[cur_LL,~,~,qo_prate] = getLLGLM_lexp(qo,cur_xv_stimemb,cur_xv_spkbns,'none'); cur_LL

%%
% for cell 5
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 100;
defmod.lambda_L1x = 20; %2
defmod.lambda_dT = 100;
defmod.lambda_d2XT = 60;
defmod.lambda_d2X = 0;
defmod.lambda_L2x = 0;
defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
defmod.kscale = 1;
basis = 'pix';

nmods = [12];

npos = 12; nneg = 12;
cur_basis  = evecs(:,[1:nneg,length(evs)-npos+1:end]);
cur_basis  = cur_basis(:,end:-1:1);
cur_basis = [cur_basis(:,1:npos) rstcs(:,1:nneg)];

n_bvs = size(cur_basis,2);

for n = length(nmods)
    % n = 4;
    fprintf('fitting model with %d subunits\n',nmods(n));
    STCcf_0 = randn(n_bvs,nmods(n));
    %normalize
    for i = 1:nmods(n); STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    init_kerns = cur_basis*STCcf_0;
    init_signs = [1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1];
    init_betas = 2*ones(nmods(n),1);
    init_thetas = zeros(nmods(n),1);
    basis = 'pix';
%     glm_stcb = createGLM_lexp(init_kerns,init_signs,init_betas,init_thetas,defmod,basis,[],[]);
    glm_stcb = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
%     glm_stcb = createGLM_rquad(init_kerns,init_signs,defmod,basis,[],[]);
    glm_stcb.image_type = '1d';
    glm_stcb.spk_nl = 'logexp';
    glm_lexp(n) = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-4,1e-6,[]);
    xvLL_lexp(n) = getLLGLM_lexp(glm_lexp(n),cur_xv_stimemb,cur_xv_spkbns,'none');
end

%%
% fo = glm_lexp;
fo_old = fo;

% fo = add_null_filter(fo,'zero',-1,'threshlin');

fo = adjust_all_reg(fo,'lambda_L1x',50); %75
fo = adjust_all_reg(fo,'lambda_d2XT',150); 
fo = adjust_all_reg(fo,'lambda_dX',400);
fo = adjust_all_reg(fo,'lambda_dT',400);
fo = fitGLM_lexp(fo,cur_tr_stimemb,cur_tr_spkbns,'tots',50,1e-4,1e-6,[]);
getLLGLM_lexp(fo,cur_xv_stimemb,cur_xv_spkbns,'none')
%%
glm_lexp2 = add_null_filter(glm_lexp,'zero',1,'threshlin');
glm_lexp2 = fitGLM_lexp(glm_lexp2,cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-4,1e-6,[]);


%%
dnlx = linspace(-3,3,21);
dnly = dnlx;% dnly(dnly < 0) = 0;
fref = fo;
fref = adjust_all_reg(fref,'nltype','uncon');
fref = adjust_all_reg(fref,'nlmon',0);
fref = adjust_all_reg(fref,'lnl2',100);
fref = adjust_all_reg(fref,'nlx',dnlx);
fref = adjust_all_reg(fref,'nly',dnly);

k_mat = get_k_mat(fref);
input = cur_tr_stimemb*k_mat;
fref_u = fitNL_lexp_adj(fref,input,cur_tr_spkbns,0,0);
fref2 = fitGLM_lexp(fref,cur_tr_stimemb,cur_tr_spkbns,'tots',50,1e-4,1e-6,[]); %fit the model
getLLGLM_lexp(fref,cur_xv_stimemb,cur_xv_spkbns,'none')
[cur_LL,~,~,fref2_prate] = getLLGLM_lexp(fref2,cur_xv_stimemb,cur_xv_spkbns,'none');cur_LL

%%

% n_iter = 3;
% lexp_mod{1} = glm_lexp(n);
% for ii = 2:n_iter+1
%     [lexp_mod{ii},norm_vals] = normalizeRFs_full(lexp_mod{ii-1},cur_tr_stimemb);
%     g_mat = cur_tr_stimemb*get_k_mat(lexp_mod{ii});
%     lexp_mod{ii} = fitWeights_full(lexp_mod{ii},g_mat,cur_tr_spkbns,0);
% %     for i = 1:length(lexp_mod{ii}.mods)
% %         lexp_mod{ii}.mods(i).lambda_dX = lexp_mod{ii}.mods(i).lambda_dX*norm_vals(i)^2;
% %         lexp_mod{ii}.mods(i).lambda_L1x = lexp_mod{ii}.mods(i).lambda_L1x*norm_vals(i);
% %         lexp_mod{ii}.mods(i).lambda_dT = lexp_mod{ii}.mods(i).lambda_dT*norm_vals(i)^2;
% %     end
%     lexp_mod{ii} = fitGLM_lexp(lexp_mod{ii},cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-3,1e-5,[]);
% end
% lexp_mod{end} = fitGLM_lexp(lexp_mod{end},cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-4,1e-6,[]);

%%
% nmods = 5;
% 
% defmod.h(1:end-1) = []; %eliminate PSC
% defmod.lnl = 0;
% defmod.lh = 0;
% defmod.lnl2 = 0;
% defmod.lh2 = 0;
% defmod.nlcon = 0;
% defmod.nlmon = 0;
% defmod.locLambda = 0;
% defmod.lambda_dX = 50;
% defmod.lambda_L1x = 15; %2
% defmod.lambda_dT = 50;
% defmod.lambda_L2x = 0;
% defmod.SDIM = sdim;
% defmod.fsdim = sdim;
% defmod.pids = 1:sdim;
% STCcf_0 = randn(n_bvs,nmods);
% init_kerns = cur_basis*STCcf_0;
% init_betas = [ones(1,nmods)*4];
% init_thetas = [ones(1,nmods)];
% init_signs = [1*ones(1,nmods)];
% glm_stcb = createGLM_lexp(init_kerns,init_signs,init_betas,init_thetas,defmod,basis,[],[]);
% glm_stcb.image_type = '1d';
% glm_stcb.spk_nl = 'logexp';
% glm_lexp = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots');
% xvLL_lexp = getLLGLM_lexp(glm_lexp,cur_xv_stimemb,cur_xv_spkbns,'none');
% 
%%
% defmod.h(1:end-1) = []; %eliminate PSC
% defmod.lnl = 0;
% defmod.lh = 0;
% defmod.lnl2 = 0;
% defmod.lh2 = 0;
% defmod.nlcon = 0;
% defmod.nlmon = 0;
% defmod.locLambda = 0;
% % % %for cell 5
% % defmod.lambda_dX = 200;
% % defmod.lambda_L1x = 10; %30
% % defmod.lambda_dT = 200;
% % 
% defmod.lambda_dX = 500;
% defmod.lambda_L1x = 15; %30
% defmod.lambda_dT = 500;
% defmod.lambda_L2x = 0;
% 
% 
% % defmod.lambda_dX = 350;
% % defmod.lambda_L1x = 30; %2
% % defmod.lambda_dT = 350;
% 
% %for cell 34
% % defmod.lambda_dX = 1000;
% % defmod.lambda_L1x = 10; %2
% % defmod.lambda_dT = 1000;
% % defmod.lambda_dX = 400;
% % defmod.lambda_L1x = 15; %2
% % defmod.lambda_dT = 400;
% 
% defmod.lambda_L2x = 0;
% defmod.SDIM = sdim;
% defmod.fsdim = sdim;
% defmod.pids = 1:sdim;
% basis = 'pix';
% 
% npos = 12; nneg = 12;
% % cur_basis  = evecs(:,[1:nneg,length(evs)-npos+1:end]);
% % cur_basis  = cur_basis(:,end:-1:1);
% % cur_basis = [cur_basis(:,1:npos) rstcs(:,1:nneg)];
% cur_basis = get_k_mat(glm_quad(end));
% 
% n_bvs = size(cur_basis,2);
% 
% nmods = [10 10];
% for n = 1:length(nmods);
% fprintf('fitting model with %d subunits\n',nmods(n));
% STCcf_0 = randn(n_bvs,nmods(n));
% %normalize
% for i = 1:nmods(n); STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
% init_kerns = cur_basis*STCcf_0;
% 
% 
% % init_signs = [1*ones(1,nmods(n)/2) -1*ones(1,nmods(n)/2)];
% % init_signs = [1*ones(1,nmods(n))];
% % init_signs = [1*ones(1,floor(nmods(n)/2)-1) -1*(ones(1,ceil(nmods(n)/2+1)))];
% % init_signs = [1*ones(1,2) -1*ones(1,10)];
% % init_signs = [1 1 1 1 1 1 1 1 1 1];
% init_signs = [1 1 1 1 -1 -1 -1 -1 -1 -1];
% 
% init_betas = 2*ones(nmods(n),1);
% init_thetas = zeros(nmods(n),1);
% % init_signs = [1*ones(1,nmods(n))];
% basis = 'pix';
% % glm_stcb = createGLM_rquad(init_kerns,init_signs,defmod,basis,[],[]);
% % glm_stcb.image_type = '1d';y
% % glm_stcb.spk_nl = 'logexp';
% % glm_rquad(n) = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots');
% % xvLL_rquad(n) = getLLGLM_lexp(glm_rquad(n),cur_xv_stimemb,cur_xv_spkbns,'none');
% 
% glm_stcb = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
% % glm_stcb = createGLM_lexp(init_kerns,init_signs,init_betas,init_thetas,defmod,basis,[],[]);
% glm_stcb.image_type = '1d';
% glm_stcb.spk_nl = 'logexp';
% % glm_tlin(n) = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-5,1e-7,sa_params);
% glm_tlin(n) = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-5,1e-7,[]);
% xvLL_tlin(n) = getLLGLM_lexp(glm_tlin(n),cur_xv_stimemb,cur_xv_spkbns,'tots');
% % glm_lexp(n) = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',100);
% % xvLL_lexp(n) = getLLGLM_lexp(glm_tlin(n),cur_xv_stimemb,cur_xv_spkbns,'none');
% 
%     % glm_stcb = createGLM_quad(init_kerns,init_signs,defmod,basis,[],[]);
%     % glm_stcb.image_type = '1d';
%     % glm_stcb.spk_nl = 'logexp';
%     % glm_quad(n) = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots');
%     % xvLL_quad(n) = getLLGLM_lexp(glm_quad(n),cur_xv_stimemb,cur_xv_spkbns,'none');
% 
% end
% 
% % glm_tlin(n) = normalizeRFs_full(glm_tlin(n),cur_tr_stimemb);
% % g_mat = cur_tr_stimemb*get_k_mat(glm_tlin(n));
% % glm_tlin2 = fitWeights_full(glm_tlin(n),g_mat,cur_tr_spkbns,1);
% % for i = 1:length(glm_tlin2.mods)
% %     glm_tlin2.mods(i).lambda_dT = 30;
% %     glm_tlin2.mods(i).lambda_dX = 30;
% %     glm_tlin2.mods(i).lambda_L1x = 2;
% % end
% % glm_tlin2 = fitGLM_lexp(glm_tlin2,cur_tr_stimemb,cur_tr_spkbns,'tots',100,1e-5,1e-7);
% % xvLL_tlin2 = getLLGLM_lexp(glm_tlin2,cur_xv_stimemb,cur_xv_spkbns,'tots');

%% GREEDY TLIN SEARCH
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
% %for cell 5
% defmod.lambda_dX = 200;
% defmod.lambda_L1x = 10; %30
% defmod.lambda_dT = 200;
% defmod.lambda_L2x = 0;
defmod.lambda_dX = 60;
defmod.lambda_L1x = 10; %2
defmod.lambda_dT = 60;
defmod.lambda_d2XT = 60;
defmod.lambda_d2X = 0;
defmod.lambda_L2x = 0;
defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
defmod.kscale = 1;
basis = 'pix';

npos = 12; nneg = 12;
cur_basis = get_k_mat(glm_quad(end));

n_bvs = size(cur_basis,2);

nmods = 1;
fprintf('fitting model with %d subunits\n',nmods);
STCcf_0 = randn(n_bvs,nmods);
%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
init_kerns = cur_basis*STCcf_0;
init_signs = 1;
init_betas = 2*ones(nmods,1);
init_thetas = zeros(nmods,1);
basis = 'pix';
glm_stcb = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
% glm_stcb = createGLM_lexp(init_kerns,init_signs,init_betas,init_thetas,defmod,basis,[],[]);
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
glm_tlin{1} = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-4,1e-6,[]);
xvLL_tlin(1) = getLLGLM_lexp(glm_tlin{1},cur_xv_stimemb,cur_xv_spkbns,'tots');

for nmods = 2:10
    fprintf('Attempting to fit model with %d filters\n\n',nmods);
    
    glm0 = add_null_filter(glm_tlin{nmods-1},'zero',1,'threshlin');
    pos_mod = fitGLM_lexp(glm0,cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-4,1e-6,[]);   
    pos_xvLL = getLLGLM_lexp(pos_mod,cur_xv_stimemb,cur_xv_spkbns,'tots');
    glm0 = add_null_filter(glm_tlin{nmods-1},'zero',-1,'threshlin');
    neg_mod = fitGLM_lexp(glm0,cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-4,1e-6,[]);
    neg_xvLL = getLLGLM_lexp(neg_mod,cur_xv_stimemb,cur_xv_spkbns,'tots');        
    if pos_xvLL < neg_xvLL
        glm_tlin{nmods} = pos_mod;
        xvLL_tlin(nmods) = pos_xvLL;
    else
        glm_tlin{nmods} = neg_mod;
        xvLL_tlin(nmods) = neg_xvLL;
    end
end

%%
init_signs = [1 ones(1,n_pos_used) -1*ones(1,n_neg_used)];
glm_stcb = createGLM_quad(used_stcs,init_signs,defmod,basis,[],[]);
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
stc_mod = fitWeights_lexp(glm_stcb,cur_tr_stimemb*get_k_mat(glm_stcb),cur_tr_spkbns);
xvLL_stc = getLLGLM_lexp(stc_mod,cur_xv_stimemb,cur_xv_spkbns,'tots')

%%
%     glm0 = add_null_filter(glm_tlin{5},'zero',-1,'lexp',2,0);
%  for i = 1:length(glm0.mods)
%     glm0.mods(i).lambda_dX = 100;
%     glm0.mods(i).L1x = 7;
%     glm0.mods(i).lambda_dT = 100;
% end
%    mm = fitGLM_lexp(glm0,cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-4,1e-6,[]);

    %%
%     [new_mod,norm_vals] = normalizeRFs_full(glm_tlin{6},cur_tr_stimemb);
% g_mat = cur_tr_stimemb*get_k_mat(new_mod);
% new_mod = fitWeights_full(new_mod,g_mat,cur_tr_spkbns,0);
% for i = 1:length(new_mod.mods)
%     new_mod.mods(i).lambda_dX = new_mod.mods(i).lambda_dX*norm_vals(i)^2;
%     new_mod.mods(i).lambda_L1x = new_mod.mods(i).lambda_L1x*norm_vals(i);
%     new_mod.mods(i).lambda_dT = new_mod.mods(i).lambda_dT*norm_vals(i)^2;
% end
% new_mod2 = fitGLM_lexp(new_mod,cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-4,1e-6,[]);

%%
% defmod.h(1:end-1) = []; %eliminate PSC
% defmod.lnl = 0;
% defmod.lh = 0;
% defmod.lnl2 = 0;
% defmod.lh2 = 0;
% defmod.nlcon = 0;
% defmod.nlmon = 0;
% defmod.locLambda = 0;
% % defmod.lambda_dX = 400; %was 500
% % defmod.lambda_L1x = 15; %2
% % defmod.lambda_dT = 400; %was 500
% defmod.lambda_dX = 200; %was 500
% defmod.lambda_L1x = 10; %2
% defmod.lambda_dT = 200; %was 500
% 
% defmod.lambda_L2x = 0;
% defmod.SDIM = sdim;
% defmod.fsdim = sdim;
% defmod.pids = 1:sdim;
% 
% basis = 'pix';
% base_mod = glm_quad(1);
% lin_mod = base_mod.mods(1);
% quad_mods = base_mod.mods(2:end);
% 
% nmods = 2*length(quad_mods);
% init_kerns = zeros(flen*sdim,nmods);
% init_signs = zeros(nmods,1);
% for i = 1:length(quad_mods)
%    init_kerns(:,i) = quad_mods(i).k; 
%    init_signs(i) = 1;
% end
% for i = 1:length(quad_mods)
%    init_kerns(:,i+length(quad_mods)) = -quad_mods(i).k; 
%    init_signs(i+length(quad_mods)) = 1;
% end
% % init_kerns(:,end) = lin_mod.k;
% % init_signs(end) = 1;
% init_betas = 2*ones(nmods,1);
% init_thetas = zeros(nmods,1);
% 
% % glm_stcb = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
% glm_stcb = createGLM_rquad(init_kerns,init_signs,defmod,basis,[],[]);
% % glm_stcb = createGLM_lexp(init_kerns,init_signs,init_betas,init_thetas,defmod,basis,[],[]);
% glm_stcb.image_type = '1d';
% glm_stcb.spk_nl = 'logexp';
% glm_tlin8 = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',[],1e-5,1e-7);
% xvLL_tlin8 = getLLGLM_lexp(glm_tlin8,cur_xv_stimemb,cur_xv_spkbns,'tots');

%%
% nmods = [11];
% for n = length(nmods);
%     fprintf('fitting model with %d subunits\n',nmods(n));
%     STCcf_0 = randn(n_bvs,nmods(n));
%     %normalize
%     for i = 1:nmods(n); STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
%     init_kerns = cur_basis*STCcf_0;
%     init_signs = [1*ones(1,3) -1*ones(1,8)];
%     basis = 'pix';
%     
%     glm_stcb = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
%     glm_stcb.image_type = '1d';
%     glm_stcb.spk_nl = 'logexp';
%     
%     for i = 1:2
%         glm_stcb.mods(i).nltype = 'quad';
%         glm_stcb.mods(i).nly = glm_stcb.mods(i).nlx.^2;
%         glm_stcb.mods(i).w = 1/glm_stcb.mods(i).nlx(end);
%         glm_stcb.mods(i).lambda_dX = 200;
%         glm_stcb.mods(i).lambda_L1x = 25; %2
%         glm_stcb.mods(i).lambda_dT = 200;
%     end
%     
%     glm_comb(n) = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',100);
%     xvLL_comb(n) = getLLGLM_lexp(glm_comb(n),cur_xv_stimemb,cur_xv_spkbns,'none');
% end

%% 
% cur_mod = glm_tlin(end);
% % for i = 1:length(cur_mod.mods)
% %     cur_mod.mods(i).lambda_L1x = 10;
% %     cur_mod.mods(i).lambda_dX = 200;
% %     cur_mod.mods(i).lamda_dT = 200;
% % end
% cur_mod = normalizeRFs_full(cur_mod,cur_tr_stimemb);
% cur_mod = fitGLM_lexp(cur_mod,cur_tr_stimemb,cur_tr_spkbns,'tots',50);
% cur_mod_xvLL = getLLGLM_lexp(cur_mod,cur_xv_stimemb,cur_xv_spkbns,'none');
% 
%%
% % cur_mod2 = glm_tlin(end);
% cur_mod2 = glm_quad;
% for i = 1:length(cur_mod2.mods)
%     cur_mod2.mods(i).nltype = 'uncon';
%     cur_mod2.mods(i).lnl2 = 0;
%     cur_mod2.mods(i).nlmon = 0;
% end
% cur_mod2 = normalizeRFs_full(cur_mod2,cur_tr_stimemb);
% k_mat = get_k_mat(cur_mod2);
% inputs = cur_tr_stimemb*k_mat;
% % inputs = applykfilters(cur_mod2,cur_tr_stimemb);
% cur_mod3 = fitNL_jmm(cur_mod2,inputs,cur_tr_spkbns,0);
% cur_mod3_xvLL = getLLGLM_lexp(cur_mod3,cur_xv_stimemb,cur_xv_spkbns,'none');

%%
% cur_mod = glm_rquad(end);
% cur_nmods = length(cur_mod.mods);
% k_mat = get_k_mat(cur_mod);
% k_out = cur_tr_stimemb*k_mat;
% for i = 1:cur_nmods
%     cur_nlx = prctile(k_out(:,i),linspace(0.1,99.99,11));
%     cur_mod.mods(i).nlx = cur_nlx;
% end
% nextNLfit = fitNL_jmm(cur_mod,k_out,cur_tr_spkbns,0);
% xvLL_rquad_NL = getLLGLM_lexp(nextNLfit,cur_xv_stimemb,cur_xv_spkbns,'none');
% glm_rquad_NL = fitGLM_lexp(nextNLfit,cur_tr_stimemb,cur_tr_spkbns,'tots');
% xvLL_rquad_NL = getLLGLM_lexp(glm_rquad_NL,cur_xv_stimemb,cur_xv_spkbns,'none');


%%
% starter = glm_quad(7);
% starter.mods(1) = [];
% nmods = 2*length(starter.mods);
% n = 1;
% %set every other filter to opposite of what the quad filter was
% init_kerns = zeros(sdim*flen,nmods);
% for i = 1:nmods
%     init_kerns(:,i) = starter.mods(ceil(i/2)).k;
%    if mod(i,2) == 0
%       init_kerns(:,i) = -init_kerns(:,i);
%    end
% end
% init_signs = [1*ones(1,nmods/2-1) -1*ones(1,nmods/2+1)];
% basis = 'pix';
% glm_stcb = createGLM_rquad(init_kerns,init_signs,defmod,basis,[],[]);
% glm_stcb.image_type = '1d';
% glm_stcb.spk_nl = 'logexp';
% glm_stcb.const = starter.const;
% k_mat = get_k_mat(glm_stcb);
% bad_ks = find(max(abs(k_mat)) == 0);
% glm_stcb.mods(bad_ks) = [];
% 
% % lin_mod = glm_quad(7).mods(1);
% % glm_stcb.mods = [glm_stcb.mods lin_mod];
% 
% glm_rquad = fitGLM_lexp(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots');
% xvLL_rquad = getLLGLM_lexp(glm_rquad,cur_xv_stimemb,cur_xv_spkbns,'none');

%%
% use_n = 5;
% cur_basis = get_k_mat(glm_quad);
% cur_basis(:,1) = [];
% n_bvs = size(cur_basis,2);
% 
% nmods_tlin = [8];
% for n = 1:length(nmods_tlin)
%     fprintf('%d subunits\n',nmods_tlin(n));
%     STCcf_0 = randn(n_bvs,nmods_tlin(n));
%     %normalize
%     for i = 1:nmods_tlin(n); STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
%     init_kerns = cur_basis*STCcf_0;
% %     init_signs = [1*ones(1,nmods_tlin(n)/2) -1*ones(1,nmods_tlin(n)/2)];
%     init_signs = [1*ones(1,nmods_tlin(n))];
%     
%     defmod.lnl = 0;
%     defmod.lh = 0;
%     defmod.lnl2 = 0;
%     defmod.lh2 = 0;
%     defmod.nlcon = 0;
%     defmod.nlmon = 0;
%     defmod.locLambda = 0;
%     defmod.lambda_dX = 0;
%     defmod.lambda_L1x = 0; %2
%     defmod.lambda_dT = 0;
%     defmod.SDIM = sdim;
%     defmod.fsdim = sdim;
%     defmod.pids = 1:sdim;
%     basis = 'pix';
%     
%     glm_stcb = createGLM_lexp_rot(cur_basis,STCcf_0,init_signs,ones(size(init_signs)),defmod,basis);
%     glm_stcb.image_type = '1d';
%     for i = 1:length(glm_stcb.mods)
%         glm_stcb.mods(i).nltype = 'threshlin';
%         glm_stcb.mods(i).nly = glm_stcb.mods(i).nlx;
%         glm_stcb.mods(i).nly(glm_stcb.mods(i).nlx < 0) = 0;
%         
% %         beta = 2; theta = 0;
% %         glm_stcb.mods(i).nly = 1/beta*log(1+exp(beta*(glm_stcb.mods(i).nlx-theta))) - ...
% %             1/beta*log(1+exp(beta*(glm_stcb.mods(i).nlx(1)-theta)));
%         
%     end
%     glm_stcb.spk_nl = 'logexp';
%     [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,cur_tr_stimemb*cur_basis);
%     kern_output = cur_tr_stimemb*cur_basis;
%     for i = 1:length(glm_stcb.mods)
%         glm_stcb.mods(i).locLambda = 0;
%     end
%     glm_lexp_rot(n) = fitGLM_lexp_rot(glm_stcb,kern_output,cur_tr_spkbns,'tots');
%     rotxvLL(n) = getLLGLM_lexp_rot(glm_lexp_rot(n),cur_xv_stimemb*cur_basis,cur_xv_spkbns,'none'); 
% %     
% %     loc_lexp(n) = glm_lexp_rot(n);
% %     for i = 1:length(loc_lexp(n).mods)
% %         loc_lexp(n).mods(i).locLambda = 40; %40
% %         loc_lexp(n).mods(i).locSigmaX = 3;
% %         loc_lexp(n).mods(i).locSigmaT = 3;
% %     end
% %     loc_lexp(n) = fitGLM_lexp_rot(loc_lexp(n),kern_output,cur_tr_spkbns,'tots');
% %     loc_xvLL(n) = getLLGLM_lexp_rot(loc_lexp(n),cur_xv_stimemb*cur_basis,cur_xv_spkbns,'none');    
% end




%%
cd /Users/James/James_scripts/
save rust_gnm_fig_temp_cell5_v3 xvLL* glm_quad* *_parts glm_tlin 
% save rust_gnm_fig_temp_cell1 loc_* xvLL* glm_quad* *_parts glm_lexp* rotxvLL glm_rquad

