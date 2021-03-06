clear all;

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
npos = 12; nneg = 12;

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
xv = 1; %use only one XV set
% for xv = 1:nfold

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
nmods = 24  ;
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
xvspkbns = cur_xv_spkbns(cur_xv_spkbns>foff & cur_xv_spkbns<xvstimlen)-foff+2;

dim_signs = nan(size(STCbvs,2),1);
dim_signs(posdims) = 1;
dim_signs(negdims) = -1;


%% 
n_its = 50;
for it = 1:n_its
fprintf('ITERATION: %d\n\n\n\n',it);
mod_signs = nan(nmods,1);
posmods = 1:(nmods/2); negmods = (nmods/2+1):nmods;
posmod_inds = 1:length(posmods);
if nmods <= Nstcbf
    negmod_inds = (length(posdims)+1):(length(posdims)+length(negmods));
    unused_stcs = setdiff(1:nSTCbvs,[posmod_inds negmod_inds]);
    mod_signs(posmods) = 1;
    mod_signs(negmods) = -1;
else
    extra_mods = (Nstcbf+1):nmods;
    mod_signs(posdims) = 1;
    mod_signs(negdims) = -1;
    mod_signs(extra_mods(mod(extra_mods,2)==0)) = 1;
    mod_signs(extra_mods(mod(extra_mods,2)==1)) = -1;
    unused_stcs = [];
end
    
%initialize on STC dims
% STCcf_0 = eye(nSTCbvs);
% STCcf_0(:,unused_stcs) = [];
% if nmods > nSTCbvs
%     n_extra = nmods-nSTCbvs;
%     STCcf_0 = [STCcf_0 randn(nSTCbvs,n_extra)];
% end
% random initialization
STCcf_0 = randn(nSTCbvs,nmods);

%make sure expansive subunits start out in expansive subspace, etc
STCcf_0(negdims,mod_signs==1) = 0;
STCcf_0(posdims,mod_signs==-1) = 0;

%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;

%initialize model
glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
[glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);

glm_stcb.lambdaW = 0; %sparseness on model weights

%fit model
stc_posneg_mod{it} = fitNLHI_stcb_nonlpsc(glm_stcb,cur_tr_stim,tr_cspkbs,'tots');

stc_posneg_xvLL(it) = getLLGLM_STCBF(stc_posneg_mod{it},xvkern_output,xvspkbns,'none');
    
[spatial_profiles, temporal_profiles, weights, mod_type, space_COM, temp_COM] = ...
    compute_mod_stats(stc_posneg_mod{it});
pos_inds = find(weights > 0); neg_inds = find(weights < 0);
[~,space_com_ord_pos] = sort(space_COM(pos_inds));
[~,space_com_ord_neg] = sort(space_COM(neg_inds));
used_ord = [pos_inds(space_com_ord_pos); neg_inds(space_com_ord_neg)];
stc_posneg_mod{it}.mods = stc_posneg_mod{it}.mods(used_ord);

cd ~/James_scripts/stc_sparse_test_figs/
dataname = sprintf('randitcalc_%dmods_%s',nmods,tuid);
save(dataname,'nfold','nmods','stc_posneg*')
end

%%
all_kmat = [];
all_weights = [];
all_coms = [];
all_itnums = [];
for it = 1:n_its
    [spatial_profiles, temporal_profiles, weights, mod_type, space_COM, temp_COM] = ...
        compute_mod_stats(stc_posneg_mod{it});
    all_kmat = [all_kmat get_k_mat(stc_posneg_mod{it})];
    all_weights = [all_weights; weights(:)];
    all_coms = [all_coms; space_COM(:)];
    all_itnums = [all_itnums; it*ones(nmods,1)];
end
all_pos = find(all_weights > 0);
all_neg = find(all_weights < 0);

npoints = 50; %define spatially interpolated axis for aligned profiles
[arfsE,profmeansE] = alignRFs(all_kmat(:,all_pos),sdim,npoints,0);
[arfsS,profmeansS] = alignRFs(all_kmat(:,all_neg),sdim,npoints,0);

[evecsE,pcofsE]    = findPCAs(arfsE,npoints,4,0);
[evecsS,pcofsS]    = findPCAs(arfsS,npoints,4,0);

% [arfs,profmeans] = alignRFs(kmat,sdim,npoints);
% [evecs,pcofs]    = findPCAs(arfs,npoints,4);
clear i
tempE = pcofsE(1,:) + i*pcofsE(2,:);
tempS = pcofsS(1,:) + i*pcofsS(2,:);
nphasE = angle(tempE);
nphasS = angle(tempS);

ncofsE = normvecs(pcofsE(1:2,:));  
ncofsS = normvecs(pcofsS(1:2,:));  %normalized PCA scores (first 2 components) for each filter
% nphasE = 180*asin(ncofsE(1,:))/pi;
% nphasS = 180*asin(ncofsS(1,:))/pi %compute direction of each filter in 2-PCA space

% nphasE = mod(nphasE,pi);
% nphasS = mod(nphasS,pi);

all_nphasE = zeros(size(all_coms));
all_nphasS = zeros(size(all_coms));
all_nphasE(all_pos) = nphasE; all_nphasS(all_neg) = nphasS;

all_coms = [all_coms; all_coms];
all_itnums = [all_itnums; all_itnums];
all_neg = [all_neg; (all_neg + nmods*n_its)];
all_pos = [all_pos; (all_pos + nmods*n_its)];
all_nphasE = [all_nphasE; (all_nphasE+2*pi)];
all_nphasS = [all_nphasS; (all_nphasS+2*pi)];
all_weights = [all_weights; all_weights];

%%
figure
plot1dfilterbank(evecsE(:,1:2),50); colormap(jet);
title('Excitatory subspace PCs')
figure
plot1dfilterbank(evecsS(:,1:2),50); colormap(jet);
title('Suppresive subspace PCs')

%%
f1 = figure; hold on; f2 = figure; hold on
cmap = colormap(jet(n_its));
for it = 1:n_its
    cur_pos = all_pos(all_itnums(all_pos) == it);
    cur_neg = all_neg(all_itnums(all_neg) == it);
figure(1)
scatter(all_coms(cur_pos),all_nphasE(cur_pos),all_weights(cur_pos)*50,repmat(cmap(it,:),length(cur_pos),1))
figure(2)
scatter(all_coms(cur_neg),all_nphasS(cur_neg),-all_weights(cur_neg)*50,repmat(cmap(it,:),length(cur_neg),1))
end
figure(1)
xlabel('Position','fontsize',14)
ylabel('PCA Phase (radians)','fontsize',14)
figure(2)
xlabel('Position','fontsize',14)
ylabel('PCA Phase (radians)','fontsize',14)

%%
n = 2^8;
min_XY = [9 -pi]; max_XY = [18 3*pi];
data = [all_coms(all_pos) all_nphasE(all_pos)];
[bandwidth_e,density_e,X,Y]=kde2d(data,n,min_XY,max_XY);
figure
contourf(X,Y,density_e,30)
xlabel('Position','fontsize',14)
ylabel('PCA Phase (radians)','fontsize',14)

data = [all_coms(all_neg) all_nphasS(all_neg)];
[bandwidth_i,density_i,X,Y]=kde2d(data,n,min_XY,max_XY);
figure
contourf(X,Y,density_i,30)
xlabel('Position','fontsize',14)
ylabel('PCA Phase (radians)','fontsize',14)

%%
cur_size = length(all_coms);
poss_nums = (1:cur_size)';

[~,n] = min(all_LLs);
fprintf('BEST MODEL LL: %.4f  LP: %.4f\n',all_LLs(n),all_LPs(n));
  figure
contourf(X,Y,density_e,30)
xlabel('Position','fontsize',14)
ylabel('PCA Phase (radians)','fontsize',14)
hold on

cur_set = find(all_itnums(1:cur_size)==n & ismember(poss_nums,all_pos));
scatter(all_coms(cur_set),all_nphasE(cur_set),round(100*all_weights(cur_set)),'w','linewidth',2)
fprintf('LL: %.4f  LP: %.4f\n',all_LLs(n),all_LPs(n));


for n = 1:n_its
  figure
contourf(X,Y,density_e,30)
xlabel('Position','fontsize',14)
ylabel('PCA Phase (radians)','fontsize',14)
hold on

cur_set = find(all_itnums(1:cur_size)==n & ismember(poss_nums,all_pos));
scatter(all_coms(cur_set),all_nphasE(cur_set),round(100*all_weights(cur_set)),'w','linewidth',2)
fprintf('LL: %.4f  LP: %.4f\n',all_LLs(n),all_LPs(n));
pause
close

end
%%
cur_size = length(all_coms);
poss_nums = (1:cur_size)';

[~,n] = min(all_LLs);
fprintf('BEST MODEL LL: %.4f  LP: %.4f\n',all_LLs(n),all_LPs(n));
  figure
contourf(X,Y,density_i,30)
xlabel('Position','fontsize',14)
ylabel('PCA Phase (radians)','fontsize',14)
hold on

cur_set = find(all_itnums(1:cur_size)==n & ismember(poss_nums,all_neg));
scatter(all_coms(cur_set),all_nphasS(cur_set),round(-100*all_weights(cur_set)),'w','linewidth',2)
fprintf('LL: %.4f  LP: %.4f\n',all_LLs(n),all_LPs(n));


for n = 1:n_its
  figure
contourf(X,Y,density_i,30)
xlabel('Position','fontsize',14)
ylabel('PCA Phase (radians)','fontsize',14)
hold on

cur_set = find(all_itnums(1:cur_size)==n & ismember(poss_nums,all_neg));
scatter(all_coms(cur_set),all_nphasS(cur_set),round(-100*all_weights(cur_set)),'w','linewidth',2)
fprintf('LL: %.4f  LP: %.4f\n',all_LLs(n),all_LPs(n));
pause
close

end

%%
e_com_rep = [];
e_phase_rep = [];
for i = 1:length(all_pos)
    e_com_rep = [e_com_rep; repmat(all_coms(all_pos(i)),round(100*all_weights(all_pos(i))),1)];
    e_phase_rep = [e_phase_rep; repmat(all_nphasE(all_pos(i)),round(100*all_weights(all_pos(i))),1)];
end
data = [e_com_rep e_phase_rep];
[bandwidth,density,X,Y]=kde2d(data,n,min_XY,max_XY,bandwidth_e);
figure
contourf(X,Y,density,30)
xlabel('Position','fontsize',14)
ylabel('PCA Phase (radians)','fontsize',14)

%%
i_com_rep = [];
i_phase_rep = [];
for i = 1:length(all_neg)
    i_com_rep = [i_com_rep; repmat(all_coms(all_neg(i)),round(-100*all_weights(all_neg(i))),1)];
    i_phase_rep = [i_phase_rep; repmat(all_nphasS(all_neg(i)),round(-100*all_weights(all_neg(i))),1)];
end
data = [i_com_rep i_phase_rep];
[bandwidth,density,X,Y]=kde2d(data,n,min_XY,max_XY,bandwidth_i);
figure
contourf(X,Y,density,30)
xlabel('Position','fontsize',14)
ylabel('PCA Phase (radians)','fontsize',14)

%%
% %% for STC
% cur_nmods = nSTCbvs;
% mod_signs = nan(cur_nmods,1);
% 
% posmods = 1:(cur_nmods/2); negmods = (cur_nmods/2+1):cur_nmods;
% posmod_inds = 1:length(posmods);
% if cur_nmods <= Nstcbf
%     negmod_inds = (length(posdims)+1):(length(posdims)+length(negmods));
%     unused_stcs = setdiff(1:nSTCbvs,[posmod_inds negmod_inds]);
%     mod_signs(posmods) = 1;
%     mod_signs(negmods) = -1;
% else
%     extra_mods = (Nstcbf+1):cur_nmods;
%     mod_signs(posdims) = 1;
%     mod_signs(negdims) = -1;
%     mod_signs(extra_mods(mod(extra_mods,2)==0)) = 1;
%     mod_signs(extra_mods(mod(extra_mods,2)==1)) = -1;
%     unused_stcs = [];
% end
% 
% %initialize on STC dims
% STCcf_0 = eye(nSTCbvs);
% STCcf_0(:,unused_stcs) = [];
% 
% %normalize
% for i = 1:cur_nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
% 
% %initialize model
% glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
% [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);
% 
% glm_stcb.lambdaW = 0; %sparseness on model weights
% 
% %fit model
% % stc_posneg_mod = fitNLHI_stcb_nonlpsc(glm_stcb,trstim,cspkbs,'tots');
% stc_posneg_mod_stc = fitNLHI_stcb_norefine(glm_stcb,cur_tr_stim,tr_cspkbs,'tots');
% 
% xvspkbns = cur_xv_spkbns(cur_xv_spkbns>foff & cur_xv_spkbns<xvstimlen)-foff;
% stc_posneg_xvLL_stc = getLLGLM_STCBF(stc_posneg_mod_stc,xvkern_output,xvspkbns,'none');

% end

