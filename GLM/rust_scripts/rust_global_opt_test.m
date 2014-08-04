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

%%
xv = 1; %use only one XV set
% for xv = 1:nfold

cur_tr_stim = stim(tr_inds{xv},:);
cur_tr_stimemb = makeStimRows(cur_tr_stim,flen);
cur_tr_spkbns = tr_spkbns{xv};

cur_xv_stim = stim(xv_inds{xv},:);
cur_xv_stimemb = makeStimRows(cur_xv_stim,flen);
cur_xv_spkbns = xv_spkbns{xv};

%% precompute the stimulus filtered by each STC kernel for training
%% stim
[klen,nSTCbvs] = size(STCbvs);
flen = klen/sdim;
stimlen = size(cur_tr_stim,1);
Nstcbf = size(STCbvs,2);
% kern_output = zeros(stimlen-flen+1,Nstcbf); %initialize filtered output to be (NT-N_tau+1)xNmods.
% for ikern = 1:Nstcbf;
%     %for the given NL module, store the output of the stimulus filtered by
%     %the internal kernel
%     kern_output(:,ikern)= kfilterInput(cur_tr_stim,STCbvs(:,ikern));
% end

kern_output = cur_tr_stimemb*STCbvs;

%% for XV data
[xvstimlen,sdim] = size(cur_xv_stim);

%precompute the stimulus filtered by each STC kernel
[klen,cur_Nstcbf] = size(STCbvs);
flen = klen/sdim;
% xvkern_output = zeros(xvstimlen-flen+1,1); %initialize filtered output to be (NT-N_tau+1)xNmods.
% for ikern = 1:Nstcbf;
%     %for the given NL module, store the output of the stimulus filtered by
%     %the internal kernel
%     xvkern_output(:,ikern)= kfilterInput(cur_xv_stim,STCbvs(:,ikern));
% end
xvkern_output = cur_xv_stimemb*STCbvs;

%% Initialize model
nmods = 4;
defmod.SDIM=sdim;
defmod.locLambda = 50;
% defmod.locLambda = 0;
defmod.lh = 0;
defmod.lh2 = 0;
defmod.lnl = 0;
defmod.lnl2 = 0;
defmod.hcon = 0;defmod.hmon = 0;defmod.nlcon = 0;defmod.nlmon = 0;
nSTCbvs = size(STCbvs,2);

foff = flen + length(defmod.h);
% tr_cspkbs = cur_tr_spkbns(cur_tr_spkbns>foff & cur_tr_spkbns<stimlen)-foff+2;
% xvspkbns = cur_xv_spkbns(cur_xv_spkbns>foff & cur_xv_spkbns<xvstimlen)-foff+2;
tr_cspkbs = cur_tr_spkbns(cur_tr_spkbns>foff & cur_tr_spkbns<stimlen);
xvspkbns = cur_xv_spkbns(cur_xv_spkbns>foff & cur_xv_spkbns<xvstimlen);

dim_signs = nan(size(STCbvs,2),1);
dim_signs(posdims) = 1;
dim_signs(negdims) = -1;


%%
max_reps = 300;
min_reps = 10;
min_LM_fract = 0.90;
eps = 0.002;
n_its = 200;

nmods = [2:2:20];
for nn = 7:length(nmods)
    fprintf('\n\n%d of %d Filters\n\n',nn,length(nmods));
    %determine distribution of random interpoint distances
    rand_reps = 500;
    init_vals = zeros(rand_reps,nSTCbvs,nmods(nn));
    for r = 1:rand_reps
        % compute average separation between NN initial points
        STCcf_0 = randn(nSTCbvs,nmods(nn));
        %normalize
        for i = 1:nmods(nn); STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
        init_vals(r,:,:) = STCcf_0;
    end
    cur_dists = zeros(1,rand_reps*(rand_reps-1)/2);
    for n = 1:nmods(nn)
        cur_dists = cur_dists + pdist(init_vals(:,:,n),'cosine');
    end
    rdmean = 0.75*mean(cur_dists);
    rdscale = 2*std(cur_dists);
    
    mod_signs = ones(nmods(nn),1);
    
    cur_reps = 0;
    used_cfs{nn} = [];
    LL_vals{nn} = [];
    xvLL_vals{nn} = [];
    LP_vals{nn} = [];
    smallest_dists = [];
    is_optimized{nn} = [];
    rotbv_mod{nn} = [];
    for r = 1:n_its
        fprintf('ITERATION: %d\n\n\n\n',r);
        % random initialization
        STCcf_0 = randn(nSTCbvs,nmods(nn));
        %normalize
        for i = 1:nmods(nn); STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
        
        %initialize model
        glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
        glm_stcb.image_type = '1d';
        glm_stcb.basis = 'pix';
        [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);
        glm_stcb.lambdaW = 0; %sparseness on model weights
        
        %determine LL and LP at current filter point
        glm_stcb = fitWeights_stcb_nonlpsc(glm_stcb,kern_output*STCcf_0,tr_cspkbs,1,1e-3);
        [ll0, ll0p] = getLLGLM_STCBF_nonlpsc(glm_stcb,kern_output,tr_cspkbs,'none');
        xvll = getLLGLM_STCBF(glm_stcb,xvkern_output,xvspkbns,'none');
        
        %store starting point
        used_cfs{nn} = cat(3,used_cfs{nn},STCcf_0);
        LL_vals{nn} = [LL_vals{nn}; ll0]; LP_vals{nn} = [LP_vals{nn}; ll0p];
        xvLL_vals{nn} = [xvLL_vals{nn}; xvll];
        is_optimized{nn} = [is_optimized{nn}; 0];
        
        
        %compute distance in filter space to nearest better value
        better_vals = find(LP_vals{nn} < ll0p);
        n_better_vals = length(better_vals);
        if n_better_vals > 0
            fprintf('%d better values found\n',n_better_vals);
            
            %compute pairwise distances to better points (along the filter
            %n-spheres) as sum of cosine distances
            cur_dists = zeros(n_better_vals,1);
            for n = 1:nmods(nn)
                cur_dists = cur_dists + (1-squeeze(sum(repmat(STCcf_0(:,n),[1 1 n_better_vals]).*used_cfs{nn}(:,n,better_vals))));
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
            rotbv_mod{nn} = [rotbv_mod{nn}; fitNLHI_stcb_nonlpsc(glm_stcb,cur_tr_stimemb,tr_cspkbs,'none',6,2)];
            xvLL = getLLGLM_STCBF(rotbv_mod{nn}(end),xvkern_output,xvspkbns,'none');
            
            %store final point
            STCcf_fin = get_STCcf_mat(rotbv_mod{nn}(end));
            %normalize filter weights to unit length for distance comparisons
            for i = 1:nmods(nn)
                STCcf_fin(:,i) = STCcf_fin(:,i)/norm(STCcf_fin(:,i));
                dist_trav{nn}(r,i) = 1-dot(STCcf_fin(:,i),STCcf_0(:,i));
            end
            used_cfs{nn} = cat(3,used_cfs{nn},STCcf_fin);
            LL_vals{nn} = [LL_vals{nn}; rotbv_mod{nn}(end).LL]; LP_vals{nn} = [LP_vals{nn}; rotbv_mod{nn}(end).LP];
            xvLL_vals{nn} = [xvLL_vals{nn}; xvLL];
            is_optimized{nn} = [is_optimized{nn}; 1];
            
            cur_opt_LP_vals = LP_vals{nn}(is_optimized{nn}==1);
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
    
    for nn = 8
        nn
    %now refine best rotated model
    [~,best_mod] = min(xvLL_vals{nn}(is_optimized{nn}==1));
    STCcf_0 = nan(nSTCbvs,nmods(nn));
    for i = 1:nmods(nn)
        STCcf_0(:,i) = rotbv_mod{nn}(best_mod).mods(i).STCcf;
    end
    basis = 'pix';
    flen = 14;
    %initialize model
    defmod.h(1:end-1) = []; %eliminate PSC
    defmod.lnl = 0;
    defmod.lh = 0;
    defmod.lnl2 = 100;
    defmod.lh2 = 0;
    defmod.nlcon = 1;
    defmod.nlmon = 1;
    defmod.locLambda = 500;
    defmod.locSigma = 2;
    defmod.maxLocPen = 1000;
    defmod.lambda_dX = 50; %350
    defmod.lambda_L1x = 0; %40
    defmod.lambda_dT = 10;
    defmod.pids = 1:sdim;
    defmod.SDIM = sdim;
    defmod.fsdim = sdim;
    
    for i = 1:nmods(nn); init_nls{i} = 'threshlin'; end;
    %define NL types: "uncon, lin, threshlin, quad"
    for i = 1:nmods(nn); nltypes{i} = 'uncon'; end;
    glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
    
    %copy over internal NLs from previous fits
    for i = 1:nmods(nn)
        glm_stcb.mods(i).nly = rotbv_mod{nn}(best_mod).mods(i).nly;
        glm_stcb.mods(i).w = rotbv_mod{nn}(best_mod).mods(i).w;
    end
    glm_stcb.const = rotbv_mod{nn}(best_mod).const;
    
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,cur_tr_stimemb);
    glm_stcb.image_type = '1d';
    fin_glm(nn) = fitNLHI2d_fullbf(glm_stcb,cur_tr_stimemb,tr_cspkbs,'tots',2,2);
    fin_xvLL(nn) = getLLGLM_FULL2d(fin_glm(nn),cur_xv_stimemb,xvspkbns,'none');
    fin_LL(nn) = fin_glm(nn).LL;
    fin_LP(nn) = fin_glm(nn).LP;
    end
end

%%
cd ~/Data/rust
save rust_global_opt_test fin_glm rotbv_mod xvLL* LL_vals LP_vals is_optimized fin_*

%%
for nn = 1:9
%     f1 = figure;
%     plot(LL_vals{nn}(is_optimized{nn}==0),xvLL_vals{nn}(is_optimized{nn}==0),'o')
%     hold on
%     plot(LL_vals{nn}(is_optimized{nn}==1),xvLL_vals{nn}(is_optimized{nn}==1),'ro')
%     title(sprintf('Nmods: %d',nmods(nn)));
%     xlim([1.05 1.4]); ylim([1.05 1.4]);
    
    opt_set = find(is_optimized{nn}==1);
    cur_opt_LP_vals = LP_vals{nn}(opt_set);
    cur_cluster_assignments = get_eps_ball_clusters(cur_opt_LP_vals,eps);
    w = length(unique(cur_cluster_assignments));
    
    xv_avg(nn) = mean(xvLL_vals{nn}(opt_set));
    xv_min(nn) = min(xvLL_vals{nn}(opt_set));
    ll_avg(nn) = mean(LL_vals{nn}(opt_set));
    lp_avg(nn) = mean(LP_vals{nn}(opt_set));
    ll_min(nn) = min(LL_vals{nn}(opt_set));
    lp_min(nn) = min(LP_vals{nn}(opt_set));
    
    %     for s = 1:w
    %         cur_set = find(cur_cluster_assignments== s);
    %         disp(mean(LL_vals{nn}(opt_set(cur_set))));
    %         plotfo1d_nopsc(rotbv_mod{nn}(cur_set(1)),10);
    %         pause
    %         close
    %     end
%     [~,best_mod] = min(xvLL_vals{nn}(opt_set));
%     plotfo1d_nopsc(rotbv_mod{nn}(best_mod),6);
%     
%     plotfo1d_nopsc(fin_glm(nn),6);
%     pause
%     close all
end

%% NOW FIT "STC" model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
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
mod_signs = ones(nSTCbvs,1);
mod_signs(7:end) = -1;
for i = 1:nmods(nn); init_nls{i} = 'pquad'; end;
%define NL types: "uncon, lin, threshlin, quad"
for i = 1:nmods(nn); nltypes{i} = 'uncon'; end;
glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize

glm_stcb.image_type = '1d';
glm_stcb.basis = 'pix';
[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,cur_tr_stimemb);
glm_stcb.lambdaW = 0; %sparseness on model weights
[glm_stcb,COMs] = get_filter_coms_1d(glm_stcb);


%determine LL and LP at current filter point
stc_mod = fitNLw_alt_full(glm_stcb,cur_tr_stimemb,tr_cspkbs);
[ll0, ll0p] = getLLGLM_FULL2d(stc_mod,cur_tr_stimemb,tr_cspkbs,'tots');

xvll = getLLGLM_FULL2d(stc_mod,cur_xv_stimemb,xvspkbns,'none');
