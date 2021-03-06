clear all;
close all;

addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/James_scripts/GLM/')

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

%load precomputed STC filt params
cd /Users/James/James_scripts/stc_sparse_test_figs/stac_allcells
load ./used_stcims.mat

ncells = length(uids)
for nn = 1:ncells
    
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
    
    [sta, stcs, fstim, evs] = getSTCfilters(cur_tr_stim,cur_tr_spkbns,flen,8,8);
    rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
    
    eig_diffs = diff(evs);
    npos = 7;
    nneg = 7;
    
    rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
    nposdims = n_used_stcdims(nn,1); nnegdims = n_used_stcdims(nn,2);
    STCbvs = [sta' stcs(:,1:nposdims) rstcs(:,1:nnegdims)]; %use only expansive subspace
    nSTCbvs = size(STCbvs,2);
    
    %% precompute the stimulus filtered by each STC kernel for training
    %% stim
    kern_output = cur_tr_stimemb*STCbvs;
    % for XV data
    xvkern_output = cur_xv_stimemb*STCbvs;
    
    %% Initialize model
    nmods = nSTCbvs;
    basis = 'pix';
    defmod.h(1:end-1) = []; %eliminate PSC
    defmod.lnl = 0;
    defmod.lh = 0;
    defmod.lnl2 = 100;
    defmod.lh2 = 0;
    defmod.nlcon = 0;
    defmod.nlmon = 0;
    defmod.locLambda = 40;
    defmod.lambda_dX = 20; %350
    defmod.lambda_L1x = 0; %40
    defmod.lambda_dT = 10;
    defmod.pids = 1:sdim;
    defmod.SDIM = sdim;
    defmod.fsdim = sdim;
    
    max_reps = 20;
    min_reps = 3;
    min_LM_fract = 0.90;
    eps = 0.002;
    
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
    for n = 1:nmods
        cur_dists = cur_dists + pdist(init_vals(:,:,n),'cosine');
    end
    rdmean = 0.75*mean(cur_dists);
    rdscale = 2*std(cur_dists);
    
    cur_reps = 0;
    used_cfs = [];
    LL_vals = [];
    xvLL_vals = [];
    LP_vals = [];
    smallest_dists = [];
    is_optimized = [];
    rotbv_mod = [];
    for r = 1:max_reps
        fprintf('ITERATION: %d\n\n\n\n',r);
        % random initialization
        STCcf_0 = randn(nSTCbvs,nmods);
        %normalize
        for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
        
        %initialize model
        glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,ones(nmods,1),ones(nSTCbvs,1),'test'); %initialize
        glm_stcb.image_type = '1d';
        glm_stcb.basis = 'pix';
        [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);
        glm_stcb.lambdaW = 0; %sparseness on model weights
        
        %determine LL and LP at current filter point
        %         glm_stcb = fitWeights_stcb_nonlpsc(glm_stcb,kern_output*STCcf_0,cur_tr_spkbns,1,1e-3);
        %         [ll0, ll0p] = getLLGLM_STCBF_nonlpsc(glm_stcb,kern_output,cur_tr_spkbns,'none');
        glm_stcb = fitNLw_alt(glm_stcb,kern_output,cur_tr_spkbns,1e-3);
        [ll0, ll0p] = getLLGLM_STCBF_nonlpsc(glm_stcb,kern_output,cur_tr_spkbns,'none');
        xvll = getLLGLM_STCBF(glm_stcb,xvkern_output,cur_xv_spkbns,'none');
        
        %store starting point
        used_cfs = cat(3,used_cfs,STCcf_0);
        LL_vals = [LL_vals; ll0]; LP_vals = [LP_vals; ll0p];
        xvLL_vals = [xvLL_vals; xvll];
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
            rotbv_mod = [rotbv_mod; fitNLHI_stcb_nopsc(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'none',6)];
            xvLL = getLLGLM_STCBF(rotbv_mod(end),xvkern_output,cur_xv_spkbns,'none');
            
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
            end
            used_cfs = cat(3,used_cfs,STCcf_fin);
            LL_vals = [LL_vals; rotbv_mod(end).LL]; LP_vals = [LP_vals; rotbv_mod(end).LP];
            xvLL_vals = [xvLL_vals; xvLL];
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
    
    %% now refine best rotated model
    opt_set = find(is_optimized==1);
    [~,best_mod] = min(LP_vals(is_optimized==1));
    xv_min = xvLL_vals(opt_set(best_mod));
    rotbv_mod(best_mod).xvLL = xv_min;
    ll_min = LL_vals(opt_set(best_mod));
    lp_min = LP_vals(opt_set(best_mod));
    STCcf_0 = nan(nSTCbvs,nmods);
    for i = 1:nmods
        STCcf_0(:,i) = rotbv_mod(best_mod).mods(i).STCcf;
    end
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
    
    clear init_nls nltypes
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
    fin_glm = fitNLHI2d_fullbf(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',2,2);
    fin_glm.xvLL = getLLGLM_FULL2d(fin_glm,cur_xv_stimemb,cur_xv_spkbns,'none');
    
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
    mod_signs = ones(nSTCbvs,1);
    mod_signs(nposdims+2:end) = -1;
    STCcf_0 = eye(nSTCbvs);
    clear init_nls nltypes
    init_nls{1} = 'lin';
    for i = 2:(1+nposdims); init_nls{i} = 'pquad'; end;
    for i = (nposdims+2):nSTCbvs; init_nls{i} = 'nquad'; end;
    %define NL types: "uncon, lin, threshlin, quad"
    for i = 1:nSTCbvs; nltypes{i} = 'uncon'; end;
    glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
    
    glm_stcb.image_type = '1d';
    glm_stcb.basis = 'pix';
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,cur_tr_stimemb);
    glm_stcb.lambdaW = 0; %sparseness on model weights
    [glm_stcb,COMs] = get_filter_coms_1d(glm_stcb);
    
    %determine LL and LP at current filter point
    stc_mod = fitNLw_alt_full(glm_stcb,cur_tr_stimemb,cur_tr_spkbns);
    [ll0, ll0p] = getLLGLM_FULL2d(stc_mod,cur_tr_stimemb,cur_tr_spkbns,'tots');
    stc_mod.LL = ll0; stc_mod.LP = ll0p;
    stc_mod.xvLL = getLLGLM_FULL2d(stc_mod,cur_xv_stimemb,cur_xv_spkbns,'none');
    
    %% NOW FIT "STA" model
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
    mod_signs = 1;
    STCcf_0 = 1;
    clear init_nls nltypes
    init_nls{1} = 'lin';
    %define NL types: "uncon, lin, threshlin, quad"
    nltypes{1} = 'uncon'
    glm_stcb = createGLM2d_fullbf(STCbvs(:,1),STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
    
    glm_stcb.image_type = '1d';
    glm_stcb.basis = 'pix';
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,cur_tr_stimemb);
    glm_stcb.lambdaW = 0; %sparseness on model weights
    [glm_stcb,COMs] = get_filter_coms_1d(glm_stcb);
    
    %determine LL and LP at current filter point
    sta_mod = fitNLw_alt_full(glm_stcb,cur_tr_stimemb,cur_tr_spkbns);
    [ll0, ll0p] = getLLGLM_FULL2d(stc_mod,cur_tr_stimemb,cur_tr_spkbns,'tots');
    sta_mod.LL = ll0; sta_mod.LP = ll0p;
    sta_mod.xvLL = getLLGLM_FULL2d(sta_mod,cur_xv_stimemb,cur_xv_spkbns,'none');
    
    %% SAVE DATA
    all_rotmods(nn) = rotbv_mod(best_mod);
    all_fin_mod(nn) = fin_glm;
    all_stc_mod(nn) = stc_mod;
    all_sta_mod(nn) = sta_mod;
    
    all_xv_params(nn).nparts = nparts;
    all_xv_params(nn).nfold = nfold;
    all_xv_params(nn).cur_tr_parts = cur_tr_parts;
    all_xv_params(nn).cur_xv_parts = cur_xv_parts;
    
    cd ~/James_scripts/stc_sparse_test_figs/
    save rust_allcells_modanal all_*
    
    %% CREATE PLOTS
    cd ~/James_scripts/stc_sparse_test_figs/modfits_allcells/
    
    f1 = plotfo1d_nopsc(stc_mod,6);
    set(f1,'PaperSize',[20 10]);
    fname = sprintf('%s-STCmod',uids{nn});
    print(fname,'-dpng');close
    
    f1 = plotfo1d_nopsc(rotbv_mod(best_mod),6);
    set(f1,'PaperSize',[20 10]);
    fname = sprintf('%s-rotmod',uids{nn});
    print(fname,'-dpng');close
    
    f1 = plotfo1d_nopsc(fin_glm,6);
    set(f1,'PaperSize',[20 10]);
    fname = sprintf('%s-finmod',uids{nn});
    print(fname,'-dpng');close
    
    
end