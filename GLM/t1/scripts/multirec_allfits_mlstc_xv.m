%% Load Data

clear all;
addpath('~/James_scripts/GLM')
addpath('~/James_scripts/GLM/t1')
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/Timm/tim1/functions/')

cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

% load spks7576-ns.mat; load PCScores-ns.mat; stype='ns'
% load spks7576-pn.mat; load PCScores-pn.mat; stype='pn'
% load spks7576-ps.mat; load PCScores-ps.mat; stype='ps'

%load z-scored PC data
cd ~/Data/blanche/matlabdata/
load spks7576-all.mat;
load PCScores_z-all.mat;
stype='all';

%throw away extra stimulus dimensions that we definitely won't be using
scorevars(1001:end) = [];
coefs(:,1001:end) = [];
scores(:,1001:end) = [];

pids =1:1024;

%% create XV data
stimlen = size(scores,1);
nparts = 15;
partlen = floor(stimlen/nparts);
nfold = 5;
nxvparts = nparts/nfold;

%boundaries of parts
pbounds  = [(0:nparts-1)*partlen+1;(1:nparts)*partlen]';

%create Nfold different sets of XV and TR data
for i = 1:nfold
    xv_inds{i} = [];
    tr_inds{i} = [];
    
    cur_perm = randperm(nparts);
    cur_xv_parts = sort(cur_perm(1:nxvparts));
    cur_tr_parts = setdiff(1:nparts,cur_xv_parts);
    
    xv_new_inds = nan(stimlen,1);
    for j = 1:length(cur_xv_parts)
        cur_start = pbounds(cur_xv_parts(j),1);
        cur_stop = pbounds(cur_xv_parts(j),2);
        
        cur_inds = (cur_start:cur_stop) - cur_start + length(xv_inds{i}) + 1;
        xv_new_inds(cur_start:cur_stop) = cur_inds;
        
        xv_inds{i} = [xv_inds{i} cur_start:cur_stop];
    end
    
    tr_new_inds = nan(stimlen,1);
    for j = 1:length(cur_tr_parts)
        cur_start = pbounds(cur_tr_parts(j),1);
        cur_stop = pbounds(cur_tr_parts(j),2);
        
        cur_inds = (cur_start:cur_stop) - cur_start + length(tr_inds{i}) + 1;
        tr_new_inds(cur_start:cur_stop) = cur_inds;
        
        tr_inds{i} = [tr_inds{i} cur_start:cur_stop];
    end
end

xv = 1;
ncomps = 300;
compids   = 1:ncomps;
scorevars = var(scores(tr_inds{xv},:)); %variance of PCA dimensions within the training set
whitemat = diag(1./sqrt(scorevars(compids))); %whitening (rescaling) transformation from training set
WX        = scores(tr_inds{xv},compids)*whitemat; %whitened training data
pix_conv_mat = coefs(:,compids)';
kern_conv_mat = coefs(:,compids);
X_avg = mean(WX); %average whitened training stim
NT_x = size(WX,1); %number of training samples

%% Cross-validation
WX_xv  = scores(xv_inds{xv},compids)*whitemat;

%% Compute whitened data, and STA/STCs
cd '/Users/James/James_scripts/GLM/t1/allcell_fits_mlstc/'
for t = 1:27
    
    % t = 3;
    tcell = t; 
    tsbs      = 1+floor(aselspks{tcell}/dt);
    
    %spike bins during training data
    tr_spbs = find(ismember(tr_inds{xv},tsbs)); 
    spikebins = tr_spbs(tr_spbs>flen & tr_spbs<(size(WX,1)+flen-1))-flen+1 ;

    %spike bins for XV set
    xv_spbs = find(ismember(xv_inds{xv},tsbs));
    xv_spikebins = xv_spbs(xv_spbs>flen & xv_spbs<(size(WX_xv,1)+flen-1))-flen+1 ;

    % compute STE and STA
    WS        = WX(spikebins,:);
    rsta      = mean(WS) - X_avg;
    rsta = rsta*whitemat; %correct for unequal variance bias in STA
    
%     %project STA out of STE
%     WXp = WX - (WX*rsta')*rsta/norm(rsta);
%     WSp = WXp(spikebins,:);
    
    %Compute STC
    stvcv = cov(WS); %covaraince of STE
    [evecs,evals] = eig(stvcv-eye(ncomps)); %assume exact whitening...
    evs{t} = diag(evals);
    
    npos=8; nneg=4; %retain a handful of positive and negative dimensions to check
    stcs  = evecs(:,[1:nneg,length(evs{t})-npos+1:end]); stcs  = stcs(:,end:-1:1);
    rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
    STCbvs = [stcs(:,1:npos) rstcs(:,1:nneg)]; 
    STCbvs = (STCbvs'*whitemat)'; %correct for unequal variance bias in STC dims
    STCbvs = [rsta' STCbvs];
        
    kimages = STCbvs'*pix_conv_mat;
    % f1 = figure('Name',stype); plotfilterbank(kimages',32,pids)

 %     n_reps = 500;
%     rand_sta_mag = zeros(n_reps,1);
%     rand_evals = zeros(n_reps,ncomps);
%     for r = 1:n_reps
%         r
%         rand_spkbins = mod(spikebins+ceil(rand*NT_x),NT_x)+1;
%         WS_rand        = WX(rand_spkbins,:);
%         rand_sta      = mean(WS_rand) - X_avg;
%         rand_sta = rand_sta*whitemat;
%         rand_sta_mag(r) = norm(rand_sta);
%         %WXp = WX - (WX*rand_sta')*rand_sta/rand_sta_mag(r);
%         %WSp = WXp(rand_spkbins,:);
%         %[evecs_rand,evals_rand] = eig(cov(WSp)-cov(WXp));
%         [evecs_rand,evals_rand] = eig(cov(WS_rand)-eye(ncomps));
%         rand_evals(r,:) = diag(evals_rand);
%     end    
%     
      
    %% TRY FITTING SEQUENCE OF STC MODELS WTIH INCREASING DIMENSIONALITY
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
    defmod.lambda_L1x = 0; 
    defmod.lambda_dT = 0;
    defmod.fsdim = ncomps;
    defmod.pids = 1:ncomps;
    basis = 'pix';

    %first fit sta model
    cur_ndims = 1;
    cur_basis = STCbvs(:,1); %just use STA
    STCcf_0 = eye(cur_ndims);
    init_nls{1} = 'lin';
    nltypes{1} = 'lin';
    glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
    glm_stcb.image_type = '2d';
    stc_glm{cur_ndims} = fitstc_fullbf(glm_stcb,WX,spikebins,'tots',4); %fit the model
    %pull out the Pixel-space filters
    ker_pix = get_k_mat(stc_glm{cur_ndims})'*pix_conv_mat;
    for i = 1:cur_ndims
        stc_glm{cur_ndims}.mods(i).pix = ker_pix(i,:)';
    end
    stc_glm{cur_ndims}.xvLL = getLLGLM_FULL2d(stc_glm{cur_ndims},WX_xv,xv_spikebins,'none'); %determine XVLL
    
    %now consider adding a sequence of positive eigenvectors
    delta_xvLL = -Inf;
    while delta_xvLL < 0 %while the XVLL is improving, try adding more eigenvectors
        fprintf('Fitting model with %d Expansive STC dims\n',cur_ndims);
        cur_ndims = cur_ndims + 1;
        cur_basis = STCbvs(:,1:cur_ndims);
        STCcf_0 = eye(cur_ndims);
        %these are expansive quadratic components
        for jj = 2:cur_ndims 
            init_nls{jj} = 'pquad';
            nltypes{jj} = 'quad';
        end
        glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
        [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
        glm_stcb.image_type = '2d';
        stc_glm{cur_ndims} = fitstc_fullbf(glm_stcb,WX,spikebins,'tots',4);
        ker_pix = get_k_mat(stc_glm{cur_ndims})'*pix_conv_mat;
        for i = 1:cur_ndims
            stc_glm{cur_ndims}.mods(i).pix = ker_pix(i,:)';
        end
        stc_glm{cur_ndims}.xvLL = getLLGLM_FULL2d(stc_glm{cur_ndims},WX_xv,xv_spikebins,'none');
        delta_xvLL = stc_glm{cur_ndims}.xvLL - stc_glm{cur_ndims-1}.xvLL; %check XVLL improvement
    end

    %throw away the last dim
    n_exp_dims = cur_ndims - 1;
    stc_glm(end) = []; 
    
    %now consider adding a sequence of suppressive eigenvectors
    delta_xvLL = -Inf;
    cur_ndims = n_exp_dims;
    while delta_xvLL < 0
        fprintf('Fitting model with %d Suppressive STC dims\n',cur_ndims);
        cur_ndims = cur_ndims + 1;
        cur_basis = STCbvs(:,[1:n_exp_dims (n_exp_dims+1):cur_ndims]);
        STCcf_0 = eye(cur_ndims);
        for jj = 2:n_exp_dims
            init_nls{jj} = 'pquad';
            nltypes{jj} = 'quad';
        end
        for jj = (n_exp_dims+1):cur_ndims
            init_nls{jj} = 'nquad';
            nltypes{jj} = 'quad';
        end
        
        glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
        [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
        glm_stcb.image_type = '2d';
        stc_glm{cur_ndims} = fitstc_fullbf(glm_stcb,WX,spikebins,'tots',4);
        ker_pix = get_k_mat(stc_glm{cur_ndims})'*pix_conv_mat;
        for i = 1:cur_ndims
            stc_glm{cur_ndims}.mods(i).pix = ker_pix(i,:)';
        end
        stc_glm{cur_ndims}.xvLL = getLLGLM_FULL2d(stc_glm{cur_ndims},WX_xv,xv_spikebins,'none');
        delta_xvLL = stc_glm{cur_ndims}.xvLL - stc_glm{cur_ndims-1}.xvLL;
    end
   
    
    % REFINE STC (ML STC)
    used_stc_dims = [1:6 8];
    STCbasis = STCbvs(:,used_stc_dims);
    [klen,Nstcbvs] = size(STCbasis);
    
    %initialize on STC dims
    STCcf_0 = eye(Nstcbvs);
    nmods = Nstcbvs;
    
    %normalize
    for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    
    %initialize model
    defmod.h(1:end-1) = []; %eliminate PSC
    defmod.lnl = 0;
    defmod.lh = 0;
    defmod.lnl2 = 0;
    defmod.lh2 = 0;
    defmod.nlcon = 0;
    defmod.nlmon = 0;
    defmod.locLambda = 0;
    defmod.lambda_dX = 500; %350
    defmod.lambda_L1x = 5; %40
    defmod.lambda_dT = 50;
    defmod.pids = pids;
    
    %define NL initializations: "lin, threshlin, pquad, nquad"
    clear init_nls nltypes
    init_nls{1} = 'lin';
    for i = 2:length(used_stc_dims); init_nls{i} = 'pquad'; end;
%     for i = 6:nmods; init_nls{i} = 'nquad'; end;
    nltypes{1} = 'lin';
    for i = 2:nmods; nltypes{i} = 'quad'; end;
    
    basis = 'pix';
    defmod_t = defmod;
    defmod_t.fsdim = ncomps;
    defmod_t.pids = 1:ncomps;
    defmod_t.lambda_dX = 0;
    defmod_t.lambda_dT = 0;
    defmod_t.lambda_L1x = 0;
    glm_stcb = createGLM2d_fullbf(STCbasis,STCcf_0,[],[],defmod_t,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
    glm_stcb.image_type = '2d';
    stc_glm{t} = fitstc_fullbf(glm_stcb,WX,spikebins,'tots',4);
    ker_filts = get_k_mat(stc_glm{t});
    ker_pix = ker_filts'*pix_conv_mat;
    for i = 1:nmods
        stc_glm{t}.mods(i).pix = ker_pix(i,:)';
    end
    
    % f1 = figure('Name',stype); plotfilterbank(ker_pix',32,pids)
    
    %     ker_filts_doub = [ker_filts -ker_filts];
    ker_filts_doub = ker_filts;
    ker_pix_doub = ker_filts_doub'*pix_conv_mat;
    nmods = size(ker_filts_doub,2);
    %     for i = 1:nmods; init_nls{i} = 'threshlin'; end;
    %     for i = 1:nmods; nltypes{i} = 'threshlin'; end;
    init_nls{1} = 'lin'; nltypes{1} = 'lin';
    for i = 2:nmods; init_nls{i} = 'pquad'; end;
    for i = 2:nmods; nltypes{i} = 'quad'; end;
    %
    basis = 'white';
    glm_stcb = createGLM2d_fullbf(STCbasis,ones(size(STCbasis,2),size(ker_filts_doub,2)),pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
    glm_stcb.image_type = '2d';
    for i = 1:nmods
        temp_k = ker_filts_doub(:,i);
        temp_pix = ker_pix_doub(i,:);
        glm_stcb.mods(i).k = temp_k(:);
        glm_stcb.mods(i).pix = temp_pix(:);
    end
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
    % pix_stc_glm = fitstc_fullbf(glm_stcb,WX,spikebins,'tots',2);
    % % [nstc_glm,norm_vals] = normalizeRFs_full(stc_glm,WX);
    glm_stcb.lambdaW = 50;
    full_glm2 = fitNLHI2d_fullbf(glm_stcb,WX,spikebins,'tots',2);
    w = arrayfun(@(x) x.w,full_glm2.mods);
    full_glm3 = full_glm2;
    full_glm3.mods(w==0) = [];
    
    % f1 = plot2d_mod(full_glm3);
    refined_mod{t} = full_glm3;
    
    used_mod = full_glm3;
    basis_vecs = get_k_mat(used_mod);
    n_bvs = size(basis_vecs,2);
    nmods = 6;
    mod_signs = ones(nmods,1);
    dim_signs = ones(n_bvs,1);
    unused_stcs = (nmods+1):n_bvs;
    
    %%
    flen = 6;
    %initialize model
    defmod.h(1:end-1) = []; %eliminate PSC
    defmod.lnl = 0;
    defmod.lh = 0;
    defmod.lnl2 = 100;
    defmod.lh2 = 0;
    defmod.nlcon = 0;
    defmod.nlmon = 0;
    defmod.locLambda = 50;
    defmod.lambda_dX = 200; %350
    defmod.lambda_L1x = 0; %40
    defmod.lambda_dT = 10;
    defmod.pids = 1:fsdim;
    defmod.SDIM = SDIM;
    defmod.fsdim = fsdim;
    
    %define NL initializations: "lin, threshlin, pquad, nquad"
    clear init_nls nltypes
    for i = 1:nmods; init_nls{i} = 'threshlin'; end;
    %define NL types: "uncon, lin, threshlin, quad"
    for i = 1:nmods; nltypes{i} = 'uncon'; end;
    
    Nstcbf = size(basis_vecs,2);
    klen = size(pix_conv_mat,2);
    flen = klen/fsdim;
    kern_output = WX*basis_vecs;
    
    clear init_vals all_filtproj all_initproj cur_LL cur_LP fin_vals dist_trav *_lp rotbv_mod
    max_reps = 300;
    min_reps = 10;
    min_LM_fract = 0.95;
    eps = 0.002;
    cur_reps = 0;
    used_cfs = [];
    LL_vals = [];
    LP_vals = [];
    smallest_dists = [];
    is_optimized = [];
    rotbv_mod = [];
    for r = 1:10
        r
        %points uniformly distributed on the set of Nmod p-spheres (p=n_bvs)
        STCcf_0 = randn(n_bvs,nmods);
        for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
        
        defmod.image_type = '2d';
        white_props.basis = 'white';
        white_props.pix_conv_mat = pix_conv_mat;
        white_props.kern_conv_mat = kern_conv_mat;
        glm_stcb = create2dGLM0_stcb(basis_vecs,STCcf_0,defmod,ones(nmods,1),ones(nmods,1),nltypes,init_nls,'test',white_props); %initialize
        glm_stcb.lambdaW = 50;
        
        %determine LL and LP at current filter point
        glm_stcb = fitWeights_stcb_nonlpsc(glm_stcb,kern_output*STCcf_0,spikebins,1,1e-3);
        [ll0, ll0p] = getLLGLM2d_STCBF_nonlpsc(glm_stcb,kern_output,spikebins,'none');
        
        rotbv_mod = [rotbv_mod; fitNLHI_stcb2d_nopsc(glm_stcb,WX,spikebins,'none',4,2)];
        LL_vals = [LL_vals rotbv_mod(end).LL];
        LP_vals = [LP_vals rotbv_mod(end).LP];
    end
    [~,best_mod] = min(LL_vals);
    w = arrayfun(@(x) x.w,rotbv_mod(best_mod).mods);
    rotbv_mod(best_mod).mods(w==0) = [];
    
    % f1 = plot2d_mod(rotbv_mod(best_mod));
    
    rotated_mod{t} = rotbv_mod(best_mod);
    
    %% NOW FIT "STC" model
%     used_stc_dims = [1:7];
    STCbasis = STCbvs(:,used_stc_dims);
    nSTCbvs = size(STCbasis,2);
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
    mod_signs = ones(nSTCbvs,1);
    mod_signs(nposdims+2:end) = -1;
    STCcf_0 = eye(nSTCbvs);
    
    clear init_nls nltypes
    init_nls{1} = 'lin';
    for i = 2:nSTCbvs; init_nls{i} = 'pquad'; end;
    nltypes{1} = 'lin';
    for i = 2:nSTCbvs; nltypes{i} = 'quad'; end;
    basis = 'white';
    glm_stcb = createGLM2d_fullbf(STCbasis,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
    glm_stcb.image_type = '2d';
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
    glm_stcb.lambdaW = 0; %sparseness on model weights
    %determine LL and LP at current filter point
    stc_mod{t} = fitNLw_alt_full(glm_stcb,WX,spikebins);
    
    %% NOW FIT "STA" model
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
    mod_signs = 1;
    STCcf_0 = 1;
    clear init_nls nltypes
    init_nls{1} = 'lin';
    %define NL types: "uncon, lin, threshlin, quad"
    nltypes{1} = 'lin';
    basis = 'white';
    glm_stcb = createGLM2d_fullbf(STCbasis(:,1),STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
    glm_stcb.image_type = '2d';
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
    glm_stcb.lambdaW = 0; %sparseness on model weights
    
    %determine LL and LP at current filter point
    sta_mod = fitNLw_alt_full(glm_stcb,WX,spikebins);
    
    %% Cross-validation
    WX_xv  = scores(xv_inds{xv},compids)*whitemat;
    xv_spbs = find(ismember(xv_inds{xv},tsbs));
    xv_spikebins = xv_spbs(xv_spbs>flen & xv_spbs<(size(WX_xv,1)+flen-1))-flen+1 ;
    
    stc_mod{t}.xvLL = getLLGLM_FULL2d(stc_mod{t},WX_xv,xv_spikebins,'none');
    sta_mod_xvLL(t) = getLLGLM_FULL2d(sta_mod,WX_xv,xv_spikebins,'none');
    
    % for r = 1:length(rotbv_mod)
    % rotbv_mod(r).xvLL = getLLGLM_FULL2d(rotbv_mod(r),WX_xv,xv_spikebins,'none');
    % end
    rotated_mod{t}.xvLL = getLLGLM_FULL2d(rotated_mod{t},WX_xv,xv_spikebins,'none');
    stc_glm{t}.xvLL = getLLGLM_FULL2d(stc_glm{t},WX_xv,xv_spikebins,'none');
    refined_mod{t}.xvLL = getLLGLM_FULL2d(refined_mod{t},WX_xv,xv_spikebins,'none');
    
    for r = 1:length(rotbv_mod)
        rotbv_xvLL(t,r) = getLLGLM_FULL2d(rotbv_mod(r),WX_xv,xv_spikebins,'none');
    end
    
end
%%

stcmod_perf = (sta_mod.xvLL - stc_mod.xvLL)/log(2);
stcglm_perf = (sta_mod.xvLL - stc_glm.xvLL)/log(2);
refined_perf = (sta_mod.xvLL - refined_mod{t}.xvLL)/log(2);
rotated_perf = (sta_mod.xvLL - rotated_mod{t}.xvLL)/log(2);

%%

xv = 1;
rate_sm = 0.05;
stimlen_xv = size(WX_xv,1);
t_axis = dt:dt:dt*stimlen_xv;
obs_spike_rate = hist(xv_spikebins,1:length(xv_inds{xv}));
obs_spike_rate = smooth(obs_spike_rate,round(rate_sm/dt))/dt;
pred_rate_stcmod = getPredictedRate_GLM2d(stc_mod{t},WX_xv,'white');
pred_rate_rotmod = getPredictedRate_GLM2d(rotated_mod{t},WX_xv,'white');

%% Get hartley basis predictions
cd /Users/James/Data/blanche/rec_73

npix=32; pids = 1:(npix*npix);
lutable = zeros(721,npix*npix);
for irow = 1:720;
    tor = ori(ctab(irow,2));
    tsf = 10*sfreqCycDeg(ctab(irow,3));
    tpha = phase0(ctab(irow,4));
    lutable(irow,:) = flatten(hartley3p(tor,tsf,tpha,npix));
end;

tstids = stids20;
istim = zeros(length(tstids),npix^2);
for itime = 1:length(tstids); istim(itime,:) = lutable(tstids(itime)+1,:); end;

hartley_stimemb = makeStimRows(istim,flen);

%%
hartley_rotmod_pred = getPredictedRate_GLM2d(rotated_mod{t},hartley_stimemb,'pix');
hartley_stcmod_pred = getPredictedRate_GLM2d(stc_mod,hartley_stimemb,'pix');
hartley_refmod_pred = getPredictedRate_GLM2d(refined_mod{t},hartley_stimemb,'pix');
hartley_stcglm_pred = getPredictedRate_GLM2d(stc_glm,hartley_stimemb,'pix');

%% CREATE FULL CELL RESPONSE MATRIX
n_lags = 10;
NT_h = size(hartley_stimemb,1);
cell_response_mat = zeros(720,n_lags);
for ss = 1:720
    cur_stim_set = find(tstids==ss);
    for nn = 1:n_lags
        temp = cur_stim_set + nn - 1;
        temp(temp > NT_h) = [];
        %         cell_response_mat(ss,nn) = mean(hartley_rotmod_pred(temp));
        %         cell_response_mat(ss,nn) = mean(hartley_stcmod_pred(temp));
        %         cell_response_mat(ss,nn) = mean(hartley_refmod_pred(temp));
        cell_response_mat(ss,nn) = mean(hartley_stcglm_pred(temp));
    end
end
cell_response_array = nan(length(sfreqCycDeg),length(ori),length(phase0),n_lags);
for i = 1:720
    cell_response_array(ctab(i,3),ctab(i,2),ctab(i,4),:) = cell_response_mat(i,:);
end
freq_response_mat = squeeze(mean(mean(cell_response_array,2),3));
ori_response_mat = squeeze(mean(mean(cell_response_array,1),3));
phase_response_mat = squeeze(mean(mean(cell_response_array,1),2));
emb_ori = [fliplr(-ori(2:end)) ori];
emb_ori_response_mat = cat(2,flipdim(ori_response_mat(:,2:end,:),2),ori_response_mat);

% FIND RELAVENT TIME LAGS (BEST ARE 1 and 2-sample lags)
freq_mod_index_mat = squeeze((max(freq_response_mat) - min(freq_response_mat))./(max(freq_response_mat) + min(freq_response_mat)));
%imagesc(freq_mod_index);
[~,best_lag] = max(freq_mod_index_mat);
[~,max_floc] = max(freq_response_mat(:,best_lag));
freq_mod_index = freq_mod_index_mat(best_lag);
preferred_freqs = sfreqCycDeg(max_floc);
% COMPUTE ORIENTATION TUNING
%pull out response at best frequency
overall_orientation_tuning = squeeze(cell_response_array(max_floc,:,:,:));

ph_inv = squeeze(mean(overall_orientation_tuning,2)); %average over phases

ori_mod_index_mat = squeeze((max(ph_inv) - min(ph_inv))./(max(ph_inv) + min(ph_inv)));
[~,best_lag] = max(ori_mod_index_mat);
orientation_tuning = ph_inv(:,best_lag);
ori_mod_index = ori_mod_index_mat(best_lag);

[~,best_orientation] = max(orientation_tuning);
preferred_orientation = ori(best_orientation);
% ALIGN ORIENTATION TUNING
orientation_modulation = (max(orientation_tuning) - min(orientation_tuning))./(max(orientation_tuning) + min(orientation_tuning));
aligned_orientation_tuning = circshift(orientation_tuning,[0 -best_orientation+9]);
phase_dep = squeeze(overall_orientation_tuning(best_orientation,:,best_lag));
phase_mod_index = (max(phase_dep) - min(phase_dep))./(max(phase_dep) + min(phase_dep));
[~,preferred_phase] = max(phase_dep);
preferred_phase = phase0(preferred_phase);

%% GENERATE DRIFTING GRATING SEQUENCE
Nsamps = 500;
dphase = 10;
phase_vals = mod(dphase:dphase:Nsamps*dphase,360);
grating_stim = zeros(Nsamps,npix*npix);
for i = 1:Nsamps
    temp = hartley3p(preferred_orientation,10*preferred_freqs,phase_vals(i),npix);
    grating_stim(i,:) = temp(:);
end
grating_stim_emb = makeStimRows(grating_stim,flen);

grating_rotmod_pred = getPredictedRate_GLM2d(rotated_mod{t},grating_stim_emb,'pix');
grating_stcmod_pred = getPredictedRate_GLM2d(stc_mod,grating_stim_emb,'pix');
grating_refmod_pred = getPredictedRate_GLM2d(refined_mod{t},grating_stim_emb,'pix');

