%% Load Data

clear all;
addpath('~/James_scripts/GLM')
addpath('~/James_scripts/GLM/t1')
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/Timm/tim1/functions/')
addpath(genpath('~/James_Scripts'));

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
ncomps = 400;
compids   = 1:ncomps;
scorevars = var(scores(tr_inds{xv},:)); %variance of PCA dimensions within the training set
whitemat = diag(1./sqrt(scorevars(compids))); %whitening (rescaling) transformation from training set
WX        = scores(tr_inds{xv},compids); %whitened training data
pix_conv_mat = coefs(:,compids)';
kern_conv_mat = coefs(:,compids);
X_avg = mean(WX); %average whitened training stim
NT_x = size(WX,1); %number of training samples

%% Cross-validation
WX_xv  = scores(xv_inds{xv},compids);

%% Compute whitened data, and STA/STCs
cd '/Users/James/James_scripts/GLM/t1/allcells_fits_mlstc_v4/'  
for t = 1:20
    
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
    
    %project STA out of STE
    proj_mat = rsta'*inv(rsta*rsta')*rsta;
    WXp = WX - WX*proj_mat;
    WSp = WXp(spikebins,:);
    
    %Compute STC
    stvcv = cov(WSp); %covaraince of STE
    %[evecs,evals] = eig(stvcv-eye(ncomps)); %assume exact whitening...
    [evecs,evals] = eig(stvcv-cov(WXp));
    evs{t} = diag(evals);
    
    npos=8; nneg=8; %retain a handful of positive and negative dimensions to check
    stcs  = evecs(:,[1:nneg,length(evs{t})-npos+1:end]); stcs  = stcs(:,end:-1:1);
    rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
    STCbvs = [stcs(:,1:npos) rstcs(:,1:nneg)]; 
    STCbvs = (STCbvs'*whitemat)'; %correct for unequal variance bias in STC dims
    STCbvs = [rsta' STCbvs];
        
    kimages = STCbvs'*pix_conv_mat;
    f1 = figure('Name',stype); plotfilterbank(kimages',32,pids)
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30 70]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = sprintf('Cell%d_Init_STC',t);
    print(f1,'-dpng',fname);close all    

      
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
    defmod.lambda_dX = 1000; 
    defmod.lambda_L1x = 10; 
    defmod.lambda_dT = 100;
    defmod.fsdim = fsdim;
    defmod.pids = 1:fsdim;
    basis = 'white';

    clear stc_glm init_nls nltypes
    %first fit sta model
    cur_ndims = 1;
    cur_basis = STCbvs(:,1); %just use STA
    STCcf_0 = eye(cur_ndims);
    init_nls{1} = 'lin';
    nltypes{1} = 'lin';
    glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
         glm_stcb.spk_nl = 'exp';
   glm_stcb.image_type = '2d';
    stc_glm{cur_ndims} = fitstc_fullbf(glm_stcb,WX,spikebins,'tots',4); %fit the model
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
        glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
        [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
        glm_stcb.image_type = '2d';
        glm_stcb.spk_nl = 'exp';
        glm_stcb.lambdaW = 50;
        stc_glm{cur_ndims} = fitstc_fullbf(glm_stcb,WX,spikebins,'tots',4);
        stc_glm{cur_ndims}.xvLL = getLLGLM_FULL2d(stc_glm{cur_ndims},WX_xv,xv_spikebins,'none');
        delta_xvLL = stc_glm{cur_ndims}.xvLL - stc_glm{cur_ndims-1}.xvLL; %check XVLL improvement
    end

    %throw away the last dim
    n_exp_dims(t) = cur_ndims - 1;
    cur_ndims = cur_ndims - 1;
    stc_glm(end) = []; 
        

    %now consider adding a sequence of suppressive eigenvectors
    delta_xvLL = -Inf;
    n_sup_dims(t) = 1;
    while delta_xvLL < 0
        fprintf('Fitting model with %d Suppressive STC dims\n',n_sup_dims(t));
        cur_ndims = cur_ndims + 1;
        cur_basis = STCbvs(:,[1:n_exp_dims(t) (npos+2):(npos+1+n_sup_dims(t))]);
        STCcf_0 = eye(cur_ndims);
        for jj = 2:n_exp_dims(t)
            init_nls{jj} = 'pquad';
            nltypes{jj} = 'quad';
        end
        for jj = (n_exp_dims(t)+1):cur_ndims
            init_nls{jj} = 'nquad';
            nltypes{jj} = 'quad';
        end
        
        glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
        [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
        glm_stcb.image_type = '2d';
        glm_stcb.spk_nl = 'exp';
        glm_stcb.lambdaW = 50;
        stc_glm{cur_ndims} = fitstc_fullbf(glm_stcb,WX,spikebins,'tots',4);
        stc_glm{cur_ndims}.xvLL = getLLGLM_FULL2d(stc_glm{cur_ndims},WX_xv,xv_spikebins,'none');
        delta_xvLL = stc_glm{cur_ndims}.xvLL - stc_glm{cur_ndims-1}.xvLL;
        n_sup_dims(t) = n_sup_dims(t) + 1;
    end

    n_sup_dims(t) = n_sup_dims(t) - 1;
    stc_glm(end) = [];
    
    refined_mod = stc_glm{length(stc_glm)}; %store used STC GLM MOD
     w = arrayfun(@(x) x.w,refined_mod.mods);
    refined_mod.mods(w==0) = [];
               

%     %% NOW REFINE STC ML FIT
%     nmods = size(ker_filts,2);
%     clear init_nls nltypes
%     init_nls{1} = 'lin'; nltypes{1} = 'lin';
%     for i = 2:n_exp_dims(t); init_nls{i} = 'pquad'; end;
%     for i = (n_exp_dims(t)+1):nmods; init_nls{i} = 'nquad'; end;
%     for i = 2:nmods; nltypes{i} = 'quad'; end;
% 
%     defmod.h(1:end-1) = []; %eliminate PSC
%     defmod.lnl = 0;
%     defmod.lh = 0;
%     defmod.lnl2 = 0;
%     defmod.lh2 = 0;
%     defmod.nlcon = 0;
%     defmod.nlmon = 0;
%     defmod.locLambda = 0;
%     defmod.lambda_dX = 500;
%     defmod.lambda_L1x = 3;
%     defmod.lambda_dT = 50;
%     defmod.pids = 1:fsdim;
%     defmod.SDIM = SDIM;
%     defmod.fsdim = fsdim;
%     basis = 'white';
%     glm_stcb = createGLM2d_fullbf(STCbvs(:,1:nmods),eye(nmods),pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
%     glm_stcb.image_type = '2d';
%     for i = 1:nmods
%         temp_k = ker_filts(:,i);
%         temp_pix = ker_pix(i,:);
%         glm_stcb.mods(i).k = temp_k(:);
%         glm_stcb.mods(i).pix = temp_pix(:);
%     end
%             glm_stcb.lambdaW = 50;
% 
%     [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
%     refined_mod = fitNLHI2d_fullbf(glm_stcb,WX,spikebins,'tots',2);
%     w = arrayfun(@(x) x.w,refined_mod.mods);
%     refined_mod.mods(w==0) = [];
    
    f1 = plot2d_mod(refined_mod)
        set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30 5*length(refined_mod.mods)]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = sprintf('Cell%d_REF_STC',t);
    print(f1,'-dpng',fname);close all    

    %%
    used_mod = refined_mod;
    basis_vecs = get_k_mat(used_mod);
    n_bvs = size(basis_vecs,2);
    
    nmods = n_bvs;
    
    %initialize model
    defmod.h(1:end-1) = []; %eliminate PSC
    defmod.lnl = 0;
    defmod.lh = 0;
    defmod.lnl2 = 100;
    defmod.lh2 = 0;
    defmod.nlcon = 0;
    defmod.nlmon = 0;
    defmod.locLambda = 30;
    defmod.lambda_dX = 500; %350
    defmod.lambda_L1x = 0; %40
    defmod.lambda_dT = 50;
    defmod.pids = 1:fsdim;
    defmod.SDIM = SDIM;
    defmod.fsdim = fsdim;
    defmod.image_type = '2d';
    white_props.basis = 'white';
    white_props.pix_conv_mat = pix_conv_mat;
    white_props.kern_conv_mat = kern_conv_mat;
    
    %define NL initializations: "lin, threshlin, pquad, nquad"
    clear init_nls nltypes
    for i = 1:nmods; init_nls{i} = 'threshlin'; end;
    %define NL types: "uncon, lin, threshlin, quad"
    for i = 1:nmods; nltypes{i} = 'uncon'; end;
    
%     kern_output = WX*basis_vecs;
    
if nmods > 1
    LL_vals = [];
    LP_vals = [];
    xv_vals = [];
    rotbv_mod = [];
    for r = 1:10
        r
        %points uniformly distributed on the set of Nmod p-spheres (p=n_bvs)
        STCcf_0 = randn(n_bvs,nmods);
        for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
        
        glm_stcb = create2dGLM0_stcb(basis_vecs,STCcf_0,defmod,ones(nmods,1),ones(nmods,1),nltypes,init_nls,'test',white_props); %initialize        
        glm_stcb.lambdaW = 50;
        glm_stcb.spk_nl = 'logexp';
        rotbv_mod = [rotbv_mod; fitNLHI_stcb2d_nopsc(glm_stcb,WX,spikebins,'tots',8,2)];
        LL_vals = [LL_vals rotbv_mod(end).LL];
        LP_vals = [LP_vals rotbv_mod(end).LP];
        xv_vals = [xv_vals getLLGLM_FULL2d(rotbv_mod(end),WX_xv,xv_spikebins,'none')];
        
    end
    
    [~,best_mod] = min(xv_vals);
    w = arrayfun(@(x) x.w,rotbv_mod(best_mod).mods);
    rotbv_mod(best_mod).mods(w==0) = [];    
    rotated_mod{t} = rotbv_mod(best_mod);
    rotated_mod{t}.xvLL = xv_vals(best_mod);
    
else
    
    STCcf_0 = 1;
    glm_stcb = create2dGLM0_stcb(basis_vecs,STCcf_0,defmod,ones(nmods,1),ones(nmods,1),nltypes,init_nls,'test',white_props); %initialize
    glm_stcb.spk_nl = 'logexp';
    rotated_mod{t} = fitNLw_alt_full(glm_stcb,WX,spikebins);
    rotated_mod{t}.xvLL = getLLGLM_FULL2d(rotated_mod{t},WX_xv,xv_spikebins,'none');
end

f1 = plot2d_mod(rotated_mod{t})
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30 5*length(rotated_mod{t}.mods)]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = sprintf('Cell%d_ROT_STC',t);
    print(f1,'-dpng',fname);close all

        if nmods > 1
    [~,min_LL] = min(LL_vals);
    rot_xvLL_LLsel(t) = xv_vals(min_LL);
    else
        rot_xvLL_LLsel(t) = rot_xvLL(t);
    end

    %% Try fitting overcomplete threshlin  model
    used_mod = refined_mod;
    basis_vecs = get_k_mat(used_mod);
    n_bvs = size(basis_vecs,2);
    
    nmods = 3*n_bvs;
    
    %initialize model
    defmod.h(1:end-1) = []; %eliminate PSC
    defmod.lnl = 0;
    defmod.lh = 0;
    defmod.lnl2 = 100;
    defmod.lh2 = 0;
    defmod.nlcon = 0;
    defmod.nlmon = 0;
    defmod.locLambda = 30;
    defmod.lambda_dX = 500; %350
    defmod.lambda_L1x = 0; %40
    defmod.lambda_dT = 50;
    defmod.pids = 1:fsdim;
    defmod.SDIM = SDIM;
    defmod.fsdim = fsdim;
    defmod.image_type = '2d';
    white_props.basis = 'white';
    white_props.pix_conv_mat = pix_conv_mat;
    white_props.kern_conv_mat = kern_conv_mat;
    
    %define NL initializations: "lin, threshlin, pquad, nquad"
    clear init_nls nltypes
    for i = 1:nmods; init_nls{i} = 'threshlin'; end;
    %define NL types: "uncon, lin, threshlin, quad"
    for i = 1:nmods; nltypes{i} = 'threshlin'; end;
    
%     kern_output = WX*basis_vecs;
    
    LL_vals = [];
    LP_vals = [];
    xv_vals = [];
    rotbv_mod = [];
    for r = 1:50
        r
        %points uniformly distributed on the set of Nmod p-spheres (p=n_bvs)
        STCcf_0 = randn(n_bvs,nmods);
        for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
        
        glm_stcb = create2dGLM0_stcb(basis_vecs,STCcf_0,defmod,ones(nmods,1),ones(nmods,1),nltypes,init_nls,'test',white_props); %initialize        
        glm_stcb.lambdaW = 200;
        glm_stcb.spk_nl = 'logexp';
        rotbv_mod = [rotbv_mod; fitNLHI_stcb2d_nonlpsc(glm_stcb,WX,spikebins,'tots',8,2)];
        LL_vals = [LL_vals rotbv_mod(end).LL];
        LP_vals = [LP_vals rotbv_mod(end).LP];
        xv_vals = [xv_vals getLLGLM_FULL2d(rotbv_mod(end),WX_xv,xv_spikebins,'none')];
        
    end
    
    [~,best_mod] = min(xv_vals);
    w = arrayfun(@(x) x.w,rotbv_mod(best_mod).mods);
    rotbv_mod(best_mod).mods(w==0) = [];    
    rotated_mod_oc{t} = rotbv_mod(best_mod);
    rotated_mod_oc{t}.xvLL = xv_vals(best_mod);
    [~,min_LL] = min(LL_vals);
        rot_xvLL_oc_LLsel(t) = xv_vals(min_LL);


f1 = plot2d_mod(rotated_mod_oc{t})
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30 5*length(rotated_mod_oc{t}.mods)]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = sprintf('Cell%d_ROTOC_STC',t);
    print(f1,'-dpng',fname);close all
    
    rotated_mod_oc{t}.lambdaW = 0;
    rotated_mod_oc_ref{t} = fitNLHI_stcb2d_nopsc(rotated_mod_oc{t},WX,spikebins,'tots',4,2);
    rotated_mod_oc_ref{t}.xvLL = getLLGLM_FULL2d(rotated_mod_oc_ref{t},WX_xv,xv_spikebins,'none');
    
f1 = plot2d_mod(rotated_mod_oc_ref{t})
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30 5*length(rotated_mod_oc_ref{t}.mods)]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = sprintf('Cell%d_ROTOCREF_STC',t);
    print(f1,'-dpng',fname);close all
    
    
    %% NOW FIT "STC" model
%     cur_basis = STCbvs(:,[1:n_exp_dims(t) (npos+2):(npos+1+n_sup_dims(t))]);
    cur_basis = STCbvs(:,[1:cur_ndims]);
    nSTCbvs = size(cur_basis,2);
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
    
    clear init_nls nltypes
    init_nls{1} = 'lin';
    for i = 2:cur_ndims; init_nls{i} = 'pquad'; end;
%     for i = 2:n_exp_dims(t); init_nls{i} = 'pquad'; end;
%     for i = (n_exp_dims(t)+1):(n_exp_dims(t)+n_sup_dims(t))
%         init_nls{i} = 'nquad';
%     end
    nltypes{1} = 'lin';
    for i = 2:nSTCbvs; nltypes{i} = 'quad'; end;
    basis = 'white';
    glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
    glm_stcb.image_type = '2d';
    glm_stcb.spk_nl = 'exp';
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
    glm_stcb.lambdaW = 0; %sparseness on model weights
    %determine LL and LP at current filter point
    stc_mod = fitNLw_alt_full(glm_stcb,WX,spikebins);
    
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
    glm_stcb = createGLM2d_fullbf(STCbvs(:,1),STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
    glm_stcb.image_type = '2d';
    glm_stcb.spk_nl = 'exp';
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
    glm_stcb.lambdaW = 0; %sparseness on model weights
    
    %determine LL and LP at current filter point
    sta_mod = fitNLw_alt_full(glm_stcb,WX,spikebins);
    
    %% Cross-validation  
    %stc model
    sta_LL(t) = sta_mod.LL;
    sta_LP(t) = sta_mod.LP;
    sta_xvLL(t) = getLLGLM_FULL2d(sta_mod,WX_xv,xv_spikebins,'none');

    %stc model
    stc_LL(t) = stc_mod.LL;
    stc_LP(t) = stc_mod.LP;
    stc_xvLL(t) = getLLGLM_FULL2d(stc_mod,WX_xv,xv_spikebins,'none');
    
%     %stc GLM model
%     stc_glm_LL(t) = stc_glm_mod.LL;
%     stc_glm_LP(t) = stc_glm_mod.LP;
%     stc_glm_xvLL(t) = getLLGLM_FULL2d(stc_glm_mod,WX_xv,xv_spikebins,'none');
% 
    %refined model
    ref_LL(t) = refined_mod.LL;
    ref_LP(t) = refined_mod.LP;
    ref_xvLL(t) = getLLGLM_FULL2d(refined_mod,WX_xv,xv_spikebins,'none');
    
    rot_LL(t) = rotated_mod{t}.LL;
    rot_LP(t) = rotated_mod{t}.LP;
    rot_xvLL(t) = rotated_mod{t}.xvLL;   
    rotoc_xvLL(t) = rotated_mod_oc{t}.xvLL;
    rotocref_xvLL(t) = rotated_mod_oc_ref{t}.xvLL;

    save cur_model_dump3 rotated_mod* *_LL *_LP *_xvLL* *_inds
    
    %%
    used_mod = rotated_mod{t};
    pix_mat = get_pix_mat(used_mod);
    n_filts = length(used_mod.mods);
    
    smooth_mat = epanechikov(5)'*epanechikov(5);
    for i = 1:n_filts
        cur_filt_kmat = reshape(pix_mat(:,i),flen,sdim^2);
        cur_filt_smkmat = zeros(size(cur_filt_kmat));
        for j = 1:flen
            temp_k = reshape(cur_filt_kmat(j,:),sdim,sdim);
            cur_sm_kmat = conv2(temp_k,smooth_mat,'same');
            cur_filt_smkmat(j,:) = cur_sm_kmat(:);
        end
        pix_mat(:,i) = cur_filt_smkmat(:);
    end
        
    for i = 1:n_filts
        nly_mat(i,:) = used_mod.mods(i).nly;
        [gabor_fits(i,:),gabor_fitvals(i,:,:)] = get_gaborfits_allslices(pix_mat(:,i),6,32);
    end
    
    nlx = used_mod.mods(1).nlx;
    [~,best_slice_ids] = get2dMaxSlices(pix_mat,flen,sdim,pids);
        
    f1 = figure;
    for ii = 1:n_filts
        
        %compute space-time projections
        gabor_filtmat = squeeze(gabor_fitvals(ii,:,:));
        cur_filtmat = reshape(pix_mat(:,ii),6,fsdim);
               
        filt_projprof = project_2drf(cur_filtmat(:),gabor_fits(ii,best_slice_ids(ii)).theta,sdim);
        subplot(n_filts,3,(ii-1)*3+1)
        imagesc(reshape(cur_filtmat(best_slice_ids(ii),:),sdim,sdim));colormap(jet);
        zmax = max(abs(cur_filtmat(best_slice_ids(ii),:))); caxis([-zmax zmax]);
        title('Best time-slice')
        
        subplot(n_filts,3,(ii-1)*3+2)
        imagesc(1:sdim,-(0:5)*0.02,filt_projprof);
        xm = max(abs(filt_projprof(:)));
        caxis([-xm xm]);
        title('Space-time projection')
        
        subplot(n_filts,3,(ii-1)*3+3)
        plot(nlx,nly_mat(ii,:),'k')
        axis tight; xlim([-3.1 3.1])
        title('Module non-linearity')
        
    end
    
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30 7*length(rotated_mod{t}.mods)]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = sprintf('Cell%d_spctime',t);
    print(f1,'-dpng',fname);close all    

end
%%

