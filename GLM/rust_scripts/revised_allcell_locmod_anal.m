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
for nn = 1:ncells
    
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
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30 10]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = strcat(uids{nn},'_initSTC');
    print(f1,'-dpng',fname);close all
    
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
    defmod.lambda_dX = 50;
    defmod.lambda_L1x = 2; %2
    defmod.lambda_dT = 50;
    defmod.SDIM = sdim;
    defmod.fsdim = sdim;
    defmod.pids = 1:sdim;
    basis = 'pix';
    
    %     clear stc_glm init_nls nltypes
    %     %first fit sta model
    %     cur_ndims = 1;
    %     cur_basis = sta; %just use STA
    %     STCcf_0 = eye(cur_ndims);
    %     init_nls{1} = 'lin';
    %     nltypes{1} = 'lin';
    %     glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('Cell%d',nn)); %initialize
    %     [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,cur_tr_stimemb);
    %     glm_stcb.image_type = '1d';
    %     glm_stcb.spk_nl = 'exp';
    %     stc_glm{cur_ndims} = fitstc_fullbf(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',4); %fit the model
    %     stc_glm{cur_ndims}.xvLL = getLLGLM_FULL2d(stc_glm{cur_ndims},cur_xv_stimemb,cur_xv_spkbns,'none'); %determine XVLL
    
    %     [a,ev_order] = sort(abs(evs{nn}),'descend');
    %     ev_sign = sign(evs{nn}(ev_order));
    %     cur_loc = 1;
    %     dim_set = [];
    %     %now consider adding a sequence of positive eigenvectors
    %     delta_xvLL = -Inf;
    %     while delta_xvLL < 0 %while the XVLL is improving, try adding more eigenvectors
    %         fprintf('Fitting model with %d Expansive and %d Suppressive dims\n',sum(ev_sign(1:cur_loc)==1),sum(ev_sign(1:cur_loc)==-1));
    %         dim_set = [dim_set ev_order(cur_loc)];
    %
    %         cur_basis = [sta STCbvs(:,dim_set)];
    %         STCcf_0 = eye(cur_loc+1);
    %
    %         if ev_sign(cur_loc) == 1
    %             init_nls{1+cur_loc} = 'pquad';
    %         elseif ev_sign(cur_loc) == -1
    %             init_nls{1+cur_loc} = 'nquad';
    %         end
    %         nltypes{1+cur_loc} = 'quad';
    %
    %         glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('Cell%d',nn)); %initialize
    %         [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,cur_tr_stimemb);
    %         glm_stcb.image_type = '1d';
    %         glm_stcb.spk_nl = 'exp';
    %         stc_glm{1+cur_loc} = fitstc_fullbf(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',3);
    %         stc_glm{1+cur_loc}.xvLL = getLLGLM_FULL2d(stc_glm{1+cur_loc},cur_xv_stimemb,cur_xv_spkbns,'none');
    %         delta_xvLL = stc_glm{1+cur_loc}.xvLL - stc_glm{cur_loc}.xvLL; %check XVLL improvement
    %
    %         cur_loc = cur_loc + 1;
    %     end
    
%     delta_evs = [Inf; diff(evs{nn})];
%     n_sup_dims(nn) = find(delta_evs < 0.01,1,'first')-1;
%     n_exp_dims(nn) = find(flipud(delta_evs) <0.01,1,'first')-1;
    n_sup_dims(nn) = n_used_stcdims(nn,2);
    n_exp_dims(nn) = n_used_stcdims(nn,1);
    used_stcs = [1:n_exp_dims(nn) (length(evs{nn}):-1:length(evs{nn})-n_sup_dims(nn)+1)];
    cur_basis = [sta STCbvs(:,used_stcs)];
    STCcf_0 = eye(n_exp_dims(nn)+n_sup_dims(nn)+1);
    clear init_nls nltypes
    init_nls{1} = 'lin'; nltypes{1} = 'lin';
    for i = 1:n_exp_dims(nn)
        init_nls{i+1} = 'pquad';
        nltypes{i+1} = 'quad';
    end
    for i = 1:n_sup_dims(nn)
        init_nls{i+n_exp_dims(nn)+1} = 'nquad';
        nltypes{i+n_exp_dims(nn)+1} = 'quad';
    end
    
    glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('Cell%d',nn)); %initialize
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,cur_tr_stimemb);
    glm_stcb.image_type = '1d';
    glm_stcb.spk_nl = 'exp';
    stc_glm_exp = fitstc_fullbf(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',3);
    stc_glm_exp_xvLL(nn) = getLLGLM_FULL2d(stc_glm_exp,cur_xv_stimemb,cur_xv_spkbns,'none');
    glm_stcb.spk_nl = 'logexp';
    stc_glm_lexp = fitstc_fullbf(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'tots',3);
    stc_glm_lexp_xvLL(nn) = getLLGLM_FULL2d(stc_glm_lexp,cur_xv_stimemb,cur_xv_spkbns,'none');
    
    
    refined_mod{nn} = stc_glm_lexp; %store used STC GLM MOD
    w = arrayfun(@(x) x.w,refined_mod{nn}.mods);
    refined_mod{nn}.mods(w==0) = [];
    
    f1 = plotfo1d_nopsc(refined_mod{nn},5);
    nmods = length(refined_mod{nn}.mods);
    ncols = ceil(nmods/5);
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30*ncols 50]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = strcat(uids{nn},'_refSTC');
    print(f1,'-dpng',fname);close all
    
    %% precompute the stimulus filtered by each STC kernel for training
    %% stim
    cur_basis = get_pix_mat(refined_mod{nn});
    n_bvs = size(cur_basis,2);
    kern_output = cur_tr_stimemb*cur_basis;
    % for XV data
    xvkern_output = cur_xv_stimemb*cur_basis;
    
    %% Initialize model
    nmods = n_bvs*3; %50
    basis = 'pix';
    defmod.h(1:end-1) = []; %eliminate PSC
    defmod.lnl = 0;
    defmod.lh = 0;
    defmod.lnl2 = 100; %100-200
    defmod.lh2 = 0;
    defmod.nlcon = 0;
    defmod.nlmon = 0;
    defmod.locLambda = 10; %10
    defmod.lambda_dX = 50; %20
    defmod.lambda_L1x = 0; %0
    defmod.lambda_dT = 50; %10
    defmod.pids = 1:sdim;
    defmod.SDIM = sdim;
    defmod.fsdim = sdim;
    
    max_reps = 10;
    LL_vals = [];
    xvLL_vals = [];
    LP_vals = [];
    rotbv_mod = [];
    for r = 1:max_reps
        fprintf('ITERATION: %d\n\n\n\n',r);
        % random initialization
        STCcf_0 = randn(n_bvs,nmods);
        %normalize
        for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
        
        %initialize model
        glm_stcb = createGLM0_stcb(cur_basis,STCcf_0,defmod,ones(nmods,1),ones(n_bvs,1),'test'); %initialize
        glm_stcb.image_type = '1d';
        glm_stcb.basis = 'pix';
        glm_stcb.spk_nl = 'logexp';
        [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);
        glm_stcb.lambdaW = 200; %sparseness on model weights
        
        rand_signs = sign(rand(nmods,1)-0.5);
        for i = 1:nmods
            glm_stcb.mods(i).w = rand_signs(i);
        end
        
        %local optimization
        rotbv_mod = [rotbv_mod; fitNLHI_stcb_nonlpsc(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'none',6,2)];
        %                 rotbv_mod = [rotbv_mod; fitNLHI_stcb_nopsc(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'none',4)];
        xvLL = getLLGLM_STCBF(rotbv_mod(end),xvkern_output,cur_xv_spkbns,'none');
        
        [spatial_profiles, temporal_profiles, weights, mod_type, space_COM, temp_COM] = ...
            compute_mod_stats(rotbv_mod(end));
        pos_inds = find(weights > 0); neg_inds = find(weights < 0);
        [~,space_com_ord_pos] = sort(space_COM(pos_inds));
        [~,space_com_ord_neg] = sort(space_COM(neg_inds));
        used_ord = [pos_inds(space_com_ord_pos); neg_inds(space_com_ord_neg)];
        rotbv_mod(end).mods = rotbv_mod(end).mods(used_ord);
        
        LL_vals = [LL_vals; rotbv_mod(end).LL]; LP_vals = [LP_vals; rotbv_mod(end).LP];
        xvLL_vals = [xvLL_vals; xvLL];
        
    end
    
    [~,min_LL] = min(LL_vals);
    overcomplete_rot{nn} = rotbv_mod(min_LL);
    overcomplete_rot_xvLL(nn) = xvLL_vals(min_LL);
    
    f1 = plotfo1d_nopsc(overcomplete_rot{nn},5);
    nmods = length(overcomplete_rot{nn}.mods);
    ncols = ceil(nmods/5);
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30*ncols 50]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = strcat(uids{nn},'_ocRot');
    print(f1,'-dpng',fname);close all
    
    
    overcomplete_rot{nn}.lambdaW = 0;
    overcomplete_rot_ref{nn} = fitNLHI_stcb_nopsc(overcomplete_rot{nn},cur_tr_stimemb,cur_tr_spkbns,'none',4,2);
    overcomplete_rot_ref_xvLL(nn) = getLLGLM_STCBF(overcomplete_rot_ref{nn},xvkern_output,cur_xv_spkbns,'none');
    
    
    f1 = plotfo1d_nopsc(overcomplete_rot_ref{nn},5);
    nmods = length(overcomplete_rot_ref{nn}.mods);
    ncols = ceil(nmods/5);
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30*ncols 50]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = strcat(uids{nn},'_ocRotRef');
    print(f1,'-dpng',fname);close all
   
    %% Initialize model
    nmods = n_bvs;
    basis = 'pix';
    defmod.h(1:end-1) = []; %eliminate PSC
    defmod.lnl = 0;
    defmod.lh = 0;
    defmod.lnl2 = 20;
    defmod.lh2 = 0;
    defmod.nlcon = 0;
    defmod.nlmon = 0;
    defmod.locLambda = 10; %10
    defmod.lambda_dX = 50; %20
    defmod.lambda_L1x = 0; %0
    defmod.lambda_dT = 50; %10
    defmod.pids = 1:sdim;
    defmod.SDIM = sdim;
    defmod.fsdim = sdim;
    
    max_reps = 10;
    LL_vals = [];
    xvLL_vals = [];
    LP_vals = [];
    rotbv_mod = [];
    for r = 1:max_reps
        fprintf('ITERATION: %d\n\n\n\n',r);
        % random initialization
        STCcf_0 = randn(n_bvs,nmods);
        %normalize
        for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
        
        clear init_nls nltypes
        init_nls{1} = 'lin'; nltypes{1} = 'lin';
        cur_signs = sign(rand(nmods-1,1)-0.5);
        for i = 1:length(cur_signs)
            if cur_signs(i)==1
                init_nls{i+1} = 'pquad';
            else
                init_nls{i+1} = 'nquad';
            end
            nltypes{i+1} = 'uncon';
        end
        glm_stcb = createGLM0_stcb(cur_basis,STCcf_0,defmod,ones(nmods,1),ones(n_bvs,1),'test'); %initialize
%         glm_stcb = createGLM0_stcb_connl(cur_basis,STCcf_0,defmod,init_nls,nltypes,'test'); %initialize
        glm_stcb.basis = 'pix';
        [glm_stcb,norm_vals] = normalizeRFs_STCB(glm_stcb,kern_output);
        glm_stcb.image_type = '1d';
        glm_stcb.spk_nl = 'logexp';
        glm_stcb.lambdaW = 0; %sparseness on model weights
        
        
        %local optimization
        rotbv_mod = [rotbv_mod; fitNLHI_stcb_nopsc(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'none',6,2)];
        %         rotbv_mod = [rotbv_mod; fitNLHI_stcb_connl(glm_stcb,cur_tr_stimemb,cur_tr_spkbns,'none',4)];
        xvLL = getLLGLM_STCBF(rotbv_mod(end),xvkern_output,cur_xv_spkbns,'none');
        
        [spatial_profiles, temporal_profiles, weights, mod_type, space_COM, temp_COM] = ...
            compute_mod_stats(rotbv_mod(end));
        pos_inds = find(weights > 0); neg_inds = find(weights < 0);
        [~,space_com_ord_pos] = sort(space_COM(pos_inds));
        [~,space_com_ord_neg] = sort(space_COM(neg_inds));
        used_ord = [pos_inds(space_com_ord_pos); neg_inds(space_com_ord_neg)];
        rotbv_mod(end).mods = rotbv_mod(end).mods(used_ord);
        
        LL_vals = [LL_vals; rotbv_mod(end).LL]; LP_vals = [LP_vals; rotbv_mod(end).LP];
        xvLL_vals = [xvLL_vals; xvLL];
        
    end
    [~,min_LL] = min(LL_vals);
    ref_rotated{nn} = rotbv_mod(min_LL);
    ref_rotated_xvLL(nn) = xvLL_vals(min_LL);
    
    f1 = plotfo1d_nopsc(ref_rotated{nn},5);
    nmods = length(ref_rotated{nn}.mods);
    ncols = ceil(nmods/5);
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30*ncols 50]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    fname = strcat(uids{nn},'_RotRef');
    print(f1,'-dpng',fname);close all

    %%
    %     stimlen = size(cur_tr_stimemb,1);
    %     Robs = zeros(1,stimlen);
    %     ftable = tabulate(cur_tr_spkbns);
    %     Robs(ftable(:,1)) = ftable(:,2);
    %
    %     n_filts = 25;
    %     n_pts = 1000;
    %     alpha = 0.9;
    %     T0 = 1;
    %     T = T0*alpha.^(0:n_pts-1);
    %
    %     radius = 0.02;
    %     old_x = get_STCcf_mat(cur_mod);
    %     old_y = 100;
    %     best_y = 100;
    %     y_vec = nan(n_pts,1);
    %     acc_prob = nan(n_pts,1);
    %     x_vec = nan(n_pts,n_bvs,n_filts);
    %     for i = 1:n_pts
    %         i
    %         new_x = propose_new_x(old_x,radius,n_bvs,n_filts);
    %         new_y = STCBF_LLinternal_nonlpsc_constmin(new_x(:),Robs,kern_output,cur_mod);
    %         delta_y = new_y - old_y;
    %         acc_prob(i) = min(1,exp(-delta_y/T(i)));
    %         if rand < acc_prob(i)
    %             y_vec(i) = new_y;
    %             x_vec(i,:,:) = new_x;
    %         else
    %             y_vec(i) = y_vec(i-1);
    %             x_vec(i,:,:) = x_vec(i-1,:,:);
    %         end
    %         old_y = y_vec(i);
    %         old_x = squeeze(x_vec(i,:,:));
    %     end
    
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
    mod_signs = ones(n_bvs,1);
    STCcf_0 = eye(n_bvs);
    clear init_nls nltypes
    init_nls{1} = 'lin';
    nltypes{1} = 'lin';
    for i = 2:1+n_exp_dims(nn);
        init_nls{i} = 'pquad';
        nltypes{i} = 'quad';
    end
    for i = (1+n_exp_dims(nn)+1):(1+1+n_exp_dims(nn)+n_sup_dims(nn))
        init_nls{i} = 'nquad';
        nltypes{i} = 'quad';
    end
    glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
    glm_stcb.image_type = '1d';
    glm_stcb.basis = 'pix';
    glm_stcb.spk_nl = 'exp';
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,cur_tr_stimemb);
    glm_stcb.lambdaW = 0; %sparseness on model weights
    [glm_stcb,COMs] = get_filter_coms_1d(glm_stcb);
    
    %determine LL and LP at current filter point
    stc_mod = fitNLw_alt_full(glm_stcb,cur_tr_stimemb,cur_tr_spkbns);
    [ll0, ll0p] = getLLGLM_FULL2d(stc_mod,cur_tr_stimemb,cur_tr_spkbns,'tots');
    stc_mod_LL(nn) = ll0;
    stc_mod_LP(nn) = ll0p;
    stc_mod_xvLL(nn) = getLLGLM_FULL2d(stc_mod,cur_xv_stimemb,cur_xv_spkbns,'none');
    
    %% NOW FIT "STA" model
    mod_signs = 1;
    STCcf_0 = 1;
    clear init_nls nltypes
    init_nls{1} = 'lin';
    nltypes{1} = 'lin';
    glm_stcb = createGLM2d_fullbf(cur_basis(:,1),STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
    glm_stcb.image_type = '1d';
    glm_stcb.basis = 'pix';
    glm_stcb.spk_nl = 'exp';
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,cur_tr_stimemb);
    glm_stcb.lambdaW = 0; %sparseness on model weights
    [glm_stcb,COMs] = get_filter_coms_1d(glm_stcb);
    
    %determine LL and LP at current filter point
    sta_mod = fitNLw_alt_full(glm_stcb,cur_tr_stimemb,cur_tr_spkbns);
    [ll0, ll0p] = getLLGLM_FULL2d(stc_mod,cur_tr_stimemb,cur_tr_spkbns,'tots');
    sta_mod_LL(nn) = ll0;
    sta_mod_LP(nn) = ll0p;
    sta_mod_xvLL(nn) = getLLGLM_FULL2d(sta_mod,cur_xv_stimemb,cur_xv_spkbns,'none');
    
    %% SAVE DATA
    all_xv_params(nn).nparts = nparts;
    all_xv_params(nn).nfold = nfold;
    all_xv_params(nn).cur_tr_parts = cur_tr_parts;
    all_xv_params(nn).cur_xv_parts = cur_xv_parts;
    
    cd ~/James_scripts/stc_sparse_test_figs/
    save rust_allcells_modanal_revised_v2 all_* *_LL *_LP *_xvLL refined_mod
    
    %% CREATE PLOTS
    %     cd ~/James_scripts/stc_sparse_test_figs/modfits_allcells/
    %
    %     f1 = plotfo1d_nopsc(stc_mod,6);
    %     set(f1,'PaperSize',[20 10]);
    %     fname = sprintf('%s-STCmod',uids{nn});
    %     print(fname,'-dpng');close
    %
    %     f1 = plotfo1d_nopsc(rotbv_mod(best_mod),6);
    %     set(f1,'PaperSize',[20 10]);
    %     fname = sprintf('%s-rotmod',uids{nn});
    %     print(fname,'-dpng');close
    %
    %     f1 = plotfo1d_nopsc(fin_glm,6);
    %     set(f1,'PaperSize',[20 10]);
    %     fname = sprintf('%s-finmod',uids{nn});
    %     print(fname,'-dpng');close
    
    
end