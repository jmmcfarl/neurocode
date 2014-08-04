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

cd ~/Data/blanche/matlabdata/
load spks7576-all.mat;

load StimXmat_z-all
stype='all';

pids =1:1024;

%% create XV data
[stimlen,k_len] = size(X);
% nparts = 100;
% nfold = 10;
nparts = 100;
nfold = 5;
partlen = floor(stimlen/nparts);
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
compids   = 1:k_len;

t = 13;
tcell = t;
tsbs      = 1+floor(aselspks{tcell}/dt);

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
defmod.lambda_dX = 10000;
defmod.lambda_L1x = 0;
defmod.lambda_dT = 10000;

defmod.fsdim = fsdim;
defmod.pids = 1:fsdim;
defmod.SDIM = sdim;
basis = 'pix';
cur_ndims = 2;

init_signs = ones(cur_ndims,1);
init_kerns = randn(k_len,cur_ndims);

for xv = 1
    X_tr = X(tr_inds{xv},compids); %whitened training data
    clear X
    %spike bins during training data
    tr_spbs = find(ismember(tr_inds{xv},tsbs));
    spikebins = tr_spbs(tr_spbs>flen & tr_spbs<(size(X_tr,1)+flen-1))-flen+1 ;
    
    [glm_stcb] = createGLM_quad(init_kerns,init_signs,defmod,basis);
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,X_tr);
    g_mat = X_tr*get_k_mat(glm_stcb);
    temp = fitWeights_full(glm_stcb,g_mat,spikebins,1);
    glm_stcb.const = temp.const;
    
    % glm_stcb.spk_nl = 'logexp';
        glm_stcb.image_type = '2d';
    quad_glm(xv) = fitGLM_lexp(glm_stcb,X_tr,spikebins,'tots',50,1e-4,1e-6);
end

%%
defmod.lambda_dX = 10000;
defmod.lambda_L1x = 20;
defmod.lambda_dT = 10000;
for i = 1:length(quad_glm(xv).mods)
   quad_glm(xv).mods(i).lambda_dX = defmod.lambda_dX; 
   quad_glm(xv).mods(i).lambda_L1x = defmod.lambda_L1x; 
   quad_glm(xv).mods(i).lambda_dT = defmod.lambda_dT; 
end
quad_glm(xv) = fitGLM_lexp(quad_glm(xv),X_tr,spikebins,'tots',50,1e-4,1e-6);

%%
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
% defmod.lambda_dX = 3000;
% defmod.lambda_L1x = 25;
% defmod.lambda_dT = 3000;
% defmod.lambda_dX = 10000;
% defmod.lambda_L1x = 5;
% defmod.lambda_dT = 10000;
defmod.lambda_dX = 0;
defmod.lambda_L1x = 1e-3;
defmod.lambda_dT = 0;

defmod.fsdim = fsdim;
defmod.pids = 1:fsdim;
defmod.SDIM = sdim;
basis = 'white';

% clear stc_glm init_nls nltypes
%first fit sta model

for xv = 1
    WX = scores(tr_inds{xv},compids); %whitened training data
    WX_xv = scores(xv_inds{xv},compids);
    pix_conv_mat = coefs(:,compids)';
    kern_conv_mat = coefs(:,compids);

    %spike bins during training data
    tr_spbs = find(ismember(tr_inds{xv},tsbs));
    spikebins = tr_spbs(tr_spbs>flen & tr_spbs<(size(WX,1)+flen-1))-flen+1 ;
    
    %spike bins for XV set
    xv_spbs = find(ismember(xv_inds{xv},tsbs));
    xv_spikebins = xv_spbs(xv_spbs>flen & xv_spbs<(size(WX_xv,1)+flen-1))-flen+1 ;
    
%     cur_ndims = 2;
%     cur_basis = STCbvs(:,1:cur_ndims); %just use STA
%     STCcf_0 = eye(cur_ndims);
%     init_signs = ones(cur_ndims,1);
%     [glm_stcb] = createGLM_quad(cur_basis*STCcf_0,init_signs,defmod,basis,pix_conv_mat,kern_conv_mat)
%     
%     [glm_stcb,norm_vals1] = normalizeRFs_full(glm_stcb,WX);
%     g_mat = WX*get_k_mat(glm_stcb);
%     temp = fitWeights_full(glm_stcb,g_mat,spikebins,0);
%     glm_stcb.const = temp.const;
%     
%     glm_stcb.spk_nl = 'logexp';
%     glm_stcb.image_type = '2d';
%     %     glm_stcb.mods(1).dX = 5e6;
%     %     glm_stcb.mods(1).L1x = 50000;
%     %     glm_stcb.mods(1).dT = 2e6;
%     
%     temp = fitGLM_lexp(glm_stcb,WX,spikebins,'tots',[],1e-4,1e-6);
%     [temp,norm_vals] = normalizeRFs_full(temp,WX);
%     g_mat = WX*get_k_mat(temp);
%     temp = fitWeights_full(temp,g_mat,spikebins,0);
%     for i = 1:length(temp.mods)
%         temp.mods(i).lambda_dX = temp.mods(i).lambda_dX*norm_vals(i)^2;
%         temp.mods(i).lambda_L1x = temp.mods(i).lambda_L1x*norm_vals(i);
%         temp.mods(i).lambda_dT = temp.mods(i).lambda_dT*norm_vals(i)^2;
%     end
%     temp = fitGLM_lexp(temp,WX,spikebins,'tots',[],1e-4,1e-6);
%     
%     cur_ndims = 3;
%     glm_stcb = add_null_filter(temp,'rand',1,'quad',2,0);
  
    cur_ndims = 3;
    cur_basis = STCbvs(:,1:cur_ndims); %just use STA
    STCcf_0 = eye(cur_ndims);
    init_signs = ones(cur_ndims,1);
    [glm_stcb] = createGLM_quad(cur_basis*STCcf_0,init_signs,defmod,basis,pix_conv_mat,kern_conv_mat)
    glm_stcb.spk_nl = 'logexp';
    glm_stcb.image_type = '2d';
    g_mat = WX*get_k_mat(glm_stcb);
    temp = fitWeights_full(glm_stcb,g_mat,spikebins,0);
    glm_stcb.const = temp.const;
    temp2 = fitGLM_lexp(glm_stcb,WX,spikebins,'tots',[],1e-4,1e-6);
    [temp2,norm_vals] = normalizeRFs_full(temp2,WX);
    g_mat = WX*get_k_mat(temp2);
    temp2 = fitWeights_full(temp2,g_mat,spikebins,0);
    for i = 1:length(temp2.mods)
        temp2.mods(i).lambda_dX = temp2.mods(i).lambda_dX*norm_vals(i)^2;
        temp2.mods(i).lambda_L1x = temp2.mods(i).lambda_L1x*norm_vals(i);
        temp2.mods(i).lambda_dT = temp2.mods(i).lambda_dT*norm_vals(i)^2;
    end
    quad_glm(xv) = fitGLM_lexp(temp2,WX,spikebins,'tots',[],1e-4,1e-6);
    quad_xvLL(xv) = getLLGLM_FULL2d(quad_glm(xv),WX_xv,xv_spikebins,'tots') %determine XVLL
end

%%
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
defmod.lambda_dX = 0;
defmod.lambda_L1x = 0;
defmod.lambda_dT = 0;

defmod.fsdim = length(compids);
defmod.pids = 1:length(compids);
defmod.SDIM = length(compids);
basis = 'pix';

% clear stc_glm init_nls nltypes
%first fit sta model

for xv = 1
    
%     WX = scores(tr_inds{xv},compids)*whitemat; %whitened training data
    WX = scores(tr_inds{xv},compids); %whitened training data


    %spike bins during training data
    tr_spbs = find(ismember(tr_inds{xv},tsbs));
    spikebins = tr_spbs(tr_spbs>flen & tr_spbs<(size(WX,1)+flen-1))-flen+1 ;
      
    cur_ndims = 3;
    cur_basis = STCbvs(:,1:cur_ndims); %just use STA
    STCcf_0 = eye(cur_ndims);
    init_signs = ones(cur_ndims,1);
    [glm_stcb] = createGLM_quad(cur_basis*STCcf_0,init_signs,defmod,basis);
    glm_stcb.spk_nl = 'logexp';
    glm_stcb.image_type = '1d';
    g_mat = WX*get_k_mat(glm_stcb);
    temp = fitWeights_full(glm_stcb,g_mat,spikebins,0);
    glm_stcb.const = temp.const;
    temp2 = fitGLM_lexp(glm_stcb,WX,spikebins,'tots',[],1e-6,1e-8);
 end
k_mat = get_k_mat(temp2);
k_matw = k_mat'*whitemat;
% images = k_matw*pix_conv_mat;
images = k_mat'*pix_conv_mat;


%%

for xv = 1:10
    WX = scores(tr_inds{xv},compids)*whitemat; %whitened training data
    WX_xv = scores(xv_inds{xv},compids)*whitemat;
    %spike bins during training data
    tr_spbs = find(ismember(tr_inds{xv},tsbs));
    spikebins = tr_spbs(tr_spbs>flen & tr_spbs<(size(WX,1)+flen-1))-flen+1 ;
    
    %spike bins for XV set
    xv_spbs = find(ismember(xv_inds{xv},tsbs));
    xv_spikebins = xv_spbs(xv_spbs>flen & xv_spbs<(size(WX_xv,1)+flen-1))-flen+1 ;

    avg_rate(xv) = length(spikebins)/size(WX,1);
    
    %compute binned spike vector
    Robs = zeros(1,size(WX,1));
    ftable = tabulate(spikebins);
    Robs(ftable(:,1)) = ftable(:,2);
    
    null_xvLL(xv) = -sum(bsxfun(@times,Robs,log(avg_rate(xv))) - avg_rate(xv))/length(spikebins);
    
    X_avg = mean(WX);
    WS        = WX(spikebins,:);
    rsta      = mean(WS) - X_avg;
    
    cur_g = WX*rsta';
    n_bins = 15;
    cur_bin_edges = prctile(cur_g,linspace(5,95,n_bins+1));
    cur_bin_centers = 0.5*cur_bin_edges(1:end-1)+0.5*cur_bin_edges(2:end);
    hist_p = zeros(n_bins,1);
    for i = 1:n_bins
        cur_set = find(cur_g >= cur_bin_edges(i) & cur_g < cur_bin_edges(i+1));
        hist_p(i) = mean(Robs(cur_set));
        n_pts(i) = length(cur_set);
    end
    pred_rate = nlin_proc_stim( cur_g, hist_p, cur_bin_centers);
    ln_xvLL(xv) = -sum(Robs'.*log(pred_rate) - pred_rate)/length(spikebins);
    
    %project STA out of STE
    proj_mat = rsta'*inv(rsta*rsta')*rsta;
    WXp = WX - WX*proj_mat;
    WSp = WXp(spikebins,:);
    
    %Compute STC
    stvcv = cov(WSp); %covaraince of STE
    %[evecs,evals] = eig(stvcv-eye(ncomps)); %assume exact whitening...
    [evecs,evals] = eig(stvcv-cov(WXp));
    evs{t} = diag(evals);
    
    npos=2; nneg=0; %retain a handful of positive and negative dimensions to check
    stcs  = evecs(:,[1:nneg,length(evs{t})-npos+1:end]); stcs  = stcs(:,end:-1:1);
    rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
    STCbvs = [stcs(:,1:npos) rstcs(:,1:nneg)];
    %     STCbvs = (STCbvs); %correct for unequal variance bias in STC dims
    STCbvs = [rsta' STCbvs];
    STCbvs = (STCbvs'*whitemat)'; %correct for unequal variance bias in STC dims
  
    all_kimages{xv} = STCbvs'*pix_conv_mat;
    
%     plotfilterbank(kimages',32,pids)
    
end


