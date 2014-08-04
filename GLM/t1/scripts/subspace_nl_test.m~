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
ncomps = 500;
compids   = 1:ncomps;
scorevars = var(scores(tr_inds{xv},:)); %variance of PCA dimensions within the training set
whitemat = diag(1./sqrt(scorevars(compids))); %whitening (rescaling) transformation from training set
unwhitemat = diag(sqrt(scorevars(compids)));
WX = scores(tr_inds{xv},compids)*whitemat; %whitened training data
pix_conv_mat = unwhitemat*coefs(:,compids)';
kern_conv_mat = coefs(:,compids)*whitemat;
X_avg = mean(WX); %average whitened training stim
NT_x = size(WX,1); %number of training samples

%% Cross-validation
WX_xv  = scores(xv_inds{xv},compids)*whitemat;

%% Compute whitened data, and STA/STCs
cd '/Users/James/James_scripts/GLM/t1/allcells_fits_mlstc_v2/'  
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
    
    %project STA out of STE
    WXp = WX - (WX*rsta')*rsta/norm(rsta);
    WSp = WXp(spikebins,:);
    
    %Compute STC
    stvcv = cov(WSp); %covaraince of STE
    %[evecs,evals] = eig(stvcv-eye(ncomps)); %assume exact whitening...
    [evecs,evals] = eig(stvcv - cov(WXp));
    evs{t} = diag(evals);
    
    npos=8; nneg=4; %retain a handful of positive and negative dimensions to check
    stcs  = evecs(:,[1:nneg,length(evs{t})-npos+1:end]); stcs  = stcs(:,end:-1:1);
    rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
    STCbvs = [stcs(:,1:npos) rstcs(:,1:nneg)]; 
%     STCbvs = (STCbvs'*whitemat.^2)'; %correct for unequal variance bias in STC dims
    STCbvs = (STCbvs'*whitemat)'; %correct for unequal variance bias in STC dims
    STCbvs = [rsta' STCbvs];
        
    kimages = STCbvs'*pix_conv_mat;
    f1 = figure('Name',stype); plotfilterbank(kimages',32,pids)
      
    %% TRY FITTING SEQUENCE OF STC MODELS WTIH INCREASING DIMENSIONALITY
    %initialize model
    defmod.h(1:end-1) = []; %eliminate PSC
    defmod.lnl = 0;
    defmod.lh = 0;
    defmod.lnl2 = 0;
    defmod.lh2 = 0;
    defmod.nlcon = 0;
    defmod.nlmon = 0;
    defmod.lambda_dX = 500;
    defmod.lambda_L1x = 4;
    defmod.lambda_dT = 50;
    defmod.pids = 1:fsdim;
    defmod.SDIM = SDIM;
    defmod.fsdim = fsdim;
    defmod.locLambda = 0;
    
    basis = 'white';

    %first fit sta model
    cur_ndims = 4;
    cur_basis = STCbvs(:,1:cur_ndims); %just use STA
    STCcf_0 = eye(cur_ndims);
    %these are expansive quadratic components
    for jj = 2:cur_ndims
        init_nls{jj} = 'pquad';
        nltypes{jj} = 'quad';
    end
    init_nls{1} = 'lin';
    nltypes{1} = 'lin';
    STCcf_0 = eye(cur_ndims);
    glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
    glm_stcb.image_type = '2d';
    stc_glm{cur_ndims} = fitNLHI2d_fullbf(glm_stcb,WX,spikebins,'tots',4);
    stc_glm{cur_ndims}.xvLL = getLLGLM_FULL2d(stc_glm{cur_ndims},WX_xv,xv_spikebins,'none');
                    
    f1 = plot2d_mod(stc_glm{cur_ndims});

    %%
    used_mod = stc_glm{cur_ndims};
    basis_vecs = get_k_mat(used_mod);
    n_bvs = size(basis_vecs,2);
    
    nmods = 6;
    
    %initialize model
    defmod.h(1:end-1) = []; %eliminate PSC
    defmod.lnl = 0;
    defmod.lh = 0;
    defmod.lnl2 = 100;
    defmod.lh2 = 0;
    defmod.nlcon = 0;
    defmod.nlmon = 0;
    defmod.locLambda = 0;
    defmod.lambda_dX = 200; %350
    defmod.lambda_L1x = 0; %40
    defmod.lambda_dT = 10;
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
        
    f1 = plot2d_mod(rotated_mod{t})

    %% Characterize NL fits
 basis_stim_proj = WX*basis_vecs;
 used_dims = [3 4];
minxy = [-4 -4];
maxxy = [4 4];
bandwidth = [0.3 0.3];
[~,density_ov,X2,Y2] = kde2d(basis_stim_proj(:,used_dims),2^8,minxy,maxxy,bandwidth);
[~,density_spk,X,Y] = kde2d(basis_stim_proj(spikebins,used_dims),2^8,minxy,maxxy,bandwidth);

eps = 1e-3;
cond_dens = density_spk./density_ov;
cond_dens(density_ov < eps) = eps;
outside = find(X.^2 + Y.^2 > 9);
cond_dens(outside) = nan;

rot_vecs = get_STCcf_mat(rotated_mod{t});
for i = 1:nmods
rot_vecs(:,i) = rot_vecs(:,i)/norm(rot_vecs(:,i));
end

bv_1 = X2(:);
bv_2 = Y2(:);
gen_fun = zeros(size(bv_1));
for i = 1:nmods
    input = bv_1*rot_vecs(used_dims(1),i)+bv_2*rot_vecs(used_dims(2),i);
    gen_fun = gen_fun + rotated_mod{t}.mods(i).w*nlin_proc_stim(input,rotated_mod{t}.mods(i).nly,...
        rotated_mod{t}.mods(i).nlx);
end
gen_fun = gen_fun + rotated_mod{t}.const;
mod_rate = log(1+exp(gen_fun));
mod_rate = reshape(mod_rate,size(X2));
mod_rate(outside) = nan;
gen_fun = reshape(gen_fun,size(X2));
gen_fun(outside) = nan;

figure
contourf(X2,Y2,cond_dens,30)
line([0 rot_vecs(used_dims(1),1)]*3,[0 rot_vecs(used_dims(2),1)]*3,'color','r','linewidth',2)
line([0 rot_vecs(used_dims(1),2)]*3,[0 rot_vecs(used_dims(2),2)]*3,'color','k','linewidth',2)
line([0 rot_vecs(used_dims(1),3)]*3,[0 rot_vecs(used_dims(2),3)]*3,'color','g','linewidth',2)
line([0 rot_vecs(used_dims(1),4)]*3,[0 rot_vecs(used_dims(2),4)]*3,'color','w','linewidth',2)

figure
contourf(X2,Y2,mod_rate,30)
line([0 rot_vecs(used_dims(1),1)]*3,[0 rot_vecs(used_dims(2),1)]*3,'color','r','linewidth',2)
line([0 rot_vecs(used_dims(1),2)]*3,[0 rot_vecs(used_dims(2),2)]*3,'color','k','linewidth',2)
line([0 rot_vecs(used_dims(1),3)]*3,[0 rot_vecs(used_dims(2),3)]*3,'color','g','linewidth',2)
line([0 rot_vecs(used_dims(1),4)]*3,[0 rot_vecs(used_dims(2),4)]*3,'color','w','linewidth',2)

% figure
% contourf(X2,Y2,gen_fun,30)
% line([0 rot_vecs(1,1)]*3,[0 rot_vecs(2,1)]*3,'color','r','linewidth',2)
% line([0 rot_vecs(1,2)]*3,[0 rot_vecs(2,2)]*3,'color','k','linewidth',2)

% figure
% contourf(X2,Y2,density_spk,30)
% 
% figure
% contourf(X2,Y2,log(density_ov),30)

end
%%

