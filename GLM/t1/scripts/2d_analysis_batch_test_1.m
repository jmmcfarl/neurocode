%% Load Data

clear all;
addpath('~/James_scripts/GLM')
addpath('~/James_scripts/GLM/2d')
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


%% Compute whitened data, and STA/STCs
cd '/Users/James/James_scripts/GLM/2d/c26_fits/'
tcell = 26;
pids =1:1024;
tsbs      = 1+floor(aselspks{tcell}/dt);

% figure; ecdf(aselspks{tcell})

ncomps = [200:200:2000];
lambda_Lx1 = [0 1e-4 5e-4 1e-3 5e-3 1e-2];
cur_lx1 = 4;
for n = 1:length(ncomps)
    compids   = 1:ncomps(n);
    pix_conv_mat = diag(sqrt(scorevars(compids)))*coefs(:,compids)';
    kern_conv_mat = coefs(:,compids)*diag(1./sqrt(scorevars(compids)));
    
    WX        = scores(:,compids)*diag(1./sqrt(scorevars(compids)));
    % mWX = mean(WX*pix_conv_mat)*kern_conv_mat;
    % WX = WX - repmat(mWX,size(WX,1),1);
    spikebins = tsbs(tsbs>flen & tsbs<(size(WX,1)+flen-1))-flen+1 ;
    WS        = WX(spikebins,:);
    rsta      = mean(WS) - mean(WX);
    
    stvcv = cov(WS);  utvcv = cov(WX);
    [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
    
    npos=10; nneg=6;
    stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
    
    rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
    nposdims = 4; nnegdims = 0;
    posdims = 1:nposdims; negdims = 1:nnegdims;
    % STCbvs = [sta' stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
    STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
    kimages = [rsta',STCbvs]'*diag(sqrt(scorevars(compids)))*coefs(:,compids)';
    % f1 = figure('Name',stype); plotfilterbank(kimages',32,pids)
    
    %add in STA
    STCbvs = [rsta' STCbvs];
    
    [klen,Nstcbvs] = size(STCbvs);
    flen = 6;
    stimlen = size(WX,1);
    kern_output = WX*STCbvs;
    
    %% INITIALIZE MODEL
    addpath('~/James_scripts/GLM/2d')
    nmods = Nstcbvs;
    mod_signs = ones(nmods,1);
    dim_signs = ones(Nstcbvs,1);
    
    unused_stcs = (nmods+1):Nstcbvs;
    %initialize on STC dims
    STCcf_0 = eye(Nstcbvs);
    STCcf_0(:,unused_stcs) = [];
    if nmods > Nstcbvs
        n_extra = nmods-Nstcbvs;
        STCcf_0 = [STCcf_0 randn(Nstcbvs,n_extra)];
    end
    %random initialization
    % n_extra = nmods;
    % STCcf_0 = randn(Nstcbvs,n_extra);
    
    %make sure expansive subunits start out in expansive subspace, etc
    STCcf_0(negdims,mod_signs==1) = 0;
    STCcf_0(posdims,mod_signs==-1) = 0;
    
    %normalize
    for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    
    %initialize model
    defmod.h(1:end-1) = []; %eliminate PSC
    defmod.lnl = 0;
    defmod.lh = 0;
    defmod.lnl2 = 200;
    defmod.lh2 = 0;
    defmod.nlcon = 0;
    defmod.nlmon = 0;
    defmod.locLambda = 0;
    defmod.lambda_dX = 50;
    defmod.lambda_L1x = lambda_Lx1(cur_lx1);
    defmod.lambda_dT = 0;
    defmod.pids = pids;
    nltype = 'uncon';
    init_nls{1} = 'l'; for i = 2:nmods; init_nls{i} = 'pq'; end;
    % basis = 'pix';
    basis = 'white';
    glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltype,init_nls,basis,'test'); %initialize
    [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
    % glm_stcb.basis = 'white';
    glm_stcb.lambdaW = 0; %sparseness on model weights
    
    % figure
    % plot2d_mod(glm_stcb)
    % stc_pix = STCbvs'*pix_conv_mat;
    % figure
    % plotfilterbank(stc_pix',32,pids);
    %     f1 = plotfo1d_nopsc(stc_posneg_mod{jj,ii});
    %         set(f1,'Position',[1300 1000 500 800]);
    
    %%
    full_glm = fitNLHI2d_fullbf(glm_stcb,WX,spikebins,'tots');
    f1 = figure;
    plot2d_mod(full_glm);
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [20 20]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    pname = sprintf('%dwdims_%ddX_%dL1',ncomps(n),defmod.lambda_dX,cur_lx1);
    print('-dpng',pname); close
    
    f1 = figure;
    plot2d_mod(glm_stcb);
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [20 20]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    pname = sprintf('%dwdims_init',ncomps(n));
    print('-dpng',pname); close
    
    fname = sprintf('%dwdims_%ddX_%dL1',n_comps(n),defmod.lambda_dX,cur_lx1);
    save(fname,'glm_stcb','full_glm');
end
%% Fit model
% clear WX scores scorevars coefs
% load StimXmat_z-all
% full_glm = fitNLHI2d_fullbf(glm_stcb,X,spikebins,'tots');

