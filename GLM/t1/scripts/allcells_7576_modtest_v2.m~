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

stcmods = 4;
ncomps = 500;
scorevars(1001:end) = [];
coefs(:,1001:end) = [];
scores(:,1001:end) = [];

% scores = scores(1:20000,:); %test limiting amount of data

%% Compute whitened data, and STA/STCs
cd '/Users/James/James_scripts/GLM/t1/allcell_fits/'
for t = 1:27
% t = 27;
tcell = t;
pids =1:1024;
tsbs      = 1+floor(aselspks{tcell}/dt);

compids   = 1:ncomps;
pix_conv_mat = diag(sqrt(scorevars(compids)))*coefs(:,compids)';
kern_conv_mat = coefs(:,compids)*diag(1./sqrt(scorevars(compids)));

WX        = scores(:,compids)*diag(1./sqrt(scorevars(compids)));
spikebins = tsbs(tsbs>flen & tsbs<(size(WX,1)+flen-1))-flen+1 ;
WS        = WX(spikebins,:);
rsta      = mean(WS) - mean(WX);

stvcv = cov(WS);  utvcv = cov(WX);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);

npos=10; nneg=10;
stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
nposdims = stcmods-1; nnegdims = 0;
posdims = 1:nposdims; negdims = 1:nnegdims;
% STCbvs = [sta' stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
        kimages = [rsta',STCbvs]'*diag(sqrt(scorevars(compids)))*coefs(:,compids)';
f1 = figure('Name',stype); plotfilterbank(kimages',32,pids)

%add in STA
STCbvs = [rsta' STCbvs];
%add in STA
STCbvs = [rsta' STCbvs];

[klen,Nstcbvs] = size(STCbvs);
flen = 6;
stimlen = size(WX,1);

%% INITIALIZE MODEL
nmods = Nstcbvs;
mod_signs = ones(nmods,1);
dim_signs = ones(Nstcbvs,1);

unused_stcs = (nmods+1):Nstcbvs;
%initialize on STC dims
STCcf_0 = eye(Nstcbvs);
STCcf_0(:,unused_stcs) = [];
% if nmods > Nstcbvs
%     n_extra = nmods-Nstcbvs;
%     STCcf_0 = [STCcf_0 randn(Nstcbvs,n_extra)];
% end
%random initialization
% n_extra = nmods;
% STCcf_0 = randn(Nstcbvs,n_extra);

% STCcf_0 = [STCcf_0 -STCcf_0];

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
defmod.lambda_dX = 100;
defmod.lambda_L1x = 10;
defmod.lambda_dT = 10;
defmod.pids = pids;

%define NL initializations: "lin, threshlin, pquad, nquad"
clear init_nls nltypes
init_nls{1} = 'lin'; for i = 2:nmods; init_nls{i} = 'pquad'; end;
%define NL types: "uncon, lin, threshlin, quad"
nltypes{1} = 'lin'; for i = 2:nmods; nltypes{i} = 'quad'; end;
% %define NL initializations: "lin, threshlin, pquad, nquad"
% for i = 1:nmods; init_nls{i} = 'threshlin'; end;
% %define NL types: "uncon, lin, threshlin, quad"
% for i = 1:nmods; nltypes{i} = 'threshlin'; end;

% basis = 'pix';
basis = 'white';
glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('Cell%d',t)); %initialize
[glm_stcb,sub_medians] = median_sub_filters(glm_stcb);

[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
glm_stcb.lambdaW = 10; %sparseness on model weights
glm_stcb.image_type = '2d';

plot2d_mod(glm_stcb,0,1:10);colormap(jet);
plot2d_mod(glm_stcb,0,11:20);colormap(jet);
% set(f1,'PaperUnits','centimeters');
% set(f1, 'PaperSize', [30 6*nmods]);
% set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
% pname = sprintf('cell%d_500wdims_%dmods_init',t,nmods);
% print('-dpng',pname); close

%%
full_glm = fitNLHI2d_fullbf(glm_stcb,WX,spikebins,'tots');
f2 = plot2d_mod(full_glm,0);
% set(f2,'PaperUnits','centimeters');
% set(f2, 'PaperSize', [30 6*nmods]);
% set(f2,'PaperPosition',[0,0,(get(f2,'PaperSize'))])
% pname = sprintf('cell%d_500wdims_%dmods_%ddX_%dL1',t,nmods,defmod.lambda_dX,defmod.lambda_L1x);
% print('-dpng',pname); close
% 
% 
% fname = sprintf('cell%d_500wdims_%dmods_%ddX_%dL1',t,nmods,defmod.lambda_dX,defmod.lambda_L1x);
% save(fname,'glm_stcb','full_glm');
end
