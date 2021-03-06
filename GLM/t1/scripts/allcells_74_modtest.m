%% Load Data

clear all;
addpath('~/James_scripts/GLM')
addpath('~/James_scripts/GLM/t1/')
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/Timm/tim1/functions/')

cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

%load z-scored PC data
cd ~/Data/blanche/rec_74/matlabdata/
load ./spks74-all.mat;
load ./PCScores_z-all.mat;
stype='all';
load spksegsRec74.mat;

stcmods = 1;
nmods = 1;
ncomps = 500;
scorevars(1001:end) = [];
coefs(:,1001:end) = [];
scores(:,1001:end) = [];

% scores = scores(1:15000,:);
%% Compute whitened data, and STA/STCs
cd '/Users/James/James_scripts/GLM/t1/rec74_modfits/'
for t = 1:25
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
    
    % stvcv = cov(WS);  utvcv = cov(WX);
    % [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
    
    % npos=10; nneg=10;
    % stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
    
    % rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
    % nposdims = stcmods-1; nnegdims = 0;
    % posdims = 1:nposdims; negdims = 1:nnegdims;
    % % STCbvs = [sta' stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
    % STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
    %         kimages = [rsta',STCbvs]'*diag(sqrt(scorevars(compids)))*coefs(:,compids)';
    
    f1 = figure('Name',stype); plotfilterbank(rsta',32,pids)
    set(f1,'PaperUnits','centimeters');
    set(f1, 'PaperSize', [30 6*nmods]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))])
    pname = sprintf('cell%d_500wdims_%dmods_pink',t,ncomps,nmods);
    print('-dpng',pname); close
    
end
