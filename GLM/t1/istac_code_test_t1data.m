clear all;
addpath('~/James_scripts/GLM')
addpath('~/James_scripts/GLM/2d')
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/Timm/tim1/functions/')
addpath('~/James_scripts/code_iSTAC/')

cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

% load spks7576-ns.mat; load PCScores-ns.mat; stype='ns'
% load spks7576-pn.mat; load PCScores-pn.mat; stype='pn'
% load spks7576-ps.mat; load PCScores-ps.mat; stype='ps'

cd ~/Data/blanche/matlabdata/
load spks7576-all.mat; load PCScores-all.mat; stype='all'


%%
tcell = 26;
pids = 1:1024;
tsbs      = 1+floor(aselspks{tcell}/dt);

% figure; ecdf(aselspks{tcell})

ncomps    = 600; compids   = 1:ncomps;
WX        = scores(:,compids)*diag(1./sqrt(scorevars(compids)));
spikebins = tsbs(tsbs>flen & tsbs<(size(WX,1)+flen-1))-flen+1 ;

WS        = WX(spikebins,:);
sta      = mean(WS)';
ov_mu = mean(WX)';

stvcv = cov(WS);  
utvcv = cov(WX);
[evecs,evals] = eig(stvcv-utvcv);
evs   = diag(evals);
npos=10; nneg=10;
stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

%%
stimlen = size(WX,1);
Robs = zeros(1,stimlen);
ftable = tabulate(spikebins);
Robs(ftable(:,1)) = ftable(:,2);

ndims = 10;  % (Only need 2, but compute 10 for demonstration purposes)
[vecs, vals, DD] = compiSTAC(sta, stvcv', ov_mu, utvcv', ndims);

%%
pix_conv_mat = diag(sqrt(scorevars(compids)))*coefs(:,compids)';

stap = sta'*pix_conv_mat;
stcs_p = stcs'*pix_conv_mat;

vecs_p = vecs'*pix_conv_mat;

