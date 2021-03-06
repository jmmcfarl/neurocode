clear all; flen =6; SDIM = 14;

%generate spatial pink noise data
DIM = [SDIM SDIM];
BETA = -2;
NT = 60000;
dt = 0.02;
surr_stim = zeros(NT,196);
u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
u = repmat(u,1,DIM(2));
v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]/DIM(2);
v = repmat(v,DIM(1),1);
S_f = (u.^2 + v.^2).^(BETA/2);
S_f(S_f==inf) = 0;
for t = 1:NT
    phi = rand(DIM);
    x = ifft2(S_f.^0.5 .* (cos(2*pi*phi)+i*sin(2*pi*phi)));
    x = real(x);
    surr_stim(t,:) = x(:);
end

X = makeStimRows(surr_stim,flen,1);

%%
clear all; flen =6; SDIM = 25;

cd ~/Data/blanche/rec_75/matlabdata/
load dstimpsRec75.mat;
stf75=load('stimfiles75.mat');
cd ~/Data/blanche/rec_76/matlabdata/
load dstimpsRec76.mat;
stf76=load('stimfileNames76.mat');

allstims = [dstimps75;dstimps76];
cnames   = [stf75.stimfiles,stf76.stimfiles]; cellfun(@disp,cnames)
nconds   = length(cnames);

keys = {'af','wn','pn','pt','ps','ns'}
rids = cellfun(@(tkey)...
    find(strcmp(cellfun(@(x)x(1:2),cnames,'UniformOutput',0),tkey)),...
    keys,'UniformOutput',0);

taf        = rids{1}; %cellfun(@disp,cnames(tpn))
twn        = rids{2};
tpn        = rids{3}; %cellfun(@disp,cnames(tpn))
tpt        = rids{4};
tps        = rids{5}; %cellfun(@disp,cnames(tps))
tns        = rids{6}; %cellfun(@disp,cnames(tns))

is_nat_contrast = ~cellfun(@(x) (x(end)=='5'),cnames);
is_nat_mean = ~cellfun(@(x) (x(18)=='1'),cnames);
stim_contrast = nan(size(is_nat_contrast));
stim_contrast(~is_nat_contrast) = cellfun(@(x) str2num(x(end-2:end)),cnames(~is_nat_contrast));
stim_mean = nan(size(is_nat_mean));
stim_mean(~is_nat_mean) = cellfun(@(x) str2num(x(18:20)),cnames(~is_nat_mean));

% DIFFERENT MEAN GROUPS (POOLED ACROSS NS, PN, AND PS)
used_conds = tns(~isnan(stim_mean(tns)));

%define subpatch of image
[XX,YY] = meshgrid(1:32,1:32);
Xr = XX(:); Yr = YY(:);
used_pixs = find(Xr <= 25 & Yr <= 25); 

selstim   = [];
cur_means = [];
cur_vars = [];
for icond=1:length(used_conds)
    tcond = used_conds(icond)
    cur_stim = allstims{tcond};
    cur_stim = cur_stim(:,used_pixs);

    cur_stim = bsxfun(@minus,cur_stim,mean(cur_stim,2));
    cur_stim = bsxfun(@rdivide,cur_stim,std(cur_stim,[],2));
    selstim =[selstim;cur_stim];
    cur_means = [cur_means; mean(cur_stim(:))];
    cur_vars = [cur_vars; var(cur_stim(:))];
end

addpath('~/James_scripts/')
NT   = size(selstim,1); SDIM  = size(selstim,2); NeK   = flen*SDIM;
X = makeStimRows(selstim,flen,1);
% clear allstims
% 
X = bsxfun(@minus,X,mean(X));

%% generate surrogate filters
SDIM = 25; flen = 6; dt = 0.02;
x = repmat(1:SDIM,[flen 1 SDIM]);
y = permute(x,[1 3 2]);

x0 = permute(repmat(8,[SDIM SDIM flen]),[3 1 2]);
y0 = permute(repmat(8,[SDIM SDIM flen]),[3 1 2]);
sigmax = permute(repmat(3,[SDIM SDIM flen]),[3 1 2]);
sigmay = permute(repmat(4,[SDIM SDIM flen]),[3 1 2]);
lambda = permute(repmat(10,[SDIM SDIM flen]),[3 1 2]);
theta = permute(repmat(pi/4,[SDIM SDIM flen]),[3 1 2]);
b = repmat(linspace(2,0,flen)',[1 SDIM SDIM]);
psi = repmat(linspace(0,pi/2,flen)',[1 SDIM SDIM]);

psi2 = psi+pi/2;
x02 = permute(repmat(15,[SDIM SDIM flen]),[3 1 2]);
y02 = permute(repmat(15,[SDIM SDIM flen]),[3 1 2]);

xp = (x-x0).*cos(theta)+(y-y0).*sin(theta);
yp = -(x-x0).*sin(theta)+(y-y0).*cos(theta);
xp2 = (x-x02).*cos(theta)+(y-y02).*sin(theta);
yp2 = -(x-x02).*sin(theta)+(y-y02).*cos(theta);

filt1 = b.*(exp(-(xp.^2./2./sigmax.^2 + yp.^2./2./sigmay.^2)) .* (cos(2*pi*xp./lambda+psi)));
filt2 = b.*(exp(-(xp2.^2./2./sigmax.^2 + yp2.^2./2./sigmay.^2)) .* (cos(2*pi*xp2./lambda+psi2)));
for i = 1:flen
    cur = filt1(i,:,:);
    filt1(i,:,:) = filt1(i,:,:) - mean(cur(:));
    cur = filt2(i,:,:);
    filt2(i,:,:) = filt2(i,:,:) - mean(cur(:));
end
filt1 = filt1/norm(filt1(:));
filt2 = filt2/norm(filt2(:));

% x0 = permute(repmat(22,[SDIM SDIM flen]),[3 1 2]);
% y0 = permute(repmat(22,[SDIM SDIM flen]),[3 1 2]);
% sigmax = permute(repmat(3,[SDIM SDIM flen]),[3 1 2]);
% sigmay = permute(repmat(2,[SDIM SDIM flen]),[3 1 2]);
% lambda = permute(repmat(10,[SDIM SDIM flen]),[3 1 2]);
% theta = permute(repmat(pi/4,[SDIM SDIM flen]),[3 1 2]);
% b = repmat(linspace(2,0,flen)',[1 SDIM SDIM]);
% psi = repmat(linspace(0,pi,flen)',[1 SDIM SDIM]);
%
% psi2 = psi;
% x02 = permute(repmat(22,[SDIM SDIM flen]),[3 1 2]);
% y02 = permute(repmat(11,[SDIM SDIM flen]),[3 1 2]);
%
% xp = (x-x0).*cos(theta)+(y-y0).*sin(theta);
% yp = -(x-x0).*sin(theta)+(y-y0).*cos(theta);
% xp2 = (x-x02).*cos(theta)+(y-y02).*sin(theta);
% yp2 = -(x-x02).*sin(theta)+(y-y02).*cos(theta);
%
% filt3 = b.*(exp(-(xp.^2./2./sigmax.^2 + yp.^2./2./sigmay.^2)) .* (cos(2*pi*xp./lambda+psi)));
% filt3 = filt3/norm(filt3(:));
%
% filt4 = b.*(exp(-(xp2.^2./2./sigmax.^2 + yp2.^2./2./sigmay.^2)) .* (cos(2*pi*xp2./lambda+psi2)));
% filt4 = filt4/norm(filt4(:));

figure
for i = 1:6
    subplot(2,1,1)
imagesc(squeeze(filt1(i,:,:)))
subplot(2,1,2)
imagesc(squeeze(filt2(i,:,:)))
pause
clf
end

%%
filt_weights = [5 5];
filt_long = [filt_weights(1)*filt1(:) filt_weights(2)*filt2(:)];
filt_proj_stim = X*filt_long;

% filt_white = filt_long'*bsxfun(@rdivide,kern_conv_mat,sqrt(sum(kern_conv_mat.^2)));
% filt_proj_stim = scores(:,compids)*filt_white';

%%
filt_stim1 = filt_proj_stim(:,1);
filt_stim1(filt_stim1 < 0) = 0;

filt_stim2 = filt_proj_stim(:,2);
filt_stim2(filt_stim2 < 0) = 0;

% filt_stim2 = filt_proj_stim(:,2).^2;

% filt_stim3 = filt_proj_stim(:,3);
% filt_stim3(filt_stim3 < 0) = 0;
%
% filt_stim4 = filt_proj_stim(:,4);
% filt_stim4(filt_stim4 < 0) = 0;

% g = filt_stim1 + filt_stim2 + filt_stim3 + filt_stim4;
g = filt_stim1 + filt_stim2;
g = 5*zscore(g);

target_rate = 20; %in Hz
max_rate = 2000;
target_pspike = target_rate*dt;
% K0 = [1 0.4];
K0 = [0.6 0.2];
% Kfit = fmincon(@(K) abs(mean(logexp(g,K)) - target_pspike),K0,[],[],[],[],[-10 .01],[10 10]);
Kfit = K0;
p_spike = logexp(g,Kfit);
p_spike(p_spike > max_rate*dt) = max_rate*dt;

% figure
% [y,x] = ksdensity(g);
% plot(x,y); hold on
% yl = ylim();
% tx = -10:.02:10;
% plot(tx,logexp(tx,Kfit),'k')
% ylim(yl)

%%
spikes = poissrnd(p_spike);
rbins = (find(spikes>0.5));
nsp = spikes(rbins);
spk_vals = unique(spikes); spk_vals(spk_vals==0) = [];
spikebins = [];
for i = 1:length(spk_vals)
    cur_set = find(spikes == spk_vals(i));
    spikebins = [spikebins; repmat(cur_set(:),spk_vals(i),1)];
end
fprintf('Nspks: %d\n',length(spikebins));
% or_spikebins = spikebins;
% spikebins = spikebins - flen + 1;
%%
nneg = 5;
npos = 5;   

[coeff,score,latent] = princomp(X);

%%
ncomps = 600;
compids   = 1:ncomps;
pix_conv_mat = diag(sqrt(latent(compids)))*coeff(:,compids)';
kern_conv_mat = coeff(:,compids)*diag(1./sqrt(latent(compids)));

WX        = X*kern_conv_mat;
WS        = WX(spikebins,:);
rsta      = mean(WS) - mean(WX);

stvcv = cov(WS);  utvcv = cov(WX);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);

stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

rstcs = fliplr(stcs); %reversed STC kernels (suppressive first)
nposdims = 5; nnegdims = 5;
posdims = 1:nposdims-1; negdims = 1:nnegdims;
STCbvs = [stcs(:,posdims) rstcs(:,negdims)]; %use only expansive subspace
kimages = [rsta',STCbvs]'*diag(1./sqrt(latent(compids)))*coeff(:,compids)';
pids = 1:(SDIM^2);
figure; 
plotfilterbank(kimages',SDIM,pids)


%%
npos=5; nneg=5;
ncomps = 600; compids = 1:ncomps;
% scalmat = diag(1./sqrt(latent(1:ncomps)));
scalmat = diag(ones(1,ncomps));

WX = score(:,compids)*scalmat;

WS        = WX(spikebins,:);
sta      = mean(WS) - mean(WX);
stvcv = cov(WS);  utvcv = cov(WX);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

% kimages = [sta',stcs]'*coeff(:,compids)';

% pix_conv_mat = diag(sqrt(latent(compids)))*coeff(:,compids)';
% kern_conv_mat = coeff(:,compids)*diag(1./sqrt(latent(compids)));
pix_conv_mat = coeff(:,compids)';
kern_conv_mat = coeff(:,compids);
kimages = [sta',stcs]'*pix_conv_mat;
figure
pids = 1:(SDIM^2);
plotfilterbank(kimages',SDIM,pids(:))

%%
spike_cond_stim = X(spikebins,:);
sta      = mean(spike_cond_stim) - mean(X);
sta = sta/norm(sta);
stvcv =	cov(spike_cond_stim);  utvcv = cov(X);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

kimages = [sta',stcs]';
figure
pids = 1:(SDIM^2);
plotfilterbank(kimages',SDIM,pids(:))



%% First, refine the STC analysis by doubling and splitting st components
used_stc_dims = [1:4];
STCbvs = [sta' stcs];
STCbvs = STCbvs(:,used_stc_dims);
% STCbvs = [STCbvs -STCbvs(:,1:end)]; %make copies of STC comps
Nstcbvs = size(STCbvs,2);
nmods = Nstcbvs;
basis = 'white';
STCcf_0 = eye(Nstcbvs);
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
flen = 6; SDIM = 25; fsdim = SDIM^2;

%initialize model
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
defmod.SDIM = SDIM;
defmod.fsdim = SDIM^2;
defmod.pids = 1:fsdim;
pids = 1:fsdim;

% clear init_nls nltypes
% for i = 1:nmods; init_nls{i} = 'threshlin'; end;
% %define NL types: "uncon, lin, threshlin, quad"
% for i = 1:nmods; nltypes{i} = 'threshlin'; end;
% 
init_nls{1} = 'lin';
for i = 2:4
    init_nls{i} = 'pquad';
end
% for i = 7:nmods
%     init_nls{i} = 'nquad';
% end
nltypes{1} = 'lin';
for i = 2:nmods; nltypes{i} = 'quad';end
nltypes{1} = 'lin';
for i = 2:nmods; nltypes{i} = 'quad';end

% glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize

glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
glm_stcb.image_type = '2d';
% for i = 1:length(glm_stcb.mods)
%    cur_rand_filt = rand(ncomps,1);
%    glm_stcb.mods(i).k = cur_rand_filt;
%    glm_stcb.mods(i).pix = pix_conv_mat'*cur_rand_filt;
% end
[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
stc_glm = fitstc_fullbf(glm_stcb,WX,spikebins,'tots');
[nstc_glm,norm_vals] = normalizeRFs_full(stc_glm,WX);


glm_stcb.image_type = '2d';

% full_glm = fitNLHI2d_fullbf(glm_stcb,WX,spikebins,'tots',2);

%% First, refine the STC analysis by doubling and splitting st components
% used_stc_dims = [1:5];
% STCbvs = [sta' stcs];
% STCbvs = STCbvs(:,used_stc_dims);
% STCbvs = [STCbvs -STCbvs(:,1:end)]; %make copies of STC comps
% Nstcbvs = size(STCbvs,2);
% nmods = Nstcbvs;
% basis = 'white';
% STCcf_0 = eye(Nstcbvs);
% cd ~/Data/blanche/rec_75/matlabdata/
% load stdparsRec75.mat
% flen = 6;
% %initialize model
% defmod.h(1:end-1) = []; %eliminate PSC
% defmod.lnl = 0;
% defmod.lh = 0;
% defmod.lnl2 = 100;
% defmod.lh2 = 0;
% defmod.nlcon = 0;
% defmod.nlmon = 1;
% defmod.locLambda = 0;
% defmod.locSigma = 0;
% defmod.maxLocPen = 0;
% defmod.lambda_dX = 500; %350
% defmod.lambda_L1x = 50; %40
% defmod.lambda_dT = 10;
% defmod.SDIM = SDIM;
% defmod.fsdim = fsdim;
% defmod.pids = 1:fsdim;
%
% clear init_nls nltypes
% for i = 1:nmods; init_nls{i} = 'threshlin'; end;
% %define NL types: "uncon, lin, threshlin, quad"
% for i = 1:nmods; nltypes{i} = 'threshlin'; end;
% glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
%
% [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,scores);
% glm_stcb.image_type = '2d';
% full_glm = fitNLHI2d_fullbf(glm_stcb,scores,spikebins,'tots',2);
% % mod_filts = get_pix_mat(full_glm);
% % fin_filtproj = filt_long*mod_filts;
%
% f2 = plot2d_mod(glm_stcb);
% f2 = plot2d_mod(full_glm);
%
% %%
% % w_vec = arrayfun(@(x) x.w,full_glm.mods);
% % used_mods = find(abs(w_vec) > 0.1);
% % full_glm.mods = full_glm.mods(used_mods);
%% NOW FIND BEST OBLIQUE ROTATION WITHIN THE NEW SUBSPACE
basis_vecs = get_k_mat(stc_glm);
n_bvs = size(basis_vecs,2);
nmods = 2;
mod_signs = ones(nmods,1);
dim_signs = ones(n_bvs,1);
unused_stcs = (nmods+1):n_bvs;

flen = 6;
%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 100;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 1;
defmod.locLambda = 500;
defmod.lambda_dX = 100; %350
defmod.lambda_L1x = 50; %40
defmod.lambda_dT = 10;
defmod.pids = 1:fsdim;
defmod.SDIM = SDIM;
defmod.fsdim = fsdim;

%define NL initializations: "lin, threshlin, pquad, nquad"
clear init_nls nltypes
for i = 1:nmods; init_nls{i} = 'threshlin'; end;
%define NL types: "uncon, lin, threshlin, quad"
for i = 1:nmods; nltypes{i} = 'threshlin'; end;

Nstcbf = size(basis_vecs,2);
klen = size(pix_conv_mat,2);
flen = klen/fsdim;
kern_output = WX*basis_vecs;

%determine distribution of random interpoint distances
rand_reps = 500;
init_vals = zeros(rand_reps,n_bvs,nmods);
for r = 1:rand_reps
    % compute average separation between NN initial points
    STCcf_0 = randn(n_bvs,nmods);
    %normalize
    for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    init_vals(r,:,:) = STCcf_0;
end
cur_dists = zeros(1,rand_reps*(rand_reps-1)/2);
for nn = 1:nmods
    cur_dists = cur_dists + pdist(init_vals(:,:,nn),'cosine');
end
rdmean = 0.75*mean(cur_dists);
rdscale = 2*std(cur_dists);
% figure
% ksdensity(cur_dists)
% hold on
% tt = linspace(0,10,100);
% plot(tt,normcdf(tt,rdmean,rdscale),'k')

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
    %     white_props.pix_conv_mat = [];
    %     white_props.kern_conv_mat = [];
    glm_stcb = create2dGLM0_stcb(basis_vecs,STCcf_0,defmod,mod_signs,dim_signs,nltypes,init_nls,'test',white_props); %initialize
    glm_stcb.lambdaW = 0;
    
    %determine LL and LP at current filter point
    glm_stcb = fitWeights_stcb_nonlpsc(glm_stcb,kern_output*STCcf_0,spikebins,1,1e-3);
    [ll0, ll0p] = getLLGLM2d_STCBF_nonlpsc(glm_stcb,kern_output,spikebins,'none');
    
    rotbv_mod = [rotbv_mod; fitNLHI_stcb2d_nonlpsc(glm_stcb,WX,spikebins,'none',6,2)];
    LL_vals = [LL_vals rotbv_mod(end).LL];
    LP_vals = [LP_vals rotbv_mod(end).LP];
end
[~,best_mod] = min(LL_vals);
