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
clear all; flen =6; SDIM = 14;

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
used_pixs = find(Xr <= 14 & Yr <= 14);

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
clear allstims

X = bsxfun(@minus,X,mean(X));

%% generate surrogate filters
SDIM = 14; flen = 6; dt = 0.02;
x = repmat(1:SDIM,[flen 1 SDIM]);
y = permute(x,[1 3 2]);

x0 = permute(repmat(10,[SDIM SDIM flen]),[3 1 2]);
y0 = permute(repmat(10,[SDIM SDIM flen]),[3 1 2]);
sigmax = permute(repmat(1.5,[SDIM SDIM flen]),[3 1 2]);
sigmay = permute(repmat(2,[SDIM SDIM flen]),[3 1 2]);
lambda = permute(repmat(15,[SDIM SDIM flen]),[3 1 2]);
theta = permute(repmat(pi/4,[SDIM SDIM flen]),[3 1 2]);
b = repmat(linspace(2,0,flen)',[1 SDIM SDIM]);
psi = repmat(linspace(0,pi/2,flen)',[1 SDIM SDIM]);

psi2 = psi+pi/2;
x02 = permute(repmat(5,[SDIM SDIM flen]),[3 1 2]);
y02 = permute(repmat(5,[SDIM SDIM flen]),[3 1 2]);

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

%%
filt_stim1 = filt_proj_stim(:,1);
filt_stim1(filt_stim1 < 0) = 0;

filt_stim2 = filt_proj_stim(:,2);
filt_stim2(filt_stim2 < 0) = 0;


g = filt_stim1 + filt_stim2;
g = 2*zscore(g);

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
ncomps = 500; compids = 1:ncomps;

% scalmat = diag(1./sqrt(latent(1:ncomps)));
scalmat = diag(ones(1,ncomps));

WX = score(:,1:ncomps)*scalmat;

WS        = WX(spikebins,:);
sta      = mean(WS) - mean(WX);
stvcv = cov(WS);  utvcv = cov(WX);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);


% pix_conv_mat = diag(sqrt(latent(compids)))*coeff(:,compids)';
% kern_conv_mat = coeff(:,compids)*diag(1./sqrt(latent(compids)));
pix_conv_mat = coeff(:,compids)';
kern_conv_mat = coeff(:,compids);

STCbvs = ([sta', stcs]'*diag(1./sqrt(latent(compids))))';
STCbvs = [sta', stcs];
kimages = [sta',stcs]'*diag(sqrt(latent(compids)))*coeff(:,compids)';
% kimages = [sta',stcs]'*coeff(:,compids)';

% figure
% pids = 1:(SDIM^2);
% plotfilterbank(kimages',SDIM,pids(:))

%%
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
SDIM = 14;
flen = 6;
%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.lambda_dX = 500;
defmod.lambda_L1x = 0.1;
defmod.lambda_dT = 2;
defmod.pids = 1:SDIM^2;
defmod.SDIM = SDIM;
defmod.fsdim = SDIM^2;
defmod.locLambda = 0;
basis = 'white';

%first fit sta model
cur_ndims = 3;
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
glm_stcb = createGLM2d_fullbf(cur_basis,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('Cell')); %initialize
[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,WX);
%     [glm_stcb,sub_means] = mean_sub_filters(glm_stcb);
glm_stcb.image_type = '2d';
stc_glm{cur_ndims} = fitNLHI2d_fullbf(glm_stcb,WX,spikebins,'tots',4);

w = arrayfun(@(x) x.w,stc_glm{cur_ndims}.mods);
stc_glm{cur_ndims}.mods(w==0) = [];

f1 = plot2d_mod(stc_glm{cur_ndims});

subspace = get_pix_mat(stc_glm{cur_ndims});

proj_mat = subspace*inv(subspace'*subspace)*subspace';
filt_proj_sub = filt_long'*proj_mat;

fract_captured = sum(filt_proj_sub'.^2)./sum(filt_long.^2);
