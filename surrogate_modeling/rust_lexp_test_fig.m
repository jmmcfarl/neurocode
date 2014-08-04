cd ~/James_scripts/surrogate_modeling/
addpath(genpath('~/James_scripts'));
clear all
% close all

dt = 0.01; %in s
max_rate = Inf; %in Hz
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);

NT = 100000; SDIM = 16; flen = 14;
stim = round(2*rand(NT,SDIM))-1;
xvNT = 30000; 
xvstim = round(2*rand(xvNT,SDIM))-1;

% GENERATE A FILTER
nfilts = 8;
LAMBDA = 4;
x = repmat(1:SDIM,flen,1);
sigma1 = repmat(1.3,flen,SDIM);
lambda = repmat(LAMBDA,flen,SDIM);
amp_vec = ncx2pdf((1:flen-2)-1,4,1);
amp_vec = [0 0 amp_vec];
amp_vec = fliplr(amp_vec);
b = repmat(amp_vec',1,SDIM);
a = repmat(0,flen,SDIM);

% desired_spacing = 3*LAMBDA/4;
% desired_spacing = LAMBDA;
desired_spacing = LAMBDA/3;
beg = SDIM/2-desired_spacing*floor(nfilts/2)+1;
xovals = beg:desired_spacing:(beg+desired_spacing*(nfilts-1));
% xovals = linspace(4,13,nfilts);
% xovals = linspace(3,14,nfilts);
for i = 1:nfilts
%     psi_end = 5*pi + randn*5;
    psi_end = 5*pi;
    psi1 = repmat(linspace(0,psi_end,flen)',1,SDIM);
    psi1(1:5,:) = repmat(psi1(6,:),5,1);
    
    
    x0 = repmat(xovals(i),flen,SDIM);
    cur_psi = mod(psi1 + rand*2*pi,2*pi);
    temp = b.*exp(-((x-x0).^2./2./sigma1.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+a;
    filt(i,:,:) = temp/norm(temp(:));
    filt_mat(i,:) = temp(:)/norm(temp(:));
end
filt_mat = filt_mat';
c_cent = SDIM/2+0.5;
c_std = 6;
max_cval = 4;
c_offset = 3;
cvals = max_cval*exp(-(xovals-c_cent).^2/(2*c_std))+c_offset;

% CREATE TIME-EMBEDDED STIMULUS
stim_emb = makeStimRows(stim,flen);
xvstim_emb = makeStimRows(xvstim,flen);

% FILTER STIMULUS
for i = 1:nfilts
    temp = filt(i,:,:);
    filt_stim(i,:) = zscore(stim_emb*temp(:));
    xvfilt_stim(i,:) = zscore(xvstim_emb*temp(:));
end

for i = 1:8
subplot(8,1,i)
imagesc(squeeze(filt(i,:,:))*cvals(i));colormap(gray)
caxis([-2.5 2.5])
end
%% create spike function
%pass through internal NLs
beta = 4;
coefs = ones(1,nfilts).*cvals;
g = zeros(1,NT);
xvg = zeros(1,xvNT);
for i = 1:nfilts
    lfilt_stim(i,:) = filt_stim(i,:);
    lfilt_stim(i,filt_stim(i,:) < 0) = 0;   
    xvlfilt_stim(i,:) = xvfilt_stim(i,:);
    xvlfilt_stim(i,xvfilt_stim(i,:) < 0) = 0;
    g = g + coefs(i)*lfilt_stim(i,:);
    xvg = xvg + coefs(i)*xvlfilt_stim(i,:);
end
g = g/std(g);

target_rate = 50; %in Hz
target_pspike = target_rate*dt;
cur_theta = 3;
cur_beta = 1;
p_spike = 1/cur_beta*log(1+exp(cur_beta*(g-cur_theta)));
xvp_spike = log(1+exp((xvg-cur_theta)));
scale_f = target_rate/(mean(p_spike)/dt);
p_spike = p_spike*scale_f;
xvp_spike = xvp_spike*scale_f;

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
% spikebins = unique(spikebins); %max 1 spike per bin
fprintf('Nspks: %d\n',length(spikebins));

xvspikes = poissrnd(xvp_spike);
xvrbins = (find(xvspikes>0.5));
xvnsp = spikes(xvrbins);
xvspk_vals = unique(xvspikes); xvspk_vals(xvspk_vals==0) = [];
xvspikebins = [];
for i = 1:length(xvspk_vals)
    cur_set = find(xvspikes == xvspk_vals(i));
    xvspikebins = [xvspikebins; repmat(cur_set(:),xvspk_vals(i),1)];
end

nneg =0;
npos = 6;
spike_cond_stim = stim_emb(spikebins,:);
sta      = mean(spike_cond_stim) - mean(stim_emb);
sta = sta/norm(sta);
proj_mat = sta'/(sta*sta')*sta;
stim_proj = stim_emb - stim_emb*proj_mat;
stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
figure
subplot(2,6,1)
imagesc(reshape(sta,flen,SDIM));
for i = 1:6
    subplot(2,6,6+i)
    imagesc(reshape(stcs(:,i),flen,SDIM));
end
% for i = 1:5
%     subplot(3,5,10+i)
%     imagesc(reshape(stcs(:,5+i),flen,SDIM));
%     colormap(gray)
% end
    colormap(gray)

% stc_dims = [sta' stcs(:,1) stcs(:,5)];
stc_dims = [sta' stcs];

%% FIND SUBSPACE USING LEXP MODEL
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
% sdim = 12; flen = 10;
sdim = 16; flen = 14;
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_L2x = 0;

% defmod.lambda_dX = 0;
% defmod.lambda_L1x = 0; 
% defmod.lambda_dT = 0;

defmod.lambda_dX = 10;
defmod.lambda_L1x = 2; 
defmod.lambda_dT = 10;

defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
basis = 'pix';

nmods = 9;
cur_basis = stc_dims;
n_bvs = size(cur_basis,2);
STCcf_0 = randn(n_bvs,nmods);

%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
init_kerns = cur_basis*STCcf_0;
init_signs = ones(1,nmods);
glm_stcb = createGLM_quad(init_kerns,init_signs,defmod,basis,[],[]);
% glm_stcb.mods(1) = [];
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
glm_quad = fitGLM_lexp(glm_stcb,stim_emb,spikebins,'tots',[],1e-4,1e-6);
orig_mod = glm_quad;

% n_iter = 2;
% glm_quad_ref{1} = orig_mod;
% for ii = 2:n_iter+1
%     [glm_quad_ref{ii},norm_vals] = normalizeRFs_full(glm_quad_ref{ii-1},stim_emb);
%     g_mat = stim_emb*get_k_mat(glm_quad_ref{ii});
%     glm_quad_ref{ii} = fitWeights_lexp(glm_quad_ref{ii},g_mat,spikebins,0,[],1:nmods);
%     glm_quad_ref{ii} = fitGLM_lexp(glm_quad_ref{ii},stim_emb,spikebins,'tots',[],1e-3,1e-5,[]);
% end
% glm_quad = glm_quad_ref{end};
% [glm_quad,coms,peak_locs] = get_filter_coms_1d(glm_quad);
% [~,ord] = sort(coms);
% glm_quad.mods = glm_quad.mods(ord);
% 
xvLL_quad = getLLGLM_lexp(glm_quad,xvstim_emb,xvspikebins,'none');
% 
% g_mat = stim_emb*get_k_mat(glm_quad);
% [glm_quad_nl_ref] = fitNL_lexp(glm_quad,g_mat,spikebins,0);
%%
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
% sdim = 12; flen = 10;
sdim = 16; flen = 14;
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;

% defmod.lambda_dX = 0;
% defmod.lambda_L1x = 0; 
% defmod.lambda_dT = 0;

defmod.lambda_dX = 10;
defmod.lambda_L1x = 5; 
defmod.lambda_dT = 10;

defmod.lambda_L2x = 0;
defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
basis = 'pix';

nmods = 8;
cur_basis = stc_dims;
n_bvs = size(cur_basis,2);
STCcf_0 = randn(n_bvs,nmods);
%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
init_kerns = cur_basis*STCcf_0;
% init_signs = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1];
init_signs = ones(1,nmods);
init_betas = 2*ones(nmods,1);
init_thetas = zeros(nmods,1);
% glm_stcb = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
glm_stcb = createGLM_lexp(init_kerns,init_signs,init_betas,init_thetas,defmod,basis,[],[]);
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
glm_lexp = fitGLM_lexp(glm_stcb,stim_emb,spikebins,'tots',[],1e-4,1e-6);
orig_mod = glm_lexp;

% n_iter = 3;
% glm_lexp_ref{1} = orig_mod;
% for ii = 2:n_iter+1
%     [glm_lexp_ref{ii},norm_vals] = normalizeRFs_full(glm_lexp_ref{ii-1},stim_emb);
%     g_mat = stim_emb*get_k_mat(glm_lexp_ref{ii});
%     glm_lexp_ref{ii} = fitWeights_full(glm_lexp_ref{ii},g_mat,spikebins,0,[],1:nmods);
%     glm_lexp_ref{ii} = fitGLM_lexp(glm_lexp_ref{ii},stim_emb,spikebins,'tots',[],1e-3,1e-5,[]);
% end
% glm_lexp = glm_lexp_ref{end};
% glm_lexp = fitGLM_lexp(glm_lexp,stim_emb,spikebins,'tots',[],1e-4,1e-6);
[glm_lexp,coms,peak_locs] = get_filter_coms_1d(glm_lexp);
[~,ord] = sort(coms);
glm_lexp.mods = glm_lexp.mods(ord);
% 
xvLL_lexp = getLLGLM_lexp(glm_lexp,xvstim_emb,xvspikebins,'none');
%%
quad_k_mat = get_k_mat(glm_quad);


overlap_lexp = subspace_overlap(lexp_k_mat,filt_mat);
overlap_quad = subspace_overlap(quad_k_mat,filt_mat);
overlap_stc = subspace_overlap(stc_dims,filt_mat);

soverlap_lexp = subspace(filt_mat,lexp_k_mat);
% soverlap_lexp2 = subspace(filt_mat,lexp2_k_mat);
soverlap_quad = subspace(filt_mat,quad_k_mat);
soverlap_stc = subspace(filt_mat,stc_dims);
