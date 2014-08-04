cd ~/James_scripts/surrogate_modeling/
addpath(genpath('~/James_scripts'));
clear all
close all

dt = 0.01; %in s
max_rate = Inf; %in Hz
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);

NT = 100000; SDIM = 16; flen = 14;
stim = round(2*rand(NT,SDIM))-1;
xvNT = 30000; 
xvstim = round(2*rand(xvNT,SDIM))-1;

% GENERATE A FILTER
nfilts = 4;
LAMBDA = 4.25;
SIGMA = 1.5;
x = repmat(1:SDIM,flen,1);
sigma1 = repmat(SIGMA,flen,SDIM);
lambda = repmat(LAMBDA,flen,SDIM);
amp_vec = ncx2pdf((1:flen-2)-1,4,1);
amp_vec = [0 0 amp_vec];
amp_vec = fliplr(amp_vec);
b = repmat(amp_vec',1,SDIM);
a = repmat(0,flen,SDIM);
psi1 = repmat(linspace(0,5*pi,flen)',1,SDIM);
psi1(1:5,:) = repmat(psi1(6,:),5,1);

% xovals = linspace(3.,12,nfilts);
xovals = linspace(4.,13,nfilts);
% xovals = linspace(6,11,nfilts);
% xovals = linspace(5,12,nfilts);
for i = 1:nfilts
    x0 = repmat(xovals(i),flen,SDIM);
    cur_psi = mod(psi1 + rand*2*pi,2*pi);
    temp = b.*exp(-((x-x0).^2./2./sigma1.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+a;
    filt(i,:,:) = temp/norm(temp(:));
    filt_mat(i,:) = temp(:)/norm(temp(:));
end
filt_mat = filt_mat';
c_cent = 8.5;
c_std = 6;
max_cval = 4;
c_offset = 4;
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
%% create spike function
%pass through internal NLs
beta = 4;
coefs = ones(1,nfilts).*cvals;
g = zeros(1,NT);
xvg = zeros(1,xvNT);
for i = 1:nfilts
    lfilt_stim(i,:) = filt_stim(i,:).^2;
    xvlfilt_stim(i,:) = xvfilt_stim(i,:).^2;
    g = g + coefs(i)*lfilt_stim(i,:);
    xvg = xvg + coefs(i)*xvlfilt_stim(i,:);
end
g = g/std(g);
xvg = xvg/std(xvg);

target_rate = 50; %in Hz
target_pspike = target_rate*dt;
cur_theta = 1.5;
cur_beta = 1;
p_spike = 1/cur_beta*log(1+exp(cur_beta*(g-cur_theta)));
xvp_spike = 1/cur_beta*log(1+exp(cur_beta*(xvg-cur_theta)));
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

nneg =5;
npos = 5;
spike_cond_stim = stim_emb(spikebins,:);
sta      = mean(spike_cond_stim) - mean(stim_emb);
sta = sta/norm(sta);
proj_mat = sta'/(sta*sta')*sta;
stim_proj = stim_emb - stim_emb*proj_mat;
stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
figure
subplot(3,5,1)
imagesc(reshape(sta,flen,SDIM));
for i = 1:5
    subplot(3,5,5+i)
    imagesc(reshape(stcs(:,i),flen,SDIM));
end
for i = 1:5
    subplot(3,5,10+i)
    imagesc(reshape(stcs(:,5+i),flen,SDIM));
    colormap(gray)
end

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
defmod.lambda_dX = 20;
defmod.lambda_L1x = 10; %2
defmod.lambda_dT = 20;
defmod.lambda_L2x = 0;
defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
basis = 'pix';

nmods = 5;
cur_basis = stc_dims;
n_bvs = size(cur_basis,2);
STCcf_0 = randn(n_bvs,nmods);

%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
init_kerns = cur_basis*STCcf_0;
init_signs = [1 1 1 1 1 1];
glm_stcb = createGLM_quad(init_kerns,init_signs,defmod,basis,[],[]);
% glm_stcb.mods(1) = [];
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
glm_quad = fitGLM_lexp(glm_stcb,stim_emb,spikebins,'tots');

[glm_quad,norm_vals] = normalizeRFs_full(glm_quad,stim_emb);
g_mat = stim_emb*get_k_mat(glm_quad);
glm_quad = fitWeights_full(glm_quad,g_mat,spikebins,0);
% for i = 1:length(glm_quad.mods)
%     glm_quad.mods(i).lambda_dX = glm_quad.mods(i).lambda_dX*norm_vals(i)^2;
%     glm_quad.mods(i).lambda_L1x = glm_quad.mods(i).lambda_L1x*norm_vals(i);
%     glm_quad.mods(i).lambda_dT = glm_quad.mods(i).lambda_dT*norm_vals(i)^2;
% end
glm_quad = fitGLM_lexp(glm_quad,stim_emb,spikebins,'tots',[],1e-3,1e-5,[]);
orig_mod = glm_quad;

n_iter = 2;
glm_quad_ref{1} = orig_mod;
for ii = 2:n_iter+1
    [glm_quad_ref{ii},norm_vals] = normalizeRFs_full(glm_quad_ref{ii-1},stim_emb);
    g_mat = stim_emb*get_k_mat(glm_quad_ref{ii});
    glm_quad_ref{ii} = fitWeights_full(glm_quad_ref{ii},g_mat,spikebins,0,[],1:nmods);
%     for i = 1:length(glm_tlin_ref{ii}.mods)
%         glm_tlin_ref{ii}.mods(i).lambda_dX = glm_tlin_ref{ii}.mods(i).lambda_dX*norm_vals(i)^2;
%         glm_tlin_ref{ii}.mods(i).lambda_L1x = glm_tlin_ref{ii}.mods(i).lambda_L1x*norm_vals(i);
%         glm_tlin_ref{ii}.mods(i).lambda_dT = glm_tlin_ref{ii}.mods(i).lambda_dT*norm_vals(i)^2;
%     end
    glm_quad_ref{ii} = fitGLM_lexp(glm_quad_ref{ii},stim_emb,spikebins,'tots',[],1e-3,1e-5,[]);
end
glm_quad_ref{end} = fitGLM_lexp(glm_quad_ref{end},stim_emb,spikebins,'tots',[],1e-4,1e-6,[]);
glm_quad = glm_quad_ref{end};
[glm_quad,coms,peak_locs] = get_filter_coms_1d(glm_quad);
[~,ord] = sort(coms);
glm_quad.mods = glm_quad.mods(ord);

xvLL_quad = getLLGLM_lexp(glm_quad,xvstim_emb,xvspikebins,'none');

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
defmod.locSigma = 0;
defmod.maxLocPen = 0;
defmod.lambda_dX = 40;
defmod.lambda_L1x = 40; %6
defmod.lambda_dT = 40;
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
% glm_stcb = createGLM_rquad(init_kerns,init_signs,defmod,basis,[],[]);
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
glm_tlin = fitGLM_lexp(glm_stcb,stim_emb,spikebins,'tots',[],1e-4,1e-6);
orig_mod = glm_tlin;

% n_iter = 3;
% glm_tlin_ref{1} = orig_mod;
% for ii = 2:n_iter+1
%     [glm_tlin_ref{ii},norm_vals] = normalizeRFs_full(glm_tlin_ref{ii-1},stim_emb);
%     g_mat = stim_emb*get_k_mat(glm_tlin_ref{ii});
%     glm_tlin_ref{ii} = fitWeights_lexp(glm_tlin_ref{ii},g_mat,spikebins,0,[],1:nmods);
% %     for i = 1:length(glm_tlin_ref{ii}.mods)
% %         glm_tlin_ref{ii}.mods(i).lambda_dX = glm_tlin_ref{ii}.mods(i).lambda_dX*norm_vals(i)^2;
% %         glm_tlin_ref{ii}.mods(i).lambda_L1x = glm_tlin_ref{ii}.mods(i).lambda_L1x*norm_vals(i);
% %         glm_tlin_ref{ii}.mods(i).lambda_dT = glm_tlin_ref{ii}.mods(i).lambda_dT*norm_vals(i)^2;
% %     end
% glm_tlin_ref{ii} = fitGLM_lexp(glm_tlin_ref{ii},stim_emb,spikebins,'tots',[],1e-3,1e-5,[]);
% end
% glm_tlin_ref{end} = fitGLM_lexp(glm_tlin_ref{end},stim_emb,spikebins,'tots',[],1e-5,1e-7,[]);
% glm_tlin = glm_tlin_ref{end};
[glm_tlin,coms,peak_locs] = get_filter_coms_1d(glm_tlin);
[~,ord] = sort(coms);
glm_tlin.mods = glm_tlin.mods(ord);
% 
% xvLL_lexp = getLLGLM_lexp(glm_tlin,xvstim_emb,xvspikebins,'none');
% 
g_mat = stim_emb*get_k_mat(glm_tlin);
[glm_tlin_nlref] = fitNL_lexp(glm_tlin,g_mat,spikebins,0);

%%
