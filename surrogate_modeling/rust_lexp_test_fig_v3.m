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
nfilts = 6;
LAMBDA = 4;
x = repmat(1:SDIM,flen,1);
sigma1 = repmat(1.5,flen,SDIM);
lambda = repmat(LAMBDA,flen,SDIM);
amp_vec = ncx2pdf((1:flen-2)-1,4,1);
amp_vec = [0 0 amp_vec];
amp_vec = fliplr(amp_vec);
b = repmat(amp_vec',1,SDIM);
a = repmat(0,flen,SDIM);

% desired_spacing = 3*LAMBDA/4;
% desired_spacing = LAMBDA;
% desired_spacing = LAMBDA/6;
desired_spacing = LAMBDA/6.5;
beg = SDIM/2-desired_spacing*floor(nfilts/2)+0.75;
% xovals = beg:desired_spacing:(beg+desired_spacing*(nfilts-1));
xovals = 0:desired_spacing:(desired_spacing*(nfilts-1));
xovals = xovals - mean(xovals);
xovals = xovals + SDIM/2+0.5;
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
% c_std = 2;
c_std = 2;
max_cval = 5;
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

figure
for i = 1:nfilts
subplot(nfilts,1,i)
% imagesc(squeeze(filt(i,:,:))*cvals(i));colormap(gray)
imagesc(squeeze(filt(i,:,:)));colormap(gray)
% caxis([-3 3])
end
%% create spike function
%pass through internal NLs
beta = 3;
theta = 0;
coefs = ones(1,nfilts).*cvals;
g = zeros(1,NT);
xvg = zeros(1,xvNT);
for i = 1:nfilts
%     lfilt_stim(i,:) = filt_stim(i,:);
%     lfilt_stim(i,filt_stim(i,:) < 0) = 0;   
%     xvlfilt_stim(i,:) = xvfilt_stim(i,:);
%     xvlfilt_stim(i,xvfilt_stim(i,:) < 0) = 0;
    lfilt_stim(i,:) = 1/beta*log(1+exp(beta*(filt_stim(i,:)-theta)));
    xvlfilt_stim(i,:) = 1/beta*log(1+exp(beta*(xvfilt_stim(i,:)-theta)));
    g = g + coefs(i)*lfilt_stim(i,:);
    xvg = xvg + coefs(i)*xvlfilt_stim(i,:);
end
g = g/std(g);

target_rate = 50; %in Hz
target_pspike = target_rate*dt;
cur_theta = 1.75; %2.5
cur_beta = 4; %4
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
defmod.kscale = 1;
% defmod.lambda_dX = 0;
% defmod.lambda_L1x = 0; 
% defmod.lambda_dT = 0;

defmod.lambda_dX = 0;
defmod.lambda_L1x = 15; %was 20
defmod.lambda_dT = 0;
defmod.lambda_d2X = 60;
defmod.lambda_d2XT = 0;
% defmod.lambda_dX = 10;
% defmod.lambda_L1x = 5; %was 20
% defmod.lambda_dT = 10;


defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
basis = 'pix';


nmods = 4;
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
glm_quad = fitGLM_lexp(glm_stcb,stim_emb,spikebins,'tots',[],1e-3,1e-5);
orig_mod = glm_quad;

n_iter = 2;
clear glm_quad_ref
glm_quad_ref{1} = orig_mod;
for ii = 2:n_iter+1
    [glm_quad_ref{ii},norm_vals] = normalizeRFs_full(glm_quad_ref{ii-1},stim_emb,1);
    g_mat = stim_emb*get_k_mat(glm_quad_ref{ii});
    glm_quad_ref{ii} = fitWeights_lexp(glm_quad_ref{ii},g_mat,spikebins,0);
%     for i = 1:length(glm_quad_ref{ii}.mods)
%         glm_quad_ref{ii}.mods(i).lambda_dX = defmod.lambda_dX*glm_quad_ref{ii}.mods(i).w^2;
%         glm_quad_ref{ii}.mods(i).lambda_dT = defmod.lambda_dT*glm_quad_ref{ii}.mods(i).w^2;
%         glm_quad_ref{ii}.mods(i).lambda_L1x = defmod.lambda_L1x*glm_quad_ref{ii}.mods(i).w;
%     end
    glm_quad_ref{ii} = fitGLM_lexp(glm_quad_ref{ii},stim_emb,spikebins,'tots',[],1e-3,1e-5,[]);
end
glm_quadf = glm_quad_ref{end};
[glm_quadf,norm_vals] = normalizeRFs_full(glm_quadf,stim_emb,1);
g_mat = stim_emb*get_k_mat(glm_quadf);
glm_quadf = fitWeights_lexp(glm_quadf,g_mat,spikebins,0);
% for i = 1:length(glm_quadf.mods)
%     glm_quadf.mods(i).lambda_dX = defmod.lambda_dX*glm_quadf.mods(i).w^2;
%     glm_quadf.mods(i).lambda_dT = defmod.lambda_dT*glm_quadf.mods(i).w^2;
%     glm_quadf.mods(i).lambda_L1x = defmod.lambda_L1x*glm_quadf.mods(i).w;
% end
glm_quadf = fitGLM_lexp(glm_quadf,stim_emb,spikebins,'tots',[],1e-4,1e-6);
[glm_quadf,coms,peak_locs] = get_filter_coms_1d(glm_quadf);
[~,ord] = sort(coms);
glm_quadf.mods = glm_quadf.mods(ord);

[glm_quadfN,norm_vals] = normalizeRFs_full(glm_quadf,stim_emb);
g_mat = stim_emb*get_k_mat(glm_quadfN);
glm_quadfN = fitWeights_lexp(glm_quadfN,g_mat,spikebins,0);

% g_mat = stim_emb*get_k_mat(glm_quad);
% [glm_quad_nl_ref] = fitNL_quad(glm_quad,g_mat,spikebins,0);

xvLL_quad = getLLGLM_lexp(glm_quadf,xvstim_emb,xvspikebins,'tots')
%%
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
old_nlx = defmod.nlx;
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
defmod.nlx = linspace(-3.1,3.1,21);
defmod.kscale = 1;
% defmod.lambda_dX = 0;
% defmod.lambda_L1x = 0; 
% defmod.lambda_dT = 0;

% defmod.lambda_dX = 10;
% defmod.lambda_L1x = 10; 
% defmod.lambda_dT = 10;

% defmod.lambda_dX = 30;
% defmod.lambda_L1x =10; 
% defmod.lambda_dT = 30;

defmod.lambda_dX = 200;
defmod.lambda_L1x = 50; 
defmod.lambda_dT = 200;
defmod.lambda_d2X = 0;
defmod.lambda_d2XT = 0;
defmod.lambda_L2x = 0;
defmod.kscale = 1;
defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
basis = 'pix';

nmods = 6;
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
glm_stcb = createGLM_tlin(init_kerns,init_signs,defmod,basis,[],[]);
% glm_stcb = createGLM_lexp(init_kerns,init_signs,init_betas,init_thetas,defmod,basis,[],[]);
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
glm_lexp = fitGLM_lexp(glm_stcb,stim_emb,spikebins,'tots',[],1e-4,1e-6);
orig_mod = glm_lexp;

% n_iter = 2;
% clear glm_lexp_ref
% glm_lexp_ref{1} = orig_mod;
% for ii = 2:n_iter+1
%     [glm_lexp_ref{ii},norm_vals] = normalizeRFs_full(glm_lexp_ref{ii-1},stim_emb,1);
%     g_mat = stim_emb*get_k_mat(glm_lexp_ref{ii});
%     glm_lexp_ref{ii} = fitWeights_lexp(glm_lexp_ref{ii},g_mat,spikebins,0,1:nmods);
%     glm_lexp_ref{ii} = fitGLM_lexp(glm_lexp_ref{ii},stim_emb,spikebins,'tots',[],1e-3,1e-5,[]);
% end
% glm_lexpf = glm_lexp_ref{end};
% [glm_lexpf,norm_vals] = normalizeRFs_full(glm_lexpf,stim_emb,1);
% g_mat = stim_emb*get_k_mat(glm_lexpf);
% glm_lexpf = fitWeights_lexp(glm_lexpf,g_mat,spikebins,0,1:nmods);
% glm_lexpf = fitGLM_lexp(glm_lexpf,stim_emb,spikebins,'tots',[],1e-4,1e-6);
% [glm_lexpf,coms,peak_locs] = get_filter_coms_1d(glm_lexpf);
% [~,ord] = sort(coms);
% glm_lexpf.mods = glm_lexpf.mods(ord);
% 
% [glm_lexpfN,norm_vals] = normalizeRFs_full(glm_lexpf,stim_emb);
% g_mat = stim_emb*get_k_mat(glm_lexpfN);
% glm_lexpfN = fitWeights_lexp(glm_lexpfN,g_mat,spikebins,0);


% g_mat = stim_emb*get_k_mat(glm_lexp);
% [glm_lexp_nl_ref] = fitNL_lexp(glm_lexp,g_mat,spikebins,0);

xvLL_lexp = getLLGLM_lexp(glm_lexp,xvstim_emb,xvspikebins,'tots')

%% RUN A FEW MORE ITERATIONS TO ENSURE CONVERGENCE
n_iter = 1;
clear glm_lexp_ref
glm_lexp_ref{1} = glm_lexp;
for n = 2:n_iter+1
    glm_lexp_ref{n} = fitGLM_lexp(glm_lexp_ref{n-1},stim_emb,spikebins,'tots',[],1e-4,1e-6);
    xvLL_lexp_ref(n) = getLLGLM_lexp(glm_lexp_ref{n},xvstim_emb,xvspikebins,'tots')
end

%%
temp = glm_lexp_ref{1};
temp = adjust_all_reg(temp,'nltype','lexp');
temp = adjust_all_reg(temp,'lambda_L1x',100);
temp = adjust_all_reg(temp,'lambda_dX',100);
temp = adjust_all_reg(temp,'lambda_dT',100);
for i = 1:length(temp.mods)
    temp.mods(i).beta = 2;
    temp.mods(i).theta = 0;
end
temp = normalizeRFs_full(temp,stim_emb,1); %This normalizes filter coefs so that the filtered stimulus distribution is standard normal (not really necessary)
k_mat = get_k_mat(temp);
input = stim_emb*k_mat;
temp = fitWeights_full(temp,input,spikebins,1);
temp = fitGLM_lexp(temp,stim_emb,spikebins,'tots',[],1e-4,1e-6,[]); %fit the model
%%
old_nlx = defmod.nlx;
old_tlin_nly = old_nlx; old_tlin_nly(old_nlx < 0) = 0;
fo = glm_lexp;
fo = adjust_all_reg(fo,'nlx',old_nlx);
fo = adjust_all_reg(fo,'nly',old_tlin_nly);
fo = adjust_all_reg(fo,'nltype','uncon');
fo = adjust_all_reg(fo,'nlmon',1);
fo = adjust_all_reg(fo,'lnl2',1000);
% fo = adjust_all_reg(fo,'lambda_L1x',100);

init_LL = fo.LL;
n_iter = 2;
ref_LL = zeros(n_iter,1);
clear fref
fref{1} = fo;
for n = 2:n_iter+1
    fref{n} = normalizeRFs_full(fref{n-1},stim_emb,1); %This normalizes filter coefs so that the filtered stimulus distribution is standard normal (not really necessary)
    k_mat = get_k_mat(fref{n});
    input = stim_emb*k_mat;
    fref{n} = fitWeights_lexp(fref{n},input,spikebins,1,1:nmods);
    fref{n} = fitNL_lexp(fref{n},input,spikebins,0);
    fref{n} = fitGLM_lexp(fref{n},stim_emb,spikebins,'tots',[],1e-4,1e-6,[]); %fit the model
    ref_LL(n) = fref{n}.LL;
end

[glm_lexpf,coms,peak_locs] = get_filter_coms_1d(fref{end});
[~,ord] = sort(coms);
glm_lexpf.mods = glm_lexpf.mods(ord);

[glm_lexpfN,norm_vals] = normalizeRFs_full(glm_lexpf,stim_emb);
g_mat = stim_emb*get_k_mat(glm_lexpfN);
glm_lexpfN = fitWeights_lexp(glm_lexpfN,g_mat,spikebins,0);

%%
quad_k_mat = get_k_mat(glm_quad);
lexp_k_mat = get_k_mat(glm_lexp);

overlap_lexp = subspace_overlap2(lexp_k_mat,filt_mat);
overlap_quad = subspace_overlap2(quad_k_mat,filt_mat);
overlap_stc = subspace_overlap2(stc_dims,filt_mat);

overlap_lexpquad = subspace_overlap(lexp_k_mat,quad_k_mat);
overlap_lexp = subspace_overlap(filt_mat,lexp_k_mat);
overlap_quad = subspace_overlap(filt_mat,quad_k_mat);
overlap_stc = subspace_overlap(stc_dims,filt_mat);

soverlap_lexp = subspace(filt_mat,lexp_k_mat);
% soverlap_lexp2 = subspace(filt_mat,lexp2_k_mat);
soverlap_quad = subspace(filt_mat,quad_k_mat);
soverlap_stc = subspace(filt_mat,stc_dims);

%%
% cd ~/Documents/GNM_paper/
% save lexp_rust_examp2 glm_quad glm_lexp filt* sdim flen xvLL*