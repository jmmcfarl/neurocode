cd ~/James_scripts/surrogate_modeling/
clear all
close all

load ./rustlike_stim
dt = 0.01; %in s
max_rate = 200; %in Hz
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);

NT = 100000; SDIM = 16;
stim = round(2*rand(NT,SDIM))-1;

% NT = size(stim,1);
% ep = 60000; stim = stim(1:ep,:);
%% GENERATE A FILTER
nfilts = 4;
x = repmat(1:SDIM,flen,1);
sigma = repmat(1.5,flen,SDIM);
lambda = repmat(6,flen,SDIM);
b = repmat(linspace(2,0,flen)',1,SDIM);a = repmat(0,flen,SDIM);
psi1 = repmat(linspace(0,pi,flen)',1,SDIM);
xovals = linspace(4,12,nfilts);
for i = 1:nfilts
    x0 = repmat(xovals(i),flen,SDIM);
    temp = b.*exp(-((x-x0).^2./2./sigma.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+a;
    filt(i,:,:) = temp/norm(temp(:));
    filt_mat(i,:) = temp(:)/norm(temp(:));
end
filt_mat = filt_mat';

%% CREATE TIME-EMBEDDED STIMULUS
stim_emb = makeStimRows(stim,flen);

%% FILTER STIMULUS
for i = 1:nfilts
    temp = filt(i,:,:);
    filt_stim(i,:) = zscore(stim_emb*temp(:));
end
%% create spike function
%pass through internal NLs
beta = 4;
% coefs = [1 -1.5 1 -1.5 1 -1.5 1 -1.5]/2;
coefs = [1 -1.5 1 -1.5]/2;
g = zeros(1,NT);
for i = 1:nfilts
%     lfilt_stim(i,:) = 1/beta*log(1+exp(beta*filt_stim(i,:)));
    lfilt_stim(i,:) = filt_stim(i,:).^2;
    g = g + coefs(i)*lfilt_stim(i,:);
end


target_rate = 20; %in Hz
target_pspike = target_rate*dt;
% K0 = [1 0.4];
K0 = [1 0.5];
% Kfit = fmincon(@(K) abs(mean(logexp(g,K)) - target_pspike),K0,[],[],[],[],[-10 .01],[10 10]);
Kfit = K0;
% p_spike = logexp(g,Kfit);
p_spike = log(1+3*exp((g-2)));
% p_spike(p_spike > 1) = 1;

figure
[y,x] = ksdensity(g);
plot(x,y); hold on
yl = ylim();
tx = -10:.02:10;
plot(tx,log(1+3*exp((tx-2))),'k')
ylim(yl)
%
% % f1 = X; f1(X < 0) = 0;
% % f2 = Y; f2(Y < 0) = 0;
% f1 = logexp(X,[0 0.2]);
% f2 = logexp(Y,[0 0.2]);
% gfun = f1+1.5*f2;
% gfun(X.^2 + Y.^2 > 9) = nan;
% figure
% subplot(2,1,1)
% contourf(X,Y,gfun,30)
% subplot(2,1,2)
% contourf(X,Y,logexp(gfun,Kfit),30);

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

%%
nneg = 4;
npos = 4;
spike_cond_stim = stim_emb(spikebins,:);
sta      = mean(spike_cond_stim) - mean(stim_emb);
sta = sta/norm(sta);
proj_mat = sta'/(sta*sta')*sta;
stim_proj = stim_emb - stim_emb*proj_mat;
stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
figure
subplot(3,4,1)
imagesc(reshape(sta,flen,SDIM));
for i = 1:4
    subplot(3,4,4+i)
    imagesc(reshape(stcs(:,i),flen,SDIM));
end
for i = 1:4
    subplot(3,4,8+i)
    imagesc(reshape(stcs(:,4+i),flen,SDIM));
    colormap(gray)
end

%%
% stc_dims = [sta' stcs(:,1:3) stcs(:,5:7)];
stc_dims = [stcs(:,1:2) stcs(:,5:6)];

%%
% minxy = [-4 -4];
% maxxy = [4 4];
% stc_stim1 = stim_emb*stc_dims(:,1);
% stc_stim2 = stim_emb*stc_dims(:,2);
% stc_proj_stim = [stc_stim1 stc_stim2];
% 
% [bandwidth2_s1,density2_s1,X2,Y2]=kde2d(stc_proj_stim(:,1:2),2^8,minxy,maxxy);
% % [bandwidth2_s2,density2_s2,X2,Y2]=kde2d(filt_proj_stim(:,3:4),2^8,minxy,maxxy);
% 
% [bandwidth_s1,density_s1,X,Y]=kde2d(stc_proj_stim(spikebins,1:2),2^8,minxy,maxxy,bandwidth2_s1);
% % [bandwidth_s2,density_s2,X,Y]=kde2d(filt_proj_stim(spikebins,3:4),2^8,minxy,maxxy,bandwidth2_s2);
% 

%%
% 
% filt_long = [filt1(:) filt2(:) filt3(:) filt4(:)]';
% 
% filt_stc_proj = filt_long*stc_dims;
% 
% eps = 1e-3;
% cond_dens_s1 = density_s1./density2_s1;
% cond_dens_s1(density2_s1 < eps) = eps;
% % cond_dens_s2 = density_s2./density2_s2;
% % cond_dens_s2(density2_s2 < eps) = eps;
% outside = find(X.^2+Y.^2>9);
% cond_dens_s1(outside) = nan;
% density_s1(outside) = nan;
% % cond_dens_s2(outside) = nan;
% % density_s2(outside) = nan;
% 
% figure
% contourf(X,Y,cond_dens_s1,50)
% line(3*[0 filt_stc_proj(1,1)],3*[0 filt_stc_proj(1,2)],'color','r','linewidth',2)
% line(3*[0 filt_stc_proj(2,1)],3*[0 filt_stc_proj(2,2)],'color','w','linewidth',2)
% line(3*[0 filt_stc_proj(3,1)],3*[0 filt_stc_proj(3,2)],'color','g','linewidth',2)
% line(3*[0 filt_stc_proj(4,1)],3*[0 filt_stc_proj(4,2)],'color','k','linewidth',2)

%%
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
sdim = 16; flen = 14;
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 0;
defmod.lambda_L1x = 0; %2
defmod.lambda_dT = 0;
defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
basis = 'pix';

nmods = 8;
cur_basis = stc_dims;
% cur_basis = [filt1(:) filt2(:)];
n_bvs = size(cur_basis,2);
STCcf_0 = randn(n_bvs,nmods);
%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
init_kerns = cur_basis*STCcf_0;
init_signs = [1 -1 1 -1 1 -1 1 -1];
init_betas = 2*ones(1,nmods);

% STCcf_0 = [1 0;
%     0 1];
% glm_stcb.mods(1).w = 2; glm_stcb.mods(2).w = 2;

glm_stcb = createGLM_lexp_rot(cur_basis,STCcf_0,init_signs,init_betas,defmod,basis);
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
glm_lexp2 = fitGLM_lexp(glm_stcb,stim_emb,spikebins,'tots');

[glm_lexp2,norm_vals] = normalizeRFs_full(glm_lexp2,stim_emb);
g_mat = stim_emb*get_k_mat(glm_lexp2);
glm_lexp2 = fitWeights_lexp(glm_lexp2,g_mat,spikebins,0);
% glm_lexp2 = fitGLM_lexp(glm_lexp2,stim_emb,spikebins,'tots');

k_mat = get_k_mat(glm_lexp);

overlap = abs(det(k_mat'*filt_mat))^(1/nfilts)/(abs(det(filt_mat'*filt_mat))^(1/(2*nfilts)) * abs(det(k_mat'*k_mat))^(1/(2*nfilts)));

overlap_stc = abs(det(stc_dims'*filt_mat))^(1/nfilts)/(abs(det(filt_mat'*filt_mat))^(1/(2*nfilts)) * abs(det(stc_dims'*stc_dims))^(1/(2*nfilts)));

kprojmat = k_mat*inv(k_mat'*k_mat)*k_mat';
proj_filts = filt_mat'*kprojmat;

sprojmat = stc_dims*inv(stc_dims'*stc_dims)*stc_dims';
proj_filts_s = filt_mat'*sprojmat;

%%
nmods = 8;
% cur_basis = get_k_mat(glm_lexp2);
cur_basis = stc_dims;
n_bvs = size(cur_basis,2);
STCcf_0 = randn(n_bvs,nmods);
%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
init_kerns = cur_basis*STCcf_0;
init_signs = [1 -1 1 -1 1 -1 1 -1];
% init_signs = ones(1,nmods);
init_betas = 4*ones(1,nmods);

defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 0;
defmod.lambda_L1x = 0; %2
defmod.lambda_dT = 0;
defmod.locLambda = 0;
defmod.SDIM = SDIM;
defmod.fsdim = SDIM;
defmod.pids = 1:SDIM;
basis = 'pix';

glm_stcb = createGLM_lexp_rot(cur_basis,STCcf_0,init_signs,init_betas,defmod,basis);
glm_stcb.image_type = '1d';
% glm_stcb.spk_nl = 'logexp';
glm_stcb.spk_nl = 'exp';
kern_output = stim_emb*cur_basis;
glm_lexp = fitGLM_lexp_rot(glm_stcb,kern_output,spikebins,'tots');

%%
nmods = 4;
% cur_basis = get_k_mat(glm_lexp2);
cur_basis = stc_dims;
n_bvs = size(cur_basis,2);
STCcf_0 = randn(n_bvs,nmods);
%normalize
for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
init_kerns = cur_basis*STCcf_0;
init_signs = [1 -1 1 -1];
% init_signs = ones(1,nmods);
init_betas = 4*ones(1,nmods);

defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 0;
defmod.lambda_L1x = 0; %2
defmod.lambda_dT = 0;
defmod.locLambda = 0;
defmod.SDIM = sdim;
defmod.fsdim = sdim;
defmod.pids = 1:sdim;
basis = 'pix';

init_nls{1} = 'pquad';
init_nls{2} = 'nquad';
init_nls{3} = 'pquad';
init_nls{4} = 'nquad';
for i = 1:nmods
    nltypes{i} = 'quad';
end
Robs = zeros(1,NT);
ftable = tabulate(spikebins);
Robs(ftable(:,1)) = ftable(:,2);

glm_stcb = createGLM0_stcb_connl(cur_basis,STCcf_0,defmod,init_nls,nltypes,'test')
glm_stcb.image_type = '1d';
glm_stcb.spk_nl = 'logexp';
for i = 1:nmods
    glm_stcb.mods(i).w = init_signs(i);
end
kern_output = stim_emb*cur_basis;
glm_lexp = fitSTCBF_connl(glm_stcb,kern_output,Robs,0);
