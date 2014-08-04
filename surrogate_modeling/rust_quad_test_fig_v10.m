cd ~/James_scripts/surrogate_modeling/
addpath(genpath('~/James_scripts'));
clear all
close all

dt = 0.01; %in s
max_rate = Inf; %in Hz
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);

NT = 200000; SDIM = 20; flen = 14;
stim = round(2*rand(NT,SDIM))-1;
xvNT = 30000; 
xvstim = round(2*rand(xvNT,SDIM))-1;

% GENERATE A FILTER
nfilts = 4;
LAMBDA = 5; %4
SIGMA = 1.5; %1.5
x = repmat(1:SDIM,flen,1);
sigma1 = repmat(SIGMA,flen,SDIM);
lambda = repmat(LAMBDA,flen,SDIM);
% amp_vec = ncx2pdf((1:flen-2)-1,4,1);
% amp_vec = ncx2pdf((1:flen-2)-1,4,0);
amp_vec = gampdf((1:flen-2)-1,3,1.3);
amp_vec = [0 0 amp_vec];
amp_vec = fliplr(amp_vec);
b = repmat(amp_vec',1,SDIM);
a = repmat(0,flen,SDIM);
% psi1 = repmat(linspace(0,5*pi,flen)',1,SDIM);
% psi1(1:5,:) = repmat(psi1(6,:),5,1);
psi1 = repmat(((1:flen)'/7).^3+4,1,SDIM);

desired_spacing = 0.69*LAMBDA;
% desired_spacing = 3.25;
xovals = 0:desired_spacing:(desired_spacing*(nfilts-1));
xovals = xovals - mean(xovals);
xovals = xovals + SDIM/2+0.5;
for i = 1:nfilts
    x0 = repmat(xovals(i),flen,SDIM);
    cur_psi = mod(psi1 + rand*2*pi,2*pi);
    temp = b.*exp(-((x-x0).^2./2./sigma1.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+a;
    filt(i,:,:) = temp/norm(temp(:));
    filt_mat(i,:) = temp(:)/norm(temp(:));
end
filt_mat = filt_mat';
c_cent = 8.5;
c_std = 5;
max_cval = 0;
c_offset = 1;
cvals = max_cval*exp(-(xovals-c_cent).^2/(2*c_std))+c_offset;
cvals = [1 1.25 1.25 1];

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
beta = 4; theta = 0;
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
cur_theta = 2.5; % 1.5 2.5
cur_beta = 4; % 5 4
p_spike = 1/cur_beta*log(1+exp(cur_beta*(g-cur_theta)));
% p_spike = (g-cur_theta).^2; p_spike(p_spike<0) = 0;
% p_spike = exp((g-1.5)/1);
xvp_spike = 1/cur_beta*log(1+exp(cur_beta*(xvg-cur_theta)));
scale_f = target_rate/(mean(p_spike)/dt);
p_spike = p_spike*scale_f;
% p_spike(p_spike/dt > 500) = 500*dt;
xvp_spike = xvp_spike*scale_f;
% xvp_spike(xvp_spike/dt > 500) = 500*dt;

figure
subplot(2,1,1)
hist(g,100)
xl = xlim();
subplot(2,1,2)
xx = linspace(xl(1),xl(2),100);
plot(xx,1/cur_beta*log(1+exp(cur_beta*(xx-cur_theta))))
xlim(xl);

figure
hist(p_spike/dt,100)
% a = input('');

%%
close all
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
% stim_proj = stim_emb;
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
colormap(gray)

stc_dims = [sta' stcs(:,1:4)];

%%
figure
for i = 1:nfilts
    temp = filt(i,:,:);
    subplot(nfilts,2,(i-1)*2+1)
    imagesc(squeeze(filt(i,:,:))*cvals(i));colormap(gray)
    set(gca,'xtick',[],'ytick',[])
    cur_inputs = coefs(i).^2*stim_emb*temp(:);
    [foy,fox] = ksdensity(cur_inputs);    
    [foys,foxs] = ksdensity(cur_inputs(spikebins));
    subplot(nfilts,2,(i-1)*2+2)
    plot(foxs,foys,'r')
    hold on
    xx = linspace(-4,4,100);
    plot(fox,foy,'k--')
    axis tight
    box off
    xlim([-4.5 4.5])
    ylim([0 0.53])
    plot(xx,xx.^2/25)
end

%% QUAD MOD
sdim = SDIM; 
stim_params.spatial_dims = 1;
stim_params.sdim = sdim;
stim_params.flen = flen;

defmod.lambda_L1x = 60; 
defmod.lambda_d2XT = 0;

nmods = 4;
init_kerns = randn(flen*sdim,nmods);
init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
init_signs = ones(1,nmods);
clear kern_types
for i = 1:nmods
    kern_types{i} = 'quad';
end
quad_mod = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
quad_mod = fitGNM_filters(quad_mod,stim_emb,spikebins,'none',[],1e-4,1e-6);
[~, ~, ~, ~, g] = getLL_GNM(quad_mod,stim_emb,spikebins,'none');
quad_mod = fitGNM_spkNL(quad_mod,g,spikebins,0);

[~,coms,peak_locs] = get_filter_coms_1d(quad_mod);
[~,ord] = sort(coms);
quad_mod.mods = quad_mod.mods(ord);

quad_xvLL = getLL_GNM(quad_mod,xvstim_emb,xvspikebins,'none')

%%
% fo = glm_quad;
% fo = adjust_all_reg(fo,'lambda_dX',0);
% fo = adjust_all_reg(fo,'lambda_dT',0);
% fo = adjust_all_reg(fo,'lambda_d2XT',40);
% fo = adjust_all_reg(fo,'lambda_L1x',60);
% fo = fitGLM_lexp(fo,stim_emb,spikebins,'tots',[],1e-4,1e-6);
% xvLL_quad_fo = getLLGLM_lexp(fo,xvstim_emb,xvspikebins,'tots')
% 
%%
sdim = SDIM; flen = 14;
stim_params.spatial_dims = 1;
stim_params.sdim = sdim;
stim_params.flen = flen;

defmod.lambda_L1x = 200; 
defmod.lambda_d2XT = 0;

nmods = 8;
init_kerns = randn(flen*sdim,nmods);
init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
init_signs = ones(1,nmods);
clear kern_types
for i = 1:nmods
    kern_types{i} = 'threshlin';
end
gnm = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
gnm = fitGNM_filters(gnm,stim_emb,spikebins,'none',[],1e-4,1e-6);
[~, ~, ~, ~, g] = getLL_GNM(gnm,stim_emb,spikebins,'none');
gnm = fitGNM_spkNL(gnm,g,spikebins,0);
gnm = fitGNM_filters(gnm,stim_emb,spikebins,'none',[],1e-4,1e-6);
gnm_xvLL = getLL_GNM(gnm,xvstim_emb,xvspikebins,'none')
[~,coms,peak_locs] = get_filter_coms_1d(gnm);
[~,ord] = sort(coms);
gnm.mods = gnm.mods(ord);

%for fitting internal NLs
gnmr = adjust_all_reg(gnm,'lnl2',1000);
gnmr = adjust_all_reg(gnmr,'nltype','uncon');
gnmr = adjust_all_reg(gnmr,'nlmon',1);
gnmr = fitGNM_internal_NLs(gnmr,stim_emb,spikebins,0,0);
gnmr = fitGNM_filters(gnmr,stim_emb,spikebins,'none',[],1e-4,1e-6);
[~, ~, ~, ~, g] = getLL_GNM(gnmr,stim_emb,spikebins,'none');
gnmr = fitGNM_spkNL(gnmr,g,spikebins,0);
gnmr_xvLL = getLL_GNM(gnmr,xvstim_emb,xvspikebins,'none')
[~,coms,peak_locs] = get_filter_coms_1d(gnmr);
[~,ord] = sort(coms);
gnmr.mods = gnmr.mods(ord);

%% 
gnmr2 = adjust_all_reg(gnmr,'lnl2',200);
gnmr2 = fitGNM_internal_NLs(gnmr2,stim_emb,spikebins,0,0);
gnmr2 = fitGNM_filters(gnmr2,stim_emb,spikebins,'none',[],1e-4,1e-6);
[~, ~, ~, ~, g] = getLL_GNM(gnmr2,stim_emb,spikebins,'none');
gnmr2 = fitGNM_spkNL(gnmr2,g,spikebins,0);
gnmr2_xvLL = getLL_GNM(gnmr2,xvstim_emb,xvspikebins,'none')
[~,coms,peak_locs] = get_filter_coms_1d(gnmr2);
[~,ord] = sort(coms);
gnmr2.mods = gnmr2.mods(ord);

%%
% cd /Users/James/James_scripts/surrogate_modeling
% save quad_mod_gnm_fig gnmr* quad* stc* sta*
%%
% quad_k_mat = get_k_mat(glm_quadf);
% lexp_k_mat = get_k_mat(glm_lexpf);
% 
% overlap_lexp = subspace_overlap2(lexp_k_mat,filt_mat);
% overlap_quad = subspace_overlap2(quad_k_mat,filt_mat);
% overlap_stc = subspace_overlap2(stc_dims,filt_mat);
% 
% overlap_lexpquad = subspace_overlap(lexp_k_mat,quad_k_mat);
% overlap_lexp = subspace_overlap(filt_mat,lexp_k_mat);
% overlap_quad = subspace_overlap(filt_mat,quad_k_mat);
% overlap_stc = subspace_overlap(stc_dims,filt_mat);
% 
% soverlap_lexp = subspace(filt_mat,lexp_k_mat);
% % soverlap_lexp2 = subspace(filt_mat,lexp2_k_mat);
% soverlap_quad = subspace(filt_mat,quad_k_mat);
% soverlap_stc = subspace(filt_mat,stc_dims);

%%
% cd ~/Documents/GNM_paper/
% save lexp_rust_examp2 glm_quad glm_lexp filt* sdim flen xvLL*

%%
gnm_k_mat = get_k_mat(gnmr(1));
quad_k_mat = get_k_mat(quad_mod(1));
[overlap_gnnm,vec_var_gnnm] = subspace_overlap2(gnm_k_mat,filt_mat);
[overlap_quad,vec_var_quad] = subspace_overlap2(quad_k_mat,filt_mat);
[overlap_stc,vec_var_stc] = subspace_overlap2(stc_dims,filt_mat);
[overlapr_gnnm,vecr_var_gnnm] = subspace_overlap2(filt_mat,gnm_k_mat);
[overlapr_quad,vecr_var_quad] = subspace_overlap2(filt_mat,quad_k_mat);
[overlapr_stc,vecr_var_stc] = subspace_overlap2(filt_mat,stc_dims);

cormat_gnm = corr(filt_mat,gnm_k_mat);
max_corrs_gnm = max(cormat_gnm);
cormat_quad = corr(filt_mat,quad_k_mat);
max_corrs_quad = max(cormat_quad);
cormat_stcs = corr(filt_mat,stc_dims);
max_corrs_stcs = max(cormat_stcs);

gnm_k_norm = bsxfun(@rdivide,gnm_k_mat,sqrt(sum(gnm_k_mat.^2)));
quad_k_norm = bsxfun(@rdivide,quad_k_mat,sqrt(sum(quad_k_mat.^2)));

clear beta_gnm
for i = 1:size(gnm_k_norm,2)
    beta_gnm(i,:)= regress(gnm_k_norm(:,i),filt_mat);
end
clear beta_quad
for i = 1:size(quad_k_norm,2)
    beta_quad(i,:)= regress(quad_k_norm(:,i),filt_mat);
end
clear beta_stc
for i = 1:size(stc_dims,2)
    beta_stc(i,:)= regress(stc_dims(:,i),filt_mat);
end
figure
for i = 1:8
subplot(4,2,i)
plot(beta_gnm(i,:),'o-')
ylim([-1 1])
xlim([1 4])
end
figure
for i = 1:size(quad_k_norm,2)
subplot(size(quad_k_norm,2),1,i)
plot(beta_quad(i,:),'o-')
ylim([-1 1])
xlim([1 4])
end
figure
for i = 1:size(stc_dims,2)
subplot(size(stc_dims,2),1,i)
plot(beta_stc(i,:),'o-')
ylim([-1 1])
xlim([1 4])
end

% figure
% plot(vecr_var_gnnm,max_corrs_gnm,'o')
% hold on
% plot(vecr_var_quad,max_corrs_quad,'ro')
% plot(vecr_var_stc,max_corrs_stcs,'ko')
% overlap_lexpquad = subspace_overlap(lexp_k_mat,quad_k_mat);
% overlap_lexp = subspace_overlap(filt_mat,lexp_k_mat);
% overlap_quad = subspace_overlap(filt_mat,quad_k_mat);
% overlap_stc = subspace_overlap(stc_dims,filt_mat);
% 
% soverlap_lexp = subspace(filt_mat,lexp_k_mat);
% % soverlap_lexp2 = subspace(filt_mat,lexp2_k_mat);
% soverlap_quad = subspace(filt_mat,quad_k_mat);
% soverlap_stc = subspace(filt_mat,stc_dims);
