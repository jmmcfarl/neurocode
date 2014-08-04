cd ~/James_scripts/surrogate_modeling/
clear all
close all

load ./rustlike_stim
dt = 0.01; %in s
max_rate = 200; %in Hz
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);


NT = size(stim,1);
% ep = 30000; stim = stim(1:ep,:);
%% GENERATE A FILTER
x = repmat(1:SDIM,flen,1);
x0 = repmat(11,flen,SDIM);
sigma = repmat(1.5,flen,SDIM);
lambda = repmat(6,flen,SDIM);
b = repmat(linspace(2,0,flen)',1,SDIM);a = repmat(0,flen,SDIM);
psi1 = repmat(linspace(0,pi,flen)',1,SDIM);
psi2 = psi1 + pi/2;
filt1 = b.*exp(-((x-x0).^2./2./sigma.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+a;
filt2 = b.*exp(-((x-x0).^2./2./sigma.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi2))+a;

filt1 = filt1/norm(filt1(:));
filt2 = filt2/norm(filt2(:));
%exact orthogonalization of filt2
temp1 = filt1(:); temp2 = filt2(:);
temp2n = temp2 - (temp1'*temp2)*temp1;
filt2 = reshape(temp2n,flen,SDIM);
filt1 = filt1/norm(filt1(:));
filt2 = filt2/norm(filt2(:));

x0 = repmat(5,flen,SDIM);
filt3 = b.*exp(-((x-x0).^2./2./sigma.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+a;
filt4 = b.*exp(-((x-x0).^2./2./sigma.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi2))+a;

filt3 = filt3/norm(filt3(:));
filt4 = filt4/norm(filt4(:));
%exact orthogonalization of filt2
temp3 = filt3(:); temp4 = filt4(:);
temp4n = temp4 - (temp3'*temp4)*temp3;
filt4 = reshape(temp4n,flen,SDIM);
filt3 = filt3/norm(filt3(:));
filt4 = filt4/norm(filt4(:));

f1 = X; f1(X < 0) = 0;
f2 = Y; f2(Y < 0) = 0;
% f1 = logexp(X,[0 1]);
% f2 = logexp(Y,[0 1]);
gfun = f1-f2;
figure
contour(X,Y,gfun,30)
hold on
line([0 3*filt1(1)],[0 3*filt1(2)],'color','r')
line([0 3*filt2(1)],[0 3*filt2(2)],'color','b')

%% CREATE TIME-EMBEDDED STIMULUS
stim_emb = makeStimRows(stim,flen);

%% FILTER STIMULUS
filt_stim1 = stim_emb*filt1(:);
filt_stim2 = stim_emb*filt2(:);
filt_stim3 = stim_emb*filt3(:);
filt_stim4 = stim_emb*filt4(:);
filt_proj_stim = [filt_stim1 filt_stim2 filt_stim3 filt_stim4];
%% create spike function
%pass through internal NLs
% filt_stim1 = logexp(filt_stim1,[0 1]);
% filt_stim2 = logexp(filt_stim2,[0 1]);
% filt_stim3 = logexp(filt_stim3,[0 1]);
% filt_stim4 = logexp(filt_stim4,[0 1]);
filt_stim1(filt_stim1 < 0) = 0;
filt_stim2(filt_stim2 < 0) = 0;
filt_stim3(filt_stim3 < 0) = 0;
filt_stim4(filt_stim4 < 0) = 0;
% filt_stim1 = filt_stim1.^2;
% filt_stim2 = filt_stim2.^2;
% g = filt_stim1 + filt_stim2;
g = filt_stim1 + filt_stim2 + filt_stim3 + filt_stim4;
g = zscore(g);
% g = g - mean(g);
target_rate = 20; %in Hz
target_pspike = target_rate*dt;
% K0 = [1 0.4];
K0 = [1 0.5];
% Kfit = fmincon(@(K) abs(mean(logexp(g,K)) - target_pspike),K0,[],[],[],[],[-10 .01],[10 10]);
Kfit = K0;
p_spike = logexp(g,Kfit);
p_spike(p_spike > max_rate*dt) = max_rate*dt;

figure
[y,x] = ksdensity(g);
plot(x,y); hold on
yl = ylim();
tx = -10:.02:10;
plot(tx,logexp(tx,Kfit),'k')
ylim(yl)

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

%%
nneg = 4;
npos = 4;
spike_cond_stim = stim_emb(spikebins,:);
sta      = mean(spike_cond_stim) - mean(stim_emb);
sta = sta/norm(sta);
stvcv = cov(spike_cond_stim);  utvcv = cov(stim_emb);
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

% pause
% close all
%%
minxy = [-4 -4];
maxxy = [4 4];
[bandwidth2_s1,density2_s1,X2,Y2]=kde2d(filt_proj_stim(:,1:2),2^8,minxy,maxxy);
[bandwidth2_s2,density2_s2,X2,Y2]=kde2d(filt_proj_stim(:,3:4),2^8,minxy,maxxy);

[bandwidth_s1,density_s1,X,Y]=kde2d(filt_proj_stim(spikebins,1:2),2^8,minxy,maxxy,bandwidth2_s1);
[bandwidth_s2,density_s2,X,Y]=kde2d(filt_proj_stim(spikebins,3:4),2^8,minxy,maxxy,bandwidth2_s2);


%%
filt_long = [filt1(:) filt2(:) filt3(:) filt4(:)]';
sta_filtproj = filt_long*sta';
stc_filtproj = filt_long*stcs;

eps = 1e-3;
cond_dens_s1 = density_s1./density2_s1;
cond_dens_s1(density2_s1 < eps) = eps;
cond_dens_s2 = density_s2./density2_s2;
cond_dens_s2(density2_s2 < eps) = eps;
outside = find(X.^2+Y.^2>9);
cond_dens_s1(outside) = nan;
density_s1(outside) = nan;
cond_dens_s2(outside) = nan;
density_s2(outside) = nan;

% figure
% contourf(X,Y,cond_dens,30)
% % surf(X,Y,cond_dens)
% hold on
% circle([0 0],3,100,'k');
% line(3*[0 stc_filtproj(1,1)],3*[0 stc_filtproj(2,1)],'color','r','linewidth',2)
% line(3*[0 stc_filtproj(1,2)],3*[0 stc_filtproj(2,2)],'color','w','linewidth',2)
% line(3*[0 sta_filtproj(1)],3*[0 sta_filtproj(2)],'color','k','linewidth',2)
% vline(0,'w'); hline(0,'w')

% figure
% contourf(X,Y,density,30)
% hold on
% circle([0 0],3,100,'k');
% line(3*[0 stc_filtproj(1,1)],3*[0 stc_filtproj(2,1)],'color','r','linewidth',2)
% line(3*[0 stc_filtproj(1,2)],3*[0 stc_filtproj(2,2)],'color','w','linewidth',2)
% line(3*[0 sta_filtproj(1)],3*[0 sta_filtproj(2)],'color','k','linewidth',2)
% vline(0,'w'); hline(0,'w')

%% First, refine the STC analysis by doubling and splitting st components
used_stc_dims = [1 2 3 4 6];
STCbvs = [sta' stcs];
STCbvs = STCbvs(:,used_stc_dims);
STCbvs = [STCbvs -STCbvs(:,1:end)]; %make copies of STC comps
Nstcbvs = size(STCbvs,2);
nmods = Nstcbvs;
basis = 'pix';
STCcf_0 = eye(Nstcbvs);
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
flen = 14;
%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 100;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 1;
defmod.locLambda = 0;
defmod.lambda_dX = 100; %350
defmod.lambda_L1x = 10; %40
defmod.lambda_dT = 10;
defmod.pids = 1:SDIM;
defmod.SDIM = SDIM;
defmod.fsdim = SDIM;

clear init_nls nltypes
for i = 1:nmods; init_nls{i} = 'threshlin'; end;
%define NL types: "uncon, lin, threshlin, quad"
for i = 1:nmods; nltypes{i} = 'threshlin'; end;
glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize

[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,stim_emb);
glm_stcb.image_type = '1d';
full_glm = fitNLHI2d_fullbf(glm_stcb,stim_emb,spikebins,'tots',2);
% mod_filts = get_pix_mat(full_glm);
% fin_filtproj = filt_long*mod_filts;

plotfo1d_nopsc(glm_stcb,4)
plotfo1d_nopsc(full_glm,4)

%% NOW FIND BEST OBLIQUE ROTATION WITHIN THE NEW SUBSPACE
basis_vecs = get_pix_mat(full_glm);
% basis_vecs = [sta' stcs];
% basis_vecs = basis_vecs(:,used_stc_dims);
n_bvs = size(basis_vecs,2);
nmods = 4;
mod_signs = ones(nmods,1);
dim_signs = ones(n_bvs,1);
unused_stcs = (nmods+1):n_bvs;

clear init_vals all_filtproj all_initproj cur_LL cur_LP
for r = 1:20
    r
    STCcf_0 = randn(n_bvs,nmods);
    
    %normalize
    for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    
    cd ~/Data/blanche/rec_75/matlabdata/
    load stdparsRec75.mat
    flen = 14;
    %initialize model
    defmod.h(1:end-1) = []; %eliminate PSC
    defmod.lnl = 0;
    defmod.lh = 0;
    defmod.lnl2 = 100;
    defmod.lh2 = 0;
    defmod.nlcon = 0;
    defmod.nlmon = 1;
    defmod.locLambda = 200;
    defmod.lambda_dX = 20; %350
    defmod.lambda_L1x = 10; %40
    defmod.lambda_dT = 10;
    defmod.pids = 1:SDIM;
    defmod.SDIM = SDIM;
    defmod.fsdim = SDIM;
    
    %define NL initializations: "lin, threshlin, pquad, nquad"
    clear init_nls nltypes
    for i = 1:nmods; init_nls{i} = 'threshlin'; end;
    %define NL types: "uncon, lin, threshlin, quad"
    for i = 1:nmods; nltypes{i} = 'threshlin'; end;    
    
    glm_stcb = createGLM0_stcb(basis_vecs,STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
    rotbv_mod(r) = fitNLHI_stcb_nonlpsc(glm_stcb,stim_emb,spikebins,'tots',3);
%     rotbv_mod(r) = fitNLHI_stcb_nopsc(glm_stcb,stim_emb,spikebins,'tots',3);
    cur_LL(r) = rotbv_mod(r).LL;
    cur_LP(r) = rotbv_mod(r).LP;
    init_vals(r,:,:) = STCcf_0;
    
    % mod_filts = get_pix_mat(full_glm);
    mod_filts = get_k_mat(rotbv_mod(r));
    mod_filtproj = filt_long*mod_filts;
    all_filtproj(r,:,:) = mod_filtproj;
    
    % init_mod_filts = get_pix_mat(glm_stcb);
    init_mod_filts = get_k_mat(glm_stcb);
    init_mod_filtproj = filt_long*init_mod_filts;
    all_initproj(r,:,:) = init_mod_filtproj;  
end

%% NOW FINAL ALIGNMENT OF ROTATED MODEL
[~,best_mod] = min(cur_LP);
STCcf_0 = nan(n_bvs,nmods);
for i = 1:nmods
    STCcf_0(:,i) = rotbv_mod(best_mod).mods(i).STCcf;
end
basis = 'pix';
load stdparsRec75.mat
flen = 14;
%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 100;
defmod.lh2 = 0;
defmod.nlcon = 1;
defmod.nlmon = 1;
defmod.locLambda = 500;
defmod.locSigma = 2;
defmod.maxLocPen = 5000;
defmod.lambda_dX = 50; %350
defmod.lambda_L1x = 0; %40
defmod.lambda_dT = 10;
defmod.pids = 1:SDIM;
defmod.SDIM = SDIM;
defmod.fsdim = SDIM;

for i = 1:nmods; init_nls{i} = 'threshlin'; end;
%define NL types: "uncon, lin, threshlin, quad"
for i = 1:nmods; nltypes{i} = 'threshlin'; end;
glm_stcb = createGLM2d_fullbf(basis_vecs,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize

%copy over internal NLs from previous fits
for i = 1:nmods
    glm_stcb.mods(i).nly = rotbv_mod(best_mod).mods(i).nly;
    glm_stcb.mods(i).w = rotbv_mod(best_mod).mods(i).w;
end
glm_stcb.const = rotbv_mod(best_mod).const;

[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,stim_emb);
glm_stcb.image_type = '1d';
fin_glm = fitNLHI2d_fullbf(glm_stcb,stim_emb,spikebins,'tots',8,3);
mod_filts = get_pix_mat(fin_glm);
fin_filtproj = filt_long*mod_filts;
plotfo1d_nopsc(fin_glm,4)

%%
figure
subplot(2,1,1)
contourf(X,Y,cond_dens_s1,30)
hold on
circle([0 0],3,100,'k');
line(3*[0 fin_filtproj(1,1)],3*[0 fin_filtproj(2,1)],'color','r','linewidth',2)
line(3*[0 fin_filtproj(1,2)],3*[0 fin_filtproj(2,2)],'color','w','linewidth',2)
line(3*[0 fin_filtproj(1,3)],3*[0 fin_filtproj(2,3)],'color','g','linewidth',2)
line(3*[0 fin_filtproj(1,4)],3*[0 fin_filtproj(2,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')
subplot(2,1,2)
contourf(X,Y,cond_dens_s2,30)
hold on
circle([0 0],3,100,'k');
line(3*[0 fin_filtproj(3,1)],3*[0 fin_filtproj(4,1)],'color','r','linewidth',2)
line(3*[0 fin_filtproj(3,2)],3*[0 fin_filtproj(4,2)],'color','w','linewidth',2)
line(3*[0 fin_filtproj(3,3)],3*[0 fin_filtproj(4,3)],'color','g','linewidth',2)
line(3*[0 fin_filtproj(3,4)],3*[0 fin_filtproj(4,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')

figure
subplot(2,1,1)
contourf(X,Y,cond_dens_s1,30)
% surf(X,Y,cond_dens)
hold on
circle([0 0],3,100,'k');
line(3*[0 all_filtproj(best_mod,1,1)],3*[0 all_filtproj(best_mod,2,1)],'color','r','linewidth',2)
line(3*[0 all_filtproj(best_mod,1,2)],3*[0 all_filtproj(best_mod,2,2)],'color','w','linewidth',2)
line(3*[0 all_filtproj(best_mod,1,3)],3*[0 all_filtproj(best_mod,2,3)],'color','g','linewidth',2)
line(3*[0 all_filtproj(best_mod,1,4)],3*[0 all_filtproj(best_mod,2,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')
subplot(2,1,2)
contourf(X,Y,cond_dens_s2,30)
% surf(X,Y,cond_dens)
hold on
circle([0 0],3,100,'k');
line(3*[0 all_filtproj(best_mod,3,1)],3*[0 all_filtproj(best_mod,4,1)],'color','r','linewidth',2)
line(3*[0 all_filtproj(best_mod,3,2)],3*[0 all_filtproj(best_mod,4,2)],'color','w','linewidth',2)
line(3*[0 all_filtproj(best_mod,3,3)],3*[0 all_filtproj(best_mod,4,3)],'color','g','linewidth',2)
line(3*[0 all_filtproj(best_mod,3,4)],3*[0 all_filtproj(best_mod,4,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')

[~,worst_model] = max(cur_LL);
figure
subplot(2,1,1)
contourf(X,Y,cond_dens_s1,30)
% surf(X,Y,cond_dens)
hold on
circle([0 0],3,100,'k');
line(3*[0 all_filtproj(worst_model,1,1)],3*[0 all_filtproj(worst_model,2,1)],'color','r','linewidth',2)
line(3*[0 all_filtproj(worst_model,1,2)],3*[0 all_filtproj(worst_model,2,2)],'color','w','linewidth',2)
line(3*[0 all_filtproj(worst_model,1,3)],3*[0 all_filtproj(worst_model,2,3)],'color','g','linewidth',2)
line(3*[0 all_filtproj(worst_model,1,4)],3*[0 all_filtproj(worst_model,2,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')
subplot(2,1,2)
contourf(X,Y,cond_dens_s2,30)
% surf(X,Y,cond_dens)
hold on
circle([0 0],3,100,'k');
line(3*[0 all_filtproj(worst_model,3,1)],3*[0 all_filtproj(worst_model,4,1)],'color','r','linewidth',2)
line(3*[0 all_filtproj(worst_model,3,2)],3*[0 all_filtproj(worst_model,4,2)],'color','w','linewidth',2)
line(3*[0 all_filtproj(worst_model,3,3)],3*[0 all_filtproj(worst_model,4,3)],'color','g','linewidth',2)
line(3*[0 all_filtproj(worst_model,3,4)],3*[0 all_filtproj(worst_model,4,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')


figure
subplot(2,1,1)
contourf(X,Y,cond_dens_s1,30)
% surf(X,Y,cond_dens)
hold on
circle([0 0],3,100,'k');
line(3*[0 stc_filtproj(1,1)],3*[0 stc_filtproj(2,1)],'color','r','linewidth',2)
line(3*[0 stc_filtproj(1,2)],3*[0 stc_filtproj(2,2)],'color','w','linewidth',2)
line(3*[0 stc_filtproj(1,3)],3*[0 stc_filtproj(2,3)],'color','g','linewidth',2)
line(3*[0 stc_filtproj(1,5)],3*[0 stc_filtproj(2,5)],'color','c','linewidth',2)
line(3*[0 sta_filtproj(1)],3*[0 sta_filtproj(2)],'color','k','linewidth',2)
vline(0,'w'); hline(0,'w')
subplot(2,1,2)
contourf(X,Y,cond_dens_s2,30)
% surf(X,Y,cond_dens)
hold on
circle([0 0],3,100,'k');
line(3*[0 stc_filtproj(3,1)],3*[0 stc_filtproj(4,1)],'color','r','linewidth',2)
line(3*[0 stc_filtproj(3,2)],3*[0 stc_filtproj(4,2)],'color','w','linewidth',2)
line(3*[0 stc_filtproj(3,3)],3*[0 stc_filtproj(4,3)],'color','g','linewidth',2)
line(3*[0 stc_filtproj(3,5)],3*[0 stc_filtproj(4,5)],'color','c','linewidth',2)
line(3*[0 sta_filtproj(3)],3*[0 sta_filtproj(4)],'color','b','linewidth',2)
vline(0,'w'); hline(0,'w')


% figure
% subplot(2,2,1)
% imagesc(reshape(mod_filts(:,1),flen,SDIM))
% subplot(2,2,2)
% imagesc(reshape(mod_filts(:,2),flen,SDIM));
% subplot(2,2,3)
% imagesc(reshape(mod_filts(:,3),flen,SDIM))
% subplot(2,2,4)
% imagesc(reshape(mod_filts(:,4),flen,SDIM));
% colormap(gray)

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

figure
subplot(2,2,1)
imagesc(filt1); colormap(gray); caxis([-0.3 0.3])
subplot(2,2,2)
imagesc(filt2); colormap(gray); caxis([-0.3 0.3])
subplot(2,2,3)
imagesc(filt3); colormap(gray); caxis([-0.3 0.3])
subplot(2,2,4)
imagesc(filt4); colormap(gray); caxis([-0.3 0.3])

%%
figure
subplot(2,1,1)
contourf(X,Y,cond_dens_s1,30)
hold on
circle([0 0],3,100,'k');
line(3*[0 fin_filtproj(1,1)],3*[0 fin_filtproj(2,1)],'color','r','linewidth',2)
line(3*[0 fin_filtproj(1,2)],3*[0 fin_filtproj(2,2)],'color','w','linewidth',2)
line(3*[0 fin_filtproj(1,3)],3*[0 fin_filtproj(2,3)],'color','g','linewidth',2)
line(3*[0 fin_filtproj(1,4)],3*[0 fin_filtproj(2,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')
subplot(2,1,2)
contourf(X,Y,cond_dens_s2,30)
hold on
circle([0 0],3,100,'k');
line(3*[0 fin_filtproj(3,1)],3*[0 fin_filtproj(4,1)],'color','r','linewidth',2)
line(3*[0 fin_filtproj(3,2)],3*[0 fin_filtproj(4,2)],'color','w','linewidth',2)
line(3*[0 fin_filtproj(3,3)],3*[0 fin_filtproj(4,3)],'color','g','linewidth',2)
line(3*[0 fin_filtproj(3,4)],3*[0 fin_filtproj(4,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')

f1 = logexp(X,[0 1]);
f2 = logexp(Y,[0 1]);
gfun = f1+f2;
gfun(X.^2+Y.^2 > 9) = nan;
figure
contourf(X,Y,gfun,30)
hold on
line([0 3*filt1(1)],[0 3*filt1(2)],'color','r')
line([0 3*filt2(1)],[0 3*filt2(2)],'color','b')

figure
clear spl_mod_out spl_out
subplot(2,1,1)
for i = 1:nmods
temp = fin_glm.mods(i).w*nlin_proc_stim((fin_filtproj(1,i)*X + fin_filtproj(2,i)*Y),fin_glm.mods(i).nly,fin_glm.mods(i).nlx);
spl_mod_out(i,:,:) = temp;
end
spl_out = squeeze(sum(spl_mod_out));
spl_out(X.^2+Y.^2 > 9) = nan;
contourf(X,Y,spl_out,30), hold on
circle([0 0],3,100,'k');
subplot(2,1,2)
for i = 1:nmods
temp = fin_glm.mods(i).w*nlin_proc_stim((fin_filtproj(3,i)*X + fin_filtproj(4,i)*Y),fin_glm.mods(i).nly,fin_glm.mods(i).nlx);
spl_mod_out(i,:,:) = temp;
end
spl_out = squeeze(sum(spl_mod_out));
spl_out(X.^2+Y.^2 > 9) = nan;
contourf(X,Y,spl_out,30), hold on
circle([0 0],3,100,'k');

%%
fin_filt_mat = get_pix_mat(fin_glm);
for i = 1:nmods
    fin_filt_mat_norm(:,i) = fin_filt_mat(:,i) - mean(fin_filt_mat(:,i));
    fin_filt_mat_norm(:,i) = fin_filt_mat_norm(:,i)/norm(fin_filt_mat_norm(:,i));
    true_filt_mat_norm(:,i) = filt_long(i,:) - mean(filt_long(i,:));
    true_filt_mat_norm(:,i) = true_filt_mat_norm(:,i)/norm(true_filt_mat_norm(:,i));
end
px = filt_long'*(filt_long*filt_long')^(-1)*filt_long;

stac = [sta' stcs];
stac_filtproj = px*stac;
fin_filtproj = px*fin_filt_mat_norm;

stac_correct_var = sum(stac_filtproj.^2);
fin_correct_var = sum(fin_filtproj.^2);

fin_true_weights = fin_filt_mat_norm'*true_filt_mat_norm;
fin_best_weights = max(fin_true_weights);
fin_best_weights_true = max(fin_true_weights,[],2);
stac_true_weights = stac'*true_filt_mat_norm;
stac_best_weights = max(stac_true_weights);
stac_best_weights_true = max(stac_true_weights,[],2);

figure
plot(stac_correct_var,'ro-')
hold on
plot(fin_correct_var,'o-')
xlabel('Filter number','fontsize',14)
ylabel('Fraction in correct subspace','fontsize',14)

figure
plot(fin_best_weights,'o-')
hold on
plot(stac_best_weights,'ro-')
xlabel('Filter number','fontsize',14)
ylabel('Best filter overlap','fontsize',14)

figure
plot(fin_best_weights_true,'o-')
hold on
plot(stac_best_weights_true,'ro-')
xlabel('Filter number','fontsize',14)
ylabel('Best filter overlap','fontsize',14)

%%
for r = 1:20
subplot(2,1,1)
contourf(X,Y,cond_dens_s1,30)
% surf(X,Y,cond_dens)
hold on
circle([0 0],3,100,'k');
line(3*[0 all_filtproj(r,1,1)],3*[0 all_filtproj(r,2,1)],'color','r','linewidth',2)
line(3*[0 all_filtproj(r,1,2)],3*[0 all_filtproj(r,2,2)],'color','w','linewidth',2)
line(3*[0 all_filtproj(r,1,3)],3*[0 all_filtproj(r,2,3)],'color','g','linewidth',2)
line(3*[0 all_filtproj(r,1,4)],3*[0 all_filtproj(r,2,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')
subplot(2,1,2)
contourf(X,Y,cond_dens_s2,30)
% surf(X,Y,cond_dens)
hold on
circle([0 0],3,100,'k');
line(3*[0 all_filtproj(r,3,1)],3*[0 all_filtproj(r,4,1)],'color','r','linewidth',2)
line(3*[0 all_filtproj(r,3,2)],3*[0 all_filtproj(r,4,2)],'color','w','linewidth',2)
line(3*[0 all_filtproj(r,3,3)],3*[0 all_filtproj(r,4,3)],'color','g','linewidth',2)
line(3*[0 all_filtproj(r,3,4)],3*[0 all_filtproj(r,4,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')


pause
clf
end
% 
% STCbvs = [sta' stcs]; 
 
% STCbvs = STCbvs(:,[1 2 3 4]);
% Nstcbvs = size(STCbvs,2);
% nmods = 4;
% mod_signs = ones(nmods,1);
% dim_signs = ones(Nstcbvs,1);
% unused_stcs = (nmods+1):Nstcbvs;
% %initialize on STC dims
% STCcf_0 = eye(Nstcbvs);
% STCcf_0(:,unused_stcs) = [];
% %STCcf_0 = [0 0;STCcf_0]; %STCcf_0(end,:) = [];
% 
% for r = 1:10
%     r
%     STCcf_0 = randn(Nstcbvs,nmods);
%     
%     %normalize
%     for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
%     
%     cd ~/Data/blanche/rec_75/matlabdata/
%     load stdparsRec75.mat
%     flen = 14;
%     %initialize model
%     defmod.h(1:end-1) = []; %eliminate PSC
%     defmod.lnl = 0;
%     defmod.lh = 0;
%     defmod.lnl2 = 0;
%     defmod.lh2 = 0;
%     defmod.nlcon = 0;
%     defmod.nlmon = 0;
%     defmod.locLambda = 0;
%     defmod.lambda_dX = 100; %350
%     defmod.lambda_L1x = 10; %40
%     defmod.lambda_dT = 50;
%     defmod.pids = 1:SDIM;
%     defmod.SDIM = SDIM;
%     defmod.fsdim = SDIM;
%     
%     %define NL initializations: "lin, threshlin, pquad, nquad"
%     clear init_nls nltypes
%     for i = 1:nmods; init_nls{i} = 'threshlin'; end;
%     %define NL types: "uncon, lin, threshlin, quad"
%     for i = 1:nmods; nltypes{i} = 'threshlin'; end;
%     
%     % basis = 'pix';
%     % % basis = 'white';
%     % glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
%     %
%     % % glm_stcb.mods(1).pix = filt_long(1,:)'; glm_stcb.mods(1).k = filt_long(1,:)';
%     % % glm_stcb.mods(2).pix = filt_long(2,:)'; glm_stcb.mods(2).k = filt_long(2,:)';
%     %
%     % [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,stim_emb);
%     % glm_stcb.image_type = '1d';
%     % full_glm = fitNLHI2d_fullbf(glm_stcb,stim_emb,spikebins,'tots');
%     % cur_LL(r) = full_glm.LL;
%     
%     
%     glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
%     % glm_stcb.mods(1).k = filt_long(1,:)';
%     % glm_stcb.mods(2).k = filt_long(2,:)';
%     stc_posneg_mod(r) = fitNLHI_stcb_nonlpsc(glm_stcb,stim_emb,spikebins,'tots');
%     % full_glm = stc_posneg_mod;
%     cur_LL(r) = stc_posneg_mod(r).LL;
%     
%     % mod_filts = get_pix_mat(full_glm);
%     mod_filts = get_k_mat(stc_posneg_mod(r));
%     mod_filtproj = filt_long*mod_filts;
%     all_filtproj(r,:,:) = mod_filtproj;
%     
%     % init_mod_filts = get_pix_mat(glm_stcb);
%     init_mod_filts = get_k_mat(glm_stcb);
%     init_mod_filtproj = filt_long*init_mod_filts;
%     all_initproj(r,:,:) = init_mod_filtproj;
%     
% end
% 
% %%
% [~,best_mod] = min(cur_LL);
% STCcf_0 = nan(Nstcbvs,nmods);
% for i = 1:nmods
%     STCcf_0(:,i) = stc_posneg_mod(best_mod).mods(i).STCcf;
% end
% basis = 'pix';
% load stdparsRec75.mat
% flen = 14;
% %initialize model
% defmod.h(1:end-1) = []; %eliminate PSC
% defmod.lnl = 0;
% defmod.lh = 0;
% defmod.lnl2 = 0;
% defmod.lh2 = 0;
% defmod.nlcon = 0;
% defmod.nlmon = 0;
% defmod.locLambda = 0;
% defmod.lambda_dX = 100; %350
% defmod.lambda_L1x = 10; %40
% defmod.lambda_dT = 10;
% defmod.pids = 1:SDIM;
% defmod.SDIM = SDIM;
% defmod.fsdim = SDIM;
% 
% glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
% 
% [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,stim_emb);
% glm_stcb.image_type = '1d';
% full_glm = fitNLHI2d_fullbf(glm_stcb,stim_emb,spikebins,'tots');
% mod_filts = get_pix_mat(full_glm);
% fin_filtproj = filt_long*mod_filts;

