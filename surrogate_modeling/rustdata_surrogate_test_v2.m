cd ~/James_scripts/surrogate_modeling/
clear all
close all

load ./rustlike_stim
dt = 0.01; %in s
max_rate = 200; %in Hz
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);


NT = size(stim,1);
% ep = 50000; stim = stim(1:ep,:);
%% GENERATE A FILTER
true_nfilts = 10;
true_x0 = [4 4 6 6 8 8 10 10 12 12];
true_sigmas = 1.5*ones(1,true_nfilts);
true_lambdas = 6*ones(1,true_nfilts);
true_amps = 2*ones(1,true_nfilts);
true_psi0 = [0 pi/2 0 pi/2 0 pi/2 0 pi/2 0 pi/2];
%initialize filts
true_filts = nan(true_nfilts,flen*SDIM);
%compute filts
x = repmat(1:SDIM,flen,1);
a = repmat(0,flen,SDIM);
for n = 1:true_nfilts
    x0 = repmat(true_x0(n),flen,SDIM);
    sigma = repmat(true_sigmas(n),flen,SDIM);
    lambda = repmat(true_lambdas(n),flen,SDIM);
    b = repmat(linspace(true_amps(n),0,flen)',1,SDIM);
    psi = repmat(linspace(true_psi0(n),true_psi0(n)+pi,flen)',1,SDIM);
    temp_filt = b.*exp(-((x-x0).^2./2./sigma.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi))+a;
    true_filts(n,:) = temp_filt(:)/norm(temp_filt(:));
end

%% CREATE TIME-EMBEDDED STIMULUS
stim_emb = makeStimRows(stim,flen);

%% FILTER STIMULUS
filt_proj_stim = stim_emb*true_filts';
%% create spike function
%pass through internal NLs
g = zeros(NT,1);
for n = 1:true_nfilts
%    g = g + logexp(filt_proj_stim(:,n),[0 1]); 

     temp = filt_proj_stim(:,n);
     temp(temp < 0) = 0;
     g = g + temp;
end
g = zscore(g);
target_rate = 20; %in Hz
target_pspike = target_rate*dt;
K0 = [1 0.4];
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
nneg = 8;
npos = 8;
spike_cond_stim = stim_emb(spikebins,:);
sta      = mean(spike_cond_stim) - mean(stim_emb);
sta = sta/norm(sta);

stvcv = cov(spike_cond_stim);  utvcv = cov(stim_emb);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

figure
subplot(3,8,1)
imagesc(reshape(sta,flen,SDIM));
for i = 1:8
    subplot(3,8,8+i)
    imagesc(reshape(stcs(:,i),flen,SDIM));
end
for i = 1:8
    subplot(3,8,16+i)
    imagesc(reshape(stcs(:,8+i),flen,SDIM));
    colormap(gray)
end

%% First, refine the STC analysis by doubling and splitting st components
used_stc_dims = [1 2 3 4 5 6 7];
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
defmod.lambda_dX = 40; %350
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
nmods = 10;
mod_signs = ones(nmods,1);
dim_signs = ones(n_bvs,1);
unused_stcs = (nmods+1):n_bvs;

clear init_vals all_filtproj all_initproj cur_LL cur_LP rotbv_mod
for r = 1:10
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
    defmod.locLambda = 100;
    defmod.lambda_dX = 10; %350
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

    [spatial_profiles, temporal_profiles, weights, mod_type, space_COM, temp_COM] = ...
        compute_mod_stats(rotbv_mod(r));
    [~,COM_ord] = sort(space_COM);
    rotbv_mod(r).mods = rotbv_mod(r).mods(COM_ord);
  
    %     rotbv_mod2(r) = fitNLHI_stcb_nopsc(glm_stcb,stim_emb,spikebins,'tots');
    cur_LL(r) = rotbv_mod(r).LL;
    cur_LP(r) = rotbv_mod(r).LP;
    init_vals(r,:,:) = STCcf_0;
    
    % mod_filts = get_pix_mat(full_glm);
    mod_filts = get_k_mat(rotbv_mod(r));
    mod_filtproj = true_filts*mod_filts;
    all_filtproj(r,:,:) = mod_filtproj;
    
    % init_mod_filts = get_pix_mat(glm_stcb);
    init_mod_filts = get_k_mat(glm_stcb);
    init_mod_filtproj = true_filts*init_mod_filts;
    all_initproj(r,:,:) = init_mod_filtproj;  
end

%% NOW FINAL ALIGNMENT OF ROTATED MODEL
[~,best_mod] = min(cur_LL);
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
defmod.nlcon = 0;
defmod.nlmon = 1;
defmod.locLambda = 100;
defmod.locSigma = 2;
defmod.maxLocPen = 5000;
defmod.lambda_dX = 30; %350
defmod.lambda_L1x = 20; %40
defmod.lambda_dT = 10;
defmod.pids = 1:SDIM;
defmod.SDIM = SDIM;
defmod.fsdim = SDIM;

for i = 1:nmods; init_nls{i} = 'threshlin'; end;
%define NL types: "uncon, lin, threshlin, quad"
for i = 1:nmods; nltypes{i} = 'threshlin'; end;
glm_stcb = createGLM2d_fullbf(basis_vecs,STCcf_0,[],[],defmod,nltypes,init_nls,basis,sprintf('test')); %initialize

for i = 1:nmods
    glm_stcb.mods(i).nly = rotbv_mod(best_mod).mods(i).nly;
    glm_stcb.mods(i).w = rotbv_mod(best_mod).mods(i).w;
end
glm_stcb.const = rotbv_mod(best_mod).const;

[glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,stim_emb);
glm_stcb.image_type = '1d';
fin_glm = fitNLHI2d_fullbf(glm_stcb,stim_emb,spikebins,'tots',3,3);
mod_filts = get_pix_mat(fin_glm);
fin_filtproj = true_filts*mod_filts;

%%

figure
subplot(3,8,1)
imagesc(reshape(sta,flen,SDIM));
for i = 1:8
    subplot(3,8,8+i)
    imagesc(reshape(stcs(:,i),flen,SDIM));
end
for i = 1:8
    subplot(3,8,16+i)
    imagesc(reshape(stcs(:,8+i),flen,SDIM));
    colormap(gray)
end

ov_max = max(true_filts(:));
figure
for i = 1:size(true_filts,1)
    subplot(3,4,i)
    imagesc(reshape(true_filts(i,:),flen,SDIM)); colormap(gray); caxis([-ov_max ov_max])
end

plotfo1d_nopsc(fin_glm,4)

%%

fin_filt_mat = get_pix_mat(fin_glm);
clear fin_filt_mat_norm true_filt_mat_norm
for i = 1:nmods
    fin_filt_mat_norm(:,i) = fin_filt_mat(:,i) - mean(fin_filt_mat(:,i));
    fin_filt_mat_norm(:,i) = fin_filt_mat_norm(:,i)/norm(fin_filt_mat_norm(:,i));
end
for i = 1:size(true_filts,1)
       true_filt_mat_norm(:,i) = true_filts(i,:) - mean(true_filts(i,:));
    true_filt_mat_norm(:,i) = true_filt_mat_norm(:,i)/norm(true_filt_mat_norm(:,i)); 
end
px = true_filts'*(true_filts*true_filts')^(-1)*true_filts;

stac = [sta' stcs];
stac_filtproj = px*stac;
fin_filtproj = px*fin_filt_mat_norm;

stac_correct_var = sum(stac_filtproj.^2);
fin_correct_var = sum(fin_filtproj.^2);

fin_true_weights = fin_filt_mat_norm'*true_filt_mat_norm;
[fin_best_weights_mod,fin_best_weights_mod_loc] = max(fin_true_weights);
[fin_best_weights_true,fin_best_weights_true_loc] = max(fin_true_weights,[],2);
stac_true_weights = stac'*true_filt_mat_norm;
[stac_best_weights_mod,stac_best_weights_mod_loc] = max(stac_true_weights);
[stac_best_weights_true,stac_best_weights_true_loc] = max(stac_true_weights,[],2);

fin_mod_weights = arrayfun(@(x) x.w,fin_glm.mods);

figure
plot(1:10,fin_best_weights_mod,'o-')
hold on
plot(1:10,stac_best_weights_mod,'ro-')
xlabel('True Filter Number','fontsize',14)
ylabel('Correlation','fontsize',14)

figure
plot(fin_correct_var,'o-')
hold on
plot(stac_correct_var,'ro-')
xlabel('Model Filter Number','fontsize',14)
ylabel('Fraction in true subspace','fontsize',14)


figure
plot(1:10,fin_mod_weights(fin_best_weights_mod_loc),'o')
line([0 11],[0 0],'color','k')
xlim([0 11])
ylim([-1 2.5])
xlabel('True Filter Number','fontsize',14)
ylabel('Best-model filter weight','fontsize',14)

figure
plot(fin_best_weights_true_loc,fin_mod_weights,'o')
xlim([0 11])
line([0 11],[0 0],'color','k')
ylim([-1 2.5])
xlabel('True Filter Number','fontsize',14)
ylabel('Best-model filter weight','fontsize',14)


%%
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

