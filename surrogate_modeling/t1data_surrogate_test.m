%% T1 data mat
clear all;
addpath('~/James_scripts/GLM')
addpath('~/James_scripts/GLM/t1')
addpath('~/Timm/rust/SparseFilterSelection/')
addpath('~/Timm/MatlabRepository/')
addpath('~/Timm/tim1/functions/')

cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat

%load z-scored PC data
% cd ~/Data/blanche/rec_76/matlabdata/
% load spks7576-NS_d2;
cd ~/Data/blanche/matlabdata/
load ALL_PCScores_ds2;

load ALL_Xmat_ds2

ncomps = 1536;
compids   = 1:ncomps;
scorevars = scorevars(compids);
coefs = coefs(:,compids);
scores = scores(:,compids);
pix_conv_mat = diag(sqrt(scorevars))*coefs';
kern_conv_mat = coefs*diag(1./sqrt(scorevars));
fsdim = 256;
SDIM = 16;
flen = 6;
pids = 1:SDIM^2;

%% generate surrogate filters

x = repmat(1:SDIM,[flen 1 SDIM]);
y = permute(x,[1 3 2]);

x0 = permute(repmat(11,[SDIM SDIM flen]),[3 1 2]);
y0 = permute(repmat(11,[SDIM SDIM flen]),[3 1 2]);
sigmax = permute(repmat(1,[SDIM SDIM flen]),[3 1 2]);
sigmay = permute(repmat(1.5,[SDIM SDIM flen]),[3 1 2]);
lambda = permute(repmat(10,[SDIM SDIM flen]),[3 1 2]);
theta = permute(repmat(pi/4,[SDIM SDIM flen]),[3 1 2]);
b = repmat(linspace(2,0,flen)',[1 SDIM SDIM]);
psi = repmat(linspace(0,pi,flen)',[1 SDIM SDIM]);

psi2 = psi+pi/2;
x02 = permute(repmat(11,[SDIM SDIM flen]),[3 1 2]);
y02 = permute(repmat(5,[SDIM SDIM flen]),[3 1 2]);

xp = (x-x0).*cos(theta)+(y-y0).*sin(theta);
yp = -(x-x0).*sin(theta)+(y-y0).*cos(theta);
xp2 = (x-x02).*cos(theta)+(y-y02).*sin(theta);
yp2 = -(x-x02).*sin(theta)+(y-y02).*cos(theta);

filt1 = b.*(exp(-(xp.^2./2./sigmax.^2 + yp.^2./2./sigmay.^2)) .* (cos(2*pi*xp./lambda+psi)));
filt1 = filt1/norm(filt1(:));
filt2 = b.*(exp(-(xp2.^2./2./sigmax.^2 + yp2.^2./2./sigmay.^2)) .* (cos(2*pi*xp2./lambda+psi2)));
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


% for i = 1:6
%     subplot(2,1,1)
% imagesc(squeeze(filt1(i,:,:)))
% subplot(2,1,2)
% imagesc(squeeze(filt2(i,:,:)))
% pause
% clf
% end

%%
% filt_long = [filt1(:) filt2(:) filt3(:) filt4(:)];
filt_long = [filt1(:) filt2(:)];
% filt_white = filt_long'*kern_conv_mat;
% filt_proj_stim = scores*filt_white';
filt_proj_stim = X*filt_long;
%%
filt_stim1 = filt_proj_stim(:,1);
filt_stim1(filt_stim1 < 0) = 0;

filt_stim2 = filt_proj_stim(:,2);
filt_stim2(filt_stim2 < 0) = 0;

% filt_stim3 = filt_proj_stim(:,3);
% filt_stim3(filt_stim3 < 0) = 0;
% 
% filt_stim4 = filt_proj_stim(:,4);
% filt_stim4(filt_stim4 < 0) = 0;

% g = filt_stim1 + filt_stim2 + filt_stim3 + filt_stim4;
g = filt_stim1 + filt_stim2;
g = zscore(g);

target_rate = 20; %in Hz
max_rate = 200;
target_pspike = target_rate*dt;
% K0 = [1 0.4];
K0 = [0.4 0.2];
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

%%
nneg = 8;
npos = 8;

% ncomps    = SDIM^2; compids   = 1:ncomps; 
% WX        = scores(:,compids)*diag(1./sqrt(scorevars(compids))); 
% WS        = WX(spikebins,:); 
% rsta      = mean(WS) - mean(WX);   
% stvcv = cov(WS);  utvcv = cov(WX);
% [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
% stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

spike_cond_stim = X(spikebins,:);
sta      = mean(spike_cond_stim) - mean(X);
sta = sta/norm(sta);
stvcv =	cov(spike_cond_stim);  utvcv = cov(X);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
  
% kimages = [rsta',stcs]'*diag(sqrt(scorevars(compids)))*coefs(:,compids)'; 
% kimages = [sta',stcs]'*diag(sqrt(scorevars))*coefs';
kimages = [sta',stcs]';
figure
plotfilterbank(kimages',SDIM,pids)



%% First, refine the STC analysis by doubling and splitting st components
% used_stc_dims = [1:4];
% STCbvs = [sta' stcs];
% STCbvs = STCbvs(:,used_stc_dims);
% % STCbvs = [STCbvs -STCbvs(:,1:end)]; %make copies of STC comps
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
% defmod.nlmon = 0;
% defmod.locLambda = 0;
% defmod.locSigma = 0;
% defmod.maxLocPen = 0;
% defmod.lambda_dX = 500; %350
% defmod.lambda_L1x = 1; %40
% defmod.lambda_dT = 10;
% defmod.SDIM = SDIM;
% defmod.fsdim = fsdim;
% defmod.pids = 1:fsdim;
% 
% clear init_nls nltypes
% for i = 1:nmods; init_nls{i} = 'pquad'; end;
% %define NL types: "uncon, lin, threshlin, quad"
% for i = 1:nmods; nltypes{i} = 'uncon'; end;
% glm_stcb = createGLM2d_fullbf(STCbvs,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
% % % randomize initial filters
% % for i = 1:length(glm_stcb.mods)
% %    cur_rand_filt = rand(ncomps,1); 
% %    glm_stcb.mods(i).k = cur_rand_filt;
% %    glm_stcb.mods(i).pix = pix_conv_mat'*cur_rand_filt;
% % end
% % % 
% [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,scores);
% glm_stcb.image_type = '2d';
% full_glm = fitNLHI2d_fullbf(glm_stcb,scores,spikebins,'tots',2);
% 
% %% First, refine the STC analysis by doubling and splitting st components
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
% %% NOW FIND BEST OBLIQUE ROTATION WITHIN THE NEW SUBSPACE
% basis_vecs = get_k_mat(full_glm);
% n_bvs = size(basis_vecs,2);
% nmods = 2;
% mod_signs = ones(nmods,1);
% dim_signs = ones(n_bvs,1);
% unused_stcs = (nmods+1):n_bvs;
% 
% flen = 6;
% %initialize model
% defmod.h(1:end-1) = []; %eliminate PSC
% defmod.lnl = 0;
% defmod.lh = 0;
% defmod.lnl2 = 100;
% defmod.lh2 = 0;
% defmod.nlcon = 0;
% defmod.nlmon = 1;
% defmod.locLambda = 500;
% defmod.lambda_dX = 0; %350
% defmod.lambda_L1x = 0; %40
% defmod.lambda_dT = 10;
% defmod.pids = 1:fsdim;
% defmod.SDIM = SDIM;
% defmod.fsdim = fsdim;
% 
% %define NL initializations: "lin, threshlin, pquad, nquad"
% clear init_nls nltypes
% for i = 1:nmods; init_nls{i} = 'threshlin'; end;
% %define NL types: "uncon, lin, threshlin, quad"
% for i = 1:nmods; nltypes{i} = 'threshlin'; end;
% 
% Nstcbf = size(basis_vecs,2);
% klen = size(pix_conv_mat,2);
% flen = klen/fsdim;
% kern_output = scores*basis_vecs;
% 
% %determine distribution of random interpoint distances
% rand_reps = 500;
% init_vals = zeros(rand_reps,n_bvs,nmods);
% for r = 1:rand_reps
%     % compute average separation between NN initial points
%     STCcf_0 = randn(n_bvs,nmods);
%     %normalize
%     for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
%     init_vals(r,:,:) = STCcf_0;
% end
% cur_dists = zeros(1,rand_reps*(rand_reps-1)/2);
% for nn = 1:nmods
%     cur_dists = cur_dists + pdist(init_vals(:,:,nn),'cosine');
% end
% rdmean = 0.75*mean(cur_dists);
% rdscale = 2*std(cur_dists);
% % figure
% % ksdensity(cur_dists)
% % hold on
% % tt = linspace(0,10,100);
% % plot(tt,normcdf(tt,rdmean,rdscale),'k')
% 
% clear init_vals all_filtproj all_initproj cur_LL cur_LP fin_vals dist_trav *_lp rotbv_mod
% max_reps = 300;
% min_reps = 10;
% min_LM_fract = 0.95;
% eps = 0.002;
% cur_reps = 0;
% used_cfs = [];
% LL_vals = [];
% LP_vals = [];
% smallest_dists = [];
% is_optimized = [];
% rotbv_mod = [];
% % while cur_reps < max_reps
% for r = 1:max_reps
%     fprintf('\n\nIteration %d\n\n\n',r);
%     %points uniformly distributed on the set of Nmod p-spheres (p=n_bvs)
%     STCcf_0 = randn(n_bvs,nmods);
%     for i = 1:nmods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
%     
%     defmod.image_type = '2d';
%     white_props.basis = 'white';
%     white_props.pix_conv_mat = pix_conv_mat;
%     white_props.kern_conv_mat = kern_conv_mat;
%     glm_stcb = create2dGLM0_stcb(basis_vecs,STCcf_0,defmod,mod_signs,dim_signs,nltypes,init_nls,'test',white_props); %initialize
%     glm_stcb.lambdaW = 500;
%     
%     %determine LL and LP at current filter point
%     glm_stcb = fitWeights_stcb_nonlpsc(glm_stcb,kern_output*STCcf_0,spikebins,1,1e-3);
%     [ll0, ll0p] = getLLGLM2d_STCBF_nonlpsc(glm_stcb,kern_output,spikebins,'none');
%     
%     %store starting point
%     used_cfs = cat(3,used_cfs,STCcf_0);
%     LL_vals = [LL_vals; ll0]; LP_vals = [LP_vals; ll0p];
%     is_optimized = [is_optimized; 0];
%     
%     %compute distance in filter space to nearest better value
%     better_vals = find(LP_vals < ll0p);
%     n_better_vals = length(better_vals);
%     if n_better_vals > 0
%         fprintf('%d better values found\n',n_better_vals);
%         
%         %compute pairwise distances to better points (along the filter
%         %n-spheres) as sum of cosine distances
%         cur_dists = zeros(n_better_vals,1);
%         for nn = 1:nmods
%             cur_dists = cur_dists + (1-squeeze(sum(repmat(STCcf_0(:,nn),[1 1 n_better_vals]).*used_cfs(:,nn,better_vals))));
%         end
%         
%         %determine acceptance probaiblity based on smallest distance
%         fprintf('Smallest dist %.4f\n',min(cur_dists));
%         smallest_dists = [smallest_dists; min(cur_dists)];
%         %         acc_prob = 2./(1+exp(-smallest_dists(end)/dscale))-1;
%         acc_prob = normcdf(smallest_dists(end),rdmean,rdscale);
%         
%     else
%         acc_prob = 1;
%     end
%     
%     %use this starting point as a seed with probability acc_prob
%     fprintf('Will start with probability %.5f\n',acc_prob);
%     if rand < acc_prob
%         disp('Starting iteration');
%         %         %store starting point
%         %         used_cfs = cat(3,used_cfs,STCcf_0);
%         %         LL_vals = [LL_vals; ll0]; LP_vals = [LP_vals; ll0p];
%         % is_optimized = [is_optimized; 0];
%         
%         %local optimization
% %         rotbv_mod = [rotbv_mod; fitNLHI_stcb_nonlpsc(glm_stcb,stim_emb,spikebins,'none',6,2)];
%         rotbv_mod = [rotbv_mod; fitNLHI_stcb2d_nonlpsc(glm_stcb,scores,spikebins,'none',6,2)];
%         
% %         plot2d_mod(rotbv_mod(end),4)   
% %         pause
% %         close all
%         
%             %store final point
%         STCcf_fin = get_STCcf_mat(rotbv_mod(end));
%         %normalize filter weights to unit length for distance comparisons
%         for i = 1:nmods
%             STCcf_fin(:,i) = STCcf_fin(:,i)/norm(STCcf_fin(:,i));
%             dist_trav(r,i) = 1-dot(STCcf_fin(:,i),STCcf_0(:,i));
%         end
%         used_cfs = cat(3,used_cfs,STCcf_fin);
%         LL_vals = [LL_vals; rotbv_mod(end).LL]; LP_vals = [LP_vals; rotbv_mod(end).LP];
%         is_optimized = [is_optimized; 1];
%         
%         
%         %check stopping criteria
%         %         w = length(unique(LP_vals(is_optimized==1)/eps)); %number of eps-separated unique vals among converged set
%         
%         cur_opt_LP_vals = LP_vals(is_optimized==1);
%         cur_cluster_assignments = get_eps_ball_clusters(cur_opt_LP_vals,eps);
%         w = length(unique(cur_cluster_assignments));
%         cur_est_fract = (r-w-1)*(r+w)/(r*(r-1));
%         fprintf('Estimating: %d local minima found\n',w);
%         if cur_est_fract > min_LM_fract & r > min_reps
%             disp('Stopping criteria satisfied.  Halting.');
%             break
%         end
%     end
% end
% 
% %%
% minxy = [-4 -4];
% maxxy = [4 4];
% [bandwidth2_s1,density2_s1,X2,Y2]=kde2d(filt_proj_stim(:,1:2),2^8,minxy,maxxy);
% % [bandwidth2_s2,density2_s2,X2,Y2]=kde2d(filt_proj_stim(:,3:4),2^8,minxy,maxxy);
% 
% [bandwidth_s1,density_s1,X,Y]=kde2d(filt_proj_stim(spikebins,1:2),2^8,minxy,maxxy,bandwidth2_s1);
% % [bandwidth_s2,density_s2,X,Y]=kde2d(filt_proj_stim(spikebins,3:4),2^8,minxy,maxxy,bandwidth2_s2);
% 
% sta_pix = sta*pix_conv_mat;
% sta_pix = sta_pix/norm(sta_pix);
% stc_pix = stcs'*pix_conv_mat;
% for c = 1:size(stc_pix,1)
%     stc_pix(c,:) = stc_pix(c,:)/norm(stc_pix(c,:));
% end
% sta_filtproj = filt_long'*sta_pix';
% stc_filtproj = filt_long'*stc_pix';
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
% 
% %%
% figure; cmap = colormap(jet(size(stc_filtproj,2)));
% contourf(X,Y,cond_dens_s1,30)
% hold on
% circle([0 0],3,100,'k');
% line(3*[0 sta_filtproj(1)],3*[0 sta_filtproj(2)],'color','r','linewidth',2)
% for c = 1:size(stc_filtproj,2)
% line(3*[0 stc_filtproj(1,c)],3*[0 stc_filtproj(2,c)],'color',cmap(c,:),'linewidth',2)
% end
% 
% 
% opt_set = find(is_optimized==1);
% figure
% for r = 1:length(opt_set)
%     r
%     [LP_vals(opt_set(r))]
%     mod_filts = get_pix_mat(rotbv_mod(r));
%     for i = 1:nmods
%        mod_filts(:,1) = mod_filts(:,1)/norm(mod_filts(:,1)); 
%     end
%     mod_filtproj = filt_long'*mod_filts;
%     
% %     subplot(2,1,1)
%     contourf(X,Y,cond_dens_s1,30)
%     hold on
%     circle([0 0],3,100,'k');
%     line(3*[0 mod_filtproj(1,1)],3*[0 mod_filtproj(2,1)],'color','r','linewidth',2)
%     line(3*[0 mod_filtproj(1,2)],3*[0 mod_filtproj(2,2)],'color','w','linewidth',2)
% %     line(3*[0 mod_filtproj(1,3)],3*[0 mod_filtproj(2,3)],'color','g','linewidth',2)
% %     line(3*[0 mod_filtproj(1,4)],3*[0 mod_filtproj(2,4)],'color','c','linewidth',2)
%     vline(0,'w'); hline(0,'w')
%     
% %     subplot(2,1,2)
% %     contourf(X,Y,cond_dens_s2,30)
% %     hold on
% %     circle([0 0],3,100,'k');
% %     line(3*[0 mod_filtproj(3,1)],3*[0 mod_filtproj(4,1)],'color','r','linewidth',2)
% %     line(3*[0 mod_filtproj(3,2)],3*[0 mod_filtproj(4,2)],'color','w','linewidth',2)
% %     line(3*[0 mod_filtproj(3,3)],3*[0 mod_filtproj(4,3)],'color','g','linewidth',2)
% %     line(3*[0 mod_filtproj(3,4)],3*[0 mod_filtproj(4,4)],'color','c','linewidth',2)
% %     
% % %     line(3*[0 all_filtproj(best_mod,1,1)],3*[0 all_filtproj(best_mod,2,1)],'color','r','linewidth',2)
% % %     line(3*[0 all_filtproj(best_mod,1,2)],3*[0 all_filtproj(best_mod,2,2)],'color','w','linewidth',2)
% % %     line(3*[0 all_initproj(best_mod,1,1)],3*[0 all_initproj(best_mod,2,1)],'color','r','linewidth',1)
% % %     line(3*[0 all_initproj(best_mod,1,2)],3*[0 all_initproj(best_mod,2,2)],'color','w','linewidth',1)
% %     vline(0,'w'); hline(0,'w')
%     
%     pause
%     clf
% end
% 
% %% NOW FINAL ALIGNMENT OF ROTATED MODEL
% [~,best_mod] = min(LP_vals(is_optimized==1));
% STCcf_0 = nan(n_bvs,nmods);
% for i = 1:nmods
%     STCcf_0(:,i) = rotbv_mod(best_mod).mods(i).STCcf;
% end
% basis = 'white';
% flen = 6;
% %initialize model
% defmod.h(1:end-1) = []; %eliminate PSC
% defmod.lnl = 0;
% defmod.lh = 0;
% defmod.lnl2 = 100;
% defmod.lh2 = 0;
% defmod.nlcon = 1;
% defmod.nlmon = 1;
% defmod.locLambda = 500;
% defmod.locSigma = 5;
% defmod.maxLocPen = 500;
% defmod.lambda_dX = 50; %350
% defmod.lambda_L1x = 0; %40
% defmod.lambda_dT = 10;
% defmod.pids = 1:fsdim;
% defmod.SDIM = SDIM;
% defmod.fsdim = fsdim;
% 
% for i = 1:nmods; init_nls{i} = 'threshlin'; end;
% %define NL types: "uncon, lin, threshlin, quad"
% for i = 1:nmods; nltypes{i} = 'threshlin'; end;
% glm_stcb = createGLM2d_fullbf(basis_vecs,STCcf_0,pix_conv_mat,kern_conv_mat,defmod,nltypes,init_nls,basis,sprintf('test')); %initialize
% glm_stcb.lambdaW = 500;
% 
% %copy over internal NLs from previous fits
% for i = 1:nmods
%     glm_stcb.mods(i).nly = rotbv_mod(best_mod).mods(i).nly;
%     glm_stcb.mods(i).w = rotbv_mod(best_mod).mods(i).w;
% end
% glm_stcb.const = rotbv_mod(best_mod).const;
% 
% [glm_stcb,norm_vals] = normalizeRFs_full(glm_stcb,scores);
% glm_stcb.image_type = '2d';
% fin_glm = fitNLHI2d_fullbf(glm_stcb,scores,spikebins,'tots',6,3);
% mod_filts = get_pix_mat(fin_glm);
% fin_filtproj = filt_long'*mod_filts;
% plot2d_mod(fin_glm,4)
% 
% figure
% contourf(X,Y,cond_dens_s1,30)
% % surf(X,Y,cond_dens)
% hold on
% circle([0 0],3,100,'k');
% line(3*[0 fin_filtproj(1,1)],3*[0 fin_filtproj(2,1)],'color','r','linewidth',2)
% line(3*[0 fin_filtproj(1,2)],3*[0 fin_filtproj(2,2)],'color','w','linewidth',2)
% % line(3*[0 fin_filtproj(1,3)],3*[0 fin_filtproj(2,3)],'color','g','linewidth',2)
% % line(3*[0 fin_filtproj(1,4)],3*[0 fin_filtproj(2,4)],'color','c','linewidth',2)
% % 
% % figure
% % contourf(X,Y,cond_dens_s2,30)
% % % surf(X,Y,cond_dens)
% % hold on
% % circle([0 0],3,100,'k');
% % line(3*[0 fin_filtproj(3,1)],3*[0 fin_filtproj(4,1)],'color','r','linewidth',2)
% % line(3*[0 fin_filtproj(3,2)],3*[0 fin_filtproj(4,2)],'color','w','linewidth',2)
% % line(3*[0 fin_filtproj(3,3)],3*[0 fin_filtproj(4,3)],'color','g','linewidth',2)
% % line(3*[0 fin_filtproj(3,4)],3*[0 fin_filtproj(4,4)],'color','c','linewidth',2)
% 
% figure
% contourf(X,Y,cond_dens_s1,30)
% % surf(X,Y,cond_dens)
% hold on
% circle([0 0],3,100,'k');
% line(3*[0 sta_filtproj(1,1)],3*[0 sta_filtproj(2,1)],'color','k','linewidth',2)
% line(3*[0 stc_filtproj(1,1)],3*[0 stc_filtproj(2,1)],'color','w','linewidth',2)
% line(3*[0 stc_filtproj(1,2)],3*[0 stc_filtproj(2,2)],'color','r','linewidth',2)
% line(3*[0 stc_filtproj(1,3)],3*[0 stc_filtproj(2,3)],'color','g','linewidth',2)
% 
