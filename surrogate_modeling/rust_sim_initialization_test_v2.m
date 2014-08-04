cd ~/James_scripts/surrogate_modeling/
addpath(genpath('~/James_scripts'));
clear all
% close all

dt = 0.01; %in s
max_rate = Inf; %in Hz
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);

NT = 200000; SDIM = 24; flen = 14;
stim = round(2*rand(NT,SDIM))-1;
xvNT = 100000; 
xvstim = round(2*rand(xvNT,SDIM))-1;

% GENERATE A FILTER
nfilts = 6;
LAMBDA = 5; %5 4
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

desired_spacing = LAMBDA/5.75; %5.5/ 5 / 4.5n
beg = SDIM/2-desired_spacing*floor(nfilts/2)+0.75;
xovals = 0:desired_spacing:(desired_spacing*(nfilts-1));
xovals = xovals - mean(xovals);
xovals = xovals + SDIM/2+0.5;
for i = 1:nfilts
    psi1 = repmat(((1:flen)'/6).^3+5,1,SDIM); 

    x0 = repmat(xovals(i),flen,SDIM);
%     cur_psi = mod(psi1 + rand*pi,2*pi);
    temp = b.*exp(-((x-x0).^2./2./sigma1.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+a;
    filt(i,:,:) = temp/norm(temp(:));
    filt_mat(i,:) = temp(:)/norm(temp(:));
end
filt_mat = filt_mat';
c_cent = SDIM/2+0.5;
c_std = 1.6; %1.6
max_cval = 0.5; %1
c_offset = 1;
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
beta = 3; %2.5
theta = 0;
coefs = ones(1,nfilts).*cvals;
g = zeros(1,NT);
xvg = zeros(1,xvNT);
for i = 1:nfilts
    lfilt_stim(i,:) = 1/beta*log(1+exp(beta*(filt_stim(i,:)-theta)));
    xvlfilt_stim(i,:) = 1/beta*log(1+exp(beta*(xvfilt_stim(i,:)-theta)));
    g = g + coefs(i)*lfilt_stim(i,:);
    xvg = xvg + coefs(i)*xvlfilt_stim(i,:);
end
sg = std(g);
g = g/sg;
xvg = xvg/sg;

target_rate = 50; %in Hz
target_pspike = target_rate*dt;
cur_theta = 1.75; %1.5  2.5
cur_beta = 4; %4  4
p_spike = 1/cur_beta*log(1+exp(cur_beta*(g-cur_theta)));
% p_spike = (g-cur_theta).^2; p_spike(p_spike<0) = 0;
% p_spike = exp((g-1.5)/1);
xvp_spike = 1/cur_beta*log(1+exp(cur_beta*(xvg-cur_theta)));
scale_f = target_rate/(mean(p_spike)/dt);
p_spike = p_spike*scale_f;
% p_spike(p_spike/dt > 500) = 500*dt;
xvp_spike = xvp_spike*scale_f;
% xvp_spike(xvp_spike/dt > 500) = 500*dt;


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

stc_dims = [sta' stcs(:,1:5)];


%%
sdim = SDIM; flen = 14;
stim_params.spatial_dims = 1;
stim_params.sdim = sdim;
stim_params.flen = flen;

defmod.lambda_d2XT = 0;
defmod.lambda_L1x =50;
nmods = 6;
n_iter = 100;
for i = 1:n_iter
%     init_kerns = randn(flen*sdim,nmods);
%     
%     init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
%     init_signs = ones(1,nmods);
%     
%     clear kern_types
%     for j = 1:nmods
%         kern_types{j} = 'threshlin';
%     end
%     gnm(i) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
%     init_proj(i,:,:) = init_kerns'*filt_mat;
%     gnm(i) = fitGNM_filters(gnm(i),stim_emb,spikebins,'none',[],1e-4,1e-6);
%     gnm_xvLL(i) = getLL_GNM(gnm(i),xvstim_emb,xvspikebins,'none')
%     fin_kern = get_k_mat(gnm(i));
%     fin_proj(i,:,:) = fin_kern'*filt_mat;
    gnm_LL(i) = getLL_GNM(gnm(i),stim_emb,spikebins,'none');
end

%%
n_iter = 100;
for cur_i = 1:n_iter;
    [~, ~, ~, ~, g] = getLL_GNM(gnm(cur_i),stim_emb,spikebins,'none');
    gnmu(cur_i) = fitGNM_spkNL(gnm(cur_i),g,spikebins,0);
    gnmu_xvLL(cur_i) = getLL_GNM(gnmu(cur_i),xvstim_emb,xvspikebins,'none');
    
    tn_iter = 2;
    %for fitting internal NLs
    gnmr(cur_i) = adjust_all_reg(gnmu(cur_i),'lnl2',500); %1000
    gnmr(cur_i) = setGNM_NLBFs(gnmr(cur_i),stim_emb);
    gnmr(cur_i) = adjust_all_reg(gnmr(cur_i),'nltype','uncon');
    gnmr(cur_i) = adjust_all_reg(gnmr(cur_i),'nlmon',1);
    for j = 1:tn_iter
        gnmr(cur_i) = fitGNM_internal_NLs(gnmr(cur_i),stim_emb,spikebins,1,2);
        gnmr(cur_i) = fitGNM_filters(gnmr(cur_i),stim_emb,spikebins,'none',[],1e-4,1e-6);
        [~, ~, ~, ~, g] = getLL_GNM(gnmr(cur_i),stim_emb,spikebins,'none');
        gnmr(cur_i) = fitGNM_spkNL(gnmr(cur_i),g,spikebins,0);
    end
    
    gnmr_xvLL(cur_i) = getLL_GNM(gnmr(cur_i),xvstim_emb,xvspikebins,'none')
    gnmr_LL(cur_i) = getLL_GNM(gnmr(cur_i),stim_emb,spikebins,'none')
end
%%
save lexp_rust_init_test gnm gnm_* gnmr* init* fin* filt*

%%
load lexp_rust_init_test

% for i = 1:n_iter
%     hold on
% %     if ismember(i,set1)
%         line([squeeze(init_proj(i,1,1)) squeeze(fin_proj(i,1,1))],[squeeze(init_proj(i,1,2)) squeeze(fin_proj(i,1,2))])
%         line([squeeze(init_proj(i,2,1)) squeeze(fin_proj(i,2,1))],[squeeze(init_proj(i,2,2)) squeeze(fin_proj(i,2,2))],'color','r')
% %     else
% %         line([squeeze(init_proj(i,2,1)) squeeze(fin_proj(i,2,1))],[squeeze(init_proj(i,2,2)) squeeze(fin_proj(i,2,2))])
% %         line([squeeze(init_proj(i,1,1)) squeeze(fin_proj(i,1,1))],[squeeze(init_proj(i,1,2)) squeeze(fin_proj(i,1,2))],'color','r')
% %     end
% end
  for i = 1:n_iter
     fin_kern = get_k_mat(gnm(i));
     fin_kern = bsxfun(@rdivide,fin_kern,sqrt(sum(fin_kern.^2)));
     fin_projn(i,:,:) = fin_kern'*filt_mat;
     fin_kern = get_k_mat(gnmr(i));
     fin_kern = bsxfun(@rdivide,fin_kern,sqrt(sum(fin_kern.^2)));
     fin_projn2(i,:,:) = fin_kern'*filt_mat;
  end
true_proj = filt_mat'*filt_mat;

fin_projn(67,:,:) = []; 
gnm_LL(67) = [];

u1 = 1;
u2 = 6;
% close all
bad_set = find(gnm_LL > 1.109);
figure
for i = 1:6
plot(squeeze(fin_projn(:,i,u1)),squeeze(fin_projn(:,i,u2)),'o')
hold on
plot(squeeze(fin_projn(bad_set,i,u1)),squeeze(fin_projn(bad_set,i,u2)),'go')
% plot(squeeze(fin_projn2(:,i,u1)),squeeze(fin_projn2(:,i,u2)),'ko')
plot(true_proj(i,u1),true_proj(i,u2),'r*','markersize',8,'linewidth',2)
end

%%
all_projn = [];
for i = 1:6
    all_projn = [all_projn; squeeze(fin_projn(:,i,u1)) squeeze(fin_projn(:,i,u2))];
end
[bandwidth,density,X,Y]=kde2d(all_projn,2^8,[-1 -1],[1.2 1.2],0.025);
imagesc(density)