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
x = repmat(1:SDIM,flen,1);
x0 = repmat(11,flen,SDIM);
sigma = repmat(1.5,flen,SDIM);
lambda = repmat(6,flen,SDIM);
b = repmat(linspace(2,0,flen)',1,SDIM);a = repmat(0,flen,SDIM);
psi1 = repmat(linspace(0,pi,flen)',1,SDIM);
psi2 = psi1 + pi/2;
filt1 = b.*exp(-((x-x0).^2./2./sigma.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi1))+a;
filt2 = b.*exp(-((x-x0).^2./2./sigma.^2)) .* (cos(2*pi.*(x-x0)./lambda+psi2))+a;

% filt1 = bsxfun(@minus,filt1,mean(filt1,2));
% filt2 = bsxfun(@minus,filt2,mean(filt2,2));

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

% filt3 = bsxfun(@minus,filt3,mean(filt3,2));
% filt4 = bsxfun(@minus,filt4,mean(filt4,2));

filt3 = filt3/norm(filt3(:));
filt4 = filt4/norm(filt4(:));
%exact orthogonalization of filt2
temp3 = filt3(:); temp4 = filt4(:);
temp4n = temp4 - (temp3'*temp4)*temp3;
filt4 = reshape(temp4n,flen,SDIM);
filt3 = filt3/norm(filt3(:));
filt4 = filt4/norm(filt4(:));


figure
subplot(2,2,1)
imagesc(filt1); colormap(gray); caxis([-0.3 0.3])
subplot(2,2,2)
imagesc(filt2); colormap(gray); caxis([-0.3 0.3])
subplot(2,2,3)
imagesc(filt3); colormap(gray); caxis([-0.3 0.3])
subplot(2,2,4)
imagesc(filt4); colormap(gray); caxis([-0.3 0.3])

f1 = X; f1(X < 0) = 0;
f2 = Y; f2(Y < 0) = 0;
gfun = f1+f2;
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
target_rate = 50; %in Hz
target_pspike = target_rate*dt;
K0 = [3 1.5];
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

%%
minxy = [-4 -4];
maxxy = [4 4];
[bandwidth_s1,density_s1,X,Y]=kde2d(filt_proj_stim(spikebins,1:2),2^8,minxy,maxxy);
[bandwidth_s2,density_s2,X,Y]=kde2d(filt_proj_stim(spikebins,3:4),2^8,minxy,maxxy);

[bandwidth2_s1,density2_s1,X2,Y2]=kde2d(filt_proj_stim(:,1:2),2^8,minxy,maxxy);
[bandwidth2_s2,density2_s2,X2,Y2]=kde2d(filt_proj_stim(:,3:4),2^8,minxy,maxxy);

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

%% CREATE WAVELET FILTER BANK
clear gabor_props gabor_bank
% psi_vals = linspace(0,2*pi,4);
psi_vals = [0];

lambda_vals = [4 6 8];
x0_vals = 2:3:16;
sigma_vals = [1 1.5 2 2.5];
n_filts = length(psi_vals)*length(lambda_vals)*length(x0_vals)*length(sigma_vals);
gabor_bank = nan(n_filts,SDIM);
cnt = 1;
for i = 1:length(psi_vals)
    for j = 1:length(lambda_vals)
        for k = 1:length(x0_vals)
            for ii = 1:length(sigma_vals)
                cur_gabor = gabor_fun_1d(1:16,x0_vals(k),sigma_vals(ii),lambda_vals(j),psi_vals(i),0,1);
                gabor_bank(cnt,:) = cur_gabor;
                gabor_props(cnt).psi = psi_vals(i);
                gabor_props(cnt).lambda = lambda_vals(j);
                gabor_props(cnt).x0 = x0_vals(k);
                gabor_props(cnt).sigma = sigma_vals(ii);
                gabor_props(cnt).a = 0;
                gabor_props(cnt).b = 1;
                cnt = cnt + 1;
            end
        end
    end
end
gabor_psis = [gabor_props(:).psi];
gabor_lambda = [gabor_props(:).lambda];
gabor_x0 = [gabor_props(:).x0];
gabor_sigma = [gabor_props(:).sigma];

%%
NT = length(spikes);
spike_rate = spikes/dt;
cur_F = ones(NT,1)*(log(exp(mean(spike_rate))-1));
y_tilda = exp(cur_F)./(1+exp(cur_F)) .* (spike_rate./log(1+exp(cur_F)) - 1);

%%
close all
for ff = 1:8

for n = 1:n_filts
    fprintf('Filt set %d of %d\n',n,n_filts);
    gabor_mat = nan(flen,flen*SDIM);
    gabor_mat_qp = nan(flen,flen*SDIM);
    zero_mat = zeros(flen,SDIM);
    for i = 1:flen
        cur_gabor = gabor_fun_1d(1:16,gabor_props(n).x0,gabor_props(n).sigma,gabor_props(n).lambda,0,0,1);
        cur_gabor_qp = gabor_fun_1d(1:16,gabor_props(n).x0,gabor_props(n).sigma,gabor_props(n).lambda,pi/2,0,1);
        temp_gabor_mat = zero_mat;
        temp_gabor_mat(i,:) = cur_gabor;
        gabor_mat(i,:) = temp_gabor_mat(:);
        temp_gabor_mat = zero_mat;
        temp_gabor_mat(i,:) = cur_gabor_qp;
        gabor_mat_qp(i,:) = temp_gabor_mat(:);
    end
    
    gabor_mat_out = [stim_emb*gabor_mat' stim_emb*gabor_mat_qp'];
    
    k0 = randn(28,1)*0.1;
    opts = optimset('GradObj','on','Algorithm','active-set','Display','off','MaxIter',500,'MaxFunEvals',10000,'TolFun',1e-5);
    [bopt(:,n),r2val(n)] = fminunc(@(K) grad_boosting_filtbank(K,y_tilda,gabor_mat_out),k0,opts);
    
end

[a,b] = min(r2val);
for i = 1:flen
    cur_gabor = gabor_fun_1d(1:16,gabor_props(b).x0,gabor_props(b).sigma,gabor_props(b).lambda,0,0,1);
    cur_gabor_qp = gabor_fun_1d(1:16,gabor_props(b).x0,gabor_props(b).sigma,gabor_props(b).lambda,pi/2,0,1);
    temp_gabor_mat = zero_mat;
    temp_gabor_mat(i,:) = cur_gabor;
    gabor_mat(i,:) = temp_gabor_mat(:);
    temp_gabor_mat = zero_mat;
    temp_gabor_mat(i,:) = cur_gabor_qp;
    gabor_mat_qp(i,:) = temp_gabor_mat(:);
end
gabor_fit(ff,:) = bopt(1:14,b)'*gabor_mat + bopt(15:end,b)'*gabor_mat_qp;
cur_gabor_out = stim_emb*gabor_fit(ff,:)';
cur_gabor_out(cur_gabor_out < 0) = 0;

cur_F = cur_F + cur_gabor_out;
y_tilda = exp(cur_F)./(1+exp(cur_F)) .* (spike_rate./log(1+exp(cur_F)) - 1);

% imagesc(reshape(gabor_fit(ff,:),14,SDIM))
% 
% pause
% close all
end

figure
subplot(2,2,1)
imagesc(reshape(gabor_fit(1,:),flen,SDIM)); colormap(gray);
subplot(2,2,2)
imagesc(reshape(gabor_fit(2,:),flen,SDIM)); colormap(gray);
subplot(2,2,3)
imagesc(reshape(gabor_fit(3,:),flen,SDIM)); colormap(gray);
subplot(2,2,4)
imagesc(reshape(gabor_fit(4,:),flen,SDIM)); colormap(gray);

figure
subplot(2,2,1)
imagesc(reshape(gabor_fit(5,:),flen,SDIM)); colormap(gray);
subplot(2,2,2)
imagesc(reshape(gabor_fit(6,:),flen,SDIM)); colormap(gray);
subplot(2,2,3)
imagesc(reshape(gabor_fit(7,:),flen,SDIM)); colormap(gray);
subplot(2,2,4)
imagesc(reshape(gabor_fit(8,:),flen,SDIM)); colormap(gray);

%%
cd ~/Data/blanche/rec_75/matlabdata/
load stdparsRec75.mat
flen = 14;
%initialize model
defmod.h(1:end-1) = []; %eliminate PSC
defmod.lnl = 0;
defmod.lh = 0;
defmod.lnl2 = 0;
defmod.lh2 = 0;
defmod.nlcon = 0;
defmod.nlmon = 0;
defmod.locLambda = 0;
defmod.lambda_dX = 100; %350
defmod.lambda_L1x = 10; %40
defmod.lambda_dT = 50;
defmod.pids = 1:SDIM;
defmod.SDIM = SDIM;
defmod.fsdim = SDIM;

n_mods = 4;
n_bfdims = size(gabor_fit,1);
mod_signs = ones(n_mods,1);
dim_signs = ones(n_bfdims,1);

clear stc_posneg_mod
for r = 1:10
    STCcf_0 = randn(n_bfdims,n_mods);    
    %normalize
    for i = 1:n_mods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
glm_stcb = createGLM0_stcb(gabor_fit',STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
stc_posneg_mod(r) = fitNLHI_stcb_nonlpsc(glm_stcb,stim_emb,spikebins,'tots');
end

[a,b] = min(arrayfun(@(x) x.LL,stc_posneg_mod));
mod_filts = get_k_mat(stc_posneg_mod(b));
mod_filtproj = filt_long*mod_filts;

figure
subplot(2,1,1)
contourf(X,Y,cond_dens_s1,30)
hold on
circle([0 0],3,100,'k');
line(3*[0 mod_filtproj(1,1)],3*[0 mod_filtproj(2,1)],'color','r','linewidth',2)
line(3*[0 mod_filtproj(1,2)],3*[0 mod_filtproj(2,2)],'color','w','linewidth',2)
line(3*[0 mod_filtproj(1,3)],3*[0 mod_filtproj(2,3)],'color','g','linewidth',2)
line(3*[0 mod_filtproj(1,4)],3*[0 mod_filtproj(2,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')
subplot(2,1,2)
contourf(X,Y,cond_dens_s2,30)
hold on
circle([0 0],3,100,'k');
line(3*[0 mod_filtproj(3,1)],3*[0 mod_filtproj(4,1)],'color','r','linewidth',2)
line(3*[0 mod_filtproj(3,2)],3*[0 mod_filtproj(4,2)],'color','w','linewidth',2)
line(3*[0 mod_filtproj(3,3)],3*[0 mod_filtproj(4,3)],'color','g','linewidth',2)
line(3*[0 mod_filtproj(3,4)],3*[0 mod_filtproj(4,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')

%%
STCbvs = [sta' stcs];
STCbvs = STCbvs(:,[1 2 3 4]);
Nstcbvs = size(STCbvs,2);
nmods = 4;
mod_signs = ones(nmods,1);
dim_signs = ones(Nstcbvs,1);

for r = 1:10
    STCcf_0 = randn(Nstcbvs,n_mods);
    %normalize
    for i = 1:n_mods; STCcf_0(:,i) = STCcf_0(:,i)/norm(STCcf_0(:,i)); end;
    glm_stcb = createGLM0_stcb(STCbvs,STCcf_0,defmod,mod_signs,dim_signs,'test'); %initialize
    stc_posneg_mod_stc(r) = fitNLHI_stcb_nonlpsc(glm_stcb,stim_emb,spikebins,'tots');
end

[a,b] = min(arrayfun(@(x) x.LL,stc_posneg_mod_stc));
mod_filts = get_k_mat(stc_posneg_mod_stc(b));
mod_filtproj = filt_long*mod_filts;

figure
subplot(2,1,1)
contourf(X,Y,cond_dens_s1,30)
hold on
circle([0 0],3,100,'k');
line(3*[0 mod_filtproj(1,1)],3*[0 mod_filtproj(2,1)],'color','r','linewidth',2)
line(3*[0 mod_filtproj(1,2)],3*[0 mod_filtproj(2,2)],'color','w','linewidth',2)
line(3*[0 mod_filtproj(1,3)],3*[0 mod_filtproj(2,3)],'color','g','linewidth',2)
line(3*[0 mod_filtproj(1,4)],3*[0 mod_filtproj(2,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')
subplot(2,1,2)
contourf(X,Y,cond_dens_s2,30)
hold on
circle([0 0],3,100,'k');
line(3*[0 mod_filtproj(3,1)],3*[0 mod_filtproj(4,1)],'color','r','linewidth',2)
line(3*[0 mod_filtproj(3,2)],3*[0 mod_filtproj(4,2)],'color','w','linewidth',2)
line(3*[0 mod_filtproj(3,3)],3*[0 mod_filtproj(4,3)],'color','g','linewidth',2)
line(3*[0 mod_filtproj(3,4)],3*[0 mod_filtproj(4,4)],'color','c','linewidth',2)
vline(0,'w'); hline(0,'w')


%%
gabor_fit = bopt(1:14)'*gabor_mat + bopt(15:end)'*gabor_mat_qp;
imagesc(reshape(gabor_fit,14,SDIM))
%%
poss_lags = 0:(flen-1);
n_filts = n_filts/2;
gabor_bank_p1 = gabor_bank(n_filts+1:end,:);
gabor_bank_p2 = gabor_bank(1:n_filts,:);
for ii = 1:n_filts
    ii
pred_mat = [];
    for ll = 1:length(poss_lags)
        lagged_stim = [zeros(poss_lags(cur_lags),SDIM); stim(1:end-poss_lags(cur_lags),:)];
        cur_out_p1 = lagged_stim*gabor_bank_p1(ii,:)';
        cur_out_p2 = lagged_stim*gabor_bank_p2(ii,:)';
        
        pred_mat = [pred_mat cur_out_p1 cur_out_p2];
        
    end
    
    [b,bint,r,rint,stats] = regress(y_tilda,pred_mat);
    r2(ii) = stats(1);
end

%%
poss_lags = 6:(flen-1);
clear gabor_corrs
for cur_lags = 1:length(poss_lags)
    cur_lags
    lagged_stim = [zeros(poss_lags(cur_lags),SDIM); stim(1:end-poss_lags(cur_lags),:)];
    gabor_out = lagged_stim*gabor_bank';
    gabor_out(gabor_out < 0) = 0;
    gabor_corrs(cur_lags,:) = corr(y_tilda,gabor_out);   
end

[best_corr,best_corr_loc] = max(gabor_corrs(:));
% [best_corr,best_corr_loc] = max(abs(gabor_corrs(:)));
[best_lag_loc,best_gabor_loc] = ind2sub(size(gabor_corrs),best_corr_loc);
best_lag = poss_lags(best_lag_loc);

%%
cur_lambda = gabor_props(best_gabor_loc).lambda;
cur_x0 = gabor_props(best_gabor_loc).x0;
cur_sigma = gabor_props(best_gabor_loc).sigma;

poss_psi_vals = linspace(0,2*pi,20);
n_filts = length(poss_psi_vals);
poss_gabor_bank = nan(n_filts,SDIM);
cnt = 1;
for i = 1:length(poss_psi_vals)
    cur_gabor = gabor_fun_1d(1:16,cur_x0,cur_sigma,cur_lambda,poss_psi_vals(i),0,1);
    poss_gabor_bank(cnt,:) = cur_gabor;
    cnt = cnt + 1;
end

%%
    cur_gabor_emb = nan(flen,SDIM);
poss_lags = 0:(flen-1);
while ~isempty(poss_lags)
    poss_lags
    cur_gabor_corrs = nan(length(poss_lags),n_filts);
    for cur_lags = 1:length(poss_lags)
        lagged_stim = [zeros(poss_lags(cur_lags),SDIM); stim(1:end-poss_lags(cur_lags),:)];
        gabor_out = lagged_stim*poss_gabor_bank';
        gabor_out(gabor_out < 0) = 0;
        cur_gabor_corrs(cur_lags,:) = corr(y_tilda,gabor_out);
    end
    [best_corr,best_corr_loc] = max(cur_gabor_corrs(:));
    [cur_best_lag_loc,cur_best_gabor_loc] = ind2sub(size(cur_gabor_corrs),best_corr_loc);
    cur_best_lag = poss_lags(cur_best_lag_loc);  
    
    cur_gabor = poss_gabor_bank(cur_best_gabor_loc,:);
    cur_gabor_emb(cur_best_lag+1,:) = cur_gabor;
    
    lagged_stim = [zeros(poss_lags(cur_best_lag_loc),SDIM); stim(1:end-poss_lags(cur_best_lag_loc),:)];
    gabor_out = lagged_stim*cur_gabor';
    gabor_out(gabor_out < 0) = 0;
    cur_beta(cur_best_lag+1) = mldivide(gabor_out,y_tilda);
    cur_F = cur_F + cur_beta(cur_best_lag+1)*gabor_out;
    y_tilda = exp(cur_F)./(1+exp(cur_F)) .* (spike_rate./log(1+exp(cur_F)) - 1);   
    poss_lags(cur_best_lag_loc) = [];  
end




