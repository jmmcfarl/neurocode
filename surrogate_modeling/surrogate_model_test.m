%% 3d stim with 2-filter (Threshlin) LN model
clear all
close all

num_dims = 100;
num_samps = 1e6;
filt1 = zeros(num_dims,1);
filt2 = zeros(num_dims,1);
% filt3 = zeros(num_dims,1);
filt1(1) = 1;
filt2([1 2]) = 1/sqrt(2);
% filt3([1 2]) = -1/sqrt(2);
g_noise_var = 0;
offset = 0;
target_pspike = 0.01;
dt = 1e-3; %'1ms bins'
[X,Y] = meshgrid(-4:.02:4,-4:.02:4);
max_rate = 0.5;
 %%
stim = randn(num_samps,num_dims);

filt1_out = stim*filt1;
filt2_out = stim*filt2;
% filt3_out = stim*filt3;

filt1_out(filt1_out < 0) = 0;
filt2_out(filt2_out < 0) = 0;
% filt1_out = 0.25*filt1_out.^2;
% filt1_out(filt2_out < 0) = 0;
% filt2_out(filt2_out > 0) = 0;
% filt3_out(filt3_out < 0) = 0;
% filt3_out(filt3_out < 0) = 0;

g = filt1_out - filt2_out + g_noise_var*randn(num_samps,1);

f1 = X*filt1(1)+Y*filt1(2); f1(f1<0) = 0;
f2 = X*filt2(1)+Y*filt2(2); f2(f2<0) = 0;
% f3 = X*filt3(1)+Y*filt3(2); f3(f3<0) = 0;
gfun = f1-f2;
figure
contour(X,Y,gfun,30)
hold on
line([0 3*filt1(1)],[0 3*filt1(2)],'color','r')
line([0 3*filt2(1)],[0 3*filt2(2)],'color','b')
% line([0 3*filt3(1)],[0 3*filt3(2)],'color','r')
%%
g = g - mean(g);
% p_spike = log(1+exp(g));
% fprintf('Final mean rate %.4f\n',mean(p_spike))

% [y,x] = ksdensity(g);
% plot(x,y); hold on
% yl = ylim();
% tx = -10:.02:4;
% plot(tx,log(1+exp(tx)),'k')
% ylim(yl)
% 
% figure
% [y,x] = ksdensity(p_spike);
% plot(x,y);

%%
K0 = [-5 1];
Kfit = fmincon(@(K) abs(mean(logexp(g,K)) - target_pspike),K0,[],[],[],[],[-10 .01],[10 10]);
p_spike = logexp(g,Kfit);
p_spike(p_spike > max_rate) = max_rate;

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
% spikebins =[];
% for isp =1:length(rbins)
%     spikebins= [spikebins; repmat(rbins(isp),nsp(isp),1)];
% end
spk_vals = unique(spikes); spk_vals(spk_vals==0) = [];
spikebins = [];
for i = 1:length(spk_vals)
    cur_set = find(spikes == spk_vals(i));
    spikebins = [spikebins; repmat(cur_set(:),spk_vals(i),1)];
end
fprintf('Nspks: %d\n',length(spikebins));

%%
nneg = 1;
npos = 1;
spike_cond_stim = stim(spikebins,:);
sta      = mean(spike_cond_stim) - mean(stim);

stvcv = cov(spike_cond_stim);  utvcv = cov(stim);
[evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
stcs  = evecs(:,[1:nneg,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);

%%
minxy = [-3 -3];
maxxy = [3 3];
[bandwidth,density,X,Y]=kde2d(spike_cond_stim(:,1:2),2^8,minxy,maxxy);

[bandwidth2,density2,X2,Y2]=kde2d(stim(:,1:2),2^8,minxy,maxxy);

%%
eps = 1e-3;
cond_dens = density./density2;
cond_dens(density2 < eps) = eps;
outside = find(X.^2+Y.^2>9);
cond_dens(outside) = nan;
density(outside) = nan;
figure
contourf(X,Y,cond_dens,30)
% surf(X,Y,cond_dens)
hold on
circle([0 0],3,100,'k');
line(3*[0 stcs(1,1)],3*[0 stcs(1,2)],'color','r','linewidth',2)
line(3*[0 stcs(2,1)],3*[0 stcs(2,2)],'color','w','linewidth',2)
line(3*[0 sta(1)],3*[0 sta(2)],'color','k','linewidth',2)

figure
contourf(X,Y,density,30)
hold on
circle([0 0],3,100,'k');
line(3*[0 stcs(1,1)],3*[0 stcs(1,2)],'color','r','linewidth',2)
line(3*[0 stcs(2,1)],3*[0 stcs(2,2)],'color','w','linewidth',2)
line(3*[0 sta(1)],3*[0 sta(2)],'color','k','linewidth',2)

