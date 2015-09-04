clear all
clc

dt = 100;
dx = 1;
dy = 1;
stim_params(1).dims = [dt dx dy];
%%
NT = 1e4;
xax = linspace(-1,1,dt);
xax = xax - mean(xax);

sigma = 0.25;
lambda = 0.5;
phi = 0;
xoff = 0.25;
true_filt = exp(-(xax + xoff).^2/(2*sigma^2)).*cos(2*pi.*(xax + xoff)/lambda+phi);
xoff = -0.25;
true_filt2 = exp(-(xax + xoff).^2/(2*sigma^2)).*sin(2*pi.*(xax + xoff)/lambda+phi+pi/2);
true_filt3 = exp(-(xax + xoff).^2/(2*sigma^2)).*cos(2*pi.*(xax + xoff)/lambda+phi+pi/2);

Xstims{1} = randn(NT,dt);
filt_out =  Xstims{1}*true_filt';
filt_out2 =  Xstims{1}*true_filt2';
filt_out3 =  Xstims{1}*true_filt3';
G = filt_out + filt_out2.^2 + filt_out3.^2;
G = zscore(G);
rate = log(1+exp(G));
rObs = poissrnd(rate);

%%
close all
mod_signs = [1 1 1];
NLtypes = {'lin','quad','quad'};
% mod_signs = [1];
% NLtypes = {'lin'};

tic
clear nim
nim = NIM(stim_params,NLtypes,mod_signs);
nim = nim.set_reg_params('d2t',10);
nim = nim.fit_filters(rObs,Xstims,'silent',0);
toc

nim = nim.init_nonpar_NLs(Xstims, 'lambda_nld2',100);
nim.subunits(1).filtK = randn(size(nim.subunits(1).filtK))*.01;

nim = nim.fit_upstreamNLs(rObs, Xstims);

% nim = nim.fit_filters(rObs,Xstims,'silent',0);

%%
tic
nim_stim_params = NMMcreate_stim_params([dt dx dy]);
reg_params = NMMcreate_reg_params('lambda_d2T',10);
old_nim = NMMinitialize_model(nim_stim_params,mod_signs,NLtypes,reg_params);
old_nim = NMMfit_filters(old_nim,rObs,Xstims,[],[],0);
toc
%% 
figure; 
subplot(3,1,1);
est_filt = nim.subunits(1).filtK;
est_filt2 = old_nim.mods(1).filtK;
plot(est_filt/max(est_filt)); hold on; plot(true_filt,'r');
plot(est_filt2/max(est_filt2),'k');
subplot(3,1,2);
est_filt = nim.subunits(2).filtK;
est_filt2 = old_nim.mods(2).filtK;
plot(est_filt/max(est_filt)); hold on; plot(true_filt2,'r')
plot(est_filt2/max(est_filt2),'k');
subplot(3,1,3);
est_filt = nim.subunits(3).filtK;
est_filt2 = old_nim.mods(3).filtK;
plot(est_filt/max(est_filt)); hold on; plot(true_filt3,'r')
plot(est_filt2/max(est_filt2),'k');


%%
% f1 = figure(); hold on
% nim = NIM(stim_params,NLtypes,mod_signs);
% nim = nim.set_reg_params('d2t',0);
% nim = nim.fit_filters(rObs,Xstims);
% plot(nim.subunits(1).filtK);
% nim = NIM(stim_params,NLtypes,mod_signs);
% nim = nim.set_reg_params('d2t',10);
% nim = nim.fit_filters(rObs,Xstims);
% plot(nim.subunits(1).filtK,'r');
% nim = NIM(stim_params,NLtypes,mod_signs);
% nim = nim.set_reg_params('d2t',1000);
% nim = nim.fit_filters(rObs,Xstims);
% plot(nim.subunits(1).filtK,'k');
% nim = NIM(stim_params,NLtypes,mod_signs);
% nim = nim.set_reg_params('d2t',10000);
% nim = nim.fit_filters(rObs,Xstims);
% plot(nim.subunits(1).filtK,'c');
%%
% [penLL, penLLgrad] = nim.internal_LL_grad_filters(params,Xstims,rObs,targets,nontarg_g)