clear all
close all
% fig_dir = '/home/james/Analysis/bruce/ET_final/';
fig_dir = '/Volumes/james/Analysis/bruce/ET_final/';

% bar_width = 0.0283;
bar_width = 0.0565;
usfac = 8;
% dx = bar_width/4;
dt = 0.01;

% xax = -1:dx:1;
pix_ax = -1.2:bar_width:1.2;
xax = linspace(-1,1,length(pix_ax)*usfac);
dx = median(diff(xax));

NT = round(3.7/dt);
taxis = 0:dt:(NT-1)*dt;
interp_taxis = 0:0.001:(NT-1)*dt;

rf_mean = 0;
rf_psf = 5;
rf_std = 1/rf_psf/2.2;
rf_phase = pi/2;
rf_amp = 1;

eye_std = 0.1;

gauss = exp(-(xax-rf_mean).^2/(2*rf_std^2));
RF = gauss.*sin(2*pi*xax*rf_psf + rf_phase)*rf_amp;
RFqp = gauss.*sin(2*pi*xax*rf_psf + rf_phase+pi/2)*rf_amp;
%%
% load /home/james/Analysis/bruce/G086/ET_final/sim_ETdata2.mat
load /Volumes/james/Analysis/bruce/G086/ET_final/sim_ETdata2.mat
interp_ep = spline(sim_t_axis,sim_ep,interp_taxis);

%% GENERATE STIM
dds = 0.12;
Npix = length(pix_ax);
% X = zeros(NT,Npix);
% nzero = rand(NT,Npix) < dds;
% polarity = rand(NT,Npix) > 0.5;
% X(polarity & nzero) = 1;
% X(nzero & ~polarity) = -1;
% 
% X = reshape(repmat(X,1,usfac),[NT Npix usfac]);
% X = permute(X,[1 3 2]);
% X = reshape(X,NT,length(xax));
% 
% save /Volumes/james/Analysis/bruce/G086/ET_final/RFsim_randXstim X

load('/Volumes/james/Analysis/bruce/G086/ET_final/RFsim_randXstim.mat','X')
%%
spk_thresh = -1.5;
spk_beta = 0.25;
mean_spk_rate = 50*dt;

filt_out = X*RF';
% filtqp_out = X*RFqp';
pred_rate = log(1+exp(spk_beta*filt_out + spk_thresh));
pred_rate = pred_rate*mean_spk_rate/mean(pred_rate);
interp_filt_out = spline(taxis,filt_out,interp_taxis);
interp_pred_rate = log(1+exp(spk_beta*interp_filt_out + spk_thresh));
interp_pred_rate = interp_pred_rate*mean_spk_rate/mean(interp_pred_rate);
% pred_rate = log(1+exp(spk_beta*(filt_out.^2 + filtqp_out.^2) + spk_thresh));


% figure; 
% [y,x] = hist(filt_out,1000);
% bar(x,y/sum(y));hold on
% 
% plot(x,log(1+exp(spk_beta*x + spk_thresh)),'r');
% 

poss_shifts = -round(0.4/dx):round(0.4/dx);
interp_shift_rate = nan(length(poss_shifts),length(interp_taxis));
for ii = 1:length(poss_shifts)
    shiftX = shift_matrix_Nd(X,-poss_shifts(ii),2);
    filt_out = shiftX*RF';
    interp_filt_out = spline(taxis,filt_out,interp_taxis);
    interp_shift_rate(ii,:) = log(1+exp(spk_beta*interp_filt_out + spk_thresh));
%     shift_rate(ii,:) = log(1+exp(spk_beta*filt_out + spk_thresh));
    %     filt_out = (shiftX*RF').^2 + (shiftX*RFqp').^2;
    %     shift_rate(ii,:) = log(1+exp(spk_beta*filt_out + spk_thresh));
end


%%    
close all
er = [-0.35 0.35];
% xr = [0.4 0.6];
xr = [1.45 1.65];
full_xr = [0.6 2.6];

tslice = 1.55/.001;

f1 = figure;
imagesc(interp_taxis,poss_shifts*dx,interp_shift_rate); hold on
caxis([0 0.6]); colorbar
plot(interp_taxis,interp_ep,'w');
ylim(er);    
set(gca,'ydir','normal');
xlim(xr);
yl = ylim();
line(interp_taxis([tslice tslice]),yl,'color','w','linewidth',2);
line
xlabel('Time (s)');
ylabel('Eye position (deg)');
colormap(gray)

% f2 = figure;
% plot(interp_taxis,interp_pred_rate/dt,'k')
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Firing rate (Hz)');

f2 = figure;
plot(interp_taxis,interp_ep,'k')
xlim(full_xr);
ylim(er);
line(xr([1 1]),er,'color','r');
line(xr([end end]),er,'color','r');
line(interp_taxis([tslice tslice]),yl,'color','b','linewidth',1);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');

f3 = figure;
plot(poss_shifts*dx,interp_shift_rate(:,tslice)/dt);
xlim(er);
yl = ylim();
line(interp_ep([tslice tslice]),yl);
xlabel('Eye position (deg)');
ylabel('Firing rate (Hz)');

% f4 = figure;
% plot(poss_shifts*dx,normpdf(poss_shifts*dx,0,eye_std),'k');
% xlim(er);
% xlabel('Eye position (deg)');
% ylabel('Probability');

f5 = figure;
plot(xax,RF,'k')
xlim(er);
yl = ylim();ya = max(abs(yl))+.01; ylim([-ya ya]);
% line(xlim(),[0 0],'color','k','linestyle','--');
hold on
plot(xax-interp_ep(tslice),RF,'r');
xlabel('Position (deg)');
ylabel('RF amplitude');

% ff1 = figure();
% plot(interp_taxis,randn(length(interp_taxis),1)*eye_std,'r');
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Eye position (deg)');
% 

%%
close all

% sim_eyepos = 40;
interp_ep_rnd = round(interp_ep/dx) - min(poss_shifts)+1;
sim_rate = nan(length(interp_taxis),1);
for ii = 1:length(interp_taxis)
sim_rate(ii) = interp_shift_rate(interp_ep_rnd(ii),ii);
end
sim_Robs = poissrnd(sim_rate);
% sim_Robs(tslice) = 1;
LL = bsxfun(@times,log(interp_shift_rate),sim_Robs') - interp_shift_rate;

% tslices = ceil(length(interp_taxis)*rand(10,1));

slice_LL = exp(LL(:,tslice));
% slice_LL = exp(LL(:,tslices));
slice_LL = bsxfun(@rdivide,slice_LL,sum(slice_LL));

f6 = figure;
plot(poss_shifts*dx,slice_LL,'linewidth',1);
xlim(er);
yl = ylim();
line(interp_ep([tslice tslice]),yl);

sim_rtimes = poissrnd(sim_rate/10);
sim_rtimes(sim_rtimes > 1) = 1;
sim_rtimes = interp_taxis(convert_to_spikebins(sim_rtimes));
sim_rtimes(diff(sim_rtimes) < 2e-3) = [];
f8 = figure; hold on
for ii = 1:length(sim_rtimes)
    line(sim_rtimes([ii ii]),[0 1],'color','k');
end
xlim(full_xr);


sim_rate = interp_shift_rate(interp_ep_rnd(tslice),:);
sim_Robs = poissrnd(sim_rate);
LL = bsxfun(@times,log(interp_shift_rate),sim_Robs) - interp_shift_rate;
tinds = ceil(length(interp_taxis)*rand(200,1));
tot_LL = exp(sum(LL(:,tinds),2));
tot_LL = tot_LL/sum(tot_LL);

f7 = figure;
plot(poss_shifts*dx,tot_LL,'k','linewidth',2)
xlim(er);
yl = ylim();
line(interp_ep([tslice tslice]),yl);

tslices = ceil(length(interp_taxis)*rand(10,1));
slice_LL = exp(LL(:,tslices));
slice_LL = bsxfun(@rdivide,slice_LL,sum(slice_LL));
f9 = figure;
plot(poss_shifts*dx,slice_LL,'linewidth',1)
xlim(er);
yl = ylim();
line(interp_ep([tslice tslice]),yl);


%%
fig_width = 4.68;
rel_height = 0.6;

figufy(f1);
fname = [fig_dir 'simExamp_ratedens.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

figufy(f2);
fname = [fig_dir 'simExamp_eyepos.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);


fig_width = 3.27;
rel_height = 0.8;

figufy(f3);
fname = [fig_dir 'simExamp_rateveyepos.pdf'];
exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f3);

figufy(f4);
fname = [fig_dir 'simExamp_eyedist.pdf'];
exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f4);

figufy(f5);
fname = [fig_dir 'simExamp_RF.pdf'];
exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f5);


figufy(f6);
fname = [fig_dir 'simExamp_singleLL.pdf'];
exportfig(f6,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f6);

figufy(f7);
fname = [fig_dir 'simExamp_totLL.pdf'];
exportfig(f7,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f7);

figufy(f8);
fname = [fig_dir 'simExamp_simrast.pdf'];
exportfig(f8,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f8);

figufy(f9);
fname = [fig_dir 'simExamp_multLL.pdf'];
exportfig(f9,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f9);

% figufy(ff1);
% fname = [fig_dir 'simExamp_epsamp.pdf'];
% exportfig(ff1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(ff1);


%%
NT = round(5000/dt);
taxis = 0:dt:(NT-1)*dt;
interp_taxis = 0:0.001:(NT-1)*dt;

X = zeros(NT,Npix);
nzero = rand(NT,Npix) < dds;
polarity = rand(NT,Npix) > 0.5;
X(polarity & nzero) = 1;
X(nzero & ~polarity) = -1;

X = reshape(repmat(X,1,2),[NT Npix 2]);
X = permute(X,[1 3 2]);
X = reshape(X,NT,length(xax));
%%

filt_out = X*RF';
% filtqp_out = X*RFqp';
pred_rate = log(1+exp(spk_beta*filt_out + spk_thresh));
Robs = poissrnd(pred_rate);

sta = mean(bsxfun(@times,X,Robs));

%%
eye_pos = round(randn(NT,1)*eye_std/dx);
shiftX = X;
for ii = 1:NT
    shiftX(ii,:) = shift_matrix_Nd(shiftX(ii,:),-eye_pos(ii),2);    
end
sta_shifted = mean(bsxfun(@times,shiftX,Robs));

%%
stim_params = NMMcreate_stim_params([1 length(xax)],dt);
reg_params = NMMcreate_reg_params('lambda_d2X',2000);
mod1 = NMMinitialize_model(stim_params,1,{'lin'},reg_params);
mod2 = mod1;

mod1 = NMMfit_filters(mod1,Robs,X);
mod2 = NMMfit_filters(mod2,Robs,shiftX);
linfilt = mod1.mods(1).filtK;
linfilt_shifted = mod2.mods(1).filtK;
%%

f6 = figure;
plot(xax,linfilt/norm(linfilt),'b');
hold on
plot(xax,linfilt_shifted/norm(linfilt),'r');
% plot(xax,RF/norm(RF),'k.-');
xlim(er);
yl = ylim();ya = max(yl); ylim([-ya ya]);
line(xlim(),[0 0],'color','k','linestyle','--');
xlabel('Position (deg)');
ylabel('RF amplitude');

%%
fig_width = 3.27;
rel_height = 0.8;

figufy(f6);
fname = [fig_dir 'simExamp_RFests.pdf'];
exportfig(f6,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f6);

