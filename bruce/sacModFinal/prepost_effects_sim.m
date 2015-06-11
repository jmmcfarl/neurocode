clear all
close all
fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';

base_dt = 0.01; %simulation base time res
usfac = 2; %temporal up-sampling factor 
flen = 15; %number of stimulus time lags
lag_tax = (1/usfac/2:1/usfac:flen)*0.01; %lag time axis

%create temporal kernel
xx = linspace(0,1,flen*usfac);
base_tkern = gampdf(xx,4,0.1);

%make the kernel go to zero somewhat faster than the gamma
weight_fun = ones(size(base_tkern));
weight_fun(lag_tax > 0.1) = (1-(20*(lag_tax(lag_tax > 0.1)-0.1)).^2);
weight_fun(weight_fun < 0) = 0;
base_tkern = base_tkern.*weight_fun;

base_tkern(1) = .001; %make non-zero to avoid having nans later
%%
sac_lags = -0.4:0.01/usfac:0.4; %saccade lags
sup_mag = 0.8; %suppression magnitude
at_sac = find(sac_lags == 0); 
cent_off = 0.03; %offset of peak suppression relative to sac timing
kern_width = 0.015; %SD of gaussian kernel width

%create pre-suppression kernel
pre_sup_kern = ones(size(sac_lags));
pre_sup_kern = pre_sup_kern - sup_mag*exp(-(sac_lags-cent_off).^2/(2*kern_width^2));

%create effective temporal filter for pre-filter
stim_mags = ones(size(sac_lags));
stim_in = stim_mags.*pre_sup_kern;
sparams = NMMcreate_stim_params(flen*usfac);
stimX = create_time_embedding(stim_in',sparams);

prestimR = bsxfun(@times,stimX,base_tkern);

%%
% conv_tkern = [zeros(1,flen*usfac+1) base_tkern];
% post_conv_kern = conv(pre_sup_kern,conv_tkern,'same');
% mval = max(post_conv_kern);
% post_conv_kern = (post_conv_kern-mval)/mval+1;
% post_conv_kern(sac_lags < 0) = 1;

cent_off = 0.06;
%create post-filter kernel
post_sup_kern = ones(size(sac_lags));
post_sup_kern = post_sup_kern - sup_mag*exp(-(sac_lags-cent_off).^2/(2*kern_width^2));

%create effective post-fliter kernel
stim_mags = ones(size(sac_lags));
stim_in = stim_mags;
sparams = NMMcreate_stim_params(flen*usfac);
stimX = create_time_embedding(stim_in',sparams);
stimX = bsxfun(@times,stimX,post_sup_kern');
postStimR = bsxfun(@times,stimX,base_tkern);

%% make test stim
% NT = length(sac_lags);
% testStim = randn(NT,1);
% testStimX = create_time_embedding(testStim,sparams);
% testStimFilt = testStimX*base_tkern';

%%

% f1 = figure();
% plot((1:NT)*0.01,testStim,'r')
% hold on
% plot((1:NT)*0.01,zscore(testStimFilt),'k')
% xlabel('Time (s)');
% ylabel('Stim');

% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'ppsim_stimex.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% plot pre- and post- gain kernels
f2 = figure();
subplot(2,1,1);
plot(sac_lags,pre_sup_kern); hold on
xlabel('Time since saccade onset (s)');
ylabel('Gain');
xlim([-.05 0.15]);

subplot(2,1,2);
plot(sac_lags,post_sup_kern,'r')
xlabel('Time since saccade onset (s)');
ylabel('Gain');
xlim([-.05 0.15]);

% fig_width = 4; rel_height = 1.6;
% figufy(f2);
% fname = [fig_dir 'ppsim_supkerns.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% plot temporal filter
f3 = figure();
plot(-lag_tax,base_tkern,'-');
xlabel('Time lag (s)');
ylabel('Filter amp');
xlim([-0.2 0]);

% fig_width = 4; rel_height = 0.8;
% figufy(f3);
% fname = [fig_dir 'ppsim_filter.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

%% plot effective kernels
xr = [-0.05 0.15];
f4 = figure();
subplot(2,1,1)
imagesc(sac_lags,-lag_tax,prestimR');
xlim(xr);
ylim([-0.12 0]);
xlabel('Time since saccade (s)');
ylabel('Stimulus latency (s)');

subplot(2,1,2)
imagesc(sac_lags,-lag_tax,postStimR');
xlim(xr);
ylim([-0.12 0]);
xlabel('Time since saccade (s)');
ylabel('Stimulus latency (s)');

% fig_width = 4; rel_height = 1.6;
% figufy(f4);
% fname = [fig_dir 'ppsim_effK.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);

%% plot change in stim sensitivity
xr = [-0.05 0.15];
f5 = figure();
subplot(2,1,1)
imagesc(sac_lags,-lag_tax,bsxfun(@rdivide,prestimR,base_tkern)');
xlim(xr);
caxis([0.2 1.8]);
ylim([-0.12 0]);
xlabel('Time since saccade (s)');
ylabel('Stimulus latency (s)');

subplot(2,1,2)
imagesc(sac_lags,-lag_tax,bsxfun(@rdivide,postStimR,base_tkern)');
xlim(xr);
caxis([0.2 1.8]);
ylim([-0.12 0]);
xlabel('Time since saccade (s)');
ylabel('Stimulus latency (s)');

% fig_width = 4; rel_height = 1.6;
% figufy(f5);
% fname = [fig_dir 'ppsim_effg.pdf'];
% exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f5);
% 

%%
plot_lags = [0.045 .075];
cmap = jet(length(plot_lags));
f6 = figure();
for ii = 1:length(plot_lags)
    subplot(2,1,1); hold on
    [~,cur_ind] = min(abs(sac_lags - plot_lags(ii)));
    plot(-lag_tax,prestimR(cur_ind,:),'color',cmap(ii,:));

    subplot(2,1,2); hold on
    plot(-lag_tax,postStimR(cur_ind,:),'color',cmap(ii,:));

end
subplot(2,1,1);
    plot(-lag_tax,base_tkern,'k');
xlim([-0.2 0]);
subplot(2,1,2);
    plot(-lag_tax,base_tkern,'k');
xlim([-0.2 0]);

% fig_width = 4; rel_height = 1.6;
% figufy(f6);
% fname = [fig_dir 'ppsim_examp_effkerns.pdf'];
% exportfig(f6,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f6);


%%
% f6 = figure();
% subplot(3,1,1);
% plot(sac_lags,testStim);
% xlim(xr);
% ylim([-3 3]);
% subplot(3,1,[2 3]);
% imagesc(sac_lags,-lag_tax,testStimX');
% xlim(xr); colormap(gray)
% caxis([-2.5 2.5]);

% f7 = figure();
% plot(-lag_tax,prestimR(46,:),'r-');
% xlim([-0.15 0]);

% 
% fig_width = 4; rel_height = 2;
% figufy(f6);
% fname = [fig_dir 'ppsim_testStim.pdf'];
% exportfig(f6,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f6);

% fig_width = 4; rel_height = 0.8;
% figufy(f7);
% fname = [fig_dir 'prestim_exampekern.pdf'];
% exportfig(f7,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f6);