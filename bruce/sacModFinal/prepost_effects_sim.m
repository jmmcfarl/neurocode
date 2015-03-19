clear all
close all
fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';

flen = 15;
xx = linspace(0,1,flen);
base_tkern = gampdf(xx,4,0.1);
base_tkern(1)=0.02;
lag_tax = (1:flen)*0.01 - 0.005;
%%
sac_lags = -0.4:0.01:0.4;
sup_mag = 0.8;
at_sac = find(sac_lags == 0);
pre_sup_kern = ones(size(sac_lags));
pre_sup_kern = pre_sup_kern - sup_mag*exp(-sac_lags.^2/(2*0.01^2));

stim_mags = ones(size(sac_lags));
stim_in = stim_mags.*pre_sup_kern;
sparams = NMMcreate_stim_params(15);
stimX = create_time_embedding(stim_in',sparams);

prestimR = bsxfun(@times,stimX,base_tkern);

%%
post_sup_kern = ones(size(sac_lags));
post_sup_kern = post_sup_kern - sup_mag*exp(-(sac_lags-0.05).^2/(2*0.01^2));

stim_mags = ones(size(sac_lags));
stim_in = stim_mags;
sparams = NMMcreate_stim_params(15);
stimX = create_time_embedding(stim_in',sparams);
stimX = bsxfun(@times,stimX,post_sup_kern');
postStimR = bsxfun(@times,stimX,base_tkern);

%%
NT = 50;
testStim = randn(NT,1);
testStimX = create_time_embedding(testStim,sparams);
testStimFilt = testStimX*base_tkern';

f1 = figure();
plot((1:NT)*0.01,testStim,'r')
hold on
plot((1:NT)*0.01,zscore(testStimFilt),'k')
xlabel('Time (s)');
ylabel('Stim');

% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'ppsim_stimex.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%%
f2 = figure();
plot(sac_lags,pre_sup_kern); hold on
plot(sac_lags,post_sup_kern,'r')
xlabel('Time since saccade onset (s)');
ylabel('Gain');
xlim([-0.25 0.25]);

% fig_width = 4; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'ppsim_supkerns.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
%%
f3 = figure();
plot(-lag_tax,base_tkern,'o-');
xlabel('Time lag (s)');
ylabel('Filter amp');
xlim([-0.15 0]);

% fig_width = 4; rel_height = 0.8;
% figufy(f3);
% fname = [fig_dir 'ppsim_filter.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

%%
xr = [-0.05 0.2];
f4 = figure();
subplot(2,1,1)
imagesc(sac_lags,-lag_tax,prestimR');
xlim(xr);
xlabel('Time since saccade (s)');
ylabel('Stimulus latency (s)');

subplot(2,1,2)
imagesc(sac_lags,-lag_tax,postStimR');
xlim(xr);
xlabel('Time since saccade (s)');
ylabel('Stimulus latency (s)');

f5 = figure();
subplot(2,1,1)
imagesc(sac_lags,-lag_tax,bsxfun(@rdivide,prestimR,base_tkern)');
xlim(xr);
caxis([0.2 1.8]);
xlabel('Time since saccade (s)');
ylabel('Stimulus latency (s)');

subplot(2,1,2)
imagesc(sac_lags,-lag_tax,bsxfun(@rdivide,postStimR,base_tkern)');
xlim(xr);
caxis([0.2 1.8]);
xlabel('Time since saccade (s)');
ylabel('Stimulus latency (s)');


% fig_width = 4; rel_height = 1.6;
% figufy(f4);
% fname = [fig_dir 'ppsim_effK.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);
% 
% figufy(f5);
% fname = [fig_dir 'ppsim_effg.pdf'];
% exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f5);