clear all
close all

addpath('~/other_code/fastBSpline/');

fig_dname = 'tbt_fig_data';
fig_dir = '/home/james/Analysis/bruce/variability/figures/';

Expt_name = 'M012';
anal_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
cd(anal_dir)
load(fig_dname);

%% get basic info
nf = size(fig_data.tbt_EP,1); %number of used frames
dt = fig_data.EP_params.direct_bin_dts(1); %time res
tax = (1:nf)*dt; %make time axis

ex_su = 12; %[3 12] select example SU
ep_range = [-0.35 0.35]; %eye position range for plotting

trange = [2 3]; %range of times for plotting
used_tinds = find(tax >= trange(1) & tax <= trange(2));

binned_spks = squeeze(fig_data.binned_spks(:,:,ex_su)); %version without nans due to sacs/blinks
binned_spks_nan = squeeze(fig_data.binned_spks_nan(:,:,ex_su)); %version with nans
tbt_EP = fig_data.tbt_EP;

%get rid of trials where rptframes create any nan values in the plotting
%window
bad_trials = find(any(isnan(binned_spks(used_tinds,:))));
binned_spks(:,bad_trials) = [];
binned_spks_nan(:,bad_trials) = [];
tbt_EP(:,bad_trials) = []; 

%create a nanned version of the EP data to match the spk data
tbt_EP_nan = tbt_EP;
tbt_EP_nan(isnan(binned_spks_nan)) = nan;

n_rpts = size(binned_spks,2);

%% plot array of spiking data
f1 = figure();
imagescnan(tax,1:n_rpts,binned_spks');
xlim(trange);
caxis([0 4]);
colorbar
xlabel('Time (s)');
ylabel('Trial number');
set(gca,'ydir','normal');

% fname = [fig_dir sprintf('Direxamp_SU%d_spkdata_%s.pdf',ex_su,Expt_name)];
% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% plot array of EP data
f1 = figure();
imagescnan(tax,1:n_rpts,tbt_EP_nan');
xlim(trange);
colorbar
caxis(ep_range);
xlabel('Time (s)');
ylabel('Trial number');
set(gca,'ydir','normal');

% fname = [fig_dir sprintf('Direxamp_EPdata_%s.pdf',Expt_name)];
% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
%%
deltaX_range = [0 0.5];
cur_EP_data = fig_data.EP_data(ex_su,1);
psth_var = cur_EP_data.pair_psth_var;
ball_var = cur_EP_data.eps_ball_var(2);

f1 = figure();
plot(cur_EP_data.EP_bin_centers,cur_EP_data.var_ep_binned);
xlim(deltaX_range);
line(deltaX_range,psth_var + [0 0],'color','r');
line(deltaX_range,ball_var + [0 0],'color','k');
xlabel('Delta X (deg)');
ylabel('Rate covariance');
if ex_su == 12
    ylim([-0.025 0.25]);
end

% fname = [fig_dir sprintf('Direxamp_SU%d_varfun_%s.pdf',ex_su,Expt_name)];
% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% plot example trial spiking and EP data
ex_trials = [121 162];
ep_range = [-0.4 0.4]; %eye position range for plotting

tax_edges = [tax(1) - median(diff(tax))/2 tax(1:end-1) + median(diff(tax))/2];
f1 = figure();
subplot(2,1,1); hold on
% plot(tax,binned_spks(:,ex_trials(1)));
stairs(tax_edges,binned_spks(:,ex_trials));
xlim(trange);
xlabel('Time (s)');
ylabel('Spike count');

subplot(2,1,2); hold on
plot(tax,tbt_EP(:,ex_trials));
xlim(trange); ylim(ep_range);
xlabel('Time (s)');
ylabel('Eye position (deg)');

fname = [fig_dir sprintf('Direxamp_SU%d_extrials_%s.pdf',ex_su,Expt_name)];
fig_width = 4; rel_height = 1.6;
figufy(f1);
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%%
% close all
% for ii = 1:n_rpts
%     ii
%    subplot(2,1,1);
%    plot(tax,binned_spks(:,ii));
%    xlim(trange);
%    subplot(2,1,2)
%    plot(tax,tbt_EP(:,ii)); hold on
%    plot(tax,tbt_EP(:,162),'r'); hold on
% %    plot(tax,tbt_EP_nan(:,ii),'r');
%    
%    xlim(trange);
%    ylim(ep_range);
%    pause
%    clf
% end
