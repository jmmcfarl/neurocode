clear all
close all

addpath('~/other_code/fastBSpline/');

fig_dname = 'tbt_fig_data';
fig_dir = '/home/james/Analysis/bruce/variability/figures/';

Expt_name = 'M012';
anal_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
cd(anal_dir)
load(fig_dname);

%%
n_rpts = size(fig_data.tbt_EP,2);
nf = size(fig_data.tbt_EP,1);
dt = fig_data.EP_params.direct_bin_dts(1);
tax = (1:nf)*dt;

ex_su = 12;

trange = [1.5 3];
used_tinds = find(tax >= trange(1) & tax <= trange(2));
ep_range = [-0.35 0.35];

binned_spks = squeeze(fig_data.binned_spks(:,:,ex_su));
binned_spks_nan = squeeze(fig_data.binned_spks_nan(:,:,ex_su));
tbt_EP = fig_data.tbt_EP;

bad_trials = find(any(isnan(binned_spks(used_tinds,:))));
binned_spks(:,bad_trials) = [];
binned_spks_nan(:,bad_trials) = [];
tbt_EP(:,bad_trials) = [];

tbt_EP_nan = tbt_EP;
tbt_EP_nan(isnan(binned_spks_nan)) = nan;

n_rpts = size(binned_spks,2);
%%
f1 = figure();
imagescnan(tax,1:n_rpts,binned_spks');
xlim(trange);
colorbar

f2 = figure();
imagescnan(tax,1:n_rpts,tbt_EP_nan');
xlim(trange);
colorbar
caxis(ep_range);

%%
close all
for ii = 1:n_rpts
   subplot(2,1,1);
   plot(tax,binned_spks(:,ii));
   xlim(trange);
   subplot(2,1,2)
   plot(tax,tbt_EP(:,ii)); hold on
   plot(tax,tbt_EP_nan(:,ii),'r');
   
   xlim(trange);
   pause
   clf
end
