clear all
close all

addpath('~/other_code/fastBSpline/');

fig_dname = 'tbt_fig_data';
fig_dir = '/home/james/Analysis/bruce/variability/figures/';

Expt_name = 'M012';
bar_ori = 0;
rec_number = 1;

anal_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
cd(anal_dir)
load(fig_dname);

data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
if rec_number > 1
    data_name = strcat(data_name,sprintf('_r%d',rec_number));
end
fprintf('Loading %s\n',data_name);
load(data_name);


%% parse EP data
nf = size(fig_data.tbt_EP,1);%total number of frames
dt = fig_data.EP_params.direct_bin_dts(1); %time res
tax = (1:nf)*dt; %time axis
trange = [1 4]; %range of time points to use (sec)

tbt_EP = fig_data.tbt_EP;
used_tinds = find(tax >= trange(1) & tax <= trange(2));
tbt_EP = tbt_EP(used_tinds,:);
bad_trials = find(any(isnan(tbt_EP))); %get rid of trials where there are NANs during the used portion

tbt_EP(:,bad_trials) = [];
n_rpts = size(tbt_EP,2);

seed = 1;
rng(seed);
tbt_EP = tbt_EP(:,randperm(n_rpts));

%% parse stim params
n_uf = length(used_tinds);
tax = (1:n_uf)*dt; %time axis
npix = params.full_nPix;
usfac = fig_data.modFitParams.spatial_usfac*fig_data.modFitParams.add_usfac;
use_nPix = fig_data.modFitParams.use_nPix_us/fig_data.modFitParams.spatial_usfac;
dds = mode(expt_data.expt_dds);
flen = fig_data.modFitParams.flen;

%create RLS stimulus for repeat trials
seed = 1;
stim = generate_RLS_stim(n_uf,npix,dds,usfac,seed);
stim = repmat(stim,n_rpts,1);

buffer_pix = floor((npix - use_nPix)/2);
%repeat for up-sampled versions of the Xmatrix
[Xinds_up,~] = meshgrid(1/usfac:1/usfac:npix,1:flen);
cur_use_pix = (1/usfac:1/usfac:use_nPix) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

%% make stim Xmats
sp_dx = mode(expt_data.expt_dw)/usfac;
all_EP_shifts = round(tbt_EP(:)/sp_dx);

%make EP_shifted stim
stim_shift = stim;
for ii=1:length(all_EP_shifts) %correct repeat trial data
    stim_shift(ii,:) = shift_matrix_Nd(stim(ii,:),-all_EP_shifts(ii),2);
end

%build xmats
stim_params = NMMcreate_stim_params([fig_data.modFitParams.flen npix*fig_data.modFitParams.spatial_usfac]);
all_Xmat_shift = create_time_embedding(stim_shift,stim_params); %incorporate time embedding for whole-trials on repeats
all_Xmat_shift = all_Xmat_shift(:,use_kInds_up); %take only usable repeat-trial data

all_Xmat = create_time_embedding(stim,stim_params); %incorporate time embedding for whole-trials on repeats
all_Xmat = all_Xmat(:,use_kInds_up); %take only usable repeat-trial data

%if using tent-basis spatial-upsampling
if fig_data.modFitParams.add_usfac > 1
    all_Xmat_shift = tb_proc_stim(all_Xmat_shift,fig_data.modFitParams.add_usfac,flen);
    all_Xmat = tb_proc_stim(all_Xmat,fig_data.modFitParams.add_usfac,flen);
end

%% make SU simulated firing rates
n_SUs = size(fig_data.EP_data,1);
mod_prates_shift = nan(length(all_EP_shifts),n_SUs);
mod_prates = nan(length(all_EP_shifts),n_SUs);
for cc = 1:n_SUs
    if ~isempty(fig_data.EP_data(cc,1).useMod)
        [~,~,mod_prates_shift(:,cc)] = NMMmodel_eval(fig_data.EP_data(cc,1).useMod,[],all_Xmat_shift);
        [~,~,mod_prates(:,cc)] = NMMmodel_eval(fig_data.EP_data(cc,1).useMod,[],all_Xmat);
    end
end
mod_prates = reshape(mod_prates,n_uf,n_rpts,n_SUs);
mod_prates_shift = reshape(mod_prates_shift,n_uf,n_rpts,n_SUs);

%get rid of a short buffer at the beginning of the trial to account for the
%lack of stim-history there
mod_prates(1:flen,:,:) = [];
mod_prates_shift(1:flen,:,:) = [];
n_uf = size(mod_prates,1);
tax = (0:(n_uf-1))*dt; %time axis
tbt_EP(1:flen,:) = [];

%adjust to desired avg spk rates
meanrates = mean(reshape(mod_prates_shift,[],n_SUs));
target_mrate = 0.4;
mod_prates = bsxfun(@times,mod_prates,reshape(target_mrate./meanrates,1,1,n_SUs));
mod_prates_shift = bsxfun(@times,mod_prates_shift,reshape(target_mrate./meanrates,1,1,n_SUs));
mod_prates = mod_prates/dt;
mod_prates_shift = mod_prates_shift/dt;

%% compute tbt-stats
atvars = squeeze(nanvar(mod_prates_shift,[],2));
psths = squeeze(nanmean(mod_prates_shift,2));
PF_psth = squeeze(mod_prates(:,1,:));
psthvars = nanvar(psths);
totvars = nanvar(reshape(mod_prates_shift,[],n_SUs));
PF_psthvars = nanvar(PF_psth);
alphas = psthvars./totvars;
actual_alphas = arrayfun(@(x) x.pair_psth_var./x.eps_ball_var(2),fig_data.EP_data(:,1));

%% select example unit and plotting ranges
ex_su = 3;
ex_tind = 209;
ep_range = [-0.35 0.35];
plot_trange = [0 1];
plot_trials = [1 150];

%% plot EP fig
ex_trials = [9 83 157];
f1 = figure();
subplot(2,1,1); hold on
plot(tax,tbt_EP(:,ex_trials(1)),'r');
plot(tax,tbt_EP(:,ex_trials(2)),'b');
plot(tax,tbt_EP(:,ex_trials(3)),'k');
xlim(plot_trange); ylim(ep_range);
xlabel('Time (s)'); ylabel('Eye position (deg)');

nbins = 25;
bin_edges = linspace(ep_range(1),ep_range(2),nbins+1);
n = histc(tbt_EP(:),bin_edges);
bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);
subplot(2,1,2);
hobj = bar(bin_cents,n(1:end-1));
set(hobj,'faceColor','k','barwidth',1);
xlim(ep_range);
xlabel('Eye position (deg)');
ylabel('Relative frequency');
set(gca,'ytick',[]);

fname = [fig_dir sprintf('Illust_EPdata_%s.pdf',Expt_name)];
fig_width = 4; rel_height = 1.6;
figufy(f1);
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% plot psths
f1 = figure();

subplot(2,2,1); hold on
plot(tax,psths(:,ex_su)); 
xlim(plot_trange);
ylim([0 200]);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');

subplot(2,2,2); hold on
plot(tax,PF_psth(:,ex_su),'r');
xlim(plot_trange);
ylim([0 200]);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');

subplot(2,2,3); hold on
plot(tax,atvars(:,ex_su));
xlim(plot_trange);
ylim([0 1e4]);
xlabel('Time (s)');
ylabel('Rate variance (Hz^2)');

subplot(2,2,4); hold on
plot(tax,zeros(size(tax)),'r');
xlim(plot_trange);
ylim([0 1e4]);
xlabel('Time (s)');
ylabel('Rate variance (Hz^2)');

fname = [fig_dir sprintf('Illust_PSTHs_%s.pdf',Expt_name)];
fig_width = 6; rel_height = 0.8;
figufy(f1);
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% plot tbt firing rates

f1 = figure();
subplot(2,1,1);
imagesc(tax,1:n_rpts,squeeze(mod_prates(:,:,ex_su))');
xlim(plot_trange); ylim(plot_trials);
set(gca,'ydir','normal');
caxis([0 200]); colorbar
xlabel('Time (s)'); ylabel('Trial number');
subplot(2,1,2);
imagesc(tax,1:n_rpts,squeeze(mod_prates_shift(:,:,ex_su))');
xlim(plot_trange); ylim(plot_trials);
caxis([0 200]); colorbar
xlabel('Time (s)'); ylabel('Trial number');
set(gca,'ydir','normal');

fname = [fig_dir sprintf('Illust_tbtRates_%s.pdf',Expt_name)];
fig_width = 4; rel_height = 1.8;
figufy(f1);
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%%
ex_time = 0.7;
ex_tind = find(tax >= ex_time,1);
rate_vals = squeeze(mod_prates_shift(ex_tind,:,ex_su));

f1 = figure();
nbins = 25;
bin_edges = linspace(0,200,nbins+1);
n = histc(rate_vals,bin_edges);
bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);
hobj = bar(bin_cents,n(1:end-1));
set(hobj,'faceColor','k','barwidth',1);
xlim([0 200]);
xlabel('Firing rate (Hz)');
ylabel('Relative frequency');
set(gca,'ytick',[]);


fname = [fig_dir sprintf('Illust_ratedist_%s.pdf',Expt_name)];
fig_width = 4; rel_height = 0.8;
figufy(f1);
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% plot tbt-EP
% f1 = figure();
% imagesc(tax,1:n_rpts,tbt_EP');
% xlim(plot_trange); ylim(plot_trials); 
% caxis(ep_range);

%% plot simulated spike rasters
dt_usfac = 5; %temporal up-sampling factor
tax_up = (1:1/dt_usfac:n_uf)*dt; %new time axis
dt_up = dt/dt_usfac;

shift_prate_up = interp1(tax,squeeze(mod_prates_shift(:,:,ex_su)),tax_up);
prate_up = interp1(tax,squeeze(mod_prates(:,:,ex_su)),tax_up);

shift_binned_spks = poissrnd(shift_prate_up*dt/dt_usfac);
shift_binned_spks(shift_binned_spks > 1) = 1;
binned_spks = poissrnd(prate_up*dt/dt_usfac);
binned_spks(binned_spks > 1) = 1;

tick_height = 0.7;

f1 = figure(); 
subplot(1,2,1); hold on
for ii = 1:n_rpts
   cur_spk_times = tax_up(binned_spks(:,ii) == 1);
   cur_spk_times(cur_spk_times <= plot_trange(1) | cur_spk_times >= plot_trange(2)) = [];
    for jj = 1:length(cur_spk_times)
       line(cur_spk_times(jj) + [0 0],ii + [0 tick_height],'color','k'); 
    end
end
xlabel('Time (s)');
ylabel('Trial');
xlim(plot_trange);
ylim(plot_trials);

subplot(1,2,2); hold on
for ii = 1:n_rpts
   cur_spk_times = tax_up(shift_binned_spks(:,ii) == 1);
   cur_spk_times(cur_spk_times <= plot_trange(1) | cur_spk_times >= plot_trange(2)) = [];
    for jj = 1:length(cur_spk_times)
       line(cur_spk_times(jj) + [0 0],ii + [0 tick_height],'color','k'); 
    end
end
xlabel('Time (s)');
ylabel('Trial');
xlim(plot_trange);
ylim(plot_trials);

fname = [fig_dir sprintf('Illust_rasters_%s.pdf',Expt_name)];
fig_width = 6; rel_height = 0.5;
figufy(f1);
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

















%%
% prate_mat = squeeze(fig_data.pred_rates(:,:,ex_su)); %tbt firing rates from model
% PF_rate = fig_data.PF_prates(:,ex_su); %PSTH if perfect fixation
% 
% %% plot tbt EP and firing rate traces
% tset = [1 2 22];
% cmap = jet(length(tset));
% cmap(2,:) = [0 1 0];
% cmap(3,:) = [1 0 0];
% 
% f1 = figure();
% for ii = 1:length(tset)
%     subplot(2,1,1); hold on
%     plot(tax,tbt_EP(:,tset(ii)),'color',cmap(ii,:));
%     subplot(2,1,2); hold on
%     plot(tax,prate_mat(:,tset(ii))/dt,'color',cmap(ii,:));
% end
% subplot(2,1,1); 
% xlim(trange); ylim(ep_range);
% xlabel('Time (s)');
% ylabel('Eye position (deg)');
% line(trange,[0 0],'color','k','linestyle','--');
% 
% subplot(2,1,2);
% xlim(trange);
% xlabel('Time (s)');
% ylabel('Firing rate (Hz)');
% yl = ylim();
% line(tax([ex_tind ex_tind]),yl,'color','k','linestyle','--');
% 
% % fname = [fig_dir sprintf('Illust_examples_%s.pdf',Expt_name)];
% % fig_width = 4; rel_height = 1.6;
% % figufy(f1);
% % exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);
% 
% %%
% close all
% f1 = figure();
% subplot(2,1,1);
% imagescnan(tax,1:n_rpts,prate_mat'/dt);
% xlim(trange)
% caxis([0 100]);
% subplot(2,1,2);
% imagescnan(tax,1:n_rpts,repmat(PF_rate'/dt,n_rpts,1));
% xlim(trange);
% caxis([0 100]);
% 
% emp_PSTH = nanmean(prate_mat/dt,2);
% emp_var = nanvar(prate_mat/dt,[],2);
% f2 = figure();
% subplot(2,1,1);
% plot(tax,emp_PSTH); hold on
% plot(tax,PF_rate/dt,'r');
% xlim(trange);
% subplot(2,1,2)
% plot(tax,emp_var); hold on
% plot(tax,zeros(size(tax)),'r');
% xlim(trange)
% %% plot across-trial dist of rates at one time lag
% nbins = 20;
% bin_edges = linspace(0,max(prate_mat(ex_tind,:)/dt)+1,nbins+1);
% n = histc(prate_mat(ex_tind,:)/dt,bin_edges);
% bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);
% f1 = figure();
% hobj = bar(bin_cents,n(1:end-1));
% set(hobj,'faceColor','k','barwidth',0.9);
% xlabel('Firing rate (Hz)');
% ylabel('Counts');
% 
% fname = [fig_dir sprintf('Illust_rate_dist_%s.pdf',Expt_name)];
% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% %% plot distribution of eye positions
% nbins = 25;
% bin_edges = linspace(ep_range(1),ep_range(2),nbins+1);
% n = histc(tbt_EP(:),bin_edges);
% bin_cents = 0.5*bin_edges(1:end-1) + 0.5*bin_edges(2:end);
% 
% f1 = figure();
% hobj = bar(bin_cents,n(1:end-1));
% set(hobj,'faceColor','k','barwidth',1);
% xlim(ep_range);
% xlabel('Eye position (deg)');
% ylabel('Relative frequency');
% 
% fname = [fig_dir sprintf('Illust_EP_dist_%s.pdf',Expt_name)];
% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% %% plot psth
% cur_psth = nanmean(prate_mat,2);
% f1 = figure();
% plot(tax,cur_psth/dt,'k');
% xlim(trange);
% ylim(yl);
% xlabel('Time (s)');
% ylabel('Firing rate (Hz)');
% 
% fname = [fig_dir sprintf('Illust_PSTH_%s.pdf',Expt_name)];
% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% %% create simulated spike raster
% dt_usfac = 5; %temporal up-sampling factor
% tax_up = (1:1/dt_usfac:nf)*dt; %new time axis
% dt_up = dt/dt_usfac;
% 
% prate_mat_up = interp1(tax,prate_mat,tax_up);
% 
% binned_spks = poissrnd(prate_mat_up/dt_usfac);
% binned_spks(binned_spks > 1) = 1;
% 
% tick_height = 0.7;
% 
% f1 = figure(); hold on
% for ii = 1:n_rpts
%    cur_spk_times = tax_up(binned_spks(:,ii) == 1);
%    cur_spk_times(cur_spk_times <= trange(1) | cur_spk_times >= trange(2)) = [];
%     for jj = 1:length(cur_spk_times)
%        line(cur_spk_times(jj) + [0 0],ii + [0 tick_height],'color','k'); 
%     end
% end
% xlim(trange);
% xlabel('Time (s)');
% ylabel('Trial');
% 
% fname = [fig_dir sprintf('Illust_raster_%s.pdf',Expt_name)];
% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
