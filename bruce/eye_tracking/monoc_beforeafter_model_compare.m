clear all
close all

% expt_list = [85 86 87 88 89 93 95];
expt_num = 93;
probe_num = 94;

load ~/Analysis/bruce/summary_analysis/su_data.mat

fig_dir = '~/Analysis/bruce/summary_analysis/eyetrack_figs/';

%%
fprintf('Loading data for expt %d\n',expt_num);
anal_dir = ['~/Analysis/bruce/' sprintf('G0%d',expt_num)];
hor_data_name = [anal_dir '/monoc_eyecorr_hbar_finf.mat'];

cur_expt_id = find(su_data.expt_nums == expt_num);
su_probes = find(su_data.is_su(cur_expt_id,:));

%load in SU models for horizontal bars
if exist(hor_data_name,'file')
    load(hor_data_name);
    hor_tr_set = et_tr_set;
    hor_su_inds = find(et_tr_set > 96);
    hor_su_probeinds = find(ismember(1:length(su_probes),hor_tr_set(hor_su_inds)-96));
    
    hor_before_mods = it_mods{1}(hor_tr_set(hor_su_inds));
    hor_after_mods = dit_xv_mods{2}(hor_tr_set(hor_su_inds));
    hor_init_LL_imps = it_LLimp(1,hor_tr_set(hor_su_inds));
    hor_fin_LL_imps = dit_xv_LLimp(2,hor_tr_set(hor_su_inds));
    hor_init_true_xvLL = it_truexvLLimp(1,hor_tr_set(hor_su_inds));
    hor_fin_true_xvLL = dit_truexvLLimp(2,hor_tr_set(hor_su_inds));
    
else
    hor_tr_set = [];
    hor_su_inds = [];
    hor_su_probeinds = [];
end


%% GENERATE FULL MODEL COMPARISON FIGURE
cd('~/Analysis/bruce/summary_analysis/eyetrack_figs/');
close all
cur_ind = find(su_probes(hor_su_probeinds) == probe_num);

pre_mod = hor_before_mods(cur_ind);
post_mod = hor_after_mods(cur_ind);

% NMMdisplay_model(pre_mod,[],[],1);
% NMMdisplay_model(post_mod,[],[],1)

flen = et_params.flen;
dt = et_params.dt;
spatial_usfac = et_params.spatial_usfac;
use_nPix_us = et_params.use_nPix*spatial_usfac;
sp_dx = 0.0565/spatial_usfac;

lag_axis = (0:flen-1)*dt*1e3;
pix_axis = (1:use_nPix_us)*sp_dx - use_nPix_us/2*sp_dx;
% pix_axis = (-floor(use_nPix_us/2):floor(use_nPix_us/2))*sp_dx;

Xtargs = [post_mod.mods(:).Xtarget];
n_stimfilts = sum(Xtargs == 1);
n_squared_filts = sum(Xtargs == 1) - 1;
cor_filts = [post_mod.mods((Xtargs == 1)).filtK];
uncor_filts = [pre_mod.mods((Xtargs == 1)).filtK];
max_vals = max(abs(cor_filts));
uc_max_vals = max(abs(uncor_filts));

cor_filts = reshape(cor_filts,[flen use_nPix_us n_squared_filts+1]);
uncor_filts = reshape(uncor_filts,[flen use_nPix_us n_squared_filts+1]);

cor_temp_profiles = squeeze(var(cor_filts,[],2));
uncor_temp_profiles = squeeze(var(uncor_filts,[],2));
[~,best_lags_cor] = max(cor_temp_profiles);
[~,best_lags_uncor] = max(uncor_temp_profiles);

% best_lags_cor = [5 5 5];

for i = 1:n_stimfilts
    spatial_profiles_cor(:,i) = squeeze(cor_filts(best_lags_cor(i),:,i));
    spatial_profiles_uncor(:,i) = squeeze(uncor_filts(best_lags_cor(i),:,i));
end

disp_pix = find(pix_axis >= -0.5 & pix_axis <= 0.5);

figure
for ii = 1:n_stimfilts
    subplot(n_stimfilts,3,(ii-1)*3+1);
    imagesc(pix_axis(disp_pix),lag_axis,squeeze(cor_filts(:,disp_pix,ii))); caxis([-max_vals(ii) max_vals(ii)]);
    line(pix_axis([1 end]),lag_axis(best_lags_cor([ii ii])),'color','r');
    xlabel('Position (deg)','fontsize',12);
    ylabel('Lag (ms)','fontsize',12);
    set(gca,'xtick',-0.5:0.25:0.5,'fontsize',10,'fontname','arial');
    set(gca,'ydir','normal');
%     xlim([-0.5 0.5]);
    
    subplot(n_stimfilts,3,(ii-1)*3+2);
    imagesc(pix_axis(disp_pix),lag_axis,squeeze(uncor_filts(:,disp_pix,ii))); caxis([-uc_max_vals(ii) uc_max_vals(ii)]);
    line(pix_axis([1 end]),lag_axis(best_lags_cor([ii ii])),'color','r');
    xlabel('Position (deg)','fontsize',12);
    ylabel('Lag (ms)','fontsize',12);
    set(gca,'xtick',-0.5:0.25:0.5,'fontsize',10,'fontname','arial');
    set(gca,'ydir','normal');
    xlim([-0.5 0.5]);
    
    subplot(n_stimfilts,3,(ii-1)*3+3);
    plot(pix_axis(disp_pix),spatial_profiles_cor(disp_pix,ii));
    hold on
    plot(pix_axis(disp_pix),spatial_profiles_uncor(disp_pix,ii),'r');
    xlim(pix_axis([1 end]));
    xlabel('Position (deg)','fontsize',12);
    ylabel('Filter amp (AU)','fontsize',12);
    set(gca,'fontsize',10,'fontname','arial');
    box off
    xlim([-0.5 0.5]);
end
colormap(gray);
subplot(n_stimfilts,3,1);
title('Corrected','fontsize',8);
subplot(n_stimfilts,3,2);
title('Uncorrected','fontsize',8);
fillPage(gcf,'papersize',[8 7]);
fname = 'before_after_SU_exampmod';

%% EXAMPLE PRE-CORRECTION FILTERS
% figure
% imagesc(pix_axis,lag_axis,squeeze(uncor_filts(:,:,1)));
% xlabel('Position (deg)','fontsize',10);
% ylabel('Lag (s)','fontsize',10);
% set(gca,'ydir','normal');
% colormap(gray);
% xlim([-0.5 0.5]);
% ylim([0-dt/2 0.085+dt/2]);
% set(gca,'xtick',-0.5:0.25:0.5,'ytick',0:0.02:0.08,'fontsize',8,'fontname','arial');
% fillPage(gcf,'papersize',[4 4]);
% fname = 'examp_uncorr_linfilt';
% 
% figure
% imagesc(pix_axis,lag_axis,squeeze(uncor_filts(:,:,2)));
% xlabel('Position (deg)','fontsize',10);
% ylabel('Lag (s)','fontsize',10);
% set(gca,'ydir','normal');
% colormap(gray);
% xlim([-0.5 0.5]);
% ylim([0-dt/2 0.085+dt/2]);
% set(gca,'xtick',-0.5:0.25:0.5,'ytick',0:0.02:0.08,'fontsize',8,'fontname','arial');
% fillPage(gcf,'papersize',[4 4]);
% fname = 'examp_uncorr_sqfilt';


%%
zpad_factor = 5;
Ft = linspace(-1/dt/2,1/dt/2,zpad_factor*flen);
Fx = linspace(-1/sp_dx/2,1/sp_dx/2,zpad_factor*use_nPix_us);

% sig_Ft = 0.05; sig_Fx = 0.1;
sig_Ft = 0.1; sig_Fx = 0.2;
[FFx,FFt] = meshgrid(Fx,Ft);
gauss_kern = FFt.^2/(2*sig_Ft^2) + FFx.^2/(2*sig_Fx^2);
gauss_kern = exp(-gauss_kern);
gauss_kern = gauss_kern/sum(gauss_kern(:));

uncor_filts_zpad = zeros(zpad_factor*flen,zpad_factor*use_nPix_us,n_stimfilts);
uncor_filts_zpad((flen+1):2*flen,(use_nPix_us+1):2*use_nPix_us,:) = uncor_filts;
uncor_ffts = zeros(size(uncor_filts_zpad));
for ii = 1:n_squared_filts+1
    cur_ffts = abs(fftshift(fft2(uncor_filts_zpad(:,:,ii))));
    cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
    uncor_ffts(:,:,ii) = cur_ffts;
end
cor_filts_zpad = zeros(zpad_factor*flen,zpad_factor*use_nPix_us,n_stimfilts);
cor_filts_zpad((flen+1):2*flen,(use_nPix_us+1):2*use_nPix_us,:) = cor_filts;
cor_ffts = zeros(size(cor_filts_zpad));
for ii = 1:n_squared_filts+1
    cur_ffts = abs(fftshift(fft2(cor_filts_zpad(:,:,ii))));
    cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
    cor_ffts(:,:,ii) = cur_ffts;
end

cur_pspec = squeeze(cor_ffts(:,:,1));
[xloc,yloc] = ind2sub([length(Ft),length(Fx)],find(cur_pspec == max(cur_pspec(:))));
figure
subplot(2,1,1)
imagesc(Fx,Ft,cur_pspec);
hold on
plot(Fx(yloc),Ft(xloc),'wo');
plot(-Fx(yloc),-Ft(xloc),'wo');
xlabel('Spatial frequency (cyc/deg)','fontsize',12);
ylabel('Temporal frequency (Hz)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
title('Corrected','fontsize',10);
xlim([-8 8]);
ylim([-30 30]);
line([-12 12],[0 0],'color','w');
line([0 0],[-40 40],'color','w');
ca = caxis();
set(gca,'ydir','normal');

cur_pspec = squeeze(uncor_ffts(:,:,1));
[xloc,yloc] = ind2sub([length(Ft),length(Fx)],find(cur_pspec == max(cur_pspec(:))));
subplot(2,1,2)
imagesc(Fx,Ft,cur_pspec);
hold on
plot(Fx(yloc),Ft(xloc),'wo');
plot(-Fx(yloc),-Ft(xloc),'wo');
xlabel('Spatial frequency (cyc/deg)','fontsize',12);
ylabel('Temporal frequency (Hz)','fontsize',12);
xlim([-8 8]);
ylim([-30 30]);
set(gca,'fontsize',10,'fontname','arial');
title('Uncorrected','fontsize',10);
fillPage(gcf,'papersize',[3 5]);
colormap(gray);
line([-12 12],[0 0],'color','w');
line([0 0],[-40 40],'color','w');
caxis(ca);
set(gca,'ydir','normal');
colormap(jet)
fillPage(gcf,'papersize',[4 8]);
fname = 'exampmod_lin_powspec';


cur_pspec = squeeze(mean(cor_ffts(:,:,2:end),3));
[xloc,yloc] = ind2sub([length(Ft),length(Fx)],find(cur_pspec == max(cur_pspec(:))));
figure
subplot(2,1,1)
imagesc(Fx,Ft,cur_pspec);
hold on
plot(Fx(yloc),Ft(xloc),'wo');
plot(-Fx(yloc),-Ft(xloc),'wo');
xlabel('Spatial frequency (cyc/deg)','fontsize',10);
ylabel('Temporal frequency (Hz)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
title('Corrected','fontsize',12);
set(gca,'ydir','normal');
xlim([-8 8]);
ylim([-30 30]);
ca = caxis();
line([-12 12],[0 0],'color','w');
line([0 0],[-40 40],'color','w');

cur_pspec = squeeze(mean(uncor_ffts(:,:,2:end),3));
[xloc,yloc] = ind2sub([length(Ft),length(Fx)],find(cur_pspec == max(cur_pspec(:))));
subplot(2,1,2)
imagesc(Fx,Ft,squeeze(mean(uncor_ffts(:,:,2:end),3)));
hold on
plot(Fx(yloc),Ft(xloc),'wo');
plot(-Fx(yloc),-Ft(xloc),'wo');
xlabel('Spatial frequency (cyc/deg)','fontsize',12);
ylabel('Temporal frequency (Hz)','fontsize',12);
xlim([-8 8]);
ylim([-30 30]);
set(gca,'fontsize',10,'fontname','arial');
title('Uncorrected');
caxis(ca);
fillPage(gcf,'papersize',[3 5]);
colormap(gray);
line([-12 12],[0 0],'color','w');
line([0 0],[-40 40],'color','w');
set(gca,'ydir','normal');
colormap(jet)
fillPage(gcf,'papersize',[4 8]);
fname = 'exampmod_sq_powspec';


