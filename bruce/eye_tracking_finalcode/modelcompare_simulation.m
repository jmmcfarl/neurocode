clear all
close all

fig_dir = '/home/james/Analysis/bruce/ET_final/';

Expt_name = 'G086';

if strcmp(Expt_name,'M270')
    scale_fac = 1.72;
else
    scale_fac = 1;
end
poss_sim_gains = [0.25 0.5 1 2 3];

%%
anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
cd(anal_dir);
rpt = 1;
sg = 5;

% anal_name = [anal_dir sprintf('ETsim_it%d_gain%d',rpt,sg)];
anal_name = [anal_dir sprintf('ETsim2_it%d_gain%d',rpt,sg)];
mod_data_name = 'monoc_eyecorr_hbar_mods.mat';
base_name = 'monoc_eyecorr_hbar.mat';
load(mod_data_name,'all_mod_SU*');
base_data = load(base_name,'et_params','dit_mods');
true_mods = base_data.dit_mods{end};
et_params = base_data.et_params;

fprintf('Loading %s\n',anal_name);
load(anal_name,'it_mods','dit_mods','tr_set');

su_numbers = all_mod_SUnum(tr_set(all_mod_SUnum(tr_set) > 0));
su_inds = find(tr_set > 96);
fprintf('%d SUs\n',length(su_inds));
%%
% close all
cur_su = 7;
pre_mod = it_mods{1}(su_inds(cur_su));
post_mod = dit_mods{end}(su_inds(cur_su));
base_mod = true_mods(find(all_mod_SUnum == su_numbers(cur_su)));

flen = et_params.flen;
dt = et_params.dt;
spatial_usfac = et_params.spatial_usfac;
use_nPix_us = et_params.use_nPix*spatial_usfac;
sp_dx = 0.0565/spatial_usfac/scale_fac;

lag_axis = (0:flen-1)*dt*1e3;
pix_axis = (1:use_nPix_us)*sp_dx - use_nPix_us/2*sp_dx;

Xtargs = [post_mod.mods(:).Xtarget];
n_stimfilts = sum(Xtargs == 1);
n_squared_filts = sum(Xtargs == 1) - 1;
base_filts = [base_mod.mods((Xtargs == 1)).filtK];
cor_filts = [post_mod.mods((Xtargs == 1)).filtK];
uncor_filts = [pre_mod.mods((Xtargs == 1)).filtK];
max_vals = max(abs(cor_filts));
uc_max_vals = max(abs(uncor_filts));
base_max_vals = max(abs(base_filts));

base_filts = reshape(base_filts,[flen use_nPix_us n_squared_filts+1]);
cor_filts = reshape(cor_filts,[flen use_nPix_us n_squared_filts+1]);
uncor_filts = reshape(uncor_filts,[flen use_nPix_us n_squared_filts+1]);

base_temp_profiles = squeeze(var(base_filts,[],2));
cor_temp_profiles = squeeze(var(cor_filts,[],2));
uncor_temp_profiles = squeeze(var(uncor_filts,[],2));
[~,best_lags_base] = max(base_temp_profiles);
[~,best_lags_cor] = max(cor_temp_profiles);
[~,best_lags_uncor] = max(uncor_temp_profiles);

% best_lags_cor = [5 5 5];

for i = 1:n_stimfilts
    spatial_profiles_base(:,i) = squeeze(base_filts(best_lags_base(i),:,i));
    spatial_profiles_cor(:,i) = squeeze(cor_filts(best_lags_base(i),:,i));
    spatial_profiles_uncor(:,i) = squeeze(uncor_filts(best_lags_base(i),:,i));
end

disp_pix = find(pix_axis >= -0.5 & pix_axis <= 0.5);

n_cols = 4;

h=figure;
for ii = 1:n_stimfilts
   
    subplot(n_stimfilts,n_cols,(ii-1)*n_cols+1);
    imagesc(pix_axis(disp_pix),lag_axis,squeeze(base_filts(:,disp_pix,ii))); caxis([-base_max_vals(ii) base_max_vals(ii)]);
    line(pix_axis([1 end]),lag_axis(best_lags_base([ii ii])),'color','r');
    xlabel('Position (deg)');
    ylabel('Lag (ms)');
    set(gca,'xtick',-0.5:0.25:0.5);
    set(gca,'ydir','normal');
%     xlim([-0.5 0.5]);

subplot(n_stimfilts,n_cols,(ii-1)*n_cols+2);
    imagesc(pix_axis(disp_pix),lag_axis,squeeze(uncor_filts(:,disp_pix,ii))); caxis([-uc_max_vals(ii) uc_max_vals(ii)]);
    line(pix_axis([1 end]),lag_axis(best_lags_base([ii ii])),'color','r');
    xlabel('Position (deg)');
    ylabel('Lag (ms)');
    set(gca,'xtick',-0.5:0.25:0.5);
    set(gca,'ydir','normal');
%     xlim([-0.5 0.5]);

subplot(n_stimfilts,n_cols,(ii-1)*n_cols+3);
    imagesc(pix_axis(disp_pix),lag_axis,squeeze(cor_filts(:,disp_pix,ii))); caxis([-max_vals(ii) max_vals(ii)]);
    line(pix_axis([1 end]),lag_axis(best_lags_base([ii ii])),'color','r');
    xlabel('Position (deg)');
    ylabel('Lag (ms)');
    set(gca,'xtick',-0.5:0.25:0.5);
    set(gca,'ydir','normal');
    %     xlim([-0.5 0.5]);
        
    subplot(n_stimfilts,n_cols,(ii-1)*n_cols+4); hold on
    plot(pix_axis(disp_pix),spatial_profiles_base(disp_pix,ii),'k','linewidth',1);
    plot(pix_axis(disp_pix),spatial_profiles_cor(disp_pix,ii),'b','linewidth',1);
    plot(pix_axis(disp_pix),spatial_profiles_uncor(disp_pix,ii),'r','linewidth',1);
    xlim(pix_axis([1 end]));
    xlabel('Position (deg)');
    ylabel('Filter amp (AU)');
    xlim(pix_axis(disp_pix([1 end])));
    yl = ylim(); ym = max(abs(yl)); ylim([-ym ym]);
    line(pix_axis(disp_pix([1 end])),[0 0],'color','k','linestyle','--');
    
end
colormap(gray);
subplot(n_stimfilts,n_cols,1);
title('Original');
subplot(n_stimfilts,n_cols,2);
title('Uncorrected');
subplot(n_stimfilts,n_cols,3);
title('Corrected');

%%
fig_width = 6.83; %3.27 4.86 6.83
rel_height = 0.6;

figufy(h);
fname = [fig_dir sprintf('%s_SU%d_SIMmod.pdf',Expt_name,cur_su)];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close;
