clear all
close all

fig_dir = '/home/james/Analysis/bruce/ET_final/';

% Expt_name = 'G093';
Expt_name = 'M270';
Expt_num = str2num(Expt_name(2:end));

if strcmp(Expt_name,'M270')
    scale_fac = 1.72;
else
    scale_fac = 1;
end

%%
anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
cd(anal_dir);
if Expt_name(1) == 'G'
    anal_name = 'monoc_eyecorr_hbar2.mat';
    mod_data_name = 'monoc_eyecorr_hbar_mods.mat';
else
%     anal_name = 'monoc_eyecorr_Cprior.mat';
    anal_name = 'monoc_eyecorr2_Cprior.mat';
    mod_data_name = 'monoc_eyecorr_mods.mat';
end
load(anal_name,'it_mods','dit_mods_LOO','et_tr_set','et_params');
load(mod_data_name,'all_mod_SU*');

su_numbers = all_mod_SUnum(et_tr_set(all_mod_SUnum(et_tr_set) > 0));
% su_inds = find(ismember(all_mod_SUnum,su_numbers));
su_inds = et_tr_set;
fprintf('%d SUs\n',length(su_inds));
%%
%FOR M270
close all
cur_su = 28; 
xr = [-0.5 0.15];
pre_mod = it_mods{1}(su_inds(cur_su));
post_mod = dit_mods_LOO{cur_su,end}(su_inds(cur_su));

%for G093
% cur_su = 6; 
% cur_unit = find(all_mod_SUnum == su_numbers(cur_su));
% xr = [-0.45 0.45];
% pre_mod = it_mods{1}(cur_unit);
% post_mod = dit_mods_LOO{cur_su,end}(cur_unit);

flen = et_params.flen;
dt = et_params.dt;
spatial_usfac = et_params.spatial_usfac;
use_nPix_us = et_params.use_nPix*spatial_usfac;

if Expt_num == 270 || Expt_num >= 281
   load(anal_name,'et_sp_dx');
   sp_dx = et_sp_dx;
else
sp_dx = 0.0565/spatial_usfac/scale_fac;
end


lag_axis = (0:flen-1)*dt*1e3;
pix_axis = (1:use_nPix_us)*sp_dx - use_nPix_us/2*sp_dx;

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

h=figure;
for ii = 1:n_stimfilts
    subplot(n_stimfilts,3,(ii-1)*3+1);
    imagesc(pix_axis(disp_pix),lag_axis,squeeze(cor_filts(:,disp_pix,ii))); caxis([-max_vals(ii) max_vals(ii)]);
    line(pix_axis([1 end]),lag_axis(best_lags_cor([ii ii])),'color','r');
    xlabel('Position (deg)');
    ylabel('Lag (ms)');
    set(gca,'xtick',-0.5:0.25:0.5);
    set(gca,'ydir','normal');
    xlim(xr)
    
    subplot(n_stimfilts,3,(ii-1)*3+2);
    imagesc(pix_axis(disp_pix),lag_axis,squeeze(uncor_filts(:,disp_pix,ii))); caxis([-uc_max_vals(ii) uc_max_vals(ii)]);
    line(pix_axis([1 end]),lag_axis(best_lags_cor([ii ii])),'color','r');
    xlabel('Position (deg)');
    ylabel('Lag (ms)');
    set(gca,'xtick',-0.5:0.25:0.5);
    set(gca,'ydir','normal');
xlim(xr)
    
    subplot(n_stimfilts,3,(ii-1)*3+3);
    plot(pix_axis(disp_pix),spatial_profiles_cor(disp_pix,ii),'linewidth',1);
    hold on
    plot(pix_axis(disp_pix),spatial_profiles_uncor(disp_pix,ii),'r','linewidth',1);
    xlim(pix_axis([1 end]));
    xlabel('Position (deg)');
    ylabel('Filter amp (AU)');
    xlim(pix_axis(disp_pix([1 end])));
    yl = ylim(); ym = max(abs(yl)); ylim([-ym ym]);
    line(pix_axis(disp_pix([1 end])),[0 0],'color','k','linestyle','--');
xlim(xr)
    
end
colormap(gray);
subplot(n_stimfilts,3,1);
title('Corrected');
subplot(n_stimfilts,3,2);
title('Uncorrected');
%%
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 0.8;

figufy(h);
fname = [fig_dir sprintf('%s_SU%d_exampmod3.pdf',Expt_name,cur_su)];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close;
