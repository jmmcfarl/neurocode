clear all

Expt_num = 239;

Expt_name = sprintf('M%d',Expt_num);

fin_anal_dir = ['~/Analysis/bruce/' Expt_name '/stim_mods/'];
fin_mod_name = 'corr_mods_quick';

%%
cd(fin_anal_dir)
load(fin_mod_name);

%%
stim_dims = cor_gqm(1).mod_fit.stim_params.stim_dims;
flen = stim_dims(1);
nPix = stim_dims(2);
n_filts = length(cor_gqm(1).mod_fit.mods);

n_probes = 24;

avg_tkern = nan(n_probes,flen);
for ss = 1:n_probes
    
   cur_filts = reshape([cor_gqm(ss).mod_fit.mods(1:n_filts).filtK],[flen nPix n_filts]); 
   filt_tkern = squeeze(std(cur_filts,[],2));
   avg_tkern(ss,:) = mean(filt_tkern,2);
end

rel_tkern = bsxfun(@rdivide,avg_tkern,max(avg_tkern,[],2));
stim_dt = 0.01;

%%
save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
cd(save_dir);
sname = 'lfp_trig_avg_data';
load(sname)
lfp_lags = lags;
lfp_dt = dt;
%%
sname = 'sac_trig_avg_data';
load(sname);

%%
close all
xr = [-0.1 0.3];

h = figure();
subplot(1,3,1)
imagesc(lfp_lags*lfp_dt,1:24,lfp_data.trial_onset_csd);
line([0 0],[0 24],'color','w');
ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]*0.8);
xlim(xr);
xlabel('Time (s)');
ylabel('Probe');
title('Trial-onset CSD');

subplot(1,3,2)
imagesc(lfp_lags*lfp_dt,1:24,lfp_data.gsac_csd);
ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]*0.8);
xlim(xr);
line([0 0],[0 24],'color','w');
xlabel('Time (s)');
ylabel('Probe');
title('Saccade-triggered CSD');

subplot(1,3,3)
imagesc(lags*dt,1:24,[mua_data(:).gsac_avg]'-1);
ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]*0.9);
xlim(xr);
line([0 0],[0 24],'color','w');
xlabel('Time (s)');
ylabel('Probe');
title('Saccade-triggered MUA');

% subplot(2,3,4)
% imagesc(lfp_lags*lfp_dt,1:24,lfp_data.trial_offset_csd);
% line([0 0],[0 24],'color','w');
% ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]*0.8);
% xlim(xr);
% 
% subplot(2,3,5)
% imagesc(lfp_lags*lfp_dt,1:24,lfp_data.msac_csd);
% ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]*0.8);
% xlim(xr);
% line([0 0],[0 24],'color','w');
% 
% subplot(2,3,6)
% imagesc(lags*dt,1:24,[mua_data(:).msac_avg]'-1);
% ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]*0.9);
% xlim(xr);
% line([0 0],[0 24],'color','w');

ename = sprintf('/home/james/Desktop/lab_meeting_figs/%s_LP_CSDMUA.png',Expt_name);
fig_width = 9; rel_height = 0.3;
figufy(h);
exportfig(h,ename,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
%%
% close all
h = figure();

mua_gsac = [mua_data(:).gsac_avg]';
mua_relgsac = bsxfun(@minus,mua_gsac,min(mua_gsac,[],2));
xr = [0 0.12];

% subplot(2,2,1)
% imagesc(lags*dt,1:24,[mua_data(:).gsac_avg]');
% % ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Probe');
% title('Saccade-trig avg');

subplot(1,2,1)
imagesc(lags*dt,1:24,mua_relgsac);
caxis([0 0.2]);
xlim(xr);
xlabel('Time (s)');
ylabel('Probe');
title('Saccade-trig relative');


% mua_msac = [mua_data(:).msac_avg]';
% mua_relmsac = bsxfun(@minus,mua_msac,min(mua_msac,[],2));
% xr = [0 0.12];
% 
% subplot(3,2,3)
% imagesc(lags*dt,1:24,[mua_data(:).msac_avg]');
% % ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
% xlim(xr);
% 
% subplot(3,2,4)
% imagesc(lags*dt,1:24,mua_relmsac);
% caxis([0 0.2]);
% xlim(xr);
% 

% subplot(2,2,3)
% imagesc((1:flen)*stim_dt-stim_dt/2,1:24,avg_tkern);
% xlim(xr)
% xlabel('Time (s)');
% ylabel('Probe');
% title('Temporal kernel');

subplot(1,2,2)
imagesc((1:flen)*stim_dt-stim_dt/2,1:24,rel_tkern);
xlim(xr)
xlabel('Time (s)');
ylabel('Probe');
title('Temporal kernel (relative)');

ename = sprintf('/home/james/Desktop/lab_meeting_figs/%s_LP_tkern_compare.png',Expt_name);
fig_width = 5; rel_height = 0.5;
figufy(h);
exportfig(h,ename,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
% [~,gsac_mua_loc] = min(mua_gsac,[],2);
% [~,msac_mua_loc] = min(mua_msac,[],2);
% [~,tkern_loc] = max(avg_tkern,[],2)
% 
%  [a,b] = corr(gsac_mua_loc,tkern_loc,'type','spearman')
%  [a,b] = corr(msac_mua_loc,tkern_loc,'type','spearman')
%  [a,b] = corr(msac_mua_loc,gsac_mua_loc,'type','spearman')