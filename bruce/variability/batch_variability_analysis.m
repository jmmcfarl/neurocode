clear all
addpath('~/James_scripts/bruce/variability/');

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
% Expt_list = {'M266','M270','M275','M277','M281','M287','M294'};
% Expt_list = {'M289'};
global Expt_name
for ee = 1:length(Expt_list)
    %     Expt_num = str2num(Expt_list{ee}(2:end));
    Expt_name = Expt_list{ee};
%     variability_anal_v2;
    variability_anal_v3;
    clearvars -except ee Expt_list Expt_name
end

%%
clear all
% Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
% Expt_list = {'M266','M270','M275','M277','M281','M287','M294'};
all_expt_num = [];
for ee = 1:length(Expt_list)
    Expt_num = str2num(Expt_list{ee}(2:end));
    Expt_name = sprintf('M%d',Expt_num);
    anal_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
    anal_name = 'variability_anal5';
    cd(anal_dir)
    load(anal_name)
    all_var_data(ee) = var_data;
    all_expt_num = cat(1,all_expt_num,Expt_num*ones(length(var_data.cor_mods),1));
end

%%
R2_imp = cell2mat(arrayfun(@(x) x.R2_imp,all_var_data,'uniformoutput',0));
cor_r2 = cell2mat(arrayfun(@(x) x.cor_r2,all_var_data,'uniformoutput',0));
mod_EM_var_frac = cell2mat(arrayfun(@(x) x.mod_EM_var_frac,all_var_data,'uniformoutput',0));
mod_tot_var = cell2mat(arrayfun(@(x) x.mod_tot_var,all_var_data,'uniformoutput',0));

shift_em_var_frac = cell2mat(arrayfun(@(x) x.shift_em_varfrac,all_var_data,'uniformoutput',0));
shift_em_var = cell2mat(arrayfun(@(x) x.shift_em_var,all_var_data,'uniformoutput',0));

% spline_EM_var_frac = cell2mat(arrayfun(@(x) x.spline_EM_var_frac',all_var_data,'uniformoutput',0));
% spline_em_var = cell2mat(arrayfun(@(x) x.spline_em_var',all_var_data,'uniformoutput',0));
% spline_psth_var = cell2mat(arrayfun(@(x) x.spline_psth_var',all_var_data,'uniformoutput',0));
% spline_tot_var = cell2mat(arrayfun(@(x) x.spline_tot_var,all_var_data,'uniformoutput',0));

psth_var = cell2mat(arrayfun(@(x) x.psth_var,all_var_data,'uniformoutput',0));
psth_var_cor = cell2mat(arrayfun(@(x) x.psth_var_cor,all_var_data,'uniformoutput',0));
tot_var = cell2mat(arrayfun(@(x) x.tot_var,all_var_data,'uniformoutput',0));
probe_nums = cell2mat(arrayfun(@(x) x.probe_nums,all_var_data,'uniformoutput',0));
tot_spikes = cell2mat(arrayfun(@(x) x.tot_spikes,all_var_data(1:end)','uniformoutput',0)');
avg_rates = cell2mat(arrayfun(@(x) x.avg_rates,all_var_data(1:end)','uniformoutput',0)');

% spline_r2 = cell2mat(arrayfun(@(x) x.spline_r2,all_var_data,'uniformoutput',0));
psth_r2 = cell2mat(arrayfun(@(x) x.psth_r2,all_var_data,'uniformoutput',0));
cormod_r2 = cell2mat(arrayfun(@(x) x.cor_r2,all_var_data,'uniformoutput',0));

beta_cc_psth = arrayfun(@(x) x.beta_cc_psth,all_var_data(1:end)','uniformoutput',0);
beta_cc_psth2 = arrayfun(@(x) x.beta_cc_psth2,all_var_data(1:end)','uniformoutput',0);
beta_cc_eye = arrayfun(@(x) x.beta_cc_eye,all_var_data(1:end)','uniformoutput',0);
beta_cc_psth_sur = arrayfun(@(x) x.beta_cc_psth_sur,all_var_data(1:end)','uniformoutput',0);
beta_cc_psth_sur2 = arrayfun(@(x) x.beta_cc_psth_sur2,all_var_data(1:end)','uniformoutput',0);
beta_cc_eye_sur = arrayfun(@(x) x.beta_cc_eye_sur,all_var_data(1:end)','uniformoutput',0);
% all_st_em_var = arrayfun(@(x) x.st_mod_emvfrac,all_var_data(2:end)','uniformoutput',0);

% all_np = cell2mat(np_em_var_frac');
% all_cormod = cell2mat(cormod_EM_vfrac');
% all_psth_var = cell2mat(psth_var');
% all_st_em_var = cell2mat(all_st_em_var');

all_beta_cc_psth = cell2mat(beta_cc_psth');
all_beta_cc_psth2 = cell2mat(beta_cc_psth2');
all_beta_cc_eye = cell2mat(beta_cc_eye');
all_beta_cc_psth_sur = cell2mat(beta_cc_psth_sur');
all_beta_cc_psth_sur2 = cell2mat(beta_cc_psth_sur2');
all_beta_cc_eye_sur = cell2mat(beta_cc_eye_sur');

all_avg_cc_obs = cell2mat(arrayfun(@(x) x.avg_cc_obs,all_var_data','uniformoutput',0)');
all_avg_cc_psth = cell2mat(arrayfun(@(x) x.avg_cc_psth,all_var_data','uniformoutput',0)');
all_avg_cc_psth2 = cell2mat(arrayfun(@(x) x.avg_cc_psth2,all_var_data','uniformoutput',0)');
all_avg_cc_eye = cell2mat(arrayfun(@(x) x.avg_cc_eye,all_var_data','uniformoutput',0)');

all_avg_cc_obs_sur = cell2mat(arrayfun(@(x) x.avg_cc_obs_sur,all_var_data','uniformoutput',0)');
all_avg_cc_psth_sur = cell2mat(arrayfun(@(x) x.avg_cc_psth_sur,all_var_data','uniformoutput',0)');
all_avg_cc_psth_sur2 = cell2mat(arrayfun(@(x) x.avg_cc_psth_sur2,all_var_data','uniformoutput',0)');
all_avg_cc_eye_sur = cell2mat(arrayfun(@(x) x.avg_cc_eye_sur,all_var_data','uniformoutput',0)');

% all_tot_spks = cell2mat(tot_spikes');
% all_mean_rate = cell2mat(avg_rates');
% 
% all_expt_nums = [];
% for ii = 1:length(all_var_data)
%     all_expt_nums = [all_expt_nums; ones(length(all_var_data(ii).cor_mods),1)*ii];
% end

% uset = find(tot_spikes > 1e3 & psth_var >= 0.02);
uset = find(tot_spikes > 1e3);
exclude_expts = [266 270 281];
uset = uset(~ismember(all_expt_num(uset),exclude_expts));

% uset = uset(tot_var(uset)./avg_rates(uset) > 0.01 & avg_rates(uset) > 0.01);
% uset = uset(tot_spikes(uset) > 1e3 & psth_var(uset) >= 0.005);
uset = uset(tot_spikes(uset) > 1e3 & mod_tot_var(uset)./avg_rates(uset) >= 0.05);

su_set = find(probe_nums > 24);
% su_set(ismember(all_expt_num(su_set),exclude_expts)) = [];
su_set(~ismember(su_set,uset)) = [];


psth_noise_var = tot_var - psth_var_cor;
em_noise_var = tot_var - shift_em_var;
psth_FF = psth_noise_var./avg_rates;
em_FF = em_noise_var./avg_rates;

mod_psth_frac = 1-mod_EM_var_frac;
psth_pow = psth_var./tot_var;
psth_pow2 = psth_var_cor./tot_var;
psth_pow_cor = psth_var./tot_var./mod_psth_frac;
psth_pow_cor2 = psth_var_cor./tot_var./mod_psth_frac;
psth_pow_shift = shift_em_var./tot_var;
%%
out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'Vfrac_est_mod_vs_shift.pdf';

h = figure();
plot(1-mod_EM_var_frac(uset),1-shift_em_var_frac(uset),'k.','markersize',10);
line([0 1],[0 1],'color','k')
xlabel('Model-based alpha');
ylabel('Shift-predictor alpha');
ylim([0.1 1])
xlim([0.1 1])
fig_width = 3;
rel_height = 0.8;

figufy(h);
exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'stimlocked_pow_scatter.pdf';

h = figure();
plot(psth_pow(uset),psth_pow_shift(uset),'.','markersize',8)
hold on
plot(psth_pow(uset),psth_pow_cor(uset),'r.','markersize',8)
line([0 1],[0 1],'color','k')
xlim([0 0.7]);
ylim([0 0.7])
% axis square

xlabel('Stim-locked power');
ylabel('Corrected stim-locked power');

fig_width = 3;
rel_height = 0.8;

figufy(h);
exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%

out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'avg_corr_psth_scatter.pdf';

h = figure();
subplot(2,1,1);
plot(all_avg_cc_obs(uset),all_avg_cc_psth(uset),'r.');
hold on
plot(all_avg_cc_obs(uset),all_avg_cc_psth2(uset),'k.');
xlabel('Avg. observed correlation');
ylabel('Avg. PSTH-corrected correlation');
xl = xlim();
line(xl,[0 0],'color','k','linestyle','--');
line([0 0.3],[0 0.3],'color','k');

subplot(2,1,2);
plot(all_avg_cc_obs_sur(uset),all_avg_cc_psth_sur(uset),'r.');
hold on
plot(all_avg_cc_obs_sur(uset),all_avg_cc_psth_sur2(uset),'k.');
xlabel('Avg. observed correlation');
ylabel('Avg. PSTH-corrected correlation');
xl = xlim();
line(xl,[0 0],'color','k','linestyle','--');
line([0 0.3],[0 0.3],'color','k');

fig_width = 3;
rel_height = 1.6;

% figufy(h);
% exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'stim_power_scatter2.pdf';

h = figure();
% plot(psth_pow,psth_pow_cor,'k.')
plot(psth_pow2,psth_pow_cor2,'k.')
hold on
% plot(psth_pow(uset),psth_pow_np(uset),'r.');
line([0 0.7],[0 0.7],'color','k')
xlabel('Stim power')
ylabel('Corrected stim power');
xlim([0 0.7]);ylim([0 0.7])
fig_width = 4;
rel_height = 0.8;

% figufy(h);
% exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
for ii = 1:length(all_var_data)
    avg_cov_mat(ii,:) = squeeze(nanmean(nanmean(all_var_data(ii).Cmat_obs)));
    avg_cov_mat_psth(ii,:) = squeeze(nanmean(nanmean(all_var_data(ii).Cmat_psth_res)));
    avg_cov_mat_psth2(ii,:) = squeeze(nanmean(nanmean(all_var_data(ii).Cmat_psth_res2)));
    avg_cov_mat_eye(ii,:) = squeeze(nanmean(nanmean(all_var_data(ii).Cmat_eye_res)));
end

%%
out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'avg_xcorrs.pdf';


lags = all_var_data(1).cc_lags;
% ue = [1 5];
ue = [2 3 4 6 7 8];

h=figure; hold on
h1=shadedErrorBar(lags,nanmean(avg_cov_mat(ue,:)),nanstd(avg_cov_mat(ue,:))./sqrt(length(ue)));
h2=shadedErrorBar(lags,nanmean(avg_cov_mat_psth(ue,:)),nanstd(avg_cov_mat_psth(ue,:))./sqrt(length(ue)),{'color','r'});
h3=shadedErrorBar(lags,nanmean(avg_cov_mat_eye(ue,:)),nanstd(avg_cov_mat_eye(ue,:))./sqrt(length(ue)),{'color','b'});
h4=shadedErrorBar(lags,nanmean(avg_cov_mat_psth2(ue,:)),nanstd(avg_cov_mat_psth2(ue,:))./sqrt(length(ue)),{'color','m'});
xlabel('Time lag (s)');
ylabel('Cross-correlation');
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Observed','PSTH-cor','Eye-shift','scaled PSTH-cor'});

fig_width = 4;
rel_height = 0.8;

figufy(h);
exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'allexpt_avg_xcorrs.pdf';

h=figure; 
uset = [2 3 4 6 7 8];
for uu = 1:length(uset)
    ue = uset(uu);
% ue = 8;

acov_mat = squeeze(nanmean(all_var_data(ue).Cmat_obs));
acov_mat_psth = squeeze(nanmean(all_var_data(ue).Cmat_psth_res));
acov_mat_psthe = squeeze(nanmean(all_var_data(ue).Cmat_psth));
acov_mat_psth2 = squeeze(nanmean(all_var_data(ue).Cmat_psth_res2));
acov_mat_eye = squeeze(nanmean(all_var_data(ue).Cmat_eye_res));

subplot(2,3,uu); hold on
h1=shadedErrorBar(lags,nanmean(acov_mat),nanstd(acov_mat)./sqrt(size(acov_mat,1)));
h2=shadedErrorBar(lags,nanmean(acov_mat_psth),nanstd(acov_mat_psth)./sqrt(size(acov_mat,1)),{'color','r'});
h3=shadedErrorBar(lags,nanmean(acov_mat_eye),nanstd(acov_mat_eye)./sqrt(size(acov_mat,1)),{'color','b'});
h4=shadedErrorBar(lags,nanmean(acov_mat_psth2),nanstd(acov_mat_psth2)./sqrt(size(acov_mat,1)),{'color','m'});
% h5=shadedErrorBar(lags,nanmean(acov_mat_psthe),nanstd(acov_mat_psthe)./sqrt(size(acov_mat,1)),{'color','g'});
xlabel('Time lag (s)');
ylabel('Cross-correlation');
% legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Observed','PSTH-cor','Eye-shift','scaled PSTH-cor'});

% pause
% close all
end


fig_width = 8;
rel_height = 0.5;

figufy(h);
exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
close all

ue = 8;
acov_mat = squeeze(all_var_data(ue).Cmat_obs);
acov_mat_psth = squeeze(all_var_data(ue).Cmat_psth_res);
acov_mat_psthe = squeeze(all_var_data(ue).Cmat_psth);
acov_mat_psth2 = squeeze(all_var_data(ue).Cmat_psth_res2);
acov_mat_eye = squeeze(all_var_data(ue).Cmat_eye_res);
n_units = size(acov_mat,1);

for ii = 1:n_units;
    
    subplot(2,2,1)
    imagescnan(lags,1:n_units,squeeze(acov_mat(ii,:,:)));
    ca = caxis(); cam = 0.75*max(abs(ca)); caxis([-cam cam]);
    yl = ylim();
    title('Measured');
    xlabel('Time lag (s)');
    ylabel('Unit');
    
    subplot(2,2,2)
    imagescnan(lags,1:n_units,squeeze(acov_mat_psth(ii,:,:)));
%     ca = caxis(); cam = 0.75*max(abs(ca)); 
caxis([-cam cam]);
    yl = ylim();
    title('Measured');
    xlabel('Time lag (s)');
    ylabel('Unit');
    
    subplot(2,2,3)
    imagescnan(lags,1:n_units,squeeze(acov_mat_eye(ii,:,:)));
%     ca = caxis(); cam = 0.75*max(abs(ca)); 
    caxis([-cam cam]);
    yl = ylim();
    title('Measured');
    xlabel('Time lag (s)');
    ylabel('Unit');

        subplot(2,2,4)
    imagescnan(lags,1:n_units,squeeze(acov_mat_psth2(ii,:,:)));
%     ca = caxis(); cam = 0.75*max(abs(ca)); 
    caxis([-cam cam]);
    yl = ylim();
    title('Measured');
    xlabel('Time lag (s)');
    ylabel('Unit');
pause
    clf
end

%%
cd /home/james/Analysis/bruce/ET_final
load('unit_summary_data2');

expt_nums = [all_unit_data(:).expt_nums];
cset = find(expt_nums < 200);
expt_nums(cset) = [];

R2_imps = all_unit_data.after_R2./all_unit_data.before_R2;
R2_imps(cset) = [];

num_blocks = all_unit_data.num_blocks;
num_blocks(cset) = [];

rf_width = all_unit_data.rf_width;
rf_width(cset) = [];

prm = all_unit_data.prm;
prm(cset) = [];

rf_eccs = all_unit_data.rf_eccs;
rf_eccs(cset) = [];

uprobe_nums = all_unit_data.probe_nums;
uprobe_nums(cset) = [];

su_set = find(probe_nums > 24);
su_set(num_blocks(su_set) <= 3) = [];

%%
out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'RFecc_v_PSTHpow_SU.pdf';

% uset = find(num_blocks >= 3);
uset = su_set;

h = figure();
% plot(rf_eccs,psth_var./tot_var,'k.');
plot(rf_eccs(uset),psth_var_cor(uset)./tot_var(uset),'k.');
hold on
% r = robustfit(rf_eccs,psth_var./tot_var);
[r,stats] = robustfit(rf_eccs(uset),psth_var_cor(uset)./tot_var(uset));
ee = linspace(0,4.5,100);
plot(ee,r(1)+r(2)*ee,'r');
xlim([0.5 4.5])
ylim([0 0.35])
xlabel('RF Eccentricity (deg)');
ylabel('Stimulus locked power')

fig_width = 4;
rel_height = 0.8;

% figufy(h);
% exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% 

%%
out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'PSTHpow_dist_SU.pdf';

% uset = find(num_blocks >= 3);
uset = su_set;

xax = linspace(0,0.35,21);
n = histc(psth_var_cor(uset)./tot_var(uset),xax);
n = n(1:end-1)/sum(n);

h = figure();
hold on
bar(x,n);
xlim([0 0.35])
xlabel('RF Eccentricity (deg)');
ylabel('Stimulus locked power')

fig_width = 4;
rel_height = 0.8;

figufy(h);
exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
