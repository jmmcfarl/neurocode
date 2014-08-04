clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/variability/');
% Expt_list = {'M266','M270'};
% Expt_list = {'M275','M277'};
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
% Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','G085','G086','G087','G088','G089','G093','G095'};
% Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
% Expt_list = {'G095'};
% Expt_list = {'M275'};
global Expt_name bar_ori use_LOOXV
bar_ori = 0;
use_LOOXV = 1;
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    bar_ori = 0;
    model_based_EMvar;
    clearvars -except ee Expt_list Expt_name bar_ori use_LOOXV
end

%%
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','G085','G086','G087','G088','G089','G091','G093'};
% Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    out_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
    Expt_num = str2num(Expt_name(2:end));
%     if Expt_num < 200
    out_name = 'mod_EM_var';
%     else
%         out_name = 'mod_EM_var2';
%     end
    cd(out_dir)
    load(out_name)
    all_mod_em(ee) = mod_em_est;
end

%%
em_var = cell2mat(arrayfun(@(x) x.em_var,all_mod_em','uniformoutput',0)');
psth_var = cell2mat(arrayfun(@(x) x.psth_var,all_mod_em','uniformoutput',0)');
tot_var = cell2mat(arrayfun(@(x) x.tot_var,all_mod_em','uniformoutput',0)');
em_var_frac = em_var./tot_var;

all_em_corrs = [];
all_psth_corrs = [];
all_tot_corrs = [];

all_em_vars = [];
all_psth_vars = [];
all_tot_vars = [];
all_sc_fac = [];
all_expt_nums = [];
s_expt_nums = [];
for ii = 1:length(all_mod_em)
    all_em_corrs = cat(1,all_em_corrs,all_mod_em(ii).em_corr(:));
    all_psth_corrs = cat(1,all_psth_corrs,all_mod_em(ii).psth_corr(:));
    all_tot_corrs = cat(1,all_tot_corrs,all_mod_em(ii).tot_corr(:));
    
    cur_expt_num = str2num(Expt_list{ii}(2:end));
    all_expt_nums = cat(1,all_expt_nums,ones(size(all_mod_em(ii).tot_corr(:)))*cur_expt_num);

    cur_expt_num = str2num(Expt_list{ii}(2:end));
    s_expt_nums = cat(1,s_expt_nums,ones(size(all_mod_em(ii).em_var(:)))*cur_expt_num);

    cur_ee = all_mod_em(ii).em_var;
    cur_pp = all_mod_em(ii).psth_var;
    cur_tt = all_mod_em(ii).tot_var;

    cur_sc_fac = sqrt(cur_tt'*cur_tt);
    all_sc_fac = cat(1,all_sc_fac,cur_sc_fac(:));
    
    cur_em_vfrac = all_mod_em(ii).em_var./all_mod_em(ii).tot_var;
%     cur_ee = (cur_em_vfrac'*cur_em_vfrac);
% %     cur_ee = cur_ee'*cur_ee;
%     cur_pp = all_mod_em(ii).psth_var'*all_mod_em(ii).psth_var;
%     cur_tt = all_mod_em(ii).tot_var'*all_mod_em(ii).tot_var;
    cur_ee = bsxfun(@plus,cur_ee,cur_ee');
    cur_pp = bsxfun(@plus,cur_pp,cur_pp');
    cur_tt = bsxfun(@plus,cur_tt,cur_tt');
    all_em_vars = cat(1,all_em_vars,cur_ee(:));
    all_psth_vars = cat(1,all_psth_vars,cur_pp(:));
    all_tot_vars = cat(1,all_tot_vars,cur_tt(:));
end

%%
cd /home/james/Analysis/bruce/ET_final
load('unit_summary_data2');

cset = find(all_unit_data.expt_nums == 86);
cset = cset(1:94);
cset = [cset find(all_unit_data.expt_nums == 95)];

R2_imps = all_unit_data.after_R2./all_unit_data.before_R2;
R2_imps(cset) = [];

after_R2 = all_unit_data.after_R2;
after_R2(cset) = [];

expt_nums = all_unit_data.expt_nums;
expt_nums(cset) = [];

num_blocks = all_unit_data.num_blocks;
num_blocks(cset) = [];

rf_width = all_unit_data.rf_width;
rf_width(cset) = [];

probe_nums = all_unit_data.probe_nums;
probe_nums(cset) = [];

rf_eccs = all_unit_data.rf_eccs;
rf_eccs(cset) = [];

prm = all_unit_data.prm;
prm(cset) = [];

rf_SF = all_unit_data.rf_SF;
rf_SF(cset) = [];

mean_rate = all_unit_data.mean_rate;
mean_rate(cset) = [];

n_spikes = all_unit_data.n_spikes;
n_spikes(cset) = [];

% used_units = find(num_blocks > 3);
% used_units = find(num_blocks > 3 & after_R2 > 0.025);
% used_units = find(num_blocks > 3 & tot_var > 0.005);
used_units = find(num_blocks > 3 & tot_var./mean_rate > 0.01 & mean_rate > 0.01);

%%
su_inds = [];
su_pair_inds = [];
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','G085','G086','G087','G088','G089','G091','G093'};
for ii = 1:length(Expt_list)
    cur_EN = str2num(Expt_list{ii}(2:end));
    curset = find(expt_nums == cur_EN);
    if cur_EN > 200
        cur_su_start = find(diff(probe_nums(curset)) < 0,1)+1;
        su_inds = [su_inds curset(cur_su_start:end)];
    else
%         cur_su_start = find(diff(probe_nums(curset)) < 0,1)+1;
        su_inds = [su_inds curset];
    end
end
su_inds = su_inds(ismember(su_inds,used_units));
%%
out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'RFwidth_vs_PSTHfrac.pdf';

h = figure();
plot(rf_width(used_units),1-em_var_frac(used_units),'k.','markersize',10)
hold on
plot(rf_width(su_inds),1-em_var_frac(su_inds),'r.','markersize',10)
set(gca,'xscale','log')
xlim([0.08 0.8])
xlabel('RF width (deg)');
ylabel('alpha');

fig_width = 4;
rel_height = 0.8;
figufy(h);
exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'RFSF_vs_PSTHfrac.pdf';

h = figure();
plot(rf_SF(used_units),1-em_var_frac(used_units),'k.','markersize',10)
hold on
plot(rf_SF(su_inds),1-em_var_frac(su_inds),'r.','markersize',10)
set(gca,'xscale','log')
xlim([0.6 8])
xlabel('RF SF (cyc/deg)');
ylabel('alpha');

fig_width = 4;
rel_height = 0.8;
figufy(h);
exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'R2imp_vs_PSTHfrac2.pdf';

h = figure();
plot(1-em_var_frac(used_units),R2_imps(used_units),'k.','markersize',10)
hold on
plot(1-em_var_frac(su_inds),R2_imps(su_inds),'r.','markersize',10)
set(gca,'yscale','log')
ylabel('R2 fold-improvement');
xlabel('alpha');
ylim([1 10])

fig_width = 4;
rel_height = 0.8;
figufy(h);
exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'R2imp_vs_PSTHfrac.pdf';

h = figure();
plot(rf_eccs(used_units),1-em_var_frac(used_units),'k.')
% set(gca,'xscale','log')
xlabel('R2 fold-improvement');
ylabel('alpha');
% xlim([1 10])

fig_width = 4;
rel_height = 0.8;
% figufy(h);
% exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
% parafov_expts = [270 277 281 287 289 294];
% parafov_expts = [266 270 275 277 281 287 289 294];
parafov_expts = [270];
parafov_neurons = find(ismember(all_expt_nums,parafov_expts));
fov_neurons = find(~ismember(all_expt_nums,parafov_expts));

h = figure();
% plot(all_psth_corrs,all_em_corrs,'.');
hold on
plot(all_psth_corrs(fov_neurons),all_em_corrs(fov_neurons),'.');
plot(all_psth_corrs(parafov_neurons),all_em_corrs(parafov_neurons),'r.');
line([-0.5 1],[-0.5 1],'color','k','linewidth',2);
xlim([-0.25 0.75]);
ylim([-0.25 0.75]);
xlabel('PSTH correlation');
ylabel('Stim-conditional correlation');

figufy(h);

out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'modpred_corr_scatter.pdf';

fig_width = 4;
rel_height = 0.8;

% figufy(h);
% exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
parafov_expts = [270 277 281 287 289 294];
% parafov_expts = [266 270 275 277 281 287 289 294];
parafov_neurons = find(ismember(all_expt_nums,parafov_expts));
fov_neurons = find(~ismember(all_expt_nums,parafov_expts));

use_pairs = find(all_tot_corrs < 0.9);

fov_neurons = find(all_expt_nums < 200);
h = figure();
plot(all_tot_corrs(use_pairs),all_em_corrs(use_pairs),'.');
hold on
r = robustfit(all_tot_corrs(use_pairs),all_em_corrs(use_pairs));
xx = linspace(-0.25,1,100);
plot(xx,r(2)*xx+r(1),'r')

% plot(all_tot_corrs(fov_neurons),all_em_corrs(fov_neurons),'.');
% plot(all_tot_corrs(parafov_neurons),all_em_corrs(parafov_neurons),'r.');
line([-0.5 1],[-0.5 1],'color','k','linewidth',2);
xlim([-0.25 1]);
ylim([-0.25 1]);
xlabel('Stimulus correlation');
ylabel('EM-induced noise correlation');
axis square
out_dir = '/home/james/Analysis/bruce/variability/';
fname = 'SUpair_EMinduced_correlation.pdf';

fig_width = 4;
rel_height = 1;

% figufy(h);
% exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
