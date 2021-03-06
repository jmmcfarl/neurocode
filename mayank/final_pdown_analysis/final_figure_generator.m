clear all
addpath(genpath('~/James_scripts/figuremaker/'))
fig_dir = '/Users/james/Analysis/Mayank/final_pdown_analysis/figures/';

min_rec_dur = 500; %minimum total duration of recording (in sec)
max_med_ctx_down = 2.5; %maximum median duration of ctx LFP down states (in sec) 

%% load EC data and select usable recs
load ~/Analysis/Mayank/final_pdown_analysis/compiled_data.mat

data_ids = [data(:).id];
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur);
used_dirs(ismember(data_ids(used_dirs),no_cell) | ismember(data_ids(used_dirs),unclear_uds)) = [];
used_dirs_nl2 = used_dirs([data(used_dirs).layer] ~= 2);

data = data(used_dirs);
data_ids = data_ids(used_dirs);

load('~/Analysis/Mayank/final_pdown_analysis/fin_pdown_core_analysis_fin.mat');
if length(data) ~= length(core_data)
    error('data misalignment');
end

%only take recordings where the anesthesia is not too deep, as judged by
%the median duration of ctx down durations
med_ctx_downdur = arrayfun(@(x) nanmedian(x.lfp_down_durs),core_data);
tot_uds_dur = [core_data(:).tot_uds_dur];
use_ec_recs = find(med_ctx_downdur < max_med_ctx_down & tot_uds_dur > min_rec_dur);
use_ec_recs_nl2 = use_ec_recs([data(use_ec_recs).layer] ~= 2);

data = data(use_ec_recs);
data_ids = data_ids(use_ec_recs);
core_data = core_data(use_ec_recs);

%% classify EC cell types
mec = find(strcmp({data.loc},'MEC')); %include the original data set plus new ones
lec = find(strcmp({data.loc},'LEC')); %these are already only the clear L3 neurons used in the original pers ups paper
ctype = {data(:).ctype};
ctype{strcmp(ctype,'pyramidal')} = 'pyr';
layer = [data(:).layer];
l3mec = mec(layer(mec) == 3);
l3lec = lec(layer(lec) == 3);
l2mec = mec(layer(mec) == 2);
l2lec = lec(layer(lec) == 2);

l3mec_pyr = l3mec(strcmp(ctype(l3mec),'pyr'));
l3mec_nonpyr = l3mec(strcmp(ctype(l3mec),'nonpyr'));
l3lec_pyr = l3lec(strcmp(ctype(l3lec),'pyr'));
l3lec_nonpyr = l3lec(strcmp(ctype(l3lec),'nonpyr'));


fract_rt2_ups = [core_data(:).fract_rt2_ups];
fract_rt2_downs = [core_data(:).fract_rt2_downs];
fprintf('Pers Ups:\n');
fprintf('L3MEC: %.2f  %.2f - %.2f\n',100*prctile(fract_rt2_ups(l3mec),[50 25 75]));
fprintf('L2MEC: %.2f  %.2f - %.2f\n',100*prctile(fract_rt2_ups(l2mec),[50 25 75]));
fprintf('L3LEC: %.2f  %.2f - %.2f\n',100*prctile(fract_rt2_ups(l3lec),[50 25 75]));
fprintf('L2LEC: %.2f  %.2f - %.2f\n',100*prctile(fract_rt2_ups(l2lec),[50 25 75]));

fprintf('Pers Downs:\n');
fprintf('L3MEC: %.2f  %.2f - %.2f\n',100*prctile(fract_rt2_downs(l3mec),[50 25 75]));
fprintf('L2MEC: %.2f  %.2f - %.2f\n',100*prctile(fract_rt2_downs(l2mec),[50 25 75]));
fprintf('L3LEC: %.2f  %.2f - %.2f\n',100*prctile(fract_rt2_downs(l3lec),[50 25 75]));
fprintf('L2LEC: %.2f  %.2f - %.2f\n',100*prctile(fract_rt2_downs(l2lec),[50 25 75]));

[a,b] = corr(fract_rt2_ups(l3mec)',fract_rt2_downs(l3mec)','type','spearman');
fprintf('L3MEC pup pdown corr: %.4f  p: %.4f\n',a,b);

p = ranksum(fract_rt2_downs(l3mec),fract_rt2_downs(l3lec));
fprintf('L3MEC vs L3LEC pdowns, p: %.4f\n',p)
% g1 = layer;
% g2 = {data.loc};
% % [p,tbl,stats] = anovan(fract_rt2_ups,{g1,g2})
% [p,tbl,stats] = anovan(fract_rt2_ups,{g1,g2})

%% plot pers ups vs pers downs probabilities
mSize = 17;
fig_width = 3.5;
rel_heigh = 1;

fract_rt2_ups = [core_data(:).fract_rt2_ups];
fract_rt2_downs = [core_data(:).fract_rt2_downs];

% eps = 1e-3;
% fract_rt2_ups(fract_rt2_ups < eps) = eps;
% fract_rt2_downs(fract_rt2_downs < eps) = eps;

h = figure;hold on
plot(fract_rt2_ups(l3mec),fract_rt2_downs(l3mec),'r.','markersize',mSize,'linewidth',1.5);
plot(fract_rt2_ups(l3lec),fract_rt2_downs(l3lec),'b.','markersize',mSize,'linewidth',1.5);
legend('L3MEC','L3LEC');
xlabel('Prob. persistent Up');
ylabel('Prob. persistent Down');
xlim([0 0.6]);
ylim([0 0.6]);
line([0 0.6],[0 0.6],'color','k');
set(gca,'xtick',[0:0.1:0.6],'ytick',[0:0.1:0.6])
% xlim([eps 0.6]);
% ylim([eps 0.6]);
% set(gca,'xscale','log','yscale','log');
% 
% figufy(h);
% fname = [fig_dir 'pers_ups_vs_pers_downs.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);


%% load in cortical MP data and get regions
ctx_data = load('~/Analysis/Mayank/final_pdown_analysis/compiled_corticalMP_data.mat');
ctx_data = ctx_data.data;
ctx_used_dirs = find([ctx_data(:).ep] > min_rec_dur & [ctx_data(:).dp] > min_rec_dur); %dont use recs that are too short
ctx_data = ctx_data(ctx_used_dirs);

ctx_core_data = load('~/Analysis/Mayank/final_pdown_analysis/fin_pdown_core_analysis_corticalMP_fin.mat');
ctx_core_data = ctx_core_data.core_data;

ctx_cell_type = arrayfun(@(x) x.cell_type,ctx_data,'uniformoutput',0);
ctx_interneurons = find(strcmp(ctx_cell_type,'interneuron')); 

%get rid of ctx interneurons
ctx_data(ctx_interneurons) = [];
ctx_core_data(ctx_interneurons) = [];

%exclude sessions where the anesthesia was too deep
med_ctx_downdur = arrayfun(@(x) nanmedian(x.lfp_down_durs),ctx_core_data);
tot_uds_dur = [ctx_core_data(:).tot_uds_dur];
to_exclude = find(med_ctx_downdur > max_med_ctx_down | tot_uds_dur < min_rec_dur);
ctx_data(to_exclude) = [];
ctx_core_data(to_exclude) = [];

frontal = find(strcmp({ctx_data(:).region},'frontal'));
prefrontal = find(strcmp({ctx_data(:).region},'prefrontal'));
parietal = find(strcmp({ctx_data(:).region},'parietal'));
barrel = find(strcmp({ctx_data(:).region},'barrel'));

%% plot boxplot comparison of EC and ctx persistence probs
all_fract_ups = [ctx_core_data(:).fract_rt2_ups core_data(:).fract_rt2_ups];
all_fract_downs = [ctx_core_data(:).fract_rt2_downs core_data(:).fract_rt2_downs];

ctx_group = nan(length(ctx_core_data),1);
ctx_group(prefrontal) = 3;
ctx_group(frontal) = 4;
ctx_group(parietal) = 5;
EC_group = nan(length(core_data),1);
EC_group(l3mec) = 1;
EC_group(l3lec) = 2;
% EC_group(l2mec) = 6;
% EC_group(l2lec) = 7;

all_group = [ctx_group; EC_group];
uset = find(~isnan(all_group));
used_ctx = find(ismember(all_group,[3,4,5]));

fprintf('Pers Ups:\n');
fprintf('prefrontal: n=%d, %.2f  %.2f - %.2f\n',sum(all_group == 3),100*prctile(all_fract_ups(all_group == 3),[50 25 75]));
fprintf('frontal: n=%d, %.2f  %.2f - %.2f\n',sum(all_group == 4),100*prctile(all_fract_ups(all_group == 4),[50 25 75]));
fprintf('parietal: n=%d, %.2f  %.2f - %.2f\n',sum(all_group == 5),100*prctile(all_fract_ups(all_group == 5),[50 25 75]));

fprintf('Pers Downs:\n');
fprintf('prefrontal: n=%d, %.2f  %.2f - %.2f\n',sum(all_group == 3),100*prctile(all_fract_downs(all_group == 3),[50 25 75]));
fprintf('frontal: n=%d, %.2f  %.2f - %.2f\n',sum(all_group == 4),100*prctile(all_fract_downs(all_group == 4),[50 25 75]));
fprintf('parietal: n=%d, %.2f  %.2f - %.2f\n',sum(all_group == 5),100*prctile(all_fract_downs(all_group == 5),[50 25 75]));

close all
% names = {'L3-MEC','L3-LEC','Prefrontal','Frontal','Parietal','L2-MEC','L2-LEC'};
names = {'L3-MEC','L3-LEC','Prefrontal','Frontal','Parietal'};
h1 = figure;
% boxplot(all_fract_ups(uset)',names(all_group(uset)),'plotstyle','compact');
boxplot(all_fract_ups(uset)',names(all_group(uset)));
ylabel('Fraction persistent UP states');
ylim([0 0.6])
figufy(h1);
set(findobj(gca,'Type','text'),'FontSize',10)

h2=figure;
% boxplot(all_fract_downs(uset)',names(all_group(uset)),'plotstyle','compact');
boxplot(all_fract_downs(uset)',names(all_group(uset)));
ylabel('Fraction persistent DOWN states');
ylim([0 0.6])
figufy(h2);
set(findobj(gca,'Type','text'),'FontSize',10)

% fig_width = 6;
% rel_heigh = 0.8;
% fname = [fig_dir 'allcell_pUp_dist_all.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);
% 
% fname = [fig_dir 'allcell_pDown_dist_all.pdf'];
% exportfig(h2,fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);
% close(h2);

%% scatter plots of pers prob vs transition latency
mSize = 2;

Fsd = 2016/8; 
fract_rt2_ups = [core_data(:).fract_rt2_ups];
fract_rt2_downs = [core_data(:).fract_rt2_downs];

%get median and mean up-transition lags in sec
median_uplag = cellfun(@nanmedian,{core_data(:).mp_uplags})/Fsd;
mean_uplag = cellfun(@nanmean,{core_data(:).mp_uplags})/Fsd;

%same for down-transition lags
median_downlag = cellfun(@nanmedian,{core_data(:).mp_downlags})/Fsd;
mean_downlag = cellfun(@nanmean,{core_data(:).mp_downlags})/Fsd;

%calculate spearman rhos
[mec_med_pup_corr,mec_med_pup_p] = corr(median_downlag(l3mec)',fract_rt2_ups(l3mec)','type','spearman');
[mec_mean_pup_corr,mec_mean_pup_p] = corr(mean_downlag(l3mec)',fract_rt2_ups(l3mec)','type','spearman');
[lec_med_pup_corr,lec_med_pup_p] = corr(median_downlag(l3lec)',fract_rt2_ups(l3lec)','type','spearman');
[lec_mean_pup_corr,lec_mean_pup_p] = corr(mean_downlag(l3lec)',fract_rt2_ups(l3lec)','type','spearman');
[mec_med_pdown_corr,mec_med_pdown_p] = corr(median_uplag(l3mec)',fract_rt2_downs(l3mec)','type','spearman');
[mec_mean_pdown_corr,mec_mean_pdown_p] = corr(mean_uplag(l3mec)',fract_rt2_downs(l3mec)','type','spearman');
[lec_med_pdown_corr,lec_med_pdown_p] = corr(median_uplag(l3lec)',fract_rt2_downs(l3lec)','type','spearman');
[lec_mean_pdown_corr,lec_mean_pdown_p] = corr(mean_uplag(l3lec)',fract_rt2_downs(l3lec)','type','spearman');

%threshold pers probs before plotting on log scale
eps = 1e-3;
fract_rt2_ups(fract_rt2_ups < eps) = eps;
fract_rt2_downs(fract_rt2_downs < eps) = eps;

h1 = figure;hold on
plot(median_downlag(l3mec),fract_rt2_ups(l3mec),'ro','markersize',mSize);
plot(median_downlag(l3lec),fract_rt2_ups(l3lec),'bo','markersize',mSize);
% plot(mean_downlag(l3mec),fract_rt2_ups(l3mec),'ro','markersize',mSize);
% plot(mean_downlag(l3lec),fract_rt2_ups(l3lec),'bo','markersize',mSize);
set(gca,'yscale','log');
xlabel('Down-transition lag (s)');
ylabel('Prob. pers. Up');
xlim([-0.2 1.2]);

h2 = figure;hold on
plot(median_uplag(l3mec),fract_rt2_downs(l3mec),'ro','markersize',mSize);
plot(median_uplag(l3lec),fract_rt2_downs(l3lec),'bo','markersize',mSize);
% plot(mean_uplag(l3mec),fract_rt2_downs(l3mec),'ro','markersize',mSize);
% plot(mean_uplag(l3lec),fract_rt2_downs(l3lec),'bo','markersize',mSize);
set(gca,'yscale','log');
xlabel('Up-transition lag (s)');
ylabel('Prob. pers. Down');

% fig_width = 3.5;
% rel_heigh = 0.8;
% 
% figufy(h1);
% fname = [fig_dir 'persUp_downLag.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);
% 
% figufy(h2);
% fname = [fig_dir 'persDown_upLag.pdf'];
% exportfig(h2,fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);
% close(h2);

%% plot ctx up-trig avg ctx LFP spectra
load ~/Analysis/Mayank/final_pdown_analysis/final_trig_spec_fin_noconst_newpeaks_wch5_fin.mat
% use_ec_nl2_recs = use_ec_recs([data(:).layer] ~= 2);
spec_data = spec_data(use_ec_recs_nl2);
% if length(spec_data) ~= length(data([data(:).layer] ~= 2))
%     error('data misalignment');
% end

xr = [-0.5 1.5]; %time axis range

%for pup and pdown-conditional analysis, only use those cells that had at
%least a min number of pers states
min_states = 5;
n_pers_downs = cellfun(@(x) length(x),{core_data(:).rt2_downs});
n_pers_ups = cellfun(@(x) length(x),{core_data(:).rt2_ups});
l3mec_pups = l3mec(n_pers_ups(l3mec) >= min_states);
l3mec_pdowns = l3mec(n_pers_downs(l3mec) >= min_states);

interp_f = linspace(min(wfreqs),max(wfreqs),100)';

%get skipped and non-skipped trig avg specgrams of cortical LFP
non_sk_utrig_lf8_spec = reshape([spec_data(l3mec_pdowns).non_sk_utrig_lfp_spec],[length(lags) length(wfreqs) length(l3mec_pdowns)]);
sk_utrig_lf8_spec = reshape([spec_data(l3mec_pdowns).sk_utrig_lfp_spec],[length(lags) length(wfreqs) length(l3mec_pdowns)]);

%avg across neurons, and interpolate onto uniform freq axis
nonsk_specgram = squeeze(nanmean(non_sk_utrig_lf8_spec,3))';
interp_nonsk_specgram = interp1(wfreqs,nonsk_specgram,interp_f);
sk_specgram = squeeze(nanmean(sk_utrig_lf8_spec,3))';
interp_sk_specgram = interp1(wfreqs,sk_specgram,interp_f);

h1 = figure;
subplot(2,1,1);
% pcolor(lags/Fsd,wfreqs,squeeze(nanmean(non_sk_utrig_lf8_spec,3))');shading flat
imagesc(lags/Fsd,wfreqs,interp_nonsk_specgram);set(gca,'ydir','normal');
xlim(xr);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Non-skipped');
colorbar
ca = caxis(); %get color axis defined by this first plot and apply it to the next
subplot(2,1,2);
% pcolor(lags/Fsd,wfreqs,squeeze(nanmean(sk_utrig_lf8_spec,3))');shading flat
imagesc(lags/Fsd,wfreqs,interp_sk_specgram);set(gca,'ydir','normal');
xlim(xr);
caxis(ca);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Skipped');
colorbar

fig_width = 4.36;
rel_height = 1.5;
figufy(h1);
fname1 = [fig_dir 'utrig_ctx_spec_5s'];
exportfig(h1,fname1,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

%% plot ctx up-trig avg hpc LFP spectra
load ~/Analysis/Mayank/final_pdown_analysis/mua_classification_fin.mat
peak_hpcmua_loc = peak_hpcmua_loc(used_dirs_nl2);
peak_hpcmua_rate = peak_hpcmua_rate(used_dirs_nl2);
peak_hpcmua_loc = peak_hpcmua_loc(use_ec_recs_nl2);
peak_hpcmua_rate = peak_hpcmua_rate(use_ec_recs_nl2);
% if length(peak_hpcmua_loc) ~= length(data([data(:).layer] ~= 2))
%     error('data alignment mismatch');
% end

%get set of recs that have enough pers downs and usable hpc MUA channels
good_hpc_lfp = ~isnan(peak_hpcmua_loc);
l3mec_pdown_hpclfp = l3mec_pdowns(good_hpc_lfp(l3mec_pdowns));

%get skipped and non-skipped trig avg specgrams of hpc LFP
non_sk_utrig_hpc_spec = reshape([spec_data(l3mec_pdown_hpclfp).non_sk_utrig_hpc_spec],[length(lags) length(wfreqs) length(l3mec_pdown_hpclfp)]);
sk_utrig_hpc_spec = reshape([spec_data(l3mec_pdown_hpclfp).sk_utrig_hpc_spec],[length(lags) length(wfreqs) length(l3mec_pdown_hpclfp)]);

%avg across neurons, and interpolate onto uniform freq axis
nonsk_specgram = squeeze(nanmean(non_sk_utrig_hpc_spec,3))';
interp_nonsk_specgram = interp1(wfreqs,nonsk_specgram,interp_f);
sk_specgram = squeeze(nanmean(sk_utrig_hpc_spec,3))';
interp_sk_specgram = interp1(wfreqs,sk_specgram,interp_f);

h1 = figure;
subplot(2,1,1);
% pcolor(lags/Fsd,wfreqs,squeeze(nanmean(non_sk_utrig_hpc_spec,3))');shading flat
imagesc(lags/Fsd,wfreqs,interp_nonsk_specgram);set(gca,'ydir','normal');
xlim(xr);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Non-skipped');
colorbar
ca = caxis(); %get color axis defined by this first plot and apply it to the next
subplot(2,1,2);
% pcolor(lags/Fsd,wfreqs,squeeze(nanmean(sk_utrig_hpc_spec,3))');shading flat
imagesc(lags/Fsd,wfreqs,interp_sk_specgram);set(gca,'ydir','normal');
xlim(xr);
caxis(ca);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Skipped');
colorbar

% fig_width = 4.36;
% rel_height = 1.5;
% figufy(h1);
% fname1 = [fig_dir 'utrig_hpc_spec_5s'];
% exportfig(h1,fname1,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);

%% plot ctx DOWN-trig avg ctx LFP spectra

xr = [-0.5 1.5]; %time axis range

%get skipped and non-skipped trig avg specgrams of cortical LFP
non_sk_dtrig_lf8_spec = reshape([spec_data(l3mec_pups).non_sk_dtrig_lfp_spec],[length(lags) length(wfreqs) length(l3mec_pups)]);
sk_dtrig_lf8_spec = reshape([spec_data(l3mec_pups).sk_dtrig_lfp_spec],[length(lags) length(wfreqs) length(l3mec_pups)]);

%avg across neurons, and interpolate onto uniform freq axis
nonsk_specgram = squeeze(nanmean(non_sk_dtrig_lf8_spec,3))';
interp_nonsk_specgram = interp1(wfreqs,nonsk_specgram,interp_f);
sk_specgram = squeeze(nanmean(sk_dtrig_lf8_spec,3))';
interp_sk_specgram = interp1(wfreqs,sk_specgram,interp_f);

h1 = figure;
subplot(2,1,1);
% pcolor(lags/Fsd,wfreqs,squeeze(nanmean(non_sk_dtrig_lf8_spec,3))');shading flat
imagesc(lags/Fsd,wfreqs,interp_nonsk_specgram);set(gca,'ydir','normal');
xlim(xr);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Non-skipped');
colorbar
ca = caxis(); %get color axis defined by this first plot and apply it to the next
subplot(2,1,2);
% pcolor(lags/Fsd,wfreqs,squeeze(nanmean(sk_dtrig_lf8_spec,3))');shading flat
imagesc(lags/Fsd,wfreqs,interp_sk_specgram);set(gca,'ydir','normal');
xlim(xr);
caxis(ca);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Skipped');
colorbar

% fig_width = 4.36;
% rel_height = 1.5;
% figufy(h1);
% fname1 = [fig_dir 'dtrig_ctx_spec_5s'];
% exportfig(h1,fname1,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);

%% plot ctx DOWN-trig avg hpc LFP spectra
load ~/Analysis/Mayank/final_pdown_analysis/mua_classification_fin.mat
peak_hpcmua_loc = peak_hpcmua_loc(used_dirs_nl2);
peak_hpcmua_rate = peak_hpcmua_rate(used_dirs_nl2);
peak_hpcmua_loc = peak_hpcmua_loc(use_ec_recs_nl2);
peak_hpcmua_rate = peak_hpcmua_rate(use_ec_recs_nl2);
% if length(peak_hpcmua_loc) ~= length(data)
%     error('data alignment mismatch');
% end

%get set of recs that have enough pers downs and usable hpc MUA channels
good_hpc_lfp = ~isnan(peak_hpcmua_loc);
l3mec_pup_hpclfp = l3mec_pups(good_hpc_lfp(l3mec_pups));

%get skipped and non-skipped trig avg specgrams of hpc LFP
non_sk_dtrig_hpc_spec = reshape([spec_data(l3mec_pup_hpclfp).non_sk_dtrig_hpc_spec],[length(lags) length(wfreqs) length(l3mec_pup_hpclfp)]);
sk_dtrig_hpc_spec = reshape([spec_data(l3mec_pup_hpclfp).sk_dtrig_hpc_spec],[length(lags) length(wfreqs) length(l3mec_pup_hpclfp)]);

%avg across neurons, and interpolate onto uniform freq axis
nonsk_specgram = squeeze(nanmean(non_sk_dtrig_hpc_spec,3))';
interp_nonsk_specgram = interp1(wfreqs,nonsk_specgram,interp_f);
sk_specgram = squeeze(nanmean(sk_dtrig_hpc_spec,3))';
interp_sk_specgram = interp1(wfreqs,sk_specgram,interp_f);

h1 = figure;
subplot(2,1,1);
% pcolor(lags/Fsd,wfreqs,squeeze(nanmean(non_sk_dtrig_hpc_spec,3))');shading flat
imagesc(lags/Fsd,wfreqs,interp_nonsk_specgram);set(gca,'ydir','normal');
xlim(xr);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Non-skipped');
colorbar
ca = caxis(); %get color axis defined by this first plot and apply it to the next
subplot(2,1,2);
% pcolor(lags/Fsd,wfreqs,squeeze(nanmean(sk_dtrig_hpc_spec,3))');shading flat
imagesc(lags/Fsd,wfreqs,interp_sk_specgram);set(gca,'ydir','normal');
xlim(xr);
caxis(ca);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Skipped');
colorbar

% fig_width = 4.36;
% rel_height = 1.5;
% 
% figufy(h1);
% fname1 = [fig_dir 'dtrig_hpc_spec_5s'];
% exportfig(h1,fname1,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);

%% plot prob of persistence vs HF power of ctx LFP
addpath('~/James_scripts/mayank/final_pdown_analysis/');
load ~/Analysis/Mayank/final_pdown_analysis/final_cortical_state_data_fin_nobuff_wc5_fin.mat
cfun_data = cfun_data(use_ec_recs_nl2);
% if length(cfun_data) ~= length(data)
%     error('data alignment mismatch');
% end
%for pup and pdown-conditional analysis, only use those cells that had at
%least a min number of pers states
min_states = 5;
n_pers_downs = cellfun(@(x) length(x),{core_data(:).rt2_downs});
n_pers_ups = cellfun(@(x) length(x),{core_data(:).rt2_ups});
l3mec_pups = l3mec(n_pers_ups(l3mec) >= min_states);
l3mec_pdowns = l3mec(n_pers_downs(l3mec) >= min_states);

n_bins = 15;

used_data = l3mec_pdowns; %used recs
X_fac = {cfun_data(:).ctx_up_amps_hf};
Y_fac = {cfun_data(:).ctx_up_isskipped};
avg_upamp_dist = nan(length(used_data),length(state_amp_xx));
[avg_fun,std_fun,N_fun,avg_xvals,sem_fun] = deal(nan(length(used_data),n_bins));
[XY_corr_pdown,XY_pval_pdown] = deal(nan(length(used_data),1));
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr_pdown(ii),XY_pval_pdown(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_upamp_dist(ii,:) = cfun_data(used_data(ii)).ctx_up_amp_hf_dist; %get avg density estimate
end
n_nonnan = sum(~isnan(avg_fun)); %number of cells contributing to avg

%plot pers down prob vs HF-power relationship
h = figure(); 
subplot(2,1,2);hold on
xax = nanmean(avg_xvals); %across-cell avg of the avg x-val within each bin
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan),{'color','g'});
% plot(xax,nanmean(avg_fun),'ko')
subplot(2,1,1); hold on
n_nonnan = sum(~isnan(avg_upamp_dist)); %number of cells contributing to density estimates
shadedErrorBar(state_amp_xx,nanmean(avg_upamp_dist),nanstd(avg_upamp_dist)./sqrt(n_nonnan),{'color','g'});

%non repeat for pups
used_data = l3mec_pups;

X_fac = {cfun_data(:).ctx_down_amps_hf};
Y_fac = {cfun_data(:).ctx_down_isskipped};
avg_downamp_dist = nan(length(used_data),length(state_amp_xx));
[avg_fun,std_fun,N_fun,avg_xvals,sem_fun] = deal(nan(length(used_data),n_bins));
[XY_corr_pup,XY_pval_pup] = deal(nan(length(used_data),1));
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr_pup(ii),XY_pval_pup(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_downamp_dist(ii,:) = cfun_data(used_data(ii)).ctx_down_amp_hf_dist; %get avg density estimate
end
n_nonnan = sum(~isnan(avg_fun)); %number of cells contributing to avg

% figure(h);
subplot(2,1,2);
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan),{'color','m'});
% plot(xax,nanmean(avg_fun),'ko')
xlabel('Ctx state amplitude (z)');
ylabel('Prob. skipped');

% figure(h2); hold on
subplot(2,1,1);
n_nonnan = sum(~isnan(avg_downamp_dist));
shadedErrorBar(state_amp_xx,nanmean(avg_downamp_dist),nanstd(avg_downamp_dist)./sqrt(n_nonnan),{'color','m'});
xlabel('Ctx state amplitude (z)');
ylabel('Probability');

subplot(2,1,2);
xlim([-1.5 2]);
ylim([0 0.6])
subplot(2,1,1);
xlim([-1.5 2]);

% fig_width = 4; rel_height = 1.6;
% figufy(h);
% fname = [fig_dir 'prob_skip_vs_hfpow_fin_5s.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%% plot prob of persistence vs time since previous MP state transition
load ~/Analysis/Mayank/final_pdown_analysis/final_cortical_state_data_fin_nobuff_wc5_fin.mat
cfun_data = cfun_data(use_ec_recs_nl2);
% if length(cfun_data) ~= length(data)
%     error('data alignment mismatch');
% end
%for pup and pdown-conditional analysis, only use those cells that had at
%least a min number of pers states
min_states = 5;
n_pers_downs = cellfun(@(x) length(x),{core_data(:).rt2_downs});
n_pers_ups = cellfun(@(x) length(x),{core_data(:).rt2_ups});
l3mec_pups = l3mec(n_pers_ups(l3mec) >= min_states);
l3mec_pdowns = l3mec(n_pers_downs(l3mec) >= min_states);

n_bins = 15;
state_dur_xx = linspace(0.1,5,100);

used_data = l3mec_pdowns; %used recs
X_fac = {cfun_data(:).TSL_mp_down};
Y_fac = {cfun_data(:).ctx_up_isskipped};
avg_tsl_dist = nan(length(used_data),length(state_dur_xx));
[avg_fun,std_fun,N_fun,avg_xvals,sem_fun] = deal(nan(length(used_data),n_bins));
[XY_corr_pdown,XY_pval_pdown] = deal(nan(length(used_data),1));
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr_pdown(ii),XY_pval_pdown(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_tsl_dist(ii,:) = ksdensity(X_fac{used_data(ii)},state_dur_xx,'support','positive');
end
n_nonnan = sum(~isnan(avg_fun)); %number of cells contributing to avg

%plot pers down prob vs time-since last MP down-trans relationship
h = figure(); 
subplot(2,1,2);hold on
xax = nanmean(avg_xvals); %across-cell avg of the avg x-val within each bin
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan),{'color','g'});
% plot(xax,nanmean(avg_fun),'ko')
subplot(2,1,1); hold on
n_nonnan = sum(~isnan(avg_upamp_dist)); %number of cells contributing to density estimates
shadedErrorBar(state_dur_xx,nanmean(avg_tsl_dist),nanstd(avg_tsl_dist)./sqrt(n_nonnan),{'color','g'});

%non repeat for pups
used_data = l3mec_pups;

X_fac = {cfun_data(:).TSL_mp_up};
Y_fac = {cfun_data(:).ctx_down_isskipped};
avg_tsl_dist = nan(length(used_data),length(state_dur_xx));
[avg_fun,std_fun,N_fun,avg_xvals,sem_fun] = deal(nan(length(used_data),n_bins));
[XY_corr_pup,XY_pval_pup] = deal(nan(length(used_data),1));
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr_pup(ii),XY_pval_pup(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_tsl_dist(ii,:) = ksdensity(X_fac{used_data(ii)},state_dur_xx,'support','positive');
end
n_nonnan = sum(~isnan(avg_fun)); %number of cells contributing to avg

% figure(h);
subplot(2,1,2);
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan),{'color','m'});
% plot(xax,nanmean(avg_fun),'ko')
xlabel('Ctx state amplitude (z)');
ylabel('Prob. skipped');

% figure(h2); hold on
subplot(2,1,1);
n_nonnan = sum(~isnan(avg_downamp_dist));
shadedErrorBar(state_dur_xx,nanmean(avg_tsl_dist),nanstd(avg_tsl_dist)./sqrt(n_nonnan),{'color','m'});
xlabel('Ctx state amplitude (z)');
ylabel('Probability');

subplot(2,1,2);
xlim([0.1 5]); 
ylim([0 0.3])
set(gca,'xscale','log');
subplot(2,1,1);
xlim([0.1 5]); 
set(gca,'xscale','log');

% fig_width = 4; rel_height = 1.6;
% figufy(h);
% fname = [fig_dir 'prob_skip_vs_tsl_fin_5s.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%% boxplot of effect of EC and ctx states on hpc MUA rates
load ~/Analysis/Mayank/final_pdown_analysis/final_trig_avg_data4_nocon_nobuff_newpeaks_wch5_dcmp_noORDrej_fin.mat
trig_data = trig_data(use_ec_recs_nl2);
% if length(trig_data) ~= length(data)
%     error('data alignment mismatch');
% end

%load MUA classification data
load ~/Analysis/Mayank/final_pdown_analysis/mua_classification_fin.mat
peak_hpcmua_loc = peak_hpcmua_loc(used_dirs_nl2);
peak_hpcmua_rate = peak_hpcmua_rate(used_dirs_nl2);
peak_hpcmua_loc = peak_hpcmua_loc(use_ec_recs_nl2);
peak_hpcmua_rate = peak_hpcmua_rate(use_ec_recs_nl2);
% if length(peak_hpcmua_loc) ~= length(data)
%     error('data alignment mismatch');
% end

%get set of recs that have enough pers downs and usable hpc MUA channels
good_hpc_mua = ~isnan(peak_hpcmua_loc);

%for pup and pdown-conditional analysis, only use those cells that had at
%least a min number of pers states
min_states = 5;
n_pers_downs = cellfun(@(x) length(x),{core_data(:).rt2_downs});
n_pers_ups = cellfun(@(x) length(x),{core_data(:).rt2_ups});
l3mec_pups = l3mec(n_pers_ups(l3mec) >= min_states);
l3mec_pdowns = l3mec(n_pers_downs(l3mec) >= min_states);

%skip non-skip conditional avgs
sk_down_ctx_mua = [trig_data(:).avg_sk_down_ctx_mua];
nonsk_down_ctx_mua = [trig_data(:).avg_nsk_down_ctx_mua];
sk_down_hpc_mua = [trig_data(:).avg_sk_down_hpc_mua];
nonsk_down_hpc_mua = [trig_data(:).avg_nsk_down_hpc_mua];
sk_up_ctx_mua = [trig_data(:).avg_sk_up_ctx_mua];
nonsk_up_ctx_mua = [trig_data(:).avg_nsk_up_ctx_mua];
sk_up_hpc_mua = [trig_data(:).avg_sk_up_hpc_mua];
nonsk_up_hpc_mua = [trig_data(:).avg_nsk_up_hpc_mua];

% joint state conditional avgs with robust criteria
avg_cd_md_hpc_mua = [trig_data(:).ravg_cd_md_hpc_mua];
avg_cu_mu_hpc_mua = [trig_data(:).ravg_cu_mu_hpc_mua];
avg_cd_mu_hpc_mua = [trig_data(:).ravg_cd_mu_hpc_mua];
avg_cu_md_hpc_mua = [trig_data(:).ravg_cu_md_hpc_mua];

avg_pdown_effect = sk_up_hpc_mua-nonsk_up_hpc_mua; %'effect of a pers down'
avg_pup_effect = sk_down_hpc_mua-nonsk_down_hpc_mua; %'effect of a pers up'
avg_ctxup_effect = sk_up_hpc_mua - avg_cd_md_hpc_mua; %effect of a skipped cortical up
avg_ctxdown_effect = sk_down_hpc_mua - avg_cu_mu_hpc_mua; %effect of a skipped cortical down

l3mec_hpcmua_pdowns = l3mec_pdowns(good_hpc_mua(l3mec_pdowns)); %good hpc MUA and sufficient pdowns
l3mec_hpcmua_pups = l3mec_pups(good_hpc_mua(l3mec_pups)); %good hpc MUA and sufficient pups

all_data = [avg_pdown_effect(l3mec_hpcmua_pdowns)'; avg_pup_effect(l3mec_hpcmua_pups)'; ...
    avg_ctxup_effect(l3mec_hpcmua_pdowns)'; avg_ctxdown_effect(l3mec_hpcmua_pups)'];
groups = [ones(length(l3mec_hpcmua_pdowns),1); 2*ones(length(l3mec_hpcmua_pups),1); ...
    3*ones(length(l3mec_hpcmua_pdowns),1); 4*ones(length(l3mec_hpcmua_pups),1)];

h = figure();
boxplot(all_data,groups);
ylabel('Hpc MUA (z)');
xl = xlim();
% ylim([-1 2])
line(xl,[0 0],'color','k','linestyle','--');

% fig_width = 3.5; rel_height = 0.8;
% figufy(h);
% fname = [fig_dir 'hpcmua_boxplot_5s_np_wch5.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

%% up-trig avgs
load ~/Analysis/Mayank/final_pdown_analysis/final_trig_avg_data4_nocon_nobuff_newpeaks_wch5_dcmp_noORDrej_fin.mat
use_ec_nl2_recs = use_ec_recs([data(:).layer] ~= 2);
trig_data = trig_data(use_ec_nl2_recs);
if length(trig_data) ~= length(data([data(:).layer] ~= 2))
    error('data alignment mismatch');
end

%get L3MEC recs with usable MUA signals
min_rate = 1;
avg_hpc_rate = [trig_data(:).avg_hpc_mua_rate];
avg_ctx_rate = [trig_data(:).avg_ctx_mua_rate];
l3mec_ctxmua = l3mec(avg_ctx_rate(l3mec) >= min_rate);
l3mec_hpcmua = l3mec(avg_hpc_rate(l3mec) >= min_rate);

%for pup and pdown-conditional analysis, only use those cells that had at
%least a min number of pers states
min_states = 5;
n_pers_downs = cellfun(@(x) length(x),{core_data(:).rt2_downs});
n_pers_ups = cellfun(@(x) length(x),{core_data(:).rt2_ups});
l3mec_pups = l3mec(n_pers_ups(l3mec) >= min_states);
l3mec_pdowns = l3mec(n_pers_downs(l3mec) >= min_states);

%set of L3MEC neurons usable for both MUA and pers analysis
l3mec_ctxmua_pdowns = intersect(l3mec_pdowns,l3mec_ctxmua);
l3mec_hpcmua_pdowns = intersect(l3mec_pdowns,l3mec_hpcmua);
l3mec_ctxmua_pups = intersect(l3mec_pups,l3mec_ctxmua);
l3mec_hpcmua_pups = intersect(l3mec_pups,l3mec_hpcmua);

xr = [-1.5 3];

%LF amp signals
l3mec_non_sk_utrig_ctx = reshape([trig_data(l3mec_pdowns).non_sk_utrig_ctx_lf],[length(lags) length(l3mec_pdowns)]);
l3mec_sk_utrig_ctx = reshape([trig_data(l3mec_pdowns).sk_utrig_ctx_lf],[length(lags) length(l3mec_pdowns)]);
l3mec_non_sk_utrig_hpc = reshape([trig_data(l3mec_hpcmua_pdowns).non_sk_utrig_hpc_lf],[length(lags) length(l3mec_hpcmua_pdowns)]);
l3mec_sk_utrig_hpc = reshape([trig_data(l3mec_hpcmua_pdowns).sk_utrig_hpc_lf],[length(lags) length(l3mec_hpcmua_pdowns)]);
l3mec_non_sk_utrig_mp = reshape([trig_data(l3mec_pdowns).non_sk_utrig_mp_dc],[length(lags) length(l3mec_pdowns)]);
l3mec_sk_utrig_mp = reshape([trig_data(l3mec_pdowns).sk_utrig_mp_dc],[length(lags) length(l3mec_pdowns)]);

%hf power signals
l3mec_non_sk_utrig_ctxhf = reshape([trig_data(l3mec_pdowns).non_sk_utrig_ctx_hf],[length(lags) length(l3mec_pdowns)]);
l3mec_sk_utrig_ctxhf = reshape([trig_data(l3mec_pdowns).sk_utrig_ctx_hf],[length(lags) length(l3mec_pdowns)]);
l3mec_non_sk_utrig_hpchf = reshape([trig_data(l3mec_hpcmua_pdowns).non_sk_utrig_hpc_hf],[length(lags) length(l3mec_hpcmua_pdowns)]);
l3mec_sk_utrig_hpchf = reshape([trig_data(l3mec_hpcmua_pdowns).sk_utrig_hpc_hf],[length(lags) length(l3mec_hpcmua_pdowns)]);
l3mec_non_sk_utrig_mphf = reshape([trig_data(l3mec_pdowns).non_sk_utrig_mp_hf],[length(lags) length(l3mec_pdowns)]);
l3mec_sk_utrig_mphf = reshape([trig_data(l3mec_pdowns).sk_utrig_mp_hf],[length(lags) length(l3mec_pdowns)]);

%ctx and hpc MUA
l3mec_non_sk_utrig_ctxmua = reshape([trig_data(l3mec_ctxmua_pdowns).non_sk_utrig_ctx_mua],[length(lags) length(l3mec_ctxmua_pdowns)]);
l3mec_sk_utrig_ctxmua = reshape([trig_data(l3mec_ctxmua_pdowns).sk_utrig_ctx_mua],[length(lags) length(l3mec_ctxmua_pdowns)]);
l3mec_non_sk_utrig_hpcmua = reshape([trig_data(l3mec_hpcmua_pdowns).non_sk_utrig_hpc_mua],[length(lags) length(l3mec_hpcmua_pdowns)]);
l3mec_sk_utrig_hpcmua = reshape([trig_data(l3mec_hpcmua_pdowns).sk_utrig_hpc_mua],[length(lags) length(l3mec_hpcmua_pdowns)]);

%MP spike rates
l3mec_meanrate = [trig_data(l3mec_pdowns).avg_mp_rate];
l3mec_stdrate = [trig_data(l3mec_pdowns).std_mp_rate];
l3mec_non_sk_utrig_mprate = reshape([trig_data(l3mec_pdowns).non_sk_utrig_mp_rate],[length(lags) length(l3mec_pdowns)]);
l3mec_sk_utrig_mprate = reshape([trig_data(l3mec_pdowns).sk_utrig_mp_rate],[length(lags) length(l3mec_pdowns)]);
l3mec_non_sk_utrig_mprate = bsxfun(@times,l3mec_non_sk_utrig_mprate,l3mec_stdrate); %this removes zscore normalization
l3mec_non_sk_utrig_mprate = bsxfun(@plus,l3mec_non_sk_utrig_mprate,l3mec_meanrate);
l3mec_sk_utrig_mprate = bsxfun(@times,l3mec_sk_utrig_mprate,l3mec_stdrate);
l3mec_sk_utrig_mprate = bsxfun(@plus,l3mec_sk_utrig_mprate,l3mec_meanrate);

f1 = figure(); 
hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_ctx,2),nanstd(l3mec_non_sk_utrig_ctx,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_ctx),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_ctx,2),nanstd(l3mec_sk_utrig_ctx,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_ctx),2)),{'color','r'});
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
% title('Ctx LF');
ylabel('Amplitude (z)');
xlabel('Time (s)');

f2 = figure();
subplot(3,2,1);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_ctxhf,2),nanstd(l3mec_non_sk_utrig_ctxhf,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_ctxhf),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_ctxhf,2),nanstd(l3mec_sk_utrig_ctxhf,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_ctxhf),2)),{'color','r'});
xlim(xr);
ylim([-1 1.5])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
% title('Ctx HF');
ylabel('HF power (z)');

subplot(3,2,2);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_ctxmua,2),nanstd(l3mec_non_sk_utrig_ctxmua,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_ctxmua),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_ctxmua,2),nanstd(l3mec_sk_utrig_ctxmua,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_ctxmua),2)),{'color','r'});
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
% title('Ctx MUA');
ylabel('MUA rate (z)');

%MP
subplot(3,2,3);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_mp,2),nanstd(l3mec_non_sk_utrig_mp,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_mp),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_mp,2),nanstd(l3mec_sk_utrig_mp,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_mp),2)),{'color','r'});
xlim(xr);
ylim([-5 20])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
% title('MP amp');
ylabel('Relative amplitude (mV)');
xlabel('Time (s)');

% subplot(3,3,8);hold on
% shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_mphf,2),nanstd(l3mec_non_sk_utrig_mphf,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_mphf),2)),{'color','k'});
% shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_mphf,2),nanstd(l3mec_sk_utrig_mphf,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_mphf),2)),{'color','r'});
% xlim(xr);
% % ylim([-0.75 0.75])
% yl = ylim();
% line([0 0],yl,'color','k','linestyle','--');
% line(xr,[0 0],'color','k','linestyle','--');
% title('MP HF-power');
% ylabel('HF power (z)');
% xlabel('Time (s)');

subplot(3,2,4);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_mprate,2),nanstd(l3mec_non_sk_utrig_mprate,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_mprate),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_mprate,2),nanstd(l3mec_sk_utrig_mprate,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_mprate),2)),{'color','r'});
xlim(xr);
% ylim([-0.75 0.75])
ylim([0 7])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
% title('MP spiking');
ylabel('Spike rate (Hz)');
xlabel('Time (s)');

% subplot(3,3,4);hold on
% shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_hpc,2),nanstd(l3mec_non_sk_utrig_hpc,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_hpc),2)),{'color','k'});
% shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_hpc,2),nanstd(l3mec_sk_utrig_hpc,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_hpc),2)),{'color','r'});
% xlim(xr);
% yl = ylim();
% line([0 0],yl,'color','k','linestyle','--');
% line(xr,[0 0],'color','k','linestyle','--');
% title('Hpc');
% ylabel('Amplitude (z)');

subplot(3,2,5);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_hpchf,2),nanstd(l3mec_non_sk_utrig_hpchf,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_hpchf),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_hpchf,2),nanstd(l3mec_sk_utrig_hpchf,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_hpchf),2)),{'color','r'});
xlim(xr);
ylim([-0.75 1.25])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
% title('Hpc HF-power');
ylabel('HF-power (z)');

subplot(3,2,6);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_hpcmua,2),nanstd(l3mec_non_sk_utrig_hpcmua,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_hpcmua),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_hpcmua,2),nanstd(l3mec_sk_utrig_hpcmua,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_hpcmua),2)),{'color','r'});
xlim(xr);
ylim([-0.5 0.4])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
% title('Hpc MUA');
ylabel('MUA rate (z)');

%
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'utrig_ctx_lf.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% fig_width = 8; rel_height = 1.4;
% figufy(f2);
% fname = [fig_dir 'utrig_signal_fin_5s_np_wch5_v2.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);


%% down-trig avgs
load ~/Analysis/Mayank/final_pdown_analysis/final_trig_avg_data4_nocon_nobuff_newpeaks_wch5_dcmp_noORDrej_fin.mat
trig_data = trig_data(use_ec_recs_nl2);
% if length(trig_data) ~= length(data([data(:).layer] ~= 2))
%     error('data alignment mismatch');
% end

%get L3MEC recs with usable MUA signals
min_rate = 1;
avg_hpc_rate = [trig_data(:).avg_hpc_mua_rate];
avg_ctx_rate = [trig_data(:).avg_ctx_mua_rate];
l3mec_ctxmua = l3mec(avg_ctx_rate(l3mec) >= min_rate);
l3mec_hpcmua = l3mec(avg_hpc_rate(l3mec) >= min_rate);

%for pup and pdown-conditional analysis, only use those cells that had at
%least a min number of pers states
min_states = 5;
n_pers_downs = cellfun(@(x) length(x),{core_data(:).rt2_downs});
n_pers_ups = cellfun(@(x) length(x),{core_data(:).rt2_ups});
l3mec_pups = l3mec(n_pers_ups(l3mec) >= min_states);
l3mec_pdowns = l3mec(n_pers_downs(l3mec) >= min_states);

%set of L3MEC neurons usable for both MUA and pers analysis
l3mec_ctxmua_pdowns = intersect(l3mec_pdowns,l3mec_ctxmua);
l3mec_hpcmua_pdowns = intersect(l3mec_pdowns,l3mec_hpcmua);
l3mec_ctxmua_pups = intersect(l3mec_pups,l3mec_ctxmua);
l3mec_hpcmua_pups = intersect(l3mec_pups,l3mec_hpcmua);

xr = [-1.5 3];

%LF amp signals
l3mec_non_sk_dtrig_ctx = reshape([trig_data(l3mec_pups).non_sk_dtrig_ctx_lf],[length(lags) length(l3mec_pups)]);
l3mec_sk_dtrig_ctx = reshape([trig_data(l3mec_pups).sk_dtrig_ctx_lf],[length(lags) length(l3mec_pups)]);
l3mec_non_sk_dtrig_hpc = reshape([trig_data(l3mec_hpcmua_pups).non_sk_dtrig_hpc_lf],[length(lags) length(l3mec_hpcmua_pups)]);
l3mec_sk_dtrig_hpc = reshape([trig_data(l3mec_hpcmua_pups).sk_dtrig_hpc_lf],[length(lags) length(l3mec_hpcmua_pups)]);
l3mec_non_sk_dtrig_mp = reshape([trig_data(l3mec_pups).non_sk_dtrig_mp_dc],[length(lags) length(l3mec_pups)]);
l3mec_sk_dtrig_mp = reshape([trig_data(l3mec_pups).sk_dtrig_mp_dc],[length(lags) length(l3mec_pups)]);

%hf power signals
l3mec_non_sk_dtrig_ctxhf = reshape([trig_data(l3mec_pups).non_sk_dtrig_ctx_hf],[length(lags) length(l3mec_pups)]);
l3mec_sk_dtrig_ctxhf = reshape([trig_data(l3mec_pups).sk_dtrig_ctx_hf],[length(lags) length(l3mec_pups)]);
l3mec_non_sk_dtrig_hpchf = reshape([trig_data(l3mec_hpcmua_pups).non_sk_dtrig_hpc_hf],[length(lags) length(l3mec_hpcmua_pups)]);
l3mec_sk_dtrig_hpchf = reshape([trig_data(l3mec_hpcmua_pups).sk_dtrig_hpc_hf],[length(lags) length(l3mec_hpcmua_pups)]);
l3mec_non_sk_dtrig_mphf = reshape([trig_data(l3mec_pups).non_sk_dtrig_mp_hf],[length(lags) length(l3mec_pups)]);
l3mec_sk_dtrig_mphf = reshape([trig_data(l3mec_pups).sk_dtrig_mp_hf],[length(lags) length(l3mec_pups)]);

%ctx and hpc MUA
l3mec_non_sk_dtrig_ctxmua = reshape([trig_data(l3mec_ctxmua_pups).non_sk_dtrig_ctx_mua],[length(lags) length(l3mec_ctxmua_pups)]);
l3mec_sk_dtrig_ctxmua = reshape([trig_data(l3mec_ctxmua_pups).sk_dtrig_ctx_mua],[length(lags) length(l3mec_ctxmua_pups)]);
l3mec_non_sk_dtrig_hpcmua = reshape([trig_data(l3mec_hpcmua_pups).non_sk_dtrig_hpc_mua],[length(lags) length(l3mec_hpcmua_pups)]);
l3mec_sk_dtrig_hpcmua = reshape([trig_data(l3mec_hpcmua_pups).sk_dtrig_hpc_mua],[length(lags) length(l3mec_hpcmua_pups)]);

%MP spike rates
l3mec_meanrate = [trig_data(l3mec_pups).avg_mp_rate];
l3mec_stdrate = [trig_data(l3mec_pups).std_mp_rate];
l3mec_non_sk_dtrig_mprate = reshape([trig_data(l3mec_pups).non_sk_dtrig_mp_rate],[length(lags) length(l3mec_pups)]);
l3mec_sk_dtrig_mprate = reshape([trig_data(l3mec_pups).sk_dtrig_mp_rate],[length(lags) length(l3mec_pups)]);
l3mec_non_sk_dtrig_mprate = bsxfun(@times,l3mec_non_sk_dtrig_mprate,l3mec_stdrate); %this removes zscore normalization
l3mec_non_sk_dtrig_mprate = bsxfun(@plus,l3mec_non_sk_dtrig_mprate,l3mec_meanrate);
l3mec_sk_dtrig_mprate = bsxfun(@times,l3mec_sk_dtrig_mprate,l3mec_stdrate);
l3mec_sk_dtrig_mprate = bsxfun(@plus,l3mec_sk_dtrig_mprate,l3mec_meanrate);

f1 = figure(); 
hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_ctx,2),nanstd(l3mec_non_sk_dtrig_ctx,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_ctx),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_ctx,2),nanstd(l3mec_sk_dtrig_ctx,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_ctx),2)),{'color','r'});
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Ctx LF');
ylabel('Amplitude (z)');

f2 = figure();
subplot(3,2,1);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_ctxhf,2),nanstd(l3mec_non_sk_dtrig_ctxhf,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_ctxhf),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_ctxhf,2),nanstd(l3mec_sk_dtrig_ctxhf,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_ctxhf),2)),{'color','r'});
xlim(xr);
ylim([-1 1.5])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
% title('Ctx HF');
ylabel('HF power (z)');

subplot(3,2,2);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_ctxmua,2),nanstd(l3mec_non_sk_dtrig_ctxmua,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_ctxmua),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_ctxmua,2),nanstd(l3mec_sk_dtrig_ctxmua,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_ctxmua),2)),{'color','r'});
xlim(xr);
ylim([-1 1.5])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
% title('Ctx MUA');
ylabel('MUA rate (z)');

subplot(3,2,3);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_mp,2),nanstd(l3mec_non_sk_dtrig_mp,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_mp),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_mp,2),nanstd(l3mec_sk_dtrig_mp,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_mp),2)),{'color','r'});
xlim(xr);
ylim([0 25])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
% line(xr,[0 0],'color','k','linestyle','--');
% title('MP amplitude');
ylabel('Rel. Amp. (mV)');
xlabel('Time (s)');

% subplot(3,3,8);hold on
% shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_mphf,2),nanstd(l3mec_non_sk_dtrig_mphf,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_mphf),2)),{'color','k'});
% shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_mphf,2),nanstd(l3mec_sk_dtrig_mphf,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_mphf),2)),{'color','r'});
% xlim(xr);
% % ylim([-0.75 0.75])
% yl = ylim();
% line([0 0],yl,'color','k','linestyle','--');
% line(xr,[0 0],'color','k','linestyle','--');
% title('MP HF-power');
% ylabel('HF power (z)');
% xlabel('Time (s)');

subplot(3,2,4);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_mprate,2),nanstd(l3mec_non_sk_dtrig_mprate,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_mprate),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_mprate,2),nanstd(l3mec_sk_dtrig_mprate,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_mprate),2)),{'color','r'});
xlim(xr);
ylim([0 7])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
% line(xr,[0 0],'color','k','linestyle','--');
% title('MP spiking');
ylabel('Spike rate (Hz)');
xlabel('Time (s)');

% subplot(3,3,4);hold on
% shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_hpc,2),nanstd(l3mec_non_sk_dtrig_hpc,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_hpc),2)),{'color','k'});
% shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_hpc,2),nanstd(l3mec_sk_dtrig_hpc,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_hpc),2)),{'color','r'});
% xlim(xr);
% yl = ylim();
% line([0 0],yl,'color','k','linestyle','--');
% line(xr,[0 0],'color','k','linestyle','--');
% title('Hpc LF');
% ylabel('Amplitude (z)');

subplot(3,2,5);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_hpchf,2),nanstd(l3mec_non_sk_dtrig_hpchf,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_hpchf),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_hpchf,2),nanstd(l3mec_sk_dtrig_hpchf,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_hpchf),2)),{'color','r'});
xlim(xr);
ylim([-0.75 1.25])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
% title('Hpc HF-power');
ylabel('HF-power (z)');

subplot(3,2,6);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_hpcmua,2),nanstd(l3mec_non_sk_dtrig_hpcmua,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_hpcmua),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_hpcmua,2),nanstd(l3mec_sk_dtrig_hpcmua,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_hpcmua),2)),{'color','r'});
xlim(xr);
% ylim([-0.5 0.5])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
% title('Hpc MUA');
ylabel('MUA rate (z)');


%
fig_width = 8; rel_height = 1.4;
figufy(f2);
fname = [fig_dir 'dtrig_spikerates_fin_5s_np_wch5_v2.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);


