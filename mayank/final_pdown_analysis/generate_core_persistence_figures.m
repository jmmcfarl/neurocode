clear all
load C:/WC_Germany/final_pdown_analysis/compiled_data.mat
addpath(genpath('C:/Code/figuremaker/'));
fig_dir = 'C:\WC_Germany\final_pdown_analysis\figures\';

min_rec_dur = 500; %in sec
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur);
data_ids = [data(:).id];

% used_dirs(ismember(data_ids(used_dirs),no_cell) | ismember(data_ids(used_dirs),unclear_uds)) = [];
used_dirs(ismember(data_ids(used_dirs),unclear_uds)) = [];

data = data(used_dirs);
data_ids = data_ids(used_dirs);

cd C:\WC_Germany\final_pdown_analysis\
load('fin_pdown_core_analysis.mat')


l3mec = find(strcmp({data.loc},'MEC'));
l3lec = find(strcmp({data.loc},'LEC'));
ctype = {data(:).ctype};
%%
% l3mec = l3mec(ismember(data_ids(l3mec),clear_l3pyr) & strcmp(ctype(l3mec),'pyr'));

l3mec(~ismember(data_ids(l3mec),clear_l3)) = [];
l3mec_possnonpyr = l3mec(~ismember(l3mec,clear_l3pyr));
l3mec_nonpyr = l3mec(strcmp(ctype(l3mec),'nonpyr'));
l3mec_pyr = l3mec(ismember(l3mec,clear_l3pyr));

l3lec_nonpyr = l3lec(strcmp(ctype(l3lec),'nonpyr'));
l3lec_pyr = l3lec(strcmp(ctype(l3lec),'pyr'));
%%
mSize = 3;
fig_width = 3.5;
rel_heigh = 1;

fract_rt2_ups = [core_data(:).fract_rt2_ups];
fract_rt2_downs = [core_data(:).fract_rt2_downs];
% fract_t2_ups = [core_data(:).fract_t2_ups];
% fract_t2_downs = [core_data(:).fract_t2_downs];
% fract_rt2_ups = sqrt(fract_rt2_ups);
% fract_rt2_downs = sqrt(fract_rt2_downs);

eps = 1e-3;
fract_rt2_ups(fract_rt2_ups < eps) = eps;
fract_rt2_downs(fract_rt2_downs < eps) = eps;

h = figure;hold on
plot(fract_rt2_ups(l3mec),fract_rt2_downs(l3mec),'ro','markersize',mSize,'linewidth',1.5);
plot(fract_rt2_ups(l3lec),fract_rt2_downs(l3lec),'bo','markersize',mSize,'linewidth',1.5);
% plot(fract_t2_ups(l3mec),fract_t2_downs(l3mec),'ro','markersize',mSize,'linewidth',1.5);
% plot(fract_t2_ups(l3lec),fract_t2_downs(l3lec),'bo','markersize',mSize,'linewidth',1.5);
legend('L3MEC','L3LEC');
xlabel('Prob. persistent Up');
ylabel('Prob. persistent Down');
xlim([0 0.6]);
ylim([0 0.6]);

% xlim([eps 0.6]);
% ylim([eps 0.6]);
% set(gca,'xscale','log','yscale','log');

figufy(h);
fname = [fig_dir 'pers_ups_vs_pers_downs.pdf'];
exportfig(fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);

%%
n_bins = 20;

xx = linspace(0,0.6,n_bins+1);
% xx = logspace(log10(eps),log10(0.6),n_bins+1);
yy_l3mec = histc(fract_rt2_ups(l3mec),xx);
yy_l3mec = yy_l3mec/sum(yy_l3mec);
yy_l3lec = histc(fract_rt2_ups(l3lec),xx);
yy_l3lec = yy_l3lec/sum(yy_l3lec);

h1 = figure;hold on
stairs(xx,yy_l3mec,'r');
stairs(xx,yy_l3lec,'b');
legend('L3MEC','L3LEC');
xlabel('Prob. persistent Up');
ylabel('Relative Frequency');
xlim([0 xx(end)]);
% xlim([eps xx(end)]);
% set(gca,'xscale','log');

xx = linspace(0,0.6,n_bins+1);
% xx = logspace(log10(eps),log10(0.6),n_bins+1);
yy_l3mec = histc(fract_rt2_downs(l3mec),xx);
yy_l3mec = yy_l3mec/sum(yy_l3mec);
yy_l3lec = histc(fract_rt2_downs(l3lec),xx);
yy_l3lec = yy_l3lec/sum(yy_l3lec);

h2 = figure;hold on
stairs(xx,yy_l3mec,'r');
stairs(xx,yy_l3lec,'b');
legend('L3MEC','L3LEC');
xlabel('Prob. persistent Down');
ylabel('Relative Frequency');
xlim([0 xx(end)]);
% xlim([eps xx(end)]);
% set(gca,'xscale','log');

%%
fig_width = 3.27;
rel_heigh = 0.8;

figufy(h1);
fname = [fig_dir 'pers_up_dist.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

figufy(h2);
fname = [fig_dir 'pers_down_dist.pdf'];
exportfig(h2,fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);
close(h2);

%%
mSize = 2;

Fsd = 2016/8;
fract_rt2_ups = [core_data(:).fract_rt2_ups];
fract_rt2_downs = [core_data(:).fract_rt2_downs];

eps = 1e-3;
fract_rt2_ups(fract_rt2_ups < eps) = eps;
fract_rt2_downs(fract_rt2_downs < eps) = eps;

median_uplag = cellfun(@nanmedian,{core_data(:).mp_uplags})/Fsd;
mean_uplag = cellfun(@nanmean,{core_data(:).mp_uplags})/Fsd;

median_downlag = cellfun(@nanmedian,{core_data(:).mp_downlags})/Fsd;
mean_downlag = cellfun(@nanmean,{core_data(:).mp_downlags})/Fsd;

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

%%
fig_width = 3.5;
rel_heigh = 0.8;

figufy(h1);
fname = [fig_dir 'persUp_downLag.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

figufy(h2);
fname = [fig_dir 'persDown_upLag.pdf'];
exportfig(h2,fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);
close(h2);
