clear all
cd ~/Analysis/Mayank/sleep/
load sleep_dirs
% load sleep_dirs_old

%%
dd = 27;
cd(data(dd).dir)
load procData

% ipsi_csc = contra_csc;
%%
% mp_dsf = 5;
mp_dsf = 1;
csc_dsf = 1;
heka_dsf = 5;

% use_cscs = [2];
use_cscs = [4 6];
% use_cscs = [1:7];

%%
mp_Fs = 1/nanmedian(diff(mp_t));
mp_Fsd = mp_Fs/mp_dsf;
[b,a] = butter(2,[0.05]/(mp_Fsd/2),'high');
mp_d = decimate(mp_d,mp_dsf);
mp_d = zscore(filtfilt(b,a,mp_d));
mp_t = downsample(mp_t,mp_dsf);

if ~isempty(heka_data)
heka_Fs = 1/nanmedian(diff(heka_time));
heka_data = decimate(heka_data,heka_dsf);
heka_Fsd = heka_Fs/heka_dsf;
[b,a] = butter(2,0.05/(heka_Fsd/2),'high');
heka_data = zscore(filtfilt(b,a,heka_data));
heka_time = downsample(heka_time,heka_dsf);
end

csc_Fs = 1/nanmedian(diff(csc_time));
csc_Fsd = csc_Fs/csc_dsf;
[b,a] = butter(2,[0.15 20]/(csc_Fsd/2));
for ii = 1:length(ipsi_csc)
    ipsi_csc{ii} = decimate(ipsi_csc{ii},csc_dsf);
    ipsi_csc{ii} = filtfilt(b,a,ipsi_csc{ii});
    ipsi_csc{ii} = zscore(ipsi_csc{ii});
end
csc_time = downsample(csc_time,csc_dsf);

mua_sm_sig = round(0.025*csc_Fsd);
ipsi_binned_spks = nan(length(csc_time),length(ipsi_mua_times));
ipsi_sm_rate = nan(length(csc_time),length(ipsi_mua_times));
for ii = 1:length(ipsi_mua_times)
    ipsi_binned_spks(:,ii) = hist(ipsi_mua_times{ii},csc_time);
    ipsi_sm_rate(:,ii) = zscore(jmm_smooth_1d_cor(ipsi_binned_spks(:,ii),mua_sm_sig));
end
%%
for ii = 1:length(ipsi_csc_hf)
    ipsi_csc_hf{ii} = decimate(ipsi_csc_hf{ii},csc_dsf);
    rob_sd = robust_std_dev(ipsi_csc_hf{ii});
    ipsi_csc_hf{ii} = ipsi_csc_hf{ii}/rob_sd;
end
%
close all

fig_dir = '/Users/james/Desktop/sleep_examples/';
base_name = '2014-10-31_';
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

mp_scfac = 1.5;
win_size = 10;
spacing = 3;
cmap = jet(length(use_cscs));
if length(use_cscs) == 1
    cmap = [0 0 1];
end
n_wins = range(mp_t)/win_size + 1;
f1 = figure();
fig_width = 9; rel_height = 0.6;
fig_cnt = 1;

for ii = 1:n_wins
    cur_tbounds = win_size*[(ii-1) ii];
    cur_linds = find(csc_time >= cur_tbounds(1) & csc_time <= cur_tbounds(2));
    cur_minds = find(mp_t >= cur_tbounds(1) & mp_t <= cur_tbounds(2));
    hold on
    if ~isempty(heka_data)
    cur_hinds = find(heka_time >= cur_tbounds(1) & heka_time <= cur_tbounds(2));
     plot(heka_time(cur_hinds),heka_data(cur_hinds)*mp_scfac,'r','linewidth',0.75)
    else
        plot(mp_t(cur_minds),mp_d(cur_minds)*mp_scfac+4,'k','linewidth',0.75)
    end
%         plot(csc_time(cur_linds),ipsi_csc{2}(cur_linds),'b','linewidth',1);
%     plot(csc_time(cur_linds),ipsi_csc_hf{3}(cur_linds)-2*spacing,'r');
    for jj = 1:length(use_cscs)
%             plot(csc_time(cur_linds),ipsi_sm_rate(cur_linds,jj)-jj*spacing,'color',cmap(jj,:),'linewidth',0.75);
            plot(csc_time(cur_linds),ipsi_csc{use_cscs(jj)}(cur_linds)-jj*spacing,'color',cmap(jj,:),'linewidth',0.75);
%             plot(csc_time(cur_linds),ipsi_csc_hf{use_cscs(jj)}(cur_linds)-jj*spacing,'--','color',cmap(jj,:),'linewidth',0.75);
    end
    s = input('print ? [p]:','s');
    if s == 'p';
        fig_name = strcat(fig_dir,base_name,sprintf('EX%d',fig_cnt));
        xlabel('Time (s)');
        ylabel('Amplitude');
        axis tight
        figufy(f1);
        exportfig(f1,fig_name,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
        fig_cnt = fig_cnt + 1;
        
    end
    clf
end

%%
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))
uu = find(~isnan(mp_t));
mp_interp = interp1(mp_t(uu),mp_d(uu),csc_time);
% mp_interp = mp_d;
params.Fs = csc_Fsd;
params.tapers = [4 7];
params.fpass = [0 40];
% movingwin = [20 5];
movingwin = [50 25];
clear C

[b,a] = butter(4,[0.01 50]/(csc_Fs/2));
mp_interp(isnan(mp_interp)) = 0;
mp_interp = filtfilt(b,a,mp_interp);

% use_cscs = [2 8 12 14 16];
% use_cscs = [1:2:16];
% use_cscs = [2:2:16];
use_cscs = 2:2:16;

clear C S2 S1
for cc = 1:length(use_cscs)
    cur_lfp = filtfilt(b,a,ipsi_csc{use_cscs(cc)});
[C{cc},phi,S12,S1,S2{cc},t,f]=cohgramc(mp_interp(:),cur_lfp(:),movingwin,params);
end

%%
xl = [1e3 2.8e3];
yl = [0.1 10];
for cc = 1:length(use_cscs)
    subplot(3,1,1)
    pcolor(t,f,log(abs(S1))');shading flat
    caxis([-5 0])
    ylim(yl);
%     xlim(xl);
set(gca,'yscale','log')

    subplot(3,1,2)
pcolor(t,f,log(abs(S2{cc}))');shading flat
caxis([-6 1])
% caxis([25 32])
set(gca,'yscale','log')
    ylim(yl);
%     xlim(xl);

    subplot(3,1,3)
pcolor(t,f,C{cc}');shading flat
set(gca,'yscale','log')
    ylim(yl);
%     xlim(xl);

    pause
    clf
end
