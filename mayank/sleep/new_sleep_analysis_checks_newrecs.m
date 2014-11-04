clear all
addpath('~/James_scripts/NeuralynxMatlabImportExport_v501/')
addpath(genpath('~/James_scripts_old/Nlx2Mat_relDec11/'))

cd /Users/james/Analysis/Mayank/sleep/2014-10-14_14-04-11
% cd /Users/james/Analysis/Mayank/sleep/2014-10-09_12-29-04
cd /Users/james/Analysis/Mayank/sleep/2014-10-07_12-30-20
cd /Users/james/Analysis/Mayank/sleep/2014-10-01_13-08-37
% cd /Users/james/Analysis/Mayank/sleep/2014-09-30_14-18-08
% cd /Users/james/Analysis/Mayank/sleep/2014-09-16_12-26-27
% cd /Users/james/Analysis/Mayank/sleep/2014-09-17_12-45-54

dsf = 32;

%% LOAD RAW DATA

FieldSelectionFlags = [1 0 1 0 1];
ExtractMode = 1;
HeaderExtractionFlag = 1;

i = 33
fprintf('Loading CSC%d\n',i)
% Filename = sprintf('CSC%d.ncs',i);
Filename = 'WC.ncs';
[CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
    Nlx2MatCSC_v3(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
mp_t = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
Fs = median(CSC_SampleFrequencies);
mp = CSC_Samples(:)*str2num(Header{15}(13:end));
clear CSC_Samples
Fsd = Fs/dsf;

mp_d = decimate(mp,dsf);
mp_t = downsample(mp_t,dsf);

[b,a] = butter(2,[0.01 50]/(Fsd/2));
mp_lf = zscore(filtfilt(b,a,-mp_d));

%%
% lfp_list = [6 10 14 18 22 26 30];
lfp_list = [1:4];
for i = 1:length(lfp_list)
    fprintf('Loading CSC%d\n',lfp_list(i))
    Filename = sprintf('CSC%d.ncs',lfp_list(i));
    [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
        Nlx2MatCSC_v3(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
    csc_time{i} = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
    Fs(i) = median(CSC_SampleFrequencies);
    csc{i} = CSC_Samples(:)*str2num(Header{15}(13:end));
    bad_samps = find(isnan(csc_time{i}));
    csc_time{i}(bad_samps) = [];
    csc{i}(bad_samps) = [];
    clear CSC_Samples
end
Fs = unique(Fs(2));

%% DOWN-SAMPLE
[b,a] = butter(2,[0.2 20]/(Fsd/2));
[b2,a2] = butter(2,[40]/(Fsd/2),'high');
hf_smooth = round(0.025*Fsd);
for i = 1:length(csc)
    i
    csc_lf{i} = decimate(csc{i},dsf);
    csc_hf{i} = zscore(filtfilt(b2,a2,csc_lf{i}));
    csc_hf{i} = zscore(sqrt(jmm_smooth_1d_cor(csc_hf{i}.^2,hf_smooth)));
    csc_lf{i} = filtfilt(b,a,csc_lf{i});
    csc_lf{i} = zscore(csc_lf{i});
end
csc_time = downsample(csc_time{1},dsf);

%%
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'));
params.Fs = Fsd;
params.tapers = [2 3];
win = [10 5];
segave = 1;
[S,t,f]=mtspecgramc(csc_lf{3},win,params);
% [S,f]=mtspectrumsegc(csc_lf{3},win,params,segave);
%%
% lfp_list = [6 10 14 18 22 26 30]+1;
lfp_list = [1:16];
for i = 1:length(lfp_list)
    fprintf('Loading CSC%d\n',lfp_list(i))
%     Filename = sprintf('CSC%d.ncs',lfp_list(i));
    Filename = sprintf('CSC%dL.ncs',lfp_list(i));
    [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
        Nlx2MatCSC_v3(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
    cur_time = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
    Fs(i) = median(CSC_SampleFrequencies);
    csc{i} = CSC_Samples(:)*str2num(Header{15}(13:end));
    bad_samps = find(isnan(cur_time));
    %     csc_time{i}(bad_samps) = [];
    csc{i}(bad_samps) = [];
    clear CSC_Samples
end
Fs = unique(Fs(2));


for i = 1:length(csc)
    csc2_lf{i} = decimate(csc{i},dsf);
    %     csc_hf{i} = zscore(filtfilt(b2,a2,csc_lf{i}));
    %     csc_hf{i} = zscore(sqrt(jmm_smooth_1d_cor(csc_hf{i}.^2,hf_smooth)));
    csc2_lf{i} = zscore(filtfilt(b,a,csc2_lf{i}));
end

%%
t_offset = csc_time(1);
csc_time = (csc_time - t_offset)/1e6;
mp_t = (mp_t - t_offset)/1e6;

%%
close all

fig_dir = '/Users/james/Desktop/sleep_examples/';
% base_name = '2014-09-16_';
% base_name = '2014-09-17_';
% base_name = '2014-10-01_';
base_name = '2014-10-09_';
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

mp_scfac = 2;
win_size = 10;
% use_cscs = [2 4 6 8 10 12 14 16];
use_cscs = [2 6 16];
% use_cscs = [3];
spacing = 3;
cmap = jet(length(use_cscs));
cmap(2,:) = [0.2 0.8 0.2];
cmap(3,:) = [1 0 0];
% cmap = [1 0 0; 0 0 1];
n_wins = range(mp_t)/win_size + 1;
f1 = figure();
fig_width = 9; rel_height = 0.6;
fig_cnt = 1;

for ii = 1:n_wins
    cur_tbounds = win_size*[(ii-1) ii];
    cur_linds = find(csc_time >= cur_tbounds(1) & csc_time <= cur_tbounds(2));
    cur_minds = find(mp_t >= cur_tbounds(1) & mp_t <= cur_tbounds(2));
    
    plot(mp_t(cur_minds),mp_lf(cur_minds)*mp_scfac+4,'k','linewidth',0.75)
    hold on
    plot(csc_time(cur_linds),csc_lf{2}(cur_linds)-spacing,'b');
%     plot(csc_time(cur_linds),csc2_lf{3}(cur_linds)-2*spacing,'r');
    plot(csc_time(cur_linds),csc_hf{2}(cur_linds)-2*spacing,'r');
    for jj = 1:length(use_cscs)
%             plot(csc_time(cur_linds),csc_lf{use_cscs(jj)}(cur_linds)-jj*spacing,'color',cmap(jj,:),'linewidth',0.75);
%             plot(csc_time(cur_linds),csc2_lf{use_cscs(jj)}(cur_linds)-jj*spacing*2,'r','linewidth',1);
%         plot(csc_time(cur_linds),csc_lf{use_cscs(jj)}(cur_linds)-jj*spacing,'b','linewidth',1);
%         plot(csc_time(cur_linds),csc_hf{use_cscs(jj)}(cur_linds)-(jj+1)*spacing,'b','linewidth',0.75);
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





