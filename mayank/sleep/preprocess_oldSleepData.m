clear all
addpath('~/James_scripts/NeuralynxMatlabImportExport_v501/')
addpath(genpath('~/James_scripts_old/Nlx2Mat_relDec11/'))
addpath('~/Analysis/Mayank/sleep/');
addpath('~/Analysis/Mayank/sven_thomas_combined/')

%MEC
data(1).dir = '/Users/james/Data/Mayank/sleep/2012-7-5_Patch/2012-7-5_17-11-29';
data(1).ipsiLFPs = [2:8];
data(1).contraLFPs = [];
data(1).MPfile = 'CSC1.ncs';
data(1).ipsi_L = 0;
data(1).contra_L = 0;
data(1).heka_dir = '/Users/james/Data/Mayank/sleep/2012-7-5_Patch/2012_07_05_1.mat';
data(1).celltype = 'MECL3';
data(1).MPloc = 'MEC';
data(1).good_bounds = [0 192];
% MP not really that stable. Maybe some decent cortical UDS but cant see clear persistence. Probably not analyzable

%MEC
data(2).dir = '/Users/james/Data/Mayank/sleep/2012-7-8_2/2012-7-8_13-9-45';
data(2).ipsiLFPs = [2:8];
data(2).contraLFPs = [];
data(2).MPfile = 'CSC1.ncs';
data(2).ipsi_L = 0;
data(2).contra_L = 0;
data(2).heka_dir = '';
data(2).celltype = 'MECL3';
data(2).MPloc = 'MEC';
data(2).good_bounds = [0 Inf];
% MP strange. LFP doesn't really show clear UDS either. Not usable


%MEC
data(3).dir = '/Users/james/Data/Mayank/sleep/2012-7-10_Sleep_WC/2012-7-10_13-49-40';
data(3).ipsiLFPs = [2:8];
data(3).contraLFPs = [];
data(3).MPfile = 'CSC1.ncs';
data(3).ipsi_L = 0;
data(3).contra_L = 0;
data(3).heka_dir = '/Users/james/Data/Mayank/sleep/2012-7-10_Sleep_WC/2012_07_10_sleep_1.mat';
data(3).celltype = 'MECL3';
data(3).MPloc = 'MEC';
data(3).good_bounds = [0 490];
%Wierd MP with 'incomplete' looking UP states. Decent cortical UDS. Maybe
%examples of pers downs? But they aren't really clear because of the wierd
%MP UDS


%MEC
data(4).dir = '/Users/james/Data/Mayank/sleep/2012-7-13#1/2012-7-13_13-31-49';
data(4).ipsiLFPs = [2:8];
data(4).contraLFPs = [];
data(4).MPfile = 'CSC1.ncs';
data(4).ipsi_L = 0;
data(4).contra_L = 0;
data(4).heka_dir = '/Users/james/Data/Mayank/sleep/2012-7-13#1/2012_07_13_sleep_3.mat';
data(4).celltype = 'MECL3';
data(4).MPloc = 'MEC';
data(4).good_bounds = [0 515];
%some pretty good cortical UDS. MP is pretty good too. Good examples of
%pers ups, and maybe a few pers downs.


%MEC
data(5).dir = '/Users/james/Data/Mayank/sleep/2012-7-13#2/2012-7-13_16-54-15';
data(5).ipsiLFPs = [2:8];
data(5).contraLFPs = [];
data(5).MPfile = 'CSC1.ncs';
data(5).ipsi_L = 0;
data(5).contra_L = 0;
data(5).heka_dir = '/Users/james/Data/Mayank/sleep/2012-7-13#2/2012_07_13_sleep_6.mat';
data(5).celltype = 'MECL3';
data(5).MPloc = 'MEC';
%not much clear UDS (a few good brief epochs). MP is pretty wierd though.
%Probably not usable. Not really clear examples of pers.


cd ~/Analysis/Mayank/sleep/
save sleep_dirs_old data

%%
% for dd = 1:length(data);
dd = 2;
cd(data(dd).dir);

dsf = 16;
init_dsf = 1;
heka_dsf = 10;

mua_amp_threshold = 30;

%%
FieldSelectionFlags = [1 0 1 0 1];
ExtractMode = 1;
HeaderExtractionFlag = 1;

i = 1
fprintf('Loading CSC%d\n',i)
Filename = sprintf('CSC%d.Ncs',i);
%     [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
%         Nlx2MatCSC(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
[CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
    Nlx2MatCSC_v3(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
mp_t = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
Fs(i) = median(CSC_SampleFrequencies);
mp = CSC_Samples(:)*str2num(Header{15}(13:end));
clear CSC_Samples
mp_d = decimate(-mp,16);
mp_t = downsample(mp_t,16);
mp_Fs = Fs/16;



%%
Fsd = Fs/dsf;
init_Fsd = Fs/init_dsf;
[b_hf,a_hf] = butter(2,[40]/(init_Fsd/2),'high');
hf_smooth = round(0.025*init_Fsd);

lfp_list = data(dd).ipsiLFPs;
if ~isempty(lfp_list)
    for i = 1:length(lfp_list)
        fprintf('Loading CSC%d\n',lfp_list(i))
        Filename = sprintf('CSC%d.Ncs',lfp_list(i));
        [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
            Nlx2MatCSC_v3(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
        csc_time = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
        Fs(i) = median(CSC_SampleFrequencies);
        csc = CSC_Samples(:)*str2num(Header{15}(13:end));
        bad_samps = find(isnan(csc_time));
        csc_time(bad_samps) = [];
        csc(bad_samps) = [];
        clear CSC_Samples
        
        ipsi_csc_hf{i} = decimate(csc,init_dsf);
        ipsi_csc_hf{i} = filtfilt(b_hf,a_hf,ipsi_csc_hf{i});
        ipsi_csc_hf{i} = zscore(sqrt(jmm_smooth_1d_cor(ipsi_csc_hf{i}.^2,hf_smooth)));
        %         ipsi_csc_hf{i} = zscore(sqrt(smooth(ipsi_csc_hf{i}.^2,hf_smooth*2)));
        ipsi_csc_hf{i} = decimate(ipsi_csc_hf{i},dsf/init_dsf);
        
        ipsi_csc{i} = csc;
        clear CSC_Samples
        
    end
end
% Fs = unique(Fs(2));

%%
amp_threshold = 25;
max_overlap = 0.5;
[mua,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua_v2(amp_threshold,max_overlap);
ipsi_mua_times = mua;
%%
t_offset = csc_time(1);
csc_time = (csc_time - t_offset)/1e6;
mp_t = (mp_t - t_offset)/1e6;
for ii = 1:length(ipsi_mua_times)
    ipsi_mua_times{ii} = (ipsi_mua_times{ii} - t_offset)/1e6;
end

%%
if ~isempty(data(dd).heka_dir)
heka_fs = 2e4;
load(data(dd).heka_dir);
heka_var = whos('V*');
heka_var = heka_var.name;
eval(sprintf('heka_data = %s.values;',heka_var));

dc_times = (1:length(heka_data))/heka_fs;
interp_heka = nanzscore(interp1(dc_times,heka_data,mp_t));
[bb,aa] = butter(2,2/(mp_Fs/2),'high');
mp_filt = zscore(filtfilt(bb,aa,mp_d));

uset = find(~isnan(mp_filt) & ~isnan(interp_heka'));
[xc,xc_lags] = xcov(mp_filt(uset),interp_heka(uset),round(mp_Fs*100),'coef');
[xc_maxval,xc_maxloc] = max(xc);
xc_offset = xc_lags(xc_maxloc)/(mp_Fs);

heka_data = decimate(heka_data,heka_dsf);
heka_time = downsample(dc_times,heka_dsf) + xc_offset;
heka_alignment.xc_maxval = xc_maxval;
heka_alignment.xc_offset = xc_offset;

else
    heka_data = [];
    heka_time = [];
    heka_alignment = [];
end

%%
    sname = 'procData';
    save(sname,'csc_time','ipsi*','data','mp_d','mp_t','heka_data','heka_time','heka_alignment');

    %%
        clearvars -except data dd

% end