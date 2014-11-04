clear all
addpath('~/James_scripts/NeuralynxMatlabImportExport_v501/')
addpath(genpath('~/James_scripts_old/Nlx2Mat_relDec11/'))
addpath('~/Analysis/Mayank/sleep/');

%cortical
data(1).dir = '/Users/james/Analysis/Mayank/sleep/2014-09-16_12-26-27';
data(1).ipsiLFPs = [];
data(1).contraLFPs = [2 6 10 14 18 22 26 30];
data(1).MPfile = 'CSC33.ncs';
data(1).ipsi_L = 0;
data(1).contra_L = 0;

%cortical
data(2).dir = '/Users/james/Analysis/Mayank/sleep/2014-09-17_12-45-54';
data(2).ipsiLFPs = [2 6 10 14 18 22 26 30];
data(2).contraLFPs = [];
data(2).MPfile = 'CSC33.ncs';
data(2).ipsi_L = 0;
data(2).contra_L = 0;
%some decent epochs of cortical UDS. The MP is a bit wierd though and
%doesnt have too clear of UDS. Not sure if usable or not...


%cortical
data(3).dir = '/Users/james/Analysis/Mayank/sleep/2014-09-30_14-18-08';
data(3).ipsiLFPs = [];
data(3).contraLFPs = [1:16];
data(3).MPfile = 'WC.ncs';
data(3).ipsi_L = 0;
data(3).contra_L = 0;

%cortical
data(4).dir = '/Users/james/Analysis/Mayank/sleep/2014-10-01_13-08-37';
data(4).ipsiLFPs = [];
data(4).contraLFPs = [1:16];
data(4).MPfile = 'WC.ncs';
data(4).ipsi_L = 0;
data(4).contra_L = 0;

%cortical
data(5).dir = '/Users/james/Analysis/Mayank/sleep/2014-10-07_12-30-20';
data(5).ipsiLFPs = [1:16];
data(5).contraLFPs = [1:16];
data(5).MPfile = 'WC.ncs';
data(5).ipsi_L = 1;
data(5).contra_L = 0;
%some pretty good UDS towards the beginning. Probably usable but not great


%cortical
data(6).dir = '/Users/james/Analysis/Mayank/sleep/2014-10-09_12-29-04';
data(6).ipsiLFPs = [1:16];
data(6).contraLFPs = [1:16];
data(6).MPfile = 'WC.ncs';
data(6).ipsi_L = 0;
data(6).contra_L = 1;
%great UDS, especially at the beginning. Both MP and LFP are great. 


%MEC
data(7).dir = '/Users/james/Analysis/Mayank/sleep/2014-10-14_14-04-11';
data(7).ipsiLFPs = [1:16];
data(7).contraLFPs = [];
data(7).MPfile = 'WC.ncs';
data(7).ipsi_L = 0;
data(7).contra_L = 0;
%not great UDS in LFP or MP. Some brief epochs with decent UDS. Might be
%able to see some examples of pers ups. Not pers downs.

%MEC
data(8).dir = '/Users/james/Analysis/Mayank/sleep/2014-10-22/2014-10-22_11-57-18';
data(8).ipsiLFPs = [1:16];
data(8).contraLFPs = [];
data(8).MPfile = 'WC.ncs';
data(8).ipsi_L = 0;
data(8).contra_L = 0;
%not really clear UDS. I would say not usable, but there might be some
%instances that could be argued for pers states.

%MEC
data(9).dir = '/Users/james/Analysis/Mayank/sleep/2014-10-31_13-43-29';
data(9).ipsiLFPs = [1:16];
data(9).contraLFPs = [];
data(9).MPfile = 'WC.ncs';
data(9).ipsi_L = 0;
data(9).contra_L = 0;

% cd ~/Analysis/Mayank/sleep/
% save sleep_dirs data

%%
for dd = 9:length(data);
%     dd = 2;
    cd(data(dd).dir);
    
    Fs = 32e3;
    dsf = 160;
    init_dsf = 10;
    mp_temp_dsf = 32;
    heka_dsf = 10;
    
    Fsd = Fs/dsf;
    init_Fsd = Fs/init_dsf;
    [b_hf,a_hf] = butter(2,[40]/(init_Fsd/2),'high');
    hf_smooth = round(0.025*init_Fsd);
    
    mua_amp_threshold = 30;
    %% LOAD MP DATA
    
    FieldSelectionFlags = [1 0 1 0 1];
    ExtractMode = 1;
    HeaderExtractionFlag = 1;
    
    [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
        Nlx2MatCSC_v3(data(dd).MPfile, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
    mp_t = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
    Fs = median(CSC_SampleFrequencies);
    mp = CSC_Samples(:)*str2num(Header{15}(13:end));
    Fsd = Fs/dsf;
    cur_time = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
    clear CSC_Samples
    
    mp_d = decimate(mp,dsf);
    mp_t = downsample(mp_t,dsf);
    mp_d = zscore(-mp_d);
    mp_up = zscore(-decimate(mp,mp_temp_dsf));
    mp_uptime = downsample(cur_time,mp_temp_dsf);
    %% FOR IPSI LFPS
    lfp_list = data(dd).ipsiLFPs;
    if ~isempty(lfp_list)
        for i = 1:length(lfp_list)
            fprintf('Loading CSC%d\n',lfp_list(i))
            Filename = sprintf('CSC%d.ncs',lfp_list(i));
            if data(dd).ipsi_L
                Filename = sprintf('CSC%dL.ncs',lfp_list(i));
            end
            [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
                Nlx2MatCSC_v3(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
            cur_time = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
            Fs(i) = median(CSC_SampleFrequencies);
            csc{i} = CSC_Samples(:)*str2num(Header{15}(13:end));
            bad_samps = find(isnan(cur_time));
            cur_time(bad_samps) = [];
            csc{i}(bad_samps) = [];
            
            ipsi_csc_hf{i} = decimate(csc{i},init_dsf);
            ipsi_csc_hf{i} = filtfilt(b_hf,a_hf,ipsi_csc_hf{i});
            ipsi_csc_hf{i} = zscore(sqrt(jmm_smooth_1d_cor(ipsi_csc_hf{i}.^2,hf_smooth)));
            %         ipsi_csc_hf{i} = zscore(sqrt(smooth(ipsi_csc_hf{i}.^2,hf_smooth*2)));
            ipsi_csc_hf{i} = decimate(ipsi_csc_hf{i},dsf/init_dsf);
            
            ipsi_csc{i} = decimate(csc{i},dsf);
            clear CSC_Samples
            
            Fname = sprintf('SE%d.nse',lfp_list(i));
            if data(dd).ipsi_L
                Filename = sprintf('SE%dL.nse',lfp_list(i));
            end
            [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua_v3(Fname,mua_amp_threshold);
            ipsi_mua_times{i} = mua_times;
        end
        Fs = unique(Fs(2));
        
        csc_time = downsample(cur_time,dsf);
    else
        ipsi_csc_hf = {};
        ipsi_csc = {};
        ipsi_mua_times = {};
    end
    
    %% FOR CONTRA LFPS
    lfp_list = data(dd).contraLFPs;
    if ~isempty(lfp_list)
        for i = 1:length(lfp_list)
            fprintf('Loading CSC%d\n',lfp_list(i))
            Filename = sprintf('CSC%d.ncs',lfp_list(i));
            if data(dd).contra_L
                Filename = sprintf('CSC%dL.ncs',lfp_list(i));
            end
            [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
                Nlx2MatCSC_v3(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
            cur_time = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
            Fs(i) = median(CSC_SampleFrequencies);
            csc{i} = CSC_Samples(:)*str2num(Header{15}(13:end));
            bad_samps = find(isnan(cur_time));
            cur_time(bad_samps) = [];
            csc{i}(bad_samps) = [];
            
            contra_csc_hf{i} = decimate(csc{i},init_dsf);
            contra_csc_hf{i} = filtfilt(b_hf,a_hf,contra_csc_hf{i});
            contra_csc_hf{i} = zscore(sqrt(jmm_smooth_1d_cor(contra_csc_hf{i}.^2,hf_smooth)));
            contra_csc_hf{i} = decimate(contra_csc_hf{i},dsf/init_dsf);
            contra_csc{i} = decimate(csc{i},dsf);
            clear CSC_Samples
            
            Fname = sprintf('SE%d.nse',lfp_list(i));
            if data(dd).contra_L
                Filename = sprintf('SE%dL.nse',lfp_list(i));
            end
            [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua_v3(Fname,mua_amp_threshold);
            contra_mua_times{i} = mua_times;
            
        end
        Fs = unique(Fs(2));
        
        csc_time = downsample(cur_time,dsf);
    else
        contra_csc_hf = {};
        contra_csc = {};
        contra_mua_times = {};
    end
    
    %%
    t_offset = csc_time(1);
    csc_time = (csc_time - t_offset)/1e6;
    mp_uptime = (mp_uptime - t_offset)/1e6;
    mp_t = (mp_t - t_offset)/1e6;
    for ii = 1:length(ipsi_mua_times)
        ipsi_mua_times{ii} = (ipsi_mua_times{ii} - t_offset)/1e6;
    end
    for ii = 1:length(contra_mua_times)
        contra_mua_times{ii} = (contra_mua_times{ii} - t_offset)/1e6;
    end
    %%
    heka_fs = 1e4;
    heka_name = dir('20*');
    if ~isempty(heka_name)
    heka_name = heka_name.name;
    load(heka_name);
    heka_var = whos('V*');
    heka_var = heka_var.name;
    eval(sprintf('heka_data = %s.values;',heka_var));
    
    dc_times = (1:length(heka_data))/heka_fs;
    interp_heka = nanzscore(interp1(dc_times,heka_data,mp_uptime));
    [bb,aa] = butter(2,2/(Fs/mp_temp_dsf/2),'high');
    mp_filt = zscore(filtfilt(bb,aa,mp_up));
    
    uset = find(~isnan(mp_filt) & ~isnan(interp_heka'));
    [xc,xc_lags] = xcov(mp_filt(uset),interp_heka(uset),(Fs/mp_temp_dsf)*100,'coef');
    [xc_maxval,xc_maxloc] = max(xc);
    xc_offset = xc_lags(xc_maxloc)/(Fs/mp_temp_dsf);
    
    heka_data = decimate(heka_data,heka_dsf);
    heka_time = downsample(dc_times,heka_dsf) + xc_offset;
    heka_alignment.xc_maxval = xc_maxval;
    heka_alignment.xc_offset = xc_offset;
    
    xc_maxval
    else
       heka_data = [];
       heka_time = [];
       heka_alignment = [];
    end
    %%
    sname = 'procData';
    save(sname,'csc_time','contra*','ipsi*','data','mp_d','mp_t','heka_data','heka_time','heka_alignment');
    
    %%
    clearvars -except data dd
end