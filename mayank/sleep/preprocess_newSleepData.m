clear all
addpath('~/James_scripts/NeuralynxMatlabImportExport_v501/')
addpath(genpath('~/James_scripts_old/Nlx2Mat_relDec11/'))
addpath('~/Analysis/Mayank/sleep/');

cd ~/Analysis/Mayank/sleep/
load sleep_dirs 

%%
Fs = 32e3; 
dsf = 160; 
mp_temp_dsf = 32; %temporary dsf for nlx MP data for alignment purposes
heka_dsf = 10; %dsf for heka data

init_dsf = 10; %initial step of down-sampling for LFPs (for HF-power comp)
lfp_hf_lcf = 40; %low-cutoff on LFP HF-power
hf_smooth_sigma = 0.025; %smoothing sigma on HF-power

mua_amp_threshold = 30; %threshold on MUA amplitude (uV)
%%
for dd = 38:length(data);
    %     dd = 2;
    cd(data(dd).dir);
            
    %% LOAD MP DATA
    
    FieldSelectionFlags = [1 0 1 0 1];
    ExtractMode = 1;
    HeaderExtractionFlag = 1;
    
    [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
        Nlx2MatCSC_v3(data(dd).MPfile, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
    mp_t = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
    Fs = median(CSC_SampleFrequencies);
    if strcmp(Header{16}(1:12),'-ADBitVolts ')
        mp = CSC_Samples(:)*str2num(Header{16}(13:end)); %mutiply by gain
    else
        error('Header format');
    end
%     Fsd = Fs/dsf;
    cur_time = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
    clear CSC_Samples
    
    mp_d = decimate(-mp,dsf);
    mp_t = downsample(mp_t,dsf);
    
    %keep a version of MP with higher sample-rate for aligning DC data
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
            if strcmp(Header{16}(1:12),'-ADBitVolts ')
                csc{i} = CSC_Samples(:)*str2num(Header{16}(13:end)); %multiply by gain
            else
                error('Header Format');
            end
            bad_samps = find(isnan(cur_time));
            cur_time(bad_samps) = [];
            csc{i}(bad_samps) = [];

            init_Fsd = Fs(i)/init_dsf;
            [b_hf,a_hf] = butter(2,lfp_hf_lcf/(init_Fsd/2),'high');
            hf_smooth = round(hf_smooth_sigma*init_Fsd);
            
            ipsi_csc_hf{i} = decimate(csc{i},init_dsf); %partial down-sampling to speed up calculation of HF-power
            ipsi_csc_hf{i} = filtfilt(b_hf,a_hf,ipsi_csc_hf{i}); %high-pass filter
            ipsi_csc_hf{i} = zscore(sqrt(jmm_smooth_1d_cor(ipsi_csc_hf{i}.^2,hf_smooth))); %compute smoothed power
            ipsi_csc_hf{i} = decimate(ipsi_csc_hf{i},dsf/init_dsf); %finish down-sampling
            
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
            if strcmp(Header{16}(1:12),'-ADBitVolts ')
                csc{i} = CSC_Samples(:)*str2num(Header{16}(13:end));
            else
                error('Header Format');
            end
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
        Fs = unique(Fs);
        
        csc_time = downsample(cur_time,dsf);
    else
        contra_csc_hf = {};
        contra_csc = {};
        contra_mua_times = {};
    end
    
    %% align time vectors to 0
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
    %% align nlx and DC data
    heka_fs = 1e4;
    nlx_hp_lcf = 2;
    max_lag = 100; %maximum latency to search for alignment
    
    heka_name = dir('20*');
    if ~isempty(heka_name)
        heka_name = heka_name.name;
        load(heka_name);
        heka_var = whos('V*');
        heka_var = heka_var.name;
        eval(sprintf('heka_data = %s.values;',heka_var));
        
        %interpolate heka data onto nlx time series
        dc_times = (1:length(heka_data))/heka_fs;
        interp_heka = nanzscore(interp1(dc_times,heka_data,mp_uptime));
        %high-pass filter nlx data
        [bb,aa] = butter(2,nlx_hp_lcf/(Fs/mp_temp_dsf/2),'high'); 
        mp_filt = zscore(filtfilt(bb,aa,mp_up));
        
        uset = find(~isnan(mp_filt) & ~isnan(interp_heka'));
        [xc,xc_lags] = xcov(mp_filt(uset),interp_heka(uset),(Fs/mp_temp_dsf)*max_lag,'coef');
        [xc_maxval,xc_maxloc] = max(xc);
        xc_offset = xc_lags(xc_maxloc)/(Fs/mp_temp_dsf);
        
        if xc_maxval < 0.5
            warning(sprintf('XC on heka alignment is %.4f\n',xc_maxval));
        end
        
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
    save(sname,'csc_time','contra*','ipsi*','data','mp_d','mp_t','heka_data','heka_time','heka_alignment');
    
    %%
%     clearvars -except data dd
end