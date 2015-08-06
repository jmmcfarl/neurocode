clear all
addpath('~/James_scripts/NeuralynxMatlabImportExport_v501/')
addpath(genpath('~/James_scripts_old/Nlx2Mat_relDec11/'))
addpath('~/Analysis/Mayank/sleep/');
addpath('~/Analysis/Mayank/sven_thomas_combined/')

cd ~/Analysis/Mayank/sleep/
load sleep_dirs_old

%%
dsf = 16; %temporal down-sample factor
lfp_dsf = 1; %initial down-sample factor
heka_dsf = 10; %down-sample factor for heka data
hf_hcf = 40; %low-cutoff freq for computing LFP HF-power
hf_smooth_sigma = 0.025; %sigma for smoothing window on HF-power (sec)

mua_amp_threshold = 30; %threshold amplitude for MUA (uV)
mua_max_overlap = 0.5;

%%
% for dd = 1:length(data);
dd = 2;
cd(data(dd).dir);


%%
FieldSelectionFlags = [1 0 1 0 1];
ExtractMode = 1;
HeaderExtractionFlag = 1;

i = 1;
fprintf('Loading CSC%d\n',i)
Filename = sprintf('CSC%d.Ncs',i);
%     [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
%         Nlx2MatCSC(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
[CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
    Nlx2MatCSC_v3(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
mp_t = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
Fs(i) = median(CSC_SampleFrequencies);
mp = CSC_Samples(:)*str2num(Header{15}(13:end)); %multiply by gain
clear CSC_Samples
mp_d = decimate(-mp,dsf); %invert and decimate
mp_t = downsample(mp_t,dsf);
mp_Fs = Fs/dsf; %sample freq of MP data



%%
%process all ipsi-lateral LFPs
lfp_list = data(dd).ipsiLFPs;
if ~isempty(lfp_list)
    Fs = nan(length(lfp_list),1);
    for i = 1:length(lfp_list)
        fprintf('Loading CSC%d\n',lfp_list(i))
        Filename = sprintf('CSC%d.Ncs',lfp_list(i));
        [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
            Nlx2MatCSC_v3(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
        csc_time = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
        Fs(i) = median(CSC_SampleFrequencies);
        init_Fsd = Fs(i)/lfp_dsf;
        [b_hf,a_hf] = butter(2,[hf_hcf]/(init_Fsd/2),'high');
        hf_smooth = round(hf_smooth_sigma*init_Fsd);
        
        csc = CSC_Samples(:)*str2num(Header{15}(13:end)); %multiply by gain
        bad_samps = find(isnan(csc_time));
        csc_time(bad_samps) = [];
        csc(bad_samps) = [];
        clear CSC_Samples
        
        %         %compute HF-power
        %         ipsi_csc_hf{i} = decimate(csc,lfp_dsf);
        %         ipsi_csc_hf{i} = filtfilt(b_hf,a_hf,ipsi_csc_hf{i});
        %         ipsi_csc_hf{i} = zscore(sqrt(jmm_smooth_1d_cor(ipsi_csc_hf{i}.^2,hf_smooth)));
        
        ipsi_csc{i} = csc;
        clear CSC_Samples
        
    end
end
% Fs = unique(Fs(2));

%% get spiking data
[mua,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua_v2(mua_amp_threshold,mua_max_overlap);
ipsi_mua_times = mua;
%% align time vectors to zero
t_offset = csc_time(1);
csc_time = (csc_time - t_offset)/1e6; 
mp_t = (mp_t - t_offset)/1e6;
for ii = 1:length(ipsi_mua_times)
    ipsi_mua_times{ii} = (ipsi_mua_times{ii} - t_offset)/1e6;
end

%% align heka and neuralynx data
if ~isempty(data(dd).heka_dir)
    heka_fs = 2e4; %sample freq of heka data
    nlx_lcf = 2; %low-cutoff freq for hpf on nlx data
    max_lag = 100; %maximum latency to search for alignment
    
    load(data(dd).heka_dir);
    heka_var = whos('V*');
    heka_var = heka_var.name;
    eval(sprintf('heka_data = %s.values;',heka_var));
    
    dc_times = (1:length(heka_data))/heka_fs;
    interp_heka = nanzscore(interp1(dc_times,heka_data,mp_t)); %interpolate heka data onto nlx time axis
    [bb,aa] = butter(2,nlx_lcf/(mp_Fs/2),'high');
    mp_filt = zscore(filtfilt(bb,aa,mp_d)); %high-pass filter the down-sampled nlx data
    
    uset = find(~isnan(mp_filt) & ~isnan(interp_heka')); %indices used for alignment
    [xc,xc_lags] = xcov(mp_filt(uset),interp_heka(uset),round(mp_Fs*max_lag),'coef');
    [xc_maxval,xc_maxloc] = max(xc);
    xc_offset = xc_lags(xc_maxloc)/(mp_Fs);
    
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
save(sname,'csc_time','ipsi*','data','mp_d','mp_t','heka_data','heka_time','heka_alignment');

%%
clearvars -except data dd %if processing in a loop only keep track of the globals

% end