clear all
addpath('~/James_scripts/NeuralynxMatlabImportExport_v501/')
addpath(genpath('~/James_scripts_old/Nlx2Mat_relDec11/'))
addpath('~/Analysis/Mayank/sleep/');

dd = 1;
%cortical
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-09-16_12-26-27';
data(dd).ipsiLFPs = [];
data(dd).contraLFPs = [2 6 10 14 18 22 26 30];
data(dd).MPfile = 'CSC33.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'Ctx';
data(dd).good_bounds = [0 Inf];

dd = dd + 1;
%d2 cortical
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-09-17_12-45-54';
data(dd).ipsiLFPs = [2 6 10 14 18 22 26 30];
data(dd).contraLFPs = [];
data(dd).MPfile = 'CSC33.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'Ctx';
data(dd).good_bounds = [0 290];
%some decent epochs of cortical UDS. The MP is a bit wierd though and
%doesnt have too clear of UDS. Not sure if usable or not...

dd = dd + 1;
%d3 cortical
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-09-30_14-18-08';
data(dd).ipsiLFPs = [];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'Ctx';
data(dd).good_bounds = [0 124];

dd = dd + 1;
%d4 cortical
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-10-01_13-08-37';
data(dd).ipsiLFPs = [];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'Ctx';
data(dd).good_bounds = [0 546];

dd = dd + 1;
%d5 cortical
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-10-07_12-30-20';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 1;
data(dd).contra_L = 0;
data(dd).MPloc = 'Ctx';
data(dd).good_bounds = [0 800];
%some pretty good UDS towards the beginning. Probably usable but not great

dd = dd + 1;
%d6 cortical
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-10-09_12-29-04';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'Ctx';
data(dd).good_bounds = [0 Inf];
%great UDS, especially at the beginning. Both MP and LFP are great.

dd = dd + 1;
%d7 MEC [probable L3MEC. can see large prolonged UPs at time with robust
%spk rates. also lower freq uds].
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-10-14_14-04-11';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 609];
%not great UDS in LFP or MP. Some brief epochs with decent UDS. Might be
%able to see some examples of pers ups. Not pers downs.

dd = dd + 1;
%d8 MEC [cell type somewhat unclear. I'd guess L2MEC based on small peak in
%theta band, and somewhat weaker US and low spk rates]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-10-22/2014-10-22_11-57-18';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 417];
%not really clear UDS. I would say not usable, but there might be some
%instances that could be argued for pers states.

dd = dd + 1;
%d9 MEC [clear L3MEC. high rate, large US, low-freq]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-10-31_13-43-29';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 0;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 Inf];
%great MEC example. Long rec. Good cortical UDS towards the end. Very clear
%MEC MP UDS, with nice pups and pdowns.

dd = dd + 1;
%d10 MEC. [probably L3MEC. high rate, large US, and low-freq]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-04_12-27-44';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 95];
%Short, not much cortical UDS. Probably not useful

dd = dd + 1;
%d11 MEC. [very likely L2 stellate. weaker 'humped' us, lower rate, and
%clear theta peak]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-04_15-31-37';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 468];
%Long enough, with some decent cortical UDS. Wierd MP, dominated by
%DS, with choppy ups. Maybe some pdowns.

dd = dd + 1;
%d12 MEC. [probably L2. Not that clear though. slight theta peak, and
%humped US]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-06_13-41-25';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 400];
%Not much clear UDS in either cortical LFP or MP

dd = dd + 1;
%d13 MEC. [not clear, but L3 if had to guess]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-06_15-20-38';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 75];
%Nice recording with clear UDS, and nice MP. But very short, so only a
%few nice examples.

dd = dd + 1;
%d14 MEC. [clear L3. large US, high rate, low-freq]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-08_14-36-03';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 510];
%Pretty good rec. Decent cortical UDS and MP. Some good pups and
%pdowns.

dd = dd + 1;
%d15 MEC. [pretty clear L2. Weak US, low rate, maybe a hint of theta]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-12_10-57-15';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 242];
%Decent rec. Some decent cortical UDS. MP is pretty DS heavy, with
%choppy UPs. Probably some pers DOWNS. Maybe some pers UPS. Kinda hard to
%tell though.

dd = dd + 1;
%d16 MEC. [unclear]. Sven says probably cortical neuron
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-13_12-48-52';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 310];
%Cortical UDS is pretty good at times. The MP is not really stable,
%and doesn't have clear UDS. So, very likely not usable.

dd = dd + 1;
%d17 MEC. [either bad MP rec or maybe an L2 cell, not clear though]. Sven
%says probably cortical neuron
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-13_13-25-58';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 540];
%Great cortical UDS. MP not really stable, or just not showing clear
%UDS. So, very likely not usable.

dd = dd + 1;
%d18 MEC. [probably L2, but not that clear. weak US, low rate, not much of a clear theta peak though]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-15_12-48-19';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 780];
%Nice cortical UDS. MP does not have very clear UDS, but occasionally
%its decently clear. Probably a few pdowns and pups here, but wont be the
%best examples. Possibly usable.

dd = dd + 1;
%d19 MEC. [likely an L2 because of weak US, somewhat sparser rate, and a
%hint of theta]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-15_13-08-31';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 1094];
%Very nice long rec. Good cortical UDS. MEC MP pretty good, with some
%nice examples of Pup and Pdown at times. The UDS are a bit strange
%(changing very often), but should be usable. MP also appears to be
%degrading slowly towards the end.

dd = dd + 1;
%d20 MEC. [almost certainly not L3. weak US, low rate, theta peak]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-20_12-54-16';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 226];
%Probably too short to use for much. Decent cortical UDS. MEC MP not
%too stable. Also has DS-heavy with choppy UPS. Probably not usable.

dd = dd + 1;
%d21 MEC. [obvious L3. huge US, high rate, low-freq]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-20_13-09-58';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC';
data(dd).good_bounds = [0 284];
%Nice recording. A bit on the short side. Cortical UDS is OK. MEC MP
%is fantastic! Very clear bimodal UDS. Nice strong PUps. Maybe some PDowns
%as well

dd = dd + 1;
%d22. [not L3. Weak US, low rate. not an obvious theta peak, so maybe
%cortical]
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-29_13-12-01';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [13];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'unsure'; %prob cortical (POR or visual)
data(dd).good_bounds = [0 882];
%Maybe MEC or cortical. Short choppy US. Some decent epochs of
%cortical UDS. Maybe some pers DS, but not the clearest examples.

dd = dd + 1;
%d23 MP has high rate and does not show clear UDS. Cortical LFP shows decent
%UDS. MIght be interpreted as lots of pers UPs, but that's not really what
%it looks like
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-12-09_15-04-21';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 1;
data(dd).contra_L = 0;
data(dd).MPloc = 'Ctx'; 
data(dd).good_bounds = [0 779];

dd = dd + 1;
%d24. Some epochs with decent UDS in MP and LFP. Seem pretty synched when
%both have clear UDS. Should have some usable epochs, though not the best
%examples.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-12-09_15-29-56';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'Ctx'; 
data(dd).good_bounds = [0 579];

dd = dd + 1;
%d25. Too short of a rec to really be usable. 
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-12-10_11-43-50';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'Ctx'; 
data(dd).good_bounds = [0 144];

dd = dd + 1;
%d26. MP not good, with choppy unstable UDS. LFP also doesn't show
%clear/sustained UDS. Probably not usable.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-12-11_11-53-25';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'Ctx'; 
data(dd).good_bounds = [0 532];

dd = dd + 1;
%d27. MP and LFP both have some OK UDS. Likely some usable epochs, but not
%very clear, so not good examples.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-12-11_12-33-16';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'Ctx'; 
data(dd).good_bounds = [0 359];

dd = dd + 1;
%d28. UDS not very clear.
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-04_13-08-16';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC'; 
data(dd).good_bounds = [0 155];

dd = dd + 1;
%d29. Pretty good UDS in LFP, and MP. Looks like some pers UPs. Not many
%(if any) clear pers downs
data(dd).dir = '/Users/james/Data/Mayank/sleep/2014-11-06_14-08-44';
data(dd).ipsiLFPs = [1:16];
data(dd).contraLFPs = [1:16];
data(dd).MPfile = 'WC.ncs';
data(dd).ipsi_L = 0;
data(dd).contra_L = 1;
data(dd).MPloc = 'MEC'; 
data(dd).good_bounds = [0 235];

cd ~/Analysis/Mayank/sleep/
save sleep_dirs data

%%
for dd = 28:length(data);
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
    if strcmp(Header{16}(1:12),'-ADBitVolts ')
        mp = CSC_Samples(:)*str2num(Header{16}(13:end));
    else
        error('Header format');
    end
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
            if strcmp(Header{16}(1:12),'-ADBitVolts ')
                csc{i} = CSC_Samples(:)*str2num(Header{16}(13:end));
            else
                error('Header Format');
            end
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