clear all
close all
FieldSelection = [1 0 1 0 1];
ExtractHeader = 1;
ExtractMode = 1;

force_recompute = 0; %overwrite existing data?

addpath(genpath('C:\Code\WC_anal\'))
addpath('C:\Code\NLx2Mat_new\');
% addpath(genpath('C:\WC_Germany\'))


cd C:/WC_Germany/persistent_downs/
% load ./overall_EC_dir
load ./new_pdown_dir
for dd = 1:length(new_pdown_dir)
    dd
    cd(new_pdown_dir{dd})
    %% Process LFP Data
%     if ~exist('./all_eeg_data.mat','file') || force_recompute
        
        temp = dir;
        for i = 1:length(temp) %cycle through data files in this dir
            if length(temp(i).name) > 3
                if strcmp(temp(i).name(end-2:end),'Ncs') | strcmp(temp(i).name(end-2:end),'ncs') %if its a nlx file
                    disp(temp(i).name)
                    [Timestamp, SampleFrequency, Samples, Header] = ...
                        Nlx2MatCSC(temp(i).name, FieldSelection, ExtractHeader, ExtractMode); %matlab conversion
                    
                    if length(temp(i).name) == 8
                        tname = temp(i).name(1:4);
                    else
                        tname = temp(i).name(1:5);
                    end
                    tname = strcat(tname,'.mat');
                    
                    conversion = str2num(Header{15}(14:end)); %get conversion factors
                    
                    eval(strcat(temp(i).name(1:4),'_TimeStamps = Timestamp;'));
                    eval(strcat(temp(i).name(1:4),'_Samples = conversion*Samples;'));
                    eval(strcat(temp(i).name(1:4),'_SampleFrequencies = SampleFrequency;'));
                end
            end
        end
        
        save all_eeg_data CSC*
        clear CSC* Samples Timestamp Header
%     end
    %% Create aligned times
    if ~exist('./sync_times.mat','file') || force_recompute
    sync_time_jmm();
    end
    %% Detect spikes
    if ~exist('./spike_time_jmm.mat','file') || force_recompute
    wcv_to_spktim();
    end
    %% process LFP data
    if ~exist('./used_data.mat','file') || force_recompute
    hip_wc_lfp_spk_shift_combined()
    end
    %% Process MUA data
    if ~exist('./mua_data5.mat','file') || force_recompute
        
        amp_threshold = 30;
        max_overlap = 0.5;
        %     [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);
        [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua_v2(amp_threshold,max_overlap);
%         save mua_data3 mua_times mua_amps mua_widths avg_waveform std_waveform
        save mua_data5 mua_times mua_amps mua_widths avg_waveform std_waveform
    end
    
end