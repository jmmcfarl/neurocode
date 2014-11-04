clear all
close all
FieldSelection = [1 0 1 0 1];
ExtractHeader = 1;
ExtractMode = 1;

force_recompute = 0;

addpath(genpath('C:\Code\WC_anal\'))
addpath('C:\Code\NLx2Mat_new\');
% addpath(genpath('C:\WC_Germany\'))

% cur_dirs{1} = 'C:\wc_data\2012_11_07\2012-11-7_11-58-44';
% cur_dirs{2} = 'C:\wc_data\2012_11_07\2012-11-7_12-17-17';
% cur_dirs{3} = 'C:\wc_data\2012_11_09\2012-11-9_12-24-11';
% cur_dirs{4} = 'C:\wc_data\2012_11_09\2012-11-9_13-22-38';
% new_pdown_dir{5} = 'C:\wc_data\2012-11-13\2012-11-13_12-29-39';
% new_pdown_dir{6} = 'C:\wc_data\2012-11-13\2012-11-13_13-36-57';
% new_pdown_dir{7} = 'C:\wc_data\2012-11-14\2012-11-14_13-3-19';
% new_pdown_dir{8} = 'C:\wc_data\2012-11-14\2012-11-14_13-29-10';

cd C:/WC_Germany/persistent_downs/
load ./overall_EC_dir
% load ./new_pdown_dir
% for dd = 13:length(new_pdown_dir)
for dd = 1:length(sess_data)
    dd
%     cd(new_pdown_dir{dd})
cd(sess_data(dd).directory)
    %% Process LFP Data
    if ~exist('./all_eeg_data.mat','file') || force_recompute
        
        temp = dir;
        for i = 1:length(temp)
            if length(temp(i).name) > 3
                if strcmp(temp(i).name(end-2:end),'Ncs') | strcmp(temp(i).name(end-2:end),'ncs')
                    disp(temp(i).name)
                    [Timestamp, SampleFrequency, Samples, Header] = ...
                        Nlx2MatCSC(temp(i).name, FieldSelection, ExtractHeader, ExtractMode);
                    
                    if length(temp(i).name) == 8
                        tname = temp(i).name(1:4);
                    else
                        tname = temp(i).name(1:5);
                    end
                    tname = strcat(tname,'.mat');
                    
                    conversion = str2num(Header{15}(14:end));
                    
                    eval(strcat(temp(i).name(1:4),'_TimeStamps = Timestamp;'));
                    eval(strcat(temp(i).name(1:4),'_Samples = conversion*Samples;'));
                    eval(strcat(temp(i).name(1:4),'_SampleFrequencies = SampleFrequency;'));
                end
            end
        end
        
        save all_eeg_data CSC*
        clear CSC* Samples Timestamp Header
    end
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
%     %% Process MUA data
%     if ~exist('./mua_data3.mat','file') || force_recompute
%         
%         amp_threshold = 30;
%         max_overlap = 0.5;
%         %     [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);
%         [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua_v2(amp_threshold,max_overlap);
%         save mua_data3 mua_times mua_amps mua_widths avg_waveform std_waveform
%     end
    
end