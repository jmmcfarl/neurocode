clear all
close all

FieldSelection = [1 0 1 0 1];
ExtractHeader = 1;
ExtractMode = 1;

%% cycle over old data
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
            
            %         csc = reshape(Samples,numel(Samples),1)*conversion;
            %         csc = reshape(Samples,numel(Samples),1);
            
            %         lfpTimestamps = interpolateTimestamps(Timestamp, 512, length(csc));
            
            %         save(tname,'lfpTimestamps','SampleFrequency','csc')
            %
            %         clear csc Timestamp SampleFrequency Header Samples
        end
    end
end

save all_eeg_data2 CSC*
clear CSC* Samples Timestamp Header

%%
clear all
sync_time_jmm;
clear all
wcv_to_spktim;
clear all
hip_wc_lfp_spk_shift_combined;

%%
clear all
amp_threshold = 30;
max_overlap = 0.5;
[mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);
save mua_data3 mua_times mua_amps mua_widths avg_waveform std_waveform
