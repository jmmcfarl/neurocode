clear all
close all
%%

polar_dir(1).dir = 'C:\wc_data\2012-07-27';
polar_dir(2).dir = 'C:\wc_data\2012-07-30';
polar_dir(3).dir = 'C:\wc_data\2012-08-01';
polar_dir(4).dir = 'C:\wc_data\2012-08-02';
polar_dir(5).dir = 'C:\wc_data\2012-08-06';

%%
for d = 1:length(polar_dir)
    cd(polar_dir(d).dir)
    tempdir = dir;
    tempdir(1:2) = [];
    tempdir([tempdir(:).isdir]==0) = [];
    for j = 1:length(tempdir)
        polar_dir(d).subdir(j).dir = tempdir(j).name;
    end
end

cd C:\WC_Germany\sven_thomas_combined
save polar_dir polar_dir

%%
FieldSelection = [1 0 1 0 1];
ExtractHeader = 1;
ExtractMode = 1;

for d = 1:length(polar_dir)
    for j = 1:length(polar_dir(d).subdir)
        cd(polar_dir(d).dir)
        cd(polar_dir(d).subdir(j).dir)
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
        
        save all_eeg_data2 CSC*
        clear CSC* Samples Timestamp Header
        
        %%
        sync_time_jmm;
        wcv_to_spktim;
        hip_wc_lfp_spk_shift_combined;
        
        %%
        amp_threshold = 30;
        max_overlap = 0.5;
        [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);
        save mua_data3 mua_times mua_amps mua_widths avg_waveform std_waveform      
        
    end
end

%%
%d2 j3
d = 5;
j = 3;
cd(polar_dir(d).dir)
cd(polar_dir(d).subdir(j).dir)

%%
% close all
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;
lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.025;
load ./used_data lf7 wcv_minus_spike wcv
tt = (1:length(wcv))/raw_Fs;
wcv = zscore(wcv);
[lf8_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
figure
hold on
plot(t_axis,wcv_lf+3)
% 0hold on
% plot(tt,wcv+3,'k')
plot(t_axis,lf8_lf,'r')


%%
for d = 1:length(polar_dir)
    for j = 1:length(polar_dir(d).subdir)
        cd(polar_dir(d).dir)
        cd(polar_dir(d).subdir(j).dir)
        
        load ./used_data lf7 wcv_minus_spike
        lf8 = lf7;
        f_name = polar_dir(d).subdir(j).dir;
     [desynch_times_ctx,desynch_inds,P_ctx,f,t] = locate_desynch_times_individual_v2(lf8);
   
    %% compute MP state sequences
    tic
    [hsmm_bbstate_seq,hsmm_state_seq,hmm_bbstate_seq,hsmm,hmm,Fs_bb,Fs] = ...
        pa_get_hsmm_uds_state_seq_new(wcv_minus_spike,raw_Fs,desynch_times_ctx,['mp_' f_name]);
    save pa_hsmm_state_seq_combined_fin_nd hsmm* hmm* Fs*
%     save pa_hsmm_state_seq_combined_fin_newdes hsmm* hmm* Fs*
    clear hsmm* hmm* fract*
    
    %% compute LF7 state sequences
    disp('computing CTX state sequences')
    tic
    [hsmm_bbstate_seq7,hsmm_state_seq7,hmm_bbstate_seq7,hsmm7,hmm7,Fs_bb,Fs] = ...
        pa_get_hsmm_uds_state_seq_new(lf8,raw_Fs,desynch_times_ctx,['lf7_' f_name]);
    save pa_hsmm_state_seq7_combined_fin_nd hsmm* hmm* Fs*
%     save pa_hsmm_state_seq7_combined_fin_newdes hsmm* hmm* Fs*
    clear hsmm* hmm* fract*

        
    end
end
