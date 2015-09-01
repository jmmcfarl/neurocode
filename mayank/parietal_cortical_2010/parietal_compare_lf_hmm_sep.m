%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'F';

addpath('F:\WC_Germany\hmm_state_detect\\')
addpath('F:\WC_Germany\hsmm_state_detection')
addpath('F:\WC_Germany\parietal_cortical_2010')
addpath('F:\WC_Germany\new_stellate_analysis\')

cd F:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load desynch_times_mp_lf8_2_24_09

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
Fs = 50.4;

% parietal = find_struct_field_vals(sess_data,'region','parietal');
% sess_data = sess_data(parietal);
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times(interneurons) = [];

% dur_hist_binedges = linspace(0,8,200);
% dur_hist_bincenters = (dur_hist_binedges(1:end-1)+dur_hist_binedges(2:end))/2;

min_state_dur = 1;
max_state_dur = round(Fs*20);
dur_range = (1:max_state_dur)/Fs;

dur_hist_bincenters = (min_state_dur:max_state_dur)/Fs;

hcf_vals = [0.5 10];
smooth_vals = [0.02 0.05 0.2];

for d = 1:length(sess_data)
    d
    direct = sess_data(d).directory;
    direct(1) = 'F';
    cd(direct)
    load used_data lf8
    for i = 1:length(hcf_vals)
        [hmm,state_durations{d,i}] = parietal_get_hmm_lf_test(lf8,raw_Fs,hcf_vals(i),desynch_times{d});
        up_dur_pmf(d,i,:) = hist(state_durations{d,i}{2},dur_hist_bincenters);
        up_dur_pmf(d,i,:) = up_dur_pmf(d,i,:)/sum(up_dur_pmf(d,i,:));
        down_dur_pmf(d,i,:) = hist(state_durations{d,i}{1},dur_hist_bincenters);
        down_dur_pmf(d,i,:) = down_dur_pmf(d,i,:)/sum(down_dur_pmf(d,i,:));
        min_up_dur(d,i) = min(state_durations{d,i}{2});
        min_down_dur(d,i) = min(state_durations{d,i}{2});
    end
end

for d = 1:length(sess_data)
    d
    direct = sess_data(d).directory;
    direct(1) = 'F';
    cd(direct)
    load used_data lf8
    for i = 1:length(smooth_vals)
        [hmm,state_durations_hf{d,i}] = parietal_get_hmm_hf_test(lf8,raw_Fs,smooth_vals(i),desynch_times{d});
        up_dur_hf_pmf(d,i,:) = hist(state_durations_hf{d,i}{2},dur_hist_bincenters);
        up_dur_hf_pmf(d,i,:) = up_dur_hf_pmf(d,i,:)/sum(up_dur_hf_pmf(d,i,:));
        down_dur_hf_pmf(d,i,:) = hist(state_durations_hf{d,i}{1},dur_hist_bincenters);
        down_dur_hf_pmf(d,i,:) = down_dur_hf_pmf(d,i,:)/sum(down_dur_hf_pmf(d,i,:));
        min_up_dur_hf(d,i) = min(state_durations_hf{d,i}{2});
        min_down_dur_hf(d,i) = min(state_durations_hf{d,i}{2});
    end
end

cd F:\WC_Germany\parietal_cortical_2010

%%
n= length(sess_data);

figure
errorbar(hcf_vals,mean(min_up_dur),std(min_up_dur)/sqrt(n))
hold on
errorbar(hcf_vals,mean(min_down_dur),std(min_down_dur)/sqrt(n),'r')

figure
e_bar
