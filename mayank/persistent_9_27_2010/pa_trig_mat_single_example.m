clear all
close all

load G:\WC_Germany\overall_EC\overall_EC_dir
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

drive_letter = 'G';

Fs = 2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);
lcf2 = 1/niqf;
hcf2 = 10/niqf;
[b2,a2] = butter(2,[lcf2 hcf2]);
lcf3 = 15/niqf;
hcf3 = 80/niqf;
[b3,a3] = butter(2,[lcf3 hcf3]);
dsf = 8;
Fsd = Fs/dsf;
pow_smooth = round(Fsd*0.05);

backlag = 4*Fsd;
forwardlag = 10*Fsd;
lags = (-backlag:forwardlag)/Fsd;

d = 35;

cdir = sess_data(d).directory;
cdir(1) = drive_letter;
disp(sprintf('session %d',d))
cd(cdir);
s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);

load ./used_data lf8 wcv_minus_spike lf3 

wcv_f = filtfilt(b,a,wcv_minus_spike);
lf8_f = filtfilt(b,a,lf8);
lf3_f = filtfilt(b,a,lf3);
wcv_f = downsample(wcv_f,dsf)/sess_data(d).gains(1);
lf8_f = downsample(lf8_f,dsf)/sess_data(d).gains(8);
lf3_f = downsample(lf3_f,dsf)/sess_data(d).gains(8);

wcv_f = zscore(wcv_f);
lf8_f = zscore(lf8_f);
lf3_f = zscore(lf3_f);

t_axis = (1:length(lf8_f))/Fsd;

%% extract up and down transition times for MP and LF8
% load pa_hsmm_state_seq
load pa_hsmm_state_seq_new2

mp_state_seq_c =  hsmm_bbstate_seq;

[new_mp_seg_inds] = round(resample_uds_seg_inds(hmm.UDS_segs,50.4,Fsd,length(t_axis)));

mp_state_seq = nan(size(t_axis));

mp_utrans = [];
mp_dtrans = [];
for n = 1:hmm.Nsegs
    mp_state_seq(new_mp_seg_inds(n,1):new_mp_seg_inds(n,2)) = mp_state_seq_c{n};
    cur_mp_utrans = new_mp_seg_inds(n,1) + find(mp_state_seq_c{n}(1:end-1) == 1 & mp_state_seq_c{n}(2:end) == 2);
    cur_mp_dtrans = new_mp_seg_inds(n,1) + find(mp_state_seq_c{n}(1:end-1) == 2 & mp_state_seq_c{n}(2:end) == 1);
    cur_mp_dtrans(cur_mp_dtrans < cur_mp_utrans(1)) = [];
    cur_mp_utrans(cur_mp_utrans > cur_mp_dtrans(end)) = [];
    mp_utrans = [mp_utrans; cur_mp_utrans];
    mp_dtrans = [mp_dtrans; cur_mp_dtrans];
end
n_mp_ups = length(mp_utrans);

%% initialize
mp_utrig_mp_mat = nan(n_mp_ups,length(lags));
mp_utrig_lf8_mat = nan(n_mp_ups,length(lags));
mp_utrig_lf3_mat = nan(n_mp_ups,length(lags));
mp_dtrig_mp_mat = nan(n_mp_ups,length(lags));
mp_dtrig_lf8_mat = nan(n_mp_ups,length(lags));
mp_dtrig_lf3_mat = nan(n_mp_ups,length(lags));

%     calculate mp utrigs
for i = 1:n_mp_ups
    if mp_utrans(i) > backlag && length(wcv_f) - mp_utrans(i) > forwardlag
        mp_utrig_mp_mat(i,:) = wcv_f(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
        mp_utrig_lf8_mat(i,:) = lf8_f(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
        mp_utrig_lf3_mat(i,:) = lf3_f(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
    end
end

%     calculate mp dtrigs
for i = 1:n_mp_ups
    if mp_dtrans(i) > backlag && length(wcv_f) - mp_dtrans(i) > forwardlag
        mp_dtrig_mp_mat(i,:) = wcv_f(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
        mp_dtrig_lf8_mat(i,:) = lf8_f(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
        mp_dtrig_lf3_mat(i,:) = lf8_f(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
    end
end
mp_updur = (mp_dtrans-mp_utrans)/Fsd;
mp_downdur = (mp_utrans(2:end)-mp_dtrans(1:end-1))/Fsd;

cd G:\WC_Germany\persistent_9_27_2010\
save pa_trigmat_cell35_new2 mp_* lags
