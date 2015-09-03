%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'G';

addpath(strcat(drive_letter,':\Code\smoothing\software'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))
addpath(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
addpath(strcat(drive_letter,':\Code\maphmmbox\'))

cd(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
load parietal_cortical_2010
load desynch_times_individual
raw_Fs = 2016;

% get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_mp(interneurons) = [];
desynch_times_lf8(interneurons) = [];
desynch_times_lf4(interneurons) = [];

% parietal = find_struct_field_vals(sess_data,'region','parietal');
% sess_data = sess_data(parietal);
% desynch_times = desynch_times(parietal);

thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
sess_data = sess_data(thomas_el);
desynch_times_lf8 = desynch_times_lf8(thomas_el);
desynch_times_lf4 = desynch_times_lf4(thomas_el);

n = length(sess_data);
for d = 1:n
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    load used_data lf8 wcv_minus_spike lf4

    t = (1:length(lf8))/raw_Fs;
    desynch_ind = zeros(size(t));
    for i = 1:size(desynch_times_lf8{d},1)
        sp = find(t > desynch_times_lf8{d}(i,1),1,'first');
        ep = find(t > desynch_times_lf8{d}(i,2),1,'first');
        desynch_ind(sp:ep) = 1;
    end
    for i = 1:size(desynch_times_lf4{d},1)
        sp = find(t > desynch_times_lf4{d}(i,1),1,'first');
        ep = find(t > desynch_times_lf4{d}(i,2),1,'first');
        desynch_ind(sp:ep) = 1;
    end
    desynch_start = find(desynch_ind(1:end-1)==0 & desynch_ind(2:end) == 1);
    desynch_stop = find(desynch_ind(1:end-1)==1 & desynch_ind(2:end) == 0);
    desynch_start = t(desynch_start); desynch_stop = t(desynch_stop);
    desynch_times = [desynch_start(:) desynch_stop(:)];

    [hsmm_bbstate_seq84,hsmm_state_seq84,hmm_bbstate_seq84,hmm84,Fs_bb,Fs] = ...
        parietal_get_hsmm_uds_state_seq_duallfp(lf8,lf4,raw_Fs,desynch_times,strcat(sess_data(d).name,'_lf4'));
    save hsmm_state_seq84_seg_lf_duallfp hsmm* hmm* Fs* 
    clear hsmm* hmm* fract*

end

cd G:\WC_Germany\parietal_cortical_2010
