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
load ./parietal_cortical_2010
load ./desynch_times_individual
raw_Fs = 2016;

% get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_mp(interneurons) = [];
desynch_times_lf8(interneurons) = [];
desynch_times_lf4(interneurons) = [];


frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);

sess_data = sess_data(thom_pfc);
desynch_times_mp = desynch_times_mp(thom_pfc);
desynch_times_lf4 = desynch_times_lf4(thom_pfc);
desynch_times_lf8 = desynch_times_lf8(thom_pfc);

% %use only data sets with frontal/prefrontal LFP
% sess_data = sess_datap(thom_el);
% desynch_times_lf4 = desynch_times_lf4(thom_el);

raw_Fs = 2016;
dsf = 8;
Fs = raw_Fs/dsf;
hmmFs = raw_Fs/40;
niqf = raw_Fs/2;
lf_lcf = 0.05; %low cut-off for the LF amp filter (Default 0.05)
lf_hcf = 2; %high cut-off for the LF amp filter
num_states = 2; %number of hidden states (must be 2 in current implementation)
min_seg_dur = 60; %minimum duration of a UDS segment

thresh_values = linspace(-3,3,100);

n = length(sess_data);
for d = 1:n
% d=1
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    load ./used_data wcv_minus_spike lf4
    
    %% initializations
    [mp_lf,t_axis,Fs] = get_lf_features(wcv_minus_spike,raw_Fs,Fs,[lf_lcf lf_hcf]);
    [lf4_lf,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fs,[lf_lcf lf_hcf]);
    
    %% compute UDS segments
    UDS_segs4 = hsmm_uds_get_uds_segments(desynch_times_lf4{d},Fs,length(lf4_lf),min_seg_dur);
    nsegs = size(UDS_segs4,1);
    UDS_samps = [];
    for i = 1:nsegs
        UDS_samps = [UDS_samps UDS_segs4(i,1):UDS_segs4(i,2)];
    end
    
    for t = 1:length(thresh_values)
        
        %set all instances where the observation is greater than threshold to the
        %up state
        mp_state_seq = nan(size(mp_lf));
        lf4_state_seq = nan(size(lf4_lf));
        mp_state_seq(UDS_samps) = 1;
        lf4_state_seq(UDS_samps) = 1;
        mp_ups = find(mp_lf(UDS_samps) > thresh_values(t));
        lf4_ups = find(lf4_lf(UDS_samps) > thresh_values(t));
        mp_state_seq(UDS_samps(mp_ups)) = 2;
        lf4_state_seq(UDS_samps(lf4_ups)) = 2;
        
        fp(d,t) = sum(lf4_state_seq == 2 & mp_state_seq == 1)/sum(mp_state_seq == 1);
        fn(d,t) = sum(lf4_state_seq == 1 & mp_state_seq == 2)/sum(mp_state_seq == 2);
        
        
    end
end

cd G:\WC_Germany\parietal_cortical_2010
save roc_data fp fn

