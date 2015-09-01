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


frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);
sess_data = sess_data(thom_par);
desynch_times_mp = desynch_times_mp(thom_par);
desynch_times_lf4 = desynch_times_lf4(thom_par);
desynch_times_lf8 = desynch_times_lf8(thom_par);

n = length(sess_data);
for d = 1:n
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    load used_data lf4


%% compute fixedmean state sequences
% disp('computing fixed mean state sequences')
% [fixmean_sm_state_seq,fixmean_state_seq,Fs,fhmm] = parietal_get_fixthresh_uds_state_seq...
%     (wcv_minus_spike,raw_Fs,desynch_times_mp{d},sess_data(d).name);
% [fixmean_sm_state_seq_hf,fixmean_state_seq_hf,Fs,fhmm_hf] = parietal_get_fixthresh_uds_state_seq_hf...
%     (wcv_minus_spike,raw_Fs,desynch_times_mp{d},sess_data(d).name);
% save fixmean_state_seqs_4_28_10_v1 fixmean* Fs fhmm*
% clear fixmean* fhmm*
% [fixmean_sm_state_seq,fixmean_state_seq,Fs,fhmm] = parietal_get_fixthresh_zm_uds_state_seq...
%     (wcv_minus_spike,raw_Fs,desynch_times_mp{d},sess_data(d).name);
% save fixmean_state_seqs_zm_4_28_10_v1 fixmean* Fs fhmm*
% clear fixmean* fhmm*

% if sess_data(d).thom_elec
%     [phase_sm_state_seq4,phase_state_seq4,Fs,fhmm4] = parietal_get_phase_uds_state_seq...
%         (lf4,raw_Fs,desynch_times_lf4{d},sess_data(d).name);
%     save phase_state_seq4_4_26_10_v2 phase* Fs fhmm*
%     clear phase* fhmm*
%     [fixmean_sm_state_seq4,fixmean_state_seq4,Fs,fhmm4] = parietal_get_fixthresh_uds_state_seq...
%         (lf4,raw_Fs,desynch_times_lf4{d},sess_data(d).name);
%     [fixmean_sm_state_seq4_hf,fixmean_state_seq4_hf,Fs,fhmm4_hf] = parietal_get_fixthresh_uds_state_seq_hf...
%         (lf4,raw_Fs,desynch_times_lf4{d},sess_data(d).name);
%     save fixmean_state_seqs4_4_28_10_v1 fixmean* Fs fhmm*
%     clear fixmean* fhmm*
%     [fixmean_sm_state_seq4,fixmean_state_seq4,Fs,fhmm4] = parietal_get_fixthresh_zm_uds_state_seq...
%         (lf4,raw_Fs,desynch_times_lf4{d},sess_data(d).name);
%     save fixmean_state_seqs4_zm_4_28_10_v1 fixmean* Fs fhmm*
%     clear fixmean* fhmm*
% end
    
%% compute MP state sequences
%     tic
%     [hsmm_bbstate_seq,hsmm_state_seq,hmm_bbstate_seq,hmm,Fs_bb,Fs] = ...
%         parietal_get_hsmm_uds_state_seq(wcv_minus_spike,raw_Fs,desynch_times_mp{d},sess_data(d).name);
%     hmm.run_dur = toc;
%     for i = 1:hmm.Nsegs
%         curT = length(hmm.UDS_segs(i,1):hmm.UDS_segs(i,2));
%         fract_low_post{i}(d) = sum(hmm.posteriors{i}(:) < 0.05)/curT;
%     end
%     save hsmm_state_seq_seg_lf_4_28_10_v1 hsmm* hmm* Fs* fract_low_post*
%     clear hsmm* hmm* fract*
    
%     [hsmm_bbstate_seq_fm,hsmm_state_seq_fm,hmm_bbstate_seq_fm,hmm_fm,Fs_bb,Fs] = ...
%         parietal_get_hsmm_uds_state_seq_fixmean(wcv_minus_spike,raw_Fs,desynch_times_mp{d},sess_data(d).name);
%     hmm.run_dur = toc;
%     save hsmm_state_seq_seg_fm_4_28_10_v1 hsmm* hmm* Fs*
%     clear hsmm* hmm* 

%     tic
%     [hsmm_bbstate_seq_hf,hsmm_state_seq_hf,hmm_state_seq_hf,hmm_hf,Fs] = parietal_get_hsmm_uds_state_seq_hf...
%         (wcv_minus_spike,raw_Fs,desynch_times_mp{d},sess_data(d).name);
%     hmm_hf.run_dur = toc;
%     save hsmm_state_seq_seg_hf_4_10_10 hsmm* hmm* Fs* 
%     clear hsmm* hmm* fract*

%% compute LF4 state sequences
%         disp('computing LF4 state sequences')
%         tic
%         [hsmm_bbstate_seq4,hsmm_state_seq4,hmm_bbstate_seq4,hmm4,Fs_bb,Fs] = ...
%             parietal_get_hsmm_uds_state_seq(lf4,raw_Fs,desynch_times_lf4{d},strcat(sess_data(d).name,'_lf4'));
%         hmm4.run_dur4 = toc;
%         for i = 1:hmm4.Nsegs
%             curT = length(hmm4.UDS_segs(i,1):hmm4.UDS_segs(i,2));
%             fract_low_post4{i}(d) = sum(hmm4.posteriors{i}(:) < 0.05)/curT;
%         end
%         save hsmm_state_seq4_seg_lf_4_28_10_v3 hsmm* hmm* Fs* fract_low_post*
%         clear hsmm* hmm* fract*
 
%         [hsmm_bbstate_seq4_fm,hsmm_state_seq4_fm,hmm_bbstate_seq4_fm,hmm4_fm,Fs_bb,Fs] = ...
%             parietal_get_hsmm_uds_state_seq_fixmean(lf4,raw_Fs,desynch_times_lf4{d},strcat(sess_data(d).name,'_lf4'));
%         save hsmm_state_seq4_seg_fm_4_28_10_v1 hsmm* hmm* Fs* 
%         clear hsmm* hmm* 
% 
%         tic
%         [hsmm_bbstate_seq4_hf,hsmm_state_seq4_hf,hmm_state_seq4_hf,hmm4_hf,Fs] = parietal_get_hsmm_uds_state_seq_hf...
%             (lf4,raw_Fs,desynch_times_lf4{d},sess_data(d).name);
%         hmm4_hf.run_dur4 = toc;
%         save hsmm_state_seq4_seg_hf_4_28_10_v3 hsmm* hmm* Fs*
%         clear hsmm* hmm* fract*
%          
        tic
        [hsmm_bbstate_seq4_lfhf,hsmm_state_seq4_lfhf,hmm_state_seq4_lfhf,hmm4_lfhf,Fs] = parietal_get_hsmm_uds_state_seq_lfhf...
            (lf4,raw_Fs,desynch_times_lf4{d},sess_data(d).name);
        hmm4_lfhf.run_dur4 = toc;
        save hsmm_state_seq4_seg_lfhf_4_28_10_v3_3 hsmm* hmm* Fs*
        clear hsmm* hmm* fract*
% 
end

cd G:\WC_Germany\parietal_cortical_2010
