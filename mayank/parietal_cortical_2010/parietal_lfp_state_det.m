%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'F';

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
desynch_times_mp = desynch_times_mp(thomas_el);
desynch_times_lf8 = desynch_times_lf8(thomas_el);
desynch_times_lf4 = desynch_times_lf4(thomas_el);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);

n = length(sess_data);
for d = 1:n
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    load used_data lf8 wcv_minus_spike lf4


%% compute fixedmean state sequences
if ismember(d,frontal)
disp('computing fixed mean state sequences')
% [fixmean_sm_state_seq,fixmean_state_seq,Fs,fhmm] = parietal_get_fixthresh_uds_state_seq...
%     (wcv_minus_spike,raw_Fs,desynch_times_mp{d},sess_data(d).name);
% [fixmean_sm_state_seq_hf,fixmean_state_seq_hf,Fs,fhmm_hf] = parietal_get_fixthresh_uds_state_seq_hf...
%     (wcv_minus_spike,raw_Fs,desynch_times_mp{d},sess_data(d).name);
% save fixmean_state_seqs_4_10_10 fixmean* Fs fhmm*
% clear fixmean* fhmm*
[fixmean_sm_state_seq,fixmean_state_seq,Fs,fhmm] = parietal_get_fixthresh_zm_uds_state_seq...
    (wcv_minus_spike,raw_Fs,desynch_times_mp{d},sess_data(d).name);
save fixmean_state_seqs_zm fixmean* Fs fhmm*
clear fixmean* fhmm*
end

% [phase_sm_state_seq8,phase_state_seq8,Fs,fhmm8] = parietal_get_phase_uds_state_seq...
%     (lf8,raw_Fs,desynch_times_lf8{d},sess_data(d).name);
% save phase_state_seq8_4_10_10 phase* Fs fhmm*
% clear phase* fhmm*
% [fixmean_sm_state_seq8,fixmean_state_seq8,Fs,fhmm8] = parietal_get_fixthresh_uds_state_seq...
%     (lf8,raw_Fs,desynch_times_lf8{d},sess_data(d).name);
% [fixmean_sm_state_seq8_hf,fixmean_state_seq8_hf,Fs,fhmm8_hf] = parietal_get_fixthresh_uds_state_seq_hf...
%     (lf8,raw_Fs,desynch_times_lf8{d},sess_data(d).name);
% save fixmean_state_seqs8_4_10_10 fixmean* Fs fhmm*
% clear fixmean* fhmm*
% 
% % [fixmean_state_seq5,Fs,fhmm5] = parietal_get_fixthresh_uds_state_seq(lf5,raw_Fs,desynch_times{d},sess_data(d).name);
% % [fixmean_state_seq5_hf,Fs,fhmm5_hf] = parietal_get_fixthresh_uds_state_seq_hf(lf5,raw_Fs,desynch_times{d},sess_data(d).name);
% 

if sess_data(d).thom_elec
%     [phase_sm_state_seq4,phase_state_seq4,Fs,fhmm4] = parietal_get_phase_uds_state_seq...
%         (lf4,raw_Fs,desynch_times_lf4{d},sess_data(d).name);
%     save phase_state_seq4_4_10_10 phase* Fs fhmm*
%     clear phase* fhmm*
%     [fixmean_sm_state_seq4,fixmean_state_seq4,Fs,fhmm4] = parietal_get_fixthresh_uds_state_seq...
%         (lf4,raw_Fs,desynch_times_lf4{d},sess_data(d).name);
%     [fixmean_sm_state_seq4_hf,fixmean_state_seq4_hf,Fs,fhmm4_hf] = parietal_get_fixthresh_uds_state_seq_hf...
%         (lf4,raw_Fs,desynch_times_lf4{d},sess_data(d).name);
%     save fixmean_state_seqs4_4_10_10 fixmean* Fs fhmm*
%     clear fixmean* fhmm*
    [fixmean_sm_state_seq4,fixmean_state_seq4,Fs,fhmm4] = parietal_get_fixthresh_zm_uds_state_seq...
        (lf4,raw_Fs,desynch_times_lf4{d},sess_data(d).name);
    save fixmean_state_seqs4_zm fixmean* Fs fhmm*
    clear fixmean* fhmm*
end
    
%% compute MP state sequences
% if ismember(d,frontal)
%     disp('computing MP state sequences')
%     tic
%     [hsmm_bbstate_seq,hsmm_state_seq,hmm_bbstate_seq,hmm,Fs_bb,Fs] = ...
%         parietal_get_hsmm_uds_state_seq(wcv_minus_spike,raw_Fs,desynch_times_mp{d},sess_data(d).name);
%     hmm.run_dur = toc;
%     for i = 1:hmm.Nsegs
%         curT = length(hmm.UDS_segs(i,1):hmm.UDS_segs(i,2));
%         fract_low_post{i}(d) = sum(hmm.posteriors{i}(:) < 0.05)/curT;
%     end
%     save hsmm_state_seq_seg_lf_4_10_10 hsmm* hmm* Fs* fract_low_post*
%     clear hsmm* hmm* fract*
%     tic
%     [hsmm_bbstate_seq_hf,hsmm_state_seq_hf,hmm_state_seq_hf,hmm_hf,Fs] = parietal_get_hsmm_uds_state_seq_hf...
%         (wcv_minus_spike,raw_Fs,desynch_times_mp{d},sess_data(d).name);
%     hmm_hf.run_dur = toc;
%     save hsmm_state_seq_seg_hf_4_10_10 hsmm* hmm* Fs* 
%     clear hsmm* hmm* fract*
% end
%% compute LF4 state sequences
%     if sess_data(d).thom_elec
%         disp('computing LF4 state sequences')
%         tic
%         [hsmm_bbstate_seq4,hsmm_state_seq4,hmm_bbstate_seq4,hmm4,Fs_bb,Fs] = ...
%             parietal_get_hsmm_uds_state_seq(lf4,raw_Fs,desynch_times_lf4{d},strcat(sess_data(d).name,'_lf4'));
%         hmm4.run_dur4 = toc;
%         for i = 1:hmm4.Nsegs
%             curT = length(hmm4.UDS_segs(i,1):hmm4.UDS_segs(i,2));
%             fract_low_post4{i}(d) = sum(hmm4.posteriors{i}(:) < 0.05)/curT;
%         end
%         save hsmm_state_seq4_seg_lf_4_10_10 hsmm* hmm* Fs* fract_low_post*
%         clear hsmm* hmm* fract*
%         tic
%         [hsmm_bbstate_seq4_hf,hsmm_state_seq4_hf,hmm_state_seq4_hf,hmm4_hf,Fs] = parietal_get_hsmm_uds_state_seq_hf...
%             (lf4,raw_Fs,desynch_times_lf4{d},sess_data(d).name);
%         hmm4_hf.run_dur4 = toc;
%         save hsmm_state_seq4_seg_hf_4_10_10 hsmm* hmm* Fs*
%         clear hsmm* hmm* fract*
%     end
%     clear hsmm* hmm* fract*

%% compute LF8 state sequence
%     disp('computing LF8 state sequences')
%     tic
%     [hsmm_bbstate_seq8,hsmm_state_seq8,hmm_bbstate_seq8,hmm8,Fs_bb,Fs] = ...
%         parietal_get_hsmm_uds_state_seq(lf8,raw_Fs,desynch_times_lf8{d},strcat(sess_data(d).name,'_lf8'));
%     hmm8.run_dur8 = toc;
%     for i = 1:hmm8.Nsegs
%         curT = length(hmm8.UDS_segs(i,1):hmm8.UDS_segs(i,2));
%         fract_low_post8{i}(d) = sum(hmm8.posteriors{i}(:) < 0.05)/curT;
%     end
%     save hsmm_state_seq8_seg_lf_4_10_10 hsmm* hmm* Fs* fract_low_post*
%     clear hsmm* hmm* fract*
%     tic
%     [hsmm_bbstate_seq8_hf,hsmm_state_seq8_hf,hmm_state_seq8_hf,hmm8_hf,Fs] = ...
%         parietal_get_hsmm_uds_state_seq_hf(lf8,raw_Fs,desynch_times_lf8{d},sess_data(d).name);
%     hmm8_hf.run_dur8 = toc;
%     save hsmm_state_seq8_seg_hf_4_10_10 hsmm* hmm* Fs* 
%     clear hsmm* hmm* fract*

%% compute LF5 state sequence
%     disp('computing LF5 state sequences')
%     tic
%     [hsmm_bbstate_seq5,hsmm_hfstate_seq5,hsmm_state_seq5,hmm_bbstate_seq5,hmm5,Fs_bb,Fs] = ...
%         parietal_get_hsmm_uds_state_seq(lf5,raw_Fs,desynch_times{d},strcat(sess_data(d).name,'_lf5'));
%     hmm5.run_dur5 = toc;
%     for i = 1:hmm5.Nsegs
%         curT = length(hmm5.UDS_segs(i,1):hmm5.UDS_segs(i,2));
%         fract_low_post5{i}(d) = sum(hmm5.posteriors{i}(:) < 0.05)/curT;
%     end
%     tic
%     [hsmm_state_seq5_hf,hmm_state_seq5_hf,hmm5_hf,Fs] = parietal_get_hsmm_uds_state_seq_hf...
%         (lf5,raw_Fs,desynch_times{d},sess_data(d).name);
%     hmm5_hf.run_dur5 = toc;
%     save hsmm_state_seq5_seg_lf_pert15 hsmm* hmm* Fs* fract_low_post*
%     clear hsmm* hmm* fract*

end

cd F:\WC_Germany\parietal_cortical_2010
