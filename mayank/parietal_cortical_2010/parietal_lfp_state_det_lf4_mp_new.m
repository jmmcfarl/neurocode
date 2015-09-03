%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'G';

addpath(strcat(drive_letter,':\Code\smoothing\software'))
% addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
% addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))
addpath(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
% addpath(strcat(drive_letter,':\Code\maphmmbox\'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_uds_code'))

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

% sess_data = sess_data(thom_pfc);
% desynch_times_mp = desynch_times_mp(thom_pfc);
% desynch_times_lf4 = desynch_times_lf4(thom_pfc);
% desynch_times_lf8 = desynch_times_lf8(thom_pfc);
% 
% % use only data sets with frontal/prefrontal LFP
sess_data = sess_data(thom_el);
desynch_times_lf4 = desynch_times_lf4(thom_el);

n = length(sess_data);
for d = 1:n
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
%     load ./used_data wcv_minus_spike lf4
    load ./used_data lf4
    
    % compute LF4 state sequences
%     disp('computing LF4 state sequence LF')
%     tic
%     [hsmm_bbstate_seq4,hsmm_state_seq4,hmm_bbstate_seq4,hsmm4,hmm4,Fs_bb,Fs] = ...
%         parietal_get_hsmm_uds_state_seq_new(lf4,raw_Fs,desynch_times_lf4{d},strcat(sess_data(d).name,'_lf4'));
%     run_dur_hmm4(d) = hmm4.run_time;
%     run_dur_hsmm4(d) = hsmm4.run_time;
%     bad_models = 0;
%     for i = 1:hsmm4.Nsegs
%         curT = length(hsmm4.UDS_segs(i,1):hsmm4.UDS_segs(i,2));
%         fract_low_post4{i}(d) = sum(hsmm4.posteriors{i}(:) < 0.05)/curT;
%         bad_models = bad_models + length(hsmm4.rob_model_inds{i});
%     end
%     fract_bad_model4(d) = bad_models/hsmm4.Fs/hsmm4.uds_dur;
%     n_hmm_it4(d) = length(hsmm4.hmm_LP);
%     n_hsmm_it4(d) = length(hsmm4.hsmm_LP);
%     nrun_dur_hmm4(d) = run_dur_hmm4(d)/hmm4.uds_dur*60; %seconds per minute of data
%     nrun_dur_hsmm4(d) = run_dur_hsmm4(d)/hmm4.uds_dur*60; %seconds per minute of data
% %     save hsmm_state_seq4_seg_lf_4_5_2011 hsmm* hmm* Fs* fract_low_post* fract_bad_model4
%     clear hsmm* hmm*
% 
%     % compute LF4 state sequences
%     disp('computing LF4 state sequence LF FM')
%     tic
%     [hsmm_bbstate_seq4,hsmm_state_seq4,hmm_bbstate_seq4,hsmm4,hmm4,Fs_bb,Fs] = ...
%         parietal_get_hsmm_uds_state_seq_fixmean_new(lf4,raw_Fs,desynch_times_lf4{d},strcat(sess_data(d).name,'_lf4'));
%     run_dur_hmm4(d) = hmm4.run_time;
%     run_dur_hsmm4(d) = hsmm4.run_time;
%     bad_models = 0;
%     for i = 1:hsmm4.Nsegs
%         curT = length(hsmm4.UDS_segs(i,1):hsmm4.UDS_segs(i,2));
% %         fract_low_post4{i}(d) = sum(hsmm4.posteriors{i}(:) < 0.05)/curT;
%         bad_models = bad_models + length(hsmm4.rob_model_inds{i});
%     end
%     fract_bad_model4(d) = bad_models/hsmm4.Fs/hsmm4.uds_dur;
%     n_hmm_it4(d) = length(hsmm4.hmm_LP);
%     n_hsmm_it4(d) = length(hsmm4.hsmm_LP);
%     save hsmm_state_seq4_seg_lffm_4_5_2011 hsmm* hmm* Fs* fract_bad_model4
%     clear hsmm* hmm*

%     % compute MP state sequences
%     disp('computing MP state sequence LF')
%     tic
%     [hsmm_bbstate_seq,hsmm_state_seq,hmm_bbstate_seq,hsmm,hmm,Fs_bb,Fs] = ...
%         parietal_get_hsmm_uds_state_seq_new(wcv_minus_spike,raw_Fs,desynch_times_lf4{d},strcat(sess_data(d).name,'_MP'));
%     run_dur_hmm(d) = hmm.run_time;
%     run_dur_hsmm(d) = hsmm.run_time;
%     bad_models = 0;
%     for i = 1:hsmm.Nsegs
%         curT = length(hsmm.UDS_segs(i,1):hsmm.UDS_segs(i,2));
%         fract_low_post{i}(d) = sum(hsmm.posteriors{i}(:) < 0.05)/curT;
%         bad_models = bad_models + length(hsmm.rob_model_inds{i});
%     end
%     fract_bad_model(d) = bad_models/hsmm.Fs/hsmm.uds_dur;
%     n_hmm_it(d) = length(hsmm.hmm_LP);
%     n_hsmm_it(d) = length(hsmm.hsmm_LP);
%    save hsmm_state_seq_seg_lf_4_5_2011 hsmm* hmm* Fs* fract_low_post* fract_bad_model
%     clear hsmm* hmm*
% 
%     % compute MP state sequences
%     disp('computing MP state sequence LF FM')
%     tic
%     [hsmm_bbstate_seq,hsmm_state_seq,hmm_bbstate_seq,hsmm,hmm,Fs_bb,Fs] = ...
%         parietal_get_hsmm_uds_state_seq_fixmean_new(wcv_minus_spike,raw_Fs,desynch_times_lf4{d},strcat(sess_data(d).name,'_MP'));
%     run_dur_hmm(d) = hmm.run_time;
%     run_dur_hsmm(d) = hsmm.run_time;
%     bad_models = 0;
%     for i = 1:hsmm.Nsegs
%         curT = length(hsmm.UDS_segs(i,1):hsmm.UDS_segs(i,2));
% %         fract_low_post{i}(d) = sum(hsmm.posteriors{i}(:) < 0.05)/curT;
%         bad_models = bad_models + length(hsmm.rob_model_inds{i});
%     end
%     fract_bad_model(d) = bad_models/hsmm.Fs/hsmm.uds_dur;
%     n_hmm_it(d) = length(hsmm.hmm_LP);
%     n_hsmm_it(d) = length(hsmm.hsmm_LP);
%    save hsmm_state_seq_seg_lffm_4_5_2011 hsmm* hmm* Fs* fract_bad_model
%     clear hsmm* hmm*
% 
%     disp('computing LF4 state sequence HF')
%     tic
%     [hsmm_bbstate_seq4_hf,hsmm_state_seq4_hf,hmm_bbstate_seq4_hf,hsmm4_hf,hmm4_hf,Fs_bb,Fs] = ...
%         parietal_get_hsmm_uds_state_seq_hf_new(lf4,raw_Fs,desynch_times_lf4{d},strcat(sess_data(d).name,'_lf4'));
%     run_dur4_hf(d) = toc;
%     hsmm4_hf.run_dur4 = run_dur4_hf(d);
%     bad_models = 0;
%     for i = 1:hsmm4_hf.Nsegs
%         curT = length(hsmm4_hf.UDS_segs(i,1):hsmm4_hf.UDS_segs(i,2));
%         fract_low_post4_hf{i}(d) = sum(hsmm4_hf.posteriors{i}(:) < 0.05)/curT;
%         bad_models = bad_models + length(hsmm4_hf.rob_model_inds{i});
%     end
%     fract_bad_model4_hf(d) = bad_models/hsmm4_hf.Fs/hsmm4_hf.uds_dur;
%     n_hmm_it4_hf(d) = length(hsmm4_hf.hmm_LP);
%     n_hsmm_it4_hf(d) = length(hsmm4_hf.hsmm_LP);
%     save hsmm_state_seq4_seg_hf_4_5_2011 hsmm* hmm* Fs* fract_low_post* fract_bad_model4_hf
%     clear hsmm* hmm*   
%  
%     disp('computing LF4 state sequence HF')
%     tic
%     [hsmm_bbstate_seq4_lfhf,hsmm_state_seq4_lfhf,hmm_bbstate_seq4_lfhf,hsmm4_lfhf,hmm4_lfhf,Fs_bb,Fs] = ...
%         parietal_get_hsmm_uds_state_seq_lfhf_new(lf4,raw_Fs,desynch_times_lf4{d},strcat(sess_data(d).name,'_lf4'));
%     run_dur4_lfhf(d) = toc;
%     hsmm4_lfhf.run_dur4 = run_dur4_lfhf(d);
%     bad_models = 0;
%     for i = 1:hsmm4_lfhf.Nsegs
%         curT = length(hsmm4_lfhf.UDS_segs(i,1):hsmm4_lfhf.UDS_segs(i,2));
%         fract_low_post4_lfhf{i}(d) = sum(hsmm4_lfhf.posteriors{i}(:) < 0.05)/curT;
%         bad_models = bad_models + length(hsmm4_lfhf.rob_model_inds{i});
%     end
%     fract_bad_model4_lfhf(d) = bad_models/hsmm4_lfhf.Fs/hsmm4_lfhf.uds_dur;
%     n_hmm_it4_lfhf(d) = length(hsmm4_lfhf.hmm_LP);
%     n_hsmm_it4_lfhf(d) = length(hsmm4_lfhf.hsmm_LP);
%     save hsmm_state_seq4_seg_lfhf_4_5_2011 hsmm* hmm* Fs* fract_low_post* fract_bad_model4_lfhf
%     clear hsmm* hmm*
%     
    % compute LF4 fixedmean state sequences
    disp('computing LF4 state sequence fm')
    tic
    [fm_state_seq4,fhmm4,Fs] = parietal_get_fixthresh_uds_state_seq_new(lf4,raw_Fs,desynch_times_lf4{d},1);
    [fm_state_seqz4,fhmmz4,Fs] = parietal_get_fixthresh_uds_state_seq_new(lf4,raw_Fs,desynch_times_lf4{d},2);
    [fm_state_seqx4,fhmmx4,Fs] = parietal_get_fixthresh_uds_state_seq_new(lf4,raw_Fs,desynch_times_lf4{d},3);
    save fm_state_seq4_lf_4_11_2011 fhmm* fm_state_seq*
    clear fhmm4* fm_state_seq*

    % compute MP fixedmean state sequences
%     disp('computing MP state sequence fm')
%     tic
%     [fm_state_seq,fhmm,Fs] = parietal_get_fixthresh_uds_state_seq_new(wcv_minus_spike,raw_Fs,desynch_times_lf4{d},1);
%     [fm_state_seqz,fhmmz,Fs] = parietal_get_fixthresh_uds_state_seq_new(wcv_minus_spike,raw_Fs,desynch_times_lf4{d},2);
%     [fm_state_seqx,fhmmx,Fs] = parietal_get_fixthresh_uds_state_seq_new(wcv_minus_spike,raw_Fs,desynch_times_lf4{d},3);
%     save fm_state_seq_lf_4_11_2011 fhmm* fm_state_seq*
%     clear fhmm* fm_state_seq*

end

cd G:\WC_Germany\parietal_cortical_2010
% save mp_state_det_data_4_5_2011 fract_bad* run_dur* n_hmm_it* n_hsmm_it*

% %%
% for d = 1:n
%         fprintf('session %d\n',d)
%     cdir = sess_data(d).directory;
%     cdir(1) = drive_letter;
%     cd(cdir)
% 
%     load ./hsmm_state_seq_seg_lf_3_31_2011
%     load ./hsmm_state_seq4_seg_lf_3_31_2011
% % load ./temp_code_test.mat
% hmm_mp_lf4_ham2(d) = compute_state_seq_seg_hamdist_varuds(hmm_bbstate_seq,hmm_bbstate_seq4,hmm.UDS_segs,hmm4.UDS_segs,252);
%     hsmm_mp_lf4_ham2(d) = compute_state_seq_seg_hamdist_varuds(hsmm_bbstate_seq,hsmm_bbstate_seq4,hsmm.UDS_segs,hsmm4.UDS_segs,252);    
% end
