%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'C';

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
% sess_data = sess_data(thom_el);
% desynch_times_lf4 = desynch_times_lf4(thom_el);
G = [0.3 0.8 0.3];
winsize = 20;

n = length(sess_data);
d=4
% for d = 1:n
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    load ./used_data wcv_minus_spike lf4
    
    load ./hsmm_state_seq4_seg_lf_3_31_2011 
    load ./hsmm_state_seq_seg_lf_3_31_2011_v2  
    load ./fm_state_seq_lf_3_31_2011_v2
    load ./fm_state_seq4_lf_3_31_2011
    wcv_d = get_lf_features(wcv_minus_spike,2016,252,[0.05 20]);
    [lf4_d,t_axis] = get_lf_features(lf4,2016,252,[0.05 20]);
    
    [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,252,length(wcv_d));
    [new_seg_inds4] = resample_uds_seg_inds(hmm4.UDS_segs,hmm4.Fs,252,length(wcv_d));
    
    fhmm.Fs = 252;
    fhmmz.Fs = 252;
    fhmm4.Fs = 252;
    fhmmz4.Fs = 252;
    
    hsmm_mp_lf4_ham(d) = hsmm_uds_hamdist(hsmm_bbstate_seq,hsmm_bbstate_seq4,hsmm,hsmm4,252,length(wcv_d));
    hmm_mp_lf4_ham(d) = hsmm_uds_hamdist(hmm_bbstate_seq,hmm_bbstate_seq4,hmm,hmm4,252,length(wcv_d));
    fm_mp_lf4_ham(d) = hsmm_uds_hamdist(fm_state_seq,fm_state_seq4,fhmm,fhmm4,252,length(wcv_d));
    fmz_mp_lf4_ham(d) = hsmm_uds_hamdist(fm_state_seqz,fm_state_seqz4,fhmmz,fhmmz4,252,length(wcv_d));
    fprintf('HSMM agreement: %.6f\n',hsmm_mp_lf4_ham(d));
    fprintf('FM agreement: %.6f\n',fm_mp_lf4_ham(d));
    fprintf('FMz agreement: %.6f\n',fmz_mp_lf4_ham(d));
  
%     [mp_uptrans,mp_downtrans,mp_upsegnums,mp_downsegnums] = compute_state_transitions_seg(new_seg_inds4,hsmm_bbstate_seq);
%     [lf4_uptrans,lf4_downtrans,lf4_upsegnums,lf4_downsegnums] = compute_state_transitions_seg(new_seg_inds4,hsmm_bbstate_seq4);
%     [corresp_lf4_upinds,corresp_lf4_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
%         mp_uptrans,mp_downtrans,lf4_uptrans,lf4_downtrans);
%     [corresp_mp_upinds,corresp_mp_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
%         lf4_uptrans,lf4_downtrans,mp_uptrans,mp_downtrans);
%     hsmm_fract_unmatched_mp_states(d) = sum(isnan(corresp_lf4_upinds))/length(mp_uptrans);
%     hsmm_fract_unmatched_lfp_states(d) = sum(isnan(corresp_mp_upinds))/length(lf4_uptrans);
%     
%     [mp_uptrans,mp_downtrans,mp_upsegnums,mp_downsegnums] = compute_state_transitions_seg(new_seg_inds4,hmm_bbstate_seq);
%     [lf4_uptrans,lf4_downtrans,lf4_upsegnums,lf4_downsegnums] = compute_state_transitions_seg(new_seg_inds4,hmm_bbstate_seq4);
%     [corresp_lf4_upinds,corresp_lf4_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
%         mp_uptrans,mp_downtrans,lf4_uptrans,lf4_downtrans);
%     [corresp_mp_upinds,corresp_mp_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
%         lf4_uptrans,lf4_downtrans,mp_uptrans,mp_downtrans);
%     hmm_fract_unmatched_mp_states(d) = sum(isnan(corresp_lf4_upinds))/length(mp_uptrans);
%     hmm_fract_unmatched_lfp_states(d) = sum(isnan(corresp_mp_upinds))/length(lf4_uptrans);
%     
%     [mp_uptrans,mp_downtrans,mp_upsegnums,mp_downsegnums] = compute_state_transitions_seg(new_seg_inds4,fm_state_seq);
%     [lf4_uptrans,lf4_downtrans,lf4_upsegnums,lf4_downsegnums] = compute_state_transitions_seg(new_seg_inds4,fm_state_seq4);
%     [corresp_lf4_upinds,corresp_lf4_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
%         mp_uptrans,mp_downtrans,lf4_uptrans,lf4_downtrans);
%     [corresp_mp_upinds,corresp_mp_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
%         lf4_uptrans,lf4_downtrans,mp_uptrans,mp_downtrans);
%     fm_fract_unmatched_mp_states(d) = sum(isnan(corresp_lf4_upinds))/length(mp_uptrans);
%     fm_fract_unmatched_lfp_states(d) = sum(isnan(corresp_mp_upinds))/length(lf4_uptrans);
% 
%     [mp_uptrans,mp_downtrans,mp_upsegnums,mp_downsegnums] = compute_state_transitions_seg(new_seg_inds4,fm_state_seqz);
%     [lf4_uptrans,lf4_downtrans,lf4_upsegnums,lf4_downsegnums] = compute_state_transitions_seg(new_seg_inds4,fm_state_seqz4);
%     [corresp_lf4_upinds,corresp_lf4_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
%         mp_uptrans,mp_downtrans,lf4_uptrans,lf4_downtrans);
%     [corresp_mp_upinds,corresp_mp_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
%         lf4_uptrans,lf4_downtrans,mp_uptrans,mp_downtrans);
%     fmz_fract_unmatched_mp_states(d) = sum(isnan(corresp_lf4_upinds))/length(mp_uptrans);
%     fmz_fract_unmatched_lfp_states(d) = sum(isnan(corresp_mp_upinds))/length(lf4_uptrans);

    figure
    plot(t_axis,wcv_d,t_axis,lf4_d,'r')
    hold on
    for i = 1:hmm.Nsegs
       plot(t_axis(new_seg_inds(i,1):new_seg_inds(i,2)),hsmm_bbstate_seq{i},'k','linewidth',2)
       plot(t_axis(new_seg_inds4(i,1):new_seg_inds4(i,2)),hsmm_bbstate_seq4{i}-1.1,'color',G,'linewidth',2)       
       plot(t_axis(fhmm.UDS_segs(i,1):fhmm.UDS_segs(i,2)),fm_state_seq{i}+0.1,'color','k','linestyle','--','linewidth',2)
       plot(t_axis(fhmm4.UDS_segs(i,1):fhmm4.UDS_segs(i,2)),fm_state_seq4{i}-1.2,'color',G,'linestyle','--','linewidth',2)       
    end
    nwins = ceil(t_axis(end)/winsize);
    for i = 1:nwins
       b = (i-1)*winsize;
       e = b+winsize;
       xlim([b e])
       pause
    end
% 
%     figure
%     plot(t_axis,wcv_d,t_axis,lf4_d,'r')
%     hold on
%     for i = 1:hmm.Nsegs
%         plot(t_axis(new_seg_inds4(i,1):new_seg_inds4(i,2)),hsmm_bbstate_seq{i},'k','linewidth',2)
%         plot(t_axis(new_seg_inds4(i,1):new_seg_inds4(i,2)),hsmm_bbstate_seq4{i}-1.1,'color',G,'linewidth',2)
%         plot(t_axis(fhmmz.UDS_segs(i,1):fhmmz.UDS_segs(i,2)),fm_state_seqz{i}+0.1,'color','k','linestyle','--','linewidth',2)
%         plot(t_axis(fhmmz4.UDS_segs(i,1):fhmmz4.UDS_segs(i,2)),fm_state_seqz4{i}-1.2,'color',G,'linestyle','--','linewidth',2)
%     end
%     nwins = ceil(t_axis(end)/winsize);
%     for i = 1:nwins
%         b = (i-1)*winsize;
%         e = b+winsize;
%         xlim([b e])
%         pause
%     end

    %     mp_state_vec = zeros(size(wcv_d));
%     lf4_state_vec = zeros(size(wcv_d));
%     for i = 1:length(mp_uptrans)
%         mp_state_vec(mp_uptrans(i):mp_downtrans(i)) = 1;
%     end
%     for i = 1:length(lf4_uptrans)
%         lf4_state_vec(lf4_uptrans(i):lf4_downtrans(i)) = 1;
%     end
%     uds_vec = zeros(size(wcv_d));
%     for i = 1:size(new_seg_inds4,1)
%         uds_vec(new_seg_inds4(i,1):new_seg_inds4(i,2)) = 1;
%     end
%     mp_state_vec(uds_vec==0) = nan;
%     lf4_state_vec(uds_vec==0) = nan;
    
%         % compute state durations
%     [lfp_state_durations] = compute_state_durations_seg(hsmm_bbstate_seq4,252);
%     [mp_state_durations] = compute_state_durations_seg(hsmm_bbstate_seq,252);

% end

cd G:\WC_Germany\parietal_cortical_2010

hsmm_state_disagreement = (hsmm_fract_unmatched_lfp_states + hsmm_fract_unmatched_mp_states)/2;
hmm_state_disagreement = (hmm_fract_unmatched_lfp_states + hmm_fract_unmatched_mp_states)/2;
fm_state_disagreement = (fm_fract_unmatched_lfp_states + fm_fract_unmatched_mp_states)/2;
fmz_state_disagreement = (fmz_fract_unmatched_lfp_states + fmz_fract_unmatched_mp_states)/2;

figure
plot(hsmm_mp_lf4_ham*100,hmm_mp_lf4_ham*100,'o','markersize',8)
hold on
plot(hsmm_mp_lf4_ham*100,fm_mp_lf4_ham*100,'r*','markersize',8)
plot(hsmm_mp_lf4_ham*100,fmz_mp_lf4_ham*100,'k.','markersize',8)
xlim([4 12]), ylim([3 20])
xlabel('HSMM Hamming distance (%)','fontsize',14)
ylabel('Alternative Hamming distance (%)','fontsize',14)
legend('HMM','Fixed-mean','Fixed-mean 2')
line([0 20],[0 20],'color','k')

figure
plot(hsmm_state_disagreement*100,hmm_state_disagreement*100,'o','markersize',8)
hold on
plot(hsmm_state_disagreement*100,fm_state_disagreement*100,'r*','markersize',8)
plot(hsmm_state_disagreement*100,fmz_state_disagreement*100,'k.','markersize',8)
xlim([0 8]), ylim([0 11])
xlabel('HSMM state disagreement (%)','fontsize',14)
ylabel('Alternative state disagreement (%)','fontsize',14)
legend('HMM','Fixed-mean','Fixed-mean 2')
line([0 12],[0 12],'color','k')

%%
for d = 1:n
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    
    load ./hsmm_state_seq_seg_lf_3_31_2011
    load ./hsmm_state_seq4_seg_lf_3_31_2011
    load ./fm_state_seq_lf_3_31_2011
    load ./fm_state_seq4_lf_3_31_2011
    load ./used_data wcv_minus_spike
    wcv_d = downsample(wcv_minus_spike,8);
    % load ./temp_code_test.mat
    fhmm.Fs = 252;
    fhmmz.Fs = 252;
    fhmm4.Fs = 252;
    fhmmz4.Fs = 252;
    hmm_mp_lf4_ham(d) = hsmm_uds_hamdist(hmm_bbstate_seq,hmm_bbstate_seq4,hmm,hmm4,252,length(wcv_d));
    hsmm_mp_lf4_ham(d) = hsmm_uds_hamdist(hsmm_bbstate_seq,hsmm_bbstate_seq4,hsmm,hsmm4,252,length(wcv_d));
    fm_mp_lf4_ham(d) = hsmm_uds_hamdist(fm_state_seq,fm_state_seq4,fhmm,fhmm4,252,length(wcv_d));
    fmz_mp_lf4_ham(d) = hsmm_uds_hamdist(fm_state_seqz,fm_state_seqz4,fhmmz,fhmmz4,252,length(wcv_d));
end
