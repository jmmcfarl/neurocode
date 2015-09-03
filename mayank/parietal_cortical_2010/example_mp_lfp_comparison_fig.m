clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_individual

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_mp(interneurons) = [];
desynch_times_lf8(interneurons) = [];
desynch_times_lf4(interneurons) = [];

raw_Fs = 2016;
dsf = 40;
Fsd = raw_Fs/dsf;
lf_lcf = 0.05;
lf_hcf = 20;
Fsd2 = raw_Fs/8;
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);

%%
d=35;  
% d=31
% d=33 %default
% xl = [780 835];
% xl = [2 42];
% xl = [550 575];
% xl = [166.5 198]; %for d1
% xl = [1005 1022]; %for 31
% xl = [150 180];
% xl = [160 195]; %for d = 33
xl = [1020 1050]; %for d=35

cdir = sess_data(d).directory;
cdir(1) = 'G';
cd(cdir)

load used_data lf4 wcv_minus_spike wcv
% load hsmm_state_seq_seg_lf_4_28_10_v3
% load hsmm_state_seq4_seg_lf_4_28_10_v3
% load fixmean_state_seqs4_4_26_10_v2
% load hsmm_state_seq4_seg_lfhf
load hsmm_state_seq_seg_lf_3_31_2011_v2
load hsmm_state_seq4_seg_lf_3_31_2011_v3
load ./fm_state_seq4_lf_3_31_2011_v3.mat
load ./fm_state_seq_lf_3_31_2011_v2.mat

lf4_state_seq = hsmm_bbstate_seq4;
mp_state_seq = hsmm_bbstate_seq;

% lf4_state_seq_lfhf = hsmm_state_seq4_lfhf;

lf_lcf = 0.05;
f_lcf = 20;
[lf4_lf,t_axis,Fs] = get_lf_features(lf4,raw_Fs,252,[lf_lcf f_lcf]);
[wcv_lf,t_axis,Fs] = get_lf_features(wcv_minus_spike,raw_Fs,252,[0.01 f_lcf]);

t = (1:length(wcv))/raw_Fs;
wcv = zscore(wcv);

figure
plot(t_axis,lf4_lf)
hold on
% plot(t_axis,wcv_lf-2.5,'k')
plot(t,wcv-2.5,'k')

[new_seg_inds] = resample_uds_seg_inds(hmm4.UDS_segs,hmm4.Fs,252,length(lf4_lf));
for i = 1:hmm4.Nsegs
%     t1 = hmm4.UDS_segs(i,1)/Fsd;
%     t2 = hmm4.UDS_segs(i,2)/Fsd;
%     tseg = linspace(t1,t2,length(lf4_state_seq{i}));
tseg = (new_seg_inds(i,1):new_seg_inds(i,2))/252;
    plot(tseg,lf4_state_seq{i},'r','linewidth',2)
%     plot(tseg,lf8_state_seq_lfhf{i}+1,'k')
% plot(tseg,fixmean_state_seq4{i}-1,'k','linewidth',2)
% tseg = (fhmm4.UDS_segs(i,1):fhmm4.UDS_segs(i,2))/252;
% plot(tseg,fm_state_seq4{i}-1,'c','linewidth',2)
end
[new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,252,length(wcv_lf));
for i = 1:hmm.Nsegs
%     t1 = hmm.UDS_segs(i,1)/Fsd;
%     t2 = hmm.UDS_segs(i,2)/Fsd;
%     tseg = linspace(t1,t2,length(mp_state_seq{i}));
tseg = (new_seg_inds(i,1):new_seg_inds(i,2))/252;
    plot(tseg,mp_state_seq{i}-3,'g','linewidth',2)
    
end
xlim(xl)

%%
% d = 10
% xl = [940 1005];
% 
% cdir = sess_data(d).directory;
% cdir(1) = 'F';
% cd(cdir)
% 
% load used_data lf8 lf4
% load hsmm_state_seq8_seg_lf_4_10_10
% load hsmm_state_seq4_seg_lf_4_10_10
% lf8_state_seq = hsmm_state_seq8;
% lf4_state_seq = hsmm_state_seq4;
% 
% lf_lcf = 0.05;
% f_lcf = 3;
% [lf8_lf,t_axis,Fs] = get_lf_features(lf8,raw_Fs,Fsd,[lf_lcf f_lcf]);
% [lf4_lf,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd,[lf_lcf f_lcf]);
% 
% figure
% plot((1:length(lf8_lf))/Fsd,lf8_lf)
% hold on
% plot((1:length(lf4_lf))/Fsd,lf4_lf-2.5,'k')
% 
% for i = 1:hmm8.Nsegs
%     tseg = (hmm8.UDS_segs(i,1):hmm8.UDS_segs(i,2))/Fsd;
%     plot(tseg,lf8_state_seq{i},'r','linewidth',2)
% end
% for i = 1:hmm4.Nsegs
%     tseg = (hmm4.UDS_segs(i,1):hmm4.UDS_segs(i,2))/Fsd;
%     plot(tseg,lf4_state_seq{i}-3,'g','linewidth',2)
% end
% xlim(xl)

