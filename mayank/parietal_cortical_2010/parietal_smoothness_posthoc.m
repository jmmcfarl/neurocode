clear all
close all

drive_letter = 'G';

addpath('G:\WC_Germany\hmm_state_detect\\')
addpath('G:\WC_Germany\hsmm_state_detection')
addpath('G:\WC_Germany\parietal_cortical_2010')
% addpath('G:\WC_Germany\new_stellate_analysis\')
addpath('G:\WC_Germany\persistent_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load desynch_times_individual

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
Fs = 50.4;

% parietal = find_struct_field_vals(sess_data,'region','parietal');
% sess_data = sess_data(parietal);
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

short_threshold = 0.2;

hcf_vals = [0.5 1 2 4 10];
smooth_vals = [0.025 0.05 0.1 0.2 0.3];

%%
% for d = thom_el
%     direct = sess_data(d).directory;
%     direct(1) = 'G';
%     cd(direct)
%     load smoothness_lf4
%     load smoothness_lf4_hf
%     for i = 1:5
%        [state_durations_lf4] = compute_state_durations_seg(hmm_bbstate_seq4{i},Fsd); 
%        [state_durations_lf4_hf] = compute_state_durations_seg(hmm_bbstate_seq4_hf{i},Fsd); 
%        tot_dur = hmm4{i}.T/hmm4{i}.Fs;
%        tot_dur_hf = hmm4_hf{i}.T/hmm4_hf{i}.Fs;
%        freq_short_ups4(d,i) = sum(state_durations_lf4{2} < short_threshold)/tot_dur;
%        freq_short_downs4(d,i) = sum(state_durations_lf4{1} < short_threshold)/tot_dur;
%        freq_short_ups4_hf(d,i) = sum(state_durations_lf4_hf{2} < short_threshold)/tot_dur;
%        freq_short_downs4_hf(d,i) = sum(state_durations_lf4_hf{1} < short_threshold)/tot_dur;
%     end
% end

%%
for d = 1:9
    
    direct = sess_data(thom_pfc(d)).directory;
    direct(1) = 'G';
    cd(direct)
    
%     load smoothness_lf4
%     load smoothness_mp
%     for j = 1:5
%         for i = 1:5
%             ham_dist_mp_lf4(d,j,i) = compute_state_seq_seg_hamdist_varuds...
%                 (hmm_bbstate_seq{j},hmm_bbstate_seq4{i},hmm{j}.UDS_segs,hmm4{i}.UDS_segs,Fsd);
%         end
%     end
    
    load smoothness_lf4
    load smoothness_hf_mp
    for j = 1:5
        for i = 1:5
            ham_dist_mp_hf_lf4(d,j,i) = compute_state_seq_seg_hamdist_varuds...
                (hmm_bbstate_seq_hf{j},hmm_bbstate_seq4{i},hmm_hf{j}.UDS_segs,hmm4{i}.UDS_segs,Fsd);
        end
    end
    
    
    load smoothness_lf4_hf
    load smoothness_hf_mp
    for j = 1:5
        for i = 1:5
            ham_dist_mp_hf_lf4_hf(d,j,i) = compute_state_seq_seg_hamdist_varuds...
                (hmm_bbstate_seq_hf{j},hmm_bbstate_seq4_hf{i},hmm_hf{j}.UDS_segs,hmm4_hf{i}.UDS_segs,Fsd);
        end
    end  
    
end

%%
figure
cmap = colormap(jet(length(hcf_vals)));
for i = 1:length(hcf_vals)
    ebar_mat(hcf_vals,ham_dist_mp_lf4(:,i,:),.2,cmap(i,:)), hold on
end

figure
cmap = colormap(jet(length(hcf_vals)));
for i = 1:length(hcf_vals)
    ebar_mat(smooth_vals,ham_dist_mp_lf4_hf(:,i,:),.1,cmap(i,:)), hold on
end

