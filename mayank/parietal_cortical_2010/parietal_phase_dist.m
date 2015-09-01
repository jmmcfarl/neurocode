%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'G';

addpath('G:\WC_Germany\hmm_state_detect\\')
addpath('G:\WC_Germany\hsmm_state_detection')
addpath('G:\WC_Germany\parietal_cortical_2010')
addpath('G:\WC_Germany\new_stellate_analysis\')
addpath('G:\Code\CircStat2009')
addpath(strcat(drive_letter,':\Code\fullBNT-1.0.4\kPMstats\'))
addpath(strcat(drive_letter,':\Code\fullBNT-1.0.4\netlab3.3\'))

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load desynch_times_individual

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
Fs = 50.4;
niqf = Fs/2;
lcf = 0.05;
hcf = 2;

% parietal = find_struct_field_vals(sess_data,'region','parietal');
% sess_data = sess_data(parietal);
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_mp(interneurons) = [];
% desynch_times_lf8(interneurons) = [];

parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);
superficial = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:length(sess_data),superficial);
sup_par = intersect(parietal,superficial);
sup_deep = intersect(parietal,deep);
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);

used_data = thom_pfc;
sess_data = sess_data(used_data);
desynch_time_mp = desynch_times_mp(used_data);

n_bins = 400;
phase_vals = linspace(-pi,pi,n_bins);

for d = 1:length(sess_data)
    d
    direct = sess_data(d).directory;
    direct(1) = 'G';
    cd(direct)
    load used_data lf4 wcv_minus_spike
    load hsmm_state_seq_seg_lf_4_28_10_v3
    mp_state_seq = hsmm_state_seq;
    
    lf4_lf = get_lf_features(lf4,raw_Fs,hmm.Fs,[lcf hcf]);
    lf4_p = angle(hilbert(lf4_lf));
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,hmm.Fs,[lcf hcf]);
    
    ustate_lf = [];
    dstate_lf = [];
    ustate_p = [];
    dstate_p = [];
    ov_p = [];
    for j = 1:hmm.Nsegs
        cur_segs = hmm.UDS_segs(j,1):hmm.UDS_segs(j,2);
        ustate_lf = [ustate_lf; lf4_lf(cur_segs(mp_state_seq{j}==2))];
        dstate_lf = [dstate_lf; lf4_lf(cur_segs(mp_state_seq{j}==1))];
        ustate_p = [ustate_p; lf4_p(cur_segs(mp_state_seq{j}==2))];
        dstate_p = [dstate_p; lf4_p(cur_segs(mp_state_seq{j}==1))];    
        ov_p = [ov_p; lf4_p(cur_segs)];
    end
    
    u_p_h = hist(ustate_p,phase_vals);
    d_p_h = hist(dstate_p,phase_vals);
    ov_p_h = hist(lf4_p,phase_vals);
    u_p_h = u_p_h./ov_p_h;
    d_p_h = d_p_h./ov_p_h;
    
    like_diff = u_p_h - d_p_h;
    
%     tlow_bound = find(like_diff > 0,1,'first');
%     zp = find(phase_vals > 0,1,'first');
%     tup_bound = find(like_diff(zp:end) < 0,1,'first');
%     tup_bound = tup_bound + zp-1;
%     low_bound(d) = phase_vals(tlow_bound);
%     up_bound(d) = phase_vals(tup_bound);
%     plot(phase_vals,like_diff), hold on
%     plot(low_bound(d),0,'ro')
%     plot(up_bound(d),0,'ro')
%     xlim([-pi pi])
%      t_names = ['F:\WC_Germany\parietal_cortical_2010\phase_distributions\like_diff_' sess_data(d).name];
%     print(t_names,'-dpng')
%     close
   
    mse = zeros(size(phase_vals));
    for j = 1:length(mse)
       mse(j) = sqrt(sum((like_diff - cos(phase_vals - phase_vals(j))).^2)); 
    end
    [dummy,tcos_phase] = min(mse);
    cos_phase(d) = phase_vals(tcos_phase);
    cos_fit = cos(phase_vals-cos_phase(d));
    plot(phase_vals,like_diff), hold on
    plot(phase_vals,cos_fit,'r')
    xlim([-pi pi])
      t_names = ['F:\WC_Germany\parietal_cortical_2010\phase_distributions\like_diff_lf4_' sess_data(d).name];
    print(t_names,'-dpng')
    close
   
   
%     for j = 1:hmm.Nsegs
%         cur_segs = hmm.UDS_segs(j,1):hmm.UDS_segs(j,2);
%         phase_state_seq{j} = ones(size(cur_segs));
%         phase_state_seq{j}(lf8_p(cur_segs) > low_bound & lf8_p(cur_segs) < up_bound) = 2;
%     end
            
end

