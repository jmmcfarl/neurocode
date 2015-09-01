clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8

% %get rid of interneurons
% interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
% sess_data(interneurons) = [];
% desynch_start_times(interneurons) = [];
% desynch_stop_times(interneurons) = [];
sess_data(42:end) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% sess_data = sess_data(parietal);
% desynch_start_times = desynch_start_times(parietal);
% desynch_stop_times = desynch_stop_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times(interneurons) = [];

raw_Fs = 2016;
Fsd = raw_Fs/8;

dt = 1/Fsd;
timing_vec = -4:dt:4;

n = length(sess_data);
for d = 1:n
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    pwd
    load hsmm_state_seq_seg_lf_pert15
    load hsmm_state_seq8_seg_lf_pert15
    load hsmm_state_seq5_seg_lf_pert15
    load hsmm_state_seq4_seg_lf_pert15

    mp_state_seq =  hsmm_bbstate_seq;
    mp_state_seqhf = hsmm_hfstate_seq;
    lf8_state_seq = hsmm_bbstate_seq8;
    lf8_state_seqhf = hsmm_hfstate_seq8;
    lf5_state_seq = hsmm_bbstate_seq5;
    lf5_state_seqhf = hsmm_hfstate_seq5;
    if sess_data(d).thom_elec
        lf4_state_seq = hsmm_bbstate_seq4;
        lf4_state_seqhf = hsmm_hfstate_seq4;
    end

    [lf8_up_state_ind,down_state_ind] = get_used_state_trans_seg(lf8_state_seq,Fsd,hmm8.Fs,hmm8.UDS_segs,2*Fsd);
    lf8_up_state_ind = down_state_ind;
    [lf8hf_up_state_ind,down_state_ind] = get_used_state_trans_seg(lf8_state_seqhf,Fsd,hmm8.Fs,hmm8.UDS_segs,2*Fsd);
    lf8hf_up_state_ind = down_state_ind;
    [lf5_up_state_ind,down_state_ind] = get_used_state_trans_seg(lf5_state_seq,Fsd,hmm5.Fs,hmm5.UDS_segs,2*Fsd);
    lf5_up_state_ind = down_state_ind;
    [lf5hf_up_state_ind,down_state_ind] = get_used_state_trans_seg(lf5_state_seqhf,Fsd,hmm5.Fs,hmm5.UDS_segs,2*Fsd);
    lf5hf_up_state_ind = down_state_ind;
    [mp_up_state_ind,down_state_ind] = get_used_state_trans_seg(mp_state_seq,Fsd,hmm8.Fs,hmm8.UDS_segs,2*Fsd);
    mp_up_state_ind = down_state_ind;
    [mphf_up_state_ind,down_state_ind] = get_used_state_trans_seg(mp_state_seqhf,Fsd,hmm8.Fs,hmm8.UDS_segs,2*Fsd);
    mphf_up_state_ind = down_state_ind;
    
if sess_data(d).thom_elec
    [lf4_up_state_ind,down_state_ind] = get_used_state_trans_seg(lf4_state_seq,Fsd,hmm4.Fs,hmm4.UDS_segs,2*Fsd);
         lf4_up_state_ind = down_state_ind;
    [lf4hf_up_state_ind,down_state_ind] = get_used_state_trans_seg(lf4_state_seqhf,Fsd,hmm4.Fs,hmm4.UDS_segs,2*Fsd);
           lf4hf_up_state_ind = down_state_ind;
end
used_ups = lf8_up_state_ind(:,1);

    n_trans = length(used_ups);
    up_lag_lf8_mp = zeros(n_trans,1);
    up_lag_lf8_lf8hf = zeros(n_trans,1);
    up_lag_lf8_lf5 = zeros(n_trans,1);
    up_lag_lf8_lf5hf = zeros(n_trans,1);
    if sess_data(d).thom_elec
        up_lag_lf8_lf4 = zeros(n_trans,1);
        up_lag_lf8_lf4hf = zeros(n_trans,1);
    end
    for i = 1:n_trans
        [dummy,nearest_mup] = min(abs(used_ups(i) - mp_up_state_ind(:,1)));
        up_lag_lf8_mp(i) = (mp_up_state_ind(nearest_mup,1)-used_ups(i))/Fsd;
        [dummy,nearest_lup] = min(abs(used_ups(i) - lf8hf_up_state_ind(:,1)));
        up_lag_lf8_lf8hf(i) = (lf8hf_up_state_ind(nearest_lup,1)-used_ups(i))/Fsd;
        [dummy,nearest_lup] = min(abs(used_ups(i) - lf5_up_state_ind(:,1)));
        up_lag_lf8_lf5(i) = (lf5_up_state_ind(nearest_lup,1)-used_ups(i))/Fsd;
        [dummy,nearest_lup] = min(abs(used_ups(i) - lf5hf_up_state_ind(:,1)));
        up_lag_lf8_lf5hf(i) = (lf5hf_up_state_ind(nearest_lup,1)-used_ups(i))/Fsd;
        if sess_data(d).thom_elec
            [dummy,nearest_lup] = min(abs(used_ups(i) - lf4_up_state_ind(:,1)));
            up_lag_lf8_lf4(i) = (lf4_up_state_ind(nearest_lup,1)-used_ups(i))/Fsd;
            [dummy,nearest_lup] = min(abs(used_ups(i) - lf4hf_up_state_ind(:,1)));
            up_lag_lf8_lf4hf(i) = (lf4hf_up_state_ind(nearest_lup,1)-used_ups(i))/Fsd;
        end
    end

    used_ups = mp_up_state_ind(:,1);

    n_trans = length(used_ups);
    up_lag_mp_lf8 = zeros(n_trans,1);
    up_lag_mp_lf8hf = zeros(n_trans,1);
    up_lag_mp_lf5 = zeros(n_trans,1);
    up_lag_mp_lf5hf = zeros(n_trans,1);
    if sess_data(d).thom_elec
        up_lag_mp_lf4 = zeros(n_trans,1);
        up_lag_mp_lf4hf = zeros(n_trans,1);
    end
    for i = 1:n_trans
        [dummy,nearest_lup] = min(abs(used_ups(i) - lf8_up_state_ind(:,1)));
        up_lag_mp_lf8(i) = (used_ups(i) - lf8_up_state_ind(nearest_lup,1))/Fsd;
        [dummy,nearest_lup] = min(abs(used_ups(i) - lf5_up_state_ind(:,1)));
        up_lag_mp_lf5(i) = (used_ups(i) - lf5_up_state_ind(nearest_lup,1))/Fsd;
        [dummy,nearest_lup] = min(abs(used_ups(i) - lf5hf_up_state_ind(:,1)));
        up_lag_mp_lf5hf(i) = (used_ups(i) - lf5hf_up_state_ind(nearest_lup,1))/Fsd;
        [dummy,nearest_lup] = min(abs(used_ups(i) - lf8hf_up_state_ind(:,1)));
        up_lag_mp_lf8hf(i) = (used_ups(i)-lf8hf_up_state_ind(nearest_lup,1))/Fsd;
        if sess_data(d).thom_elec
            [dummy,nearest_lup] = min(abs(used_ups(i) - lf4_up_state_ind(:,1)));
            up_lag_mp_lf4(i) = (used_ups(i)-lf4_up_state_ind(nearest_lup,1))/Fsd;
            [dummy,nearest_lup] = min(abs(used_ups(i) - lf4hf_up_state_ind(:,1)));
            up_lag_mp_lf4hf(i) = (used_ups(i)-lf4hf_up_state_ind(nearest_lup,1))/Fsd;
        end
    end

    used_ups = lf8hf_up_state_ind(:,1);

    n_trans = length(used_ups);
    if sess_data(d).thom_elec
        up_lag_lf8hf_lf4hf = zeros(n_trans,1);
    end
    for i = 1:n_trans
        if sess_data(d).thom_elec
            [dummy,nearest_lup] = min(abs(used_ups(i) - lf4hf_up_state_ind(:,1)));
            up_lag_lf8hf_lf4hf(i) = (used_ups(i) - lf4hf_up_state_ind(nearest_lup,1))/Fsd;
        end
    end
    
    if sess_data(d).thom_elec
        used_ups = lf4_up_state_ind(:,1);
        n_trans = length(used_ups);
        up_lag_lf4_mp = zeros(n_trans,1);
        for i = 1:n_trans
            [dummy,nearest_mup] = min(abs(used_ups(i) - mp_up_state_ind(:,1)));
            up_lag_lf4_mp(i) = (used_ups(i) - mp_up_state_ind(nearest_mup,1))/Fsd;
        end 
    end

        
    up_trans_lags_mp_lf8(d,:) = hist(up_lag_mp_lf8,timing_vec);
    up_trans_lags_mp_lf8hf(d,:) = hist(up_lag_mp_lf8hf,timing_vec);
    up_trans_lags_mp_lf5(d,:) = hist(up_lag_mp_lf5,timing_vec);
    up_trans_lags_mp_lf5hf(d,:) = hist(up_lag_mp_lf5hf,timing_vec);
    up_trans_lags_lf8_mp(d,:) = hist(up_lag_lf8_mp,timing_vec);
    up_trans_lags_lf8_lf8hf(d,:) = hist(up_lag_lf8_lf8hf,timing_vec);
    up_trans_lags_lf8_lf5(d,:) = hist(up_lag_lf8_lf5,timing_vec);
    up_trans_lags_lf8_lf5hf(d,:) = hist(up_lag_lf8_lf5hf,timing_vec);
    if sess_data(d).thom_elec
        up_trans_lags_mp_lf4(d,:) = hist(up_lag_mp_lf4,timing_vec);
        up_trans_lags_mp_lf4hf(d,:) = hist(up_lag_mp_lf4hf,timing_vec);
        up_trans_lags_lf8_lf4(d,:) = hist(up_lag_lf8_lf4,timing_vec);
        up_trans_lags_lf8_lf4hf(d,:) = hist(up_lag_lf8_lf4hf,timing_vec);
        up_trans_lags_lf8hf_lf4hf(d,:) = hist(up_lag_lf8hf_lf4hf,timing_vec);
        up_trans_lags_lf4_mp(d,:) = hist(up_lag_lf4_mp,timing_vec);
    end

    mean_mp_lf8(d) = mean(up_lag_mp_lf8);
    mean_mp_lf8hf(d) = mean(up_lag_mp_lf8hf);
    mean_mp_lf5(d) = mean(up_lag_mp_lf5);
    mean_mp_lf5hf(d) = mean(up_lag_mp_lf5hf);
    mean_lf8_mp(d) = mean(up_lag_lf8_mp);
    mean_lf8_lf8hf(d) = mean(up_lag_lf8_lf8hf);
    mean_lf8_lf5(d) = mean(up_lag_lf8_lf5);
    mean_lf8_lf5hf(d) = mean(up_lag_lf8_lf5hf);
    if sess_data(d).thom_elec
        mean_mp_lf4(d) = mean(up_lag_mp_lf4);
        mean_mp_lf4hf(d) = mean(up_lag_mp_lf4hf);
        mean_lf8_lf4(d) = mean(up_lag_lf8_lf4);
        mean_lf8_lf4hf(d) = mean(up_lag_lf8_lf4hf);
        mean_lf8hf_lf4hf(d) = mean(up_lag_lf8hf_lf4hf);
        mean_lf4_mp(d) = mean(up_lag_lf4_mp);
    end
 
     median_mp_lf8(d) = median(up_lag_mp_lf8);
    median_mp_lf8hf(d) = median(up_lag_mp_lf8hf);
     median_mp_lf5(d) = median(up_lag_mp_lf5);
    median_mp_lf5hf(d) = median(up_lag_mp_lf5hf);
    median_lf8_mp(d) = median(up_lag_lf8_mp);
    median_lf8_lf8hf(d) = median(up_lag_lf8_lf8hf);
    median_lf8_lf5(d) = median(up_lag_lf8_lf5);
    median_lf8_lf5hf(d) = median(up_lag_lf8_lf5hf);
    if sess_data(d).thom_elec
        median_mp_lf4(d) = median(up_lag_mp_lf4);
        median_mp_lf4hf(d) = median(up_lag_mp_lf4hf);
        median_lf8_lf4(d) = median(up_lag_lf8_lf4);
        median_lf8_lf4hf(d) = median(up_lag_lf8_lf4hf);
        median_lf8hf_lf4hf(d) = median(up_lag_lf8hf_lf4hf);
        median_lf4_mp(d) = median(up_lag_lf4_mp);
    end

         std_mp_lf8(d) = std(up_lag_mp_lf8);
    std_mp_lf8hf(d) = std(up_lag_mp_lf8hf);
         std_mp_lf5(d) = std(up_lag_mp_lf5);
    std_mp_lf5hf(d) = std(up_lag_mp_lf5hf);
    std_lf8_mp(d) = std(up_lag_lf8_mp);
    std_lf8_lf8hf(d) = std(up_lag_lf8_lf8hf);
    std_lf8_lf5(d) = std(up_lag_lf8_lf5);
    std_lf8_lf5hf(d) = std(up_lag_lf8_lf5hf);
    if sess_data(d).thom_elec
        std_mp_lf4(d) = std(up_lag_mp_lf4);
        std_mp_lf4hf(d) = std(up_lag_mp_lf4hf);
        std_lf8_lf4(d) = std(up_lag_lf8_lf4);
        std_lf8_lf4hf(d) = std(up_lag_lf8_lf4hf);
        std_lf8hf_lf4hf(d) = std(up_lag_lf8hf_lf4hf);
        std_lf4_mp(d) = std(up_lag_lf4_mp);
    end

             iqr_mp_lf8(d) = iqr(up_lag_mp_lf8);
    iqr_mp_lf8hf(d) = iqr(up_lag_mp_lf8hf);
             iqr_mp_lf5(d) = iqr(up_lag_mp_lf5);
    iqr_mp_lf5hf(d) = iqr(up_lag_mp_lf5hf);
    iqr_lf8_mp(d) = iqr(up_lag_lf8_mp);
    iqr_lf8_lf8hf(d) = iqr(up_lag_lf8_lf8hf);
    iqr_lf8_lf5(d) = iqr(up_lag_lf8_lf5);
    iqr_lf8_lf5hf(d) = iqr(up_lag_lf8_lf5hf);
    if sess_data(d).thom_elec
        iqr_mp_lf4(d) = iqr(up_lag_mp_lf4);
        iqr_mp_lf4hf(d) = iqr(up_lag_mp_lf4hf);
        iqr_lf8_lf4(d) = iqr(up_lag_lf8_lf4);
        iqr_lf8_lf4hf(d) = iqr(up_lag_lf8_lf4hf);
        iqr_lf8hf_lf4hf(d) = iqr(up_lag_lf8hf_lf4hf);
        iqr_lf4_mp(d) = iqr(up_lag_lf4_mp);
    end

end

up_trans_lags_mp_lf8(:,[1 end]) = [];
up_trans_lags_mp_lf8hf(:,[1 end]) = [];
up_trans_lags_mp_lf5(:,[1 end]) = [];
up_trans_lags_mp_lf5hf(:,[1 end]) = [];
up_trans_lags_lf8_mp(:,[1 end]) = [];
up_trans_lags_lf8_lf8hf(:,[1 end]) = [];
up_trans_lags_lf8_lf5(:,[1 end]) = [];
up_trans_lags_lf8_lf5hf(:,[1 end]) = [];
up_trans_lags_mp_lf4(:,[1 end]) = [];
up_trans_lags_mp_lf4hf(:,[1 end]) = [];
up_trans_lags_lf8_lf4(:,[1 end]) = [];
up_trans_lags_lf8_lf4hf(:,[1 end]) = [];
up_trans_lags_lf8hf_lf4hf(:,[1 end]) = [];
up_trans_lags_lf4_mp(:,[1 end]) = [];


timing_vec([1 end]) = [];

cd G:\WC_Germany\parietal_cortical_2010\
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);
superficial = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:length(sess_data),superficial);
sup_par = intersect(superficial,parietal);
deep_par = intersect(deep,parietal);
sup_fro = intersect(superficial,frontal);
deep_fro = intersect(deep,frontal);
deep_thom_pfc = intersect(deep_fro,thom_el);
sup_thom_pfc = intersect(sup_fro,thom_el);
%% Parietal Up Trans
figure
used_set = parietal;
h=errorbar(timing_vec,mean(up_trans_lags_mp_lf8(used_set,:)),std(up_trans_lags_mp_lf8(used_set,:))/sqrt(length(used_set))), hold on
errorbar_tick(h,0.001,'units')
h=errorbar(timing_vec,mean(up_trans_lags_mp_lf5(used_set,:)),std(up_trans_lags_mp_lf5(used_set,:))/sqrt(length(used_set)),'k'), hold on
errorbar_tick(h,0.001,'units')
h=errorbar(timing_vec,mean(up_trans_lags_mp_lf8hf(used_set,:)),std(up_trans_lags_mp_lf8hf(used_set,:))/sqrt(length(used_set)),'r'), hold on
errorbar_tick(h,0.001,'units')
h=errorbar(timing_vec,mean(up_trans_lags_mp_lf5hf(used_set,:)),std(up_trans_lags_mp_lf5hf(used_set,:))/sqrt(length(used_set)),'g'), hold on
errorbar_tick(h,0.001,'units')

xlim([-0.5 0.5])
yl = ylim;
ylim([0 yl(2)])

%%
figure
used_set = thom_el;
h=errorbar(timing_vec,mean(up_trans_lags_lf8_lf4(used_set,:)),std(up_trans_lags_lf8_lf4(used_set,:))/sqrt(length(used_set))); hold on
errorbar_tick(h,0.001,'units')
h=errorbar(timing_vec,mean(up_trans_lags_lf8hf_lf4hf(used_set,:)),std(up_trans_lags_lf8hf_lf4hf(used_set,:))/sqrt(length(used_set)),'k');
errorbar_tick(h,0.001,'units')
used_set = thom_par;
h=errorbar(timing_vec,mean(up_trans_lags_lf8_mp(used_set,:)),std(up_trans_lags_lf8_mp(used_set,:))/sqrt(length(used_set)),'r');
errorbar_tick(h,0.001,'units')
used_set = thom_pfc;
h=errorbar(timing_vec,mean(up_trans_lags_lf4_mp(used_set,:)),std(up_trans_lags_lf4_mp(used_set,:))/sqrt(length(used_set)),'g');
errorbar_tick(h,0.001,'units')

xlim([-0.5 0.5])
yl = ylim;
ylim([0 yl(2)])

%%
figure
used_set = parietal;
h=errorbar(timing_vec,mean(up_trans_lags_mp_lf8(used_set,:)),std(up_trans_lags_mp_lf8(used_set,:))/sqrt(length(used_set))); hold on
errorbar_tick(h,0.001,'units')
h=errorbar(timing_vec,mean(up_trans_lags_mp_lf8hf(used_set,:)),std(up_trans_lags_mp_lf8hf(used_set,:))/sqrt(length(used_set)),'r');
errorbar_tick(h,0.001,'units')
% used_set = thom_pfc;
% errorbar(timing_vec,mean(up_trans_lags_mp_lf4(used_set,:)),std(up_trans_lags_mp_lf4(used_set,:))/sqrt(length(used_set)),'k'), hold on
% errorbar(timing_vec,mean(up_trans_lags_mp_lf4hf(used_set,:)),std(up_trans_lags_mp_lf4hf(used_set,:))/sqrt(length(used_set)),'g'), hold on
used_set = thom_par;
h=errorbar(timing_vec,mean(up_trans_lags_mp_lf4(used_set,:)),std(up_trans_lags_mp_lf4(used_set,:))/sqrt(length(used_set)),'c');
errorbar_tick(h,0.001,'units')
h=errorbar(timing_vec,mean(up_trans_lags_mp_lf4hf(used_set,:)),std(up_trans_lags_mp_lf4hf(used_set,:))/sqrt(length(used_set)),'y');
errorbar_tick(h,0.001,'units')

xlim([-0.5 0.5])
yl = ylim;
ylim([0 yl(2)])
