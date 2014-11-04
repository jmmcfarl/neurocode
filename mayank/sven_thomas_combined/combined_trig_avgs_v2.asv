
clear all
close all

%%
addpath('C:\WC_Germany\persistent_9_27_2010\')
addpath('C:\WC_Germany\new_mec\')
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd.mat
load ./combined_core_analysis_fin_nd.mat
uset = sort([l3mec l3lec]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
combined_dir = combined_dir(uset);
hpc_mua = hpc_mua(uset);
hpc_lfp = hpc_lfp(uset);
ctx_lfp = ctx_lfp(uset);

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;
lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.025;

maxlag = round(4*Fsd);
backlag = 4*Fsd;
forwardlag = 4*Fsd;

params.Fs = raw_Fs;
params.fpass = [0 100];
params.tapers = [3 5];
win = 25;

[b,a] = butter(2,[lcf hcf]/niqf);

thresh_lf8_updur = 0.5;
thresh_lf8_downdur = 0.5;

rate_sm = round(Fsd*0.05);

min_n_states = 5;
%%

for d = 1:length(combined_dir)
    cd(combined_dir{d})
    pwd
    load ./used_data
    for i = 1:7
        [lfp_lf(i,:),t_axis] = get_lf_features(eval(sprintf('lf%d',i+1)),raw_Fs,Fsd,[lcf hcf],0);
    end
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    
    csd = lfp_lf(1:end-2,:) + lfp_lf(3:end,:) - 2*lfp_lf(2:end-1,:);
    
    for i = 1:7
        [lfp_hf(i,:),t_axis] = get_hf_features(eval(sprintf('lf%d',i+1)),raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm,0);
    end
    [desynch_times_lf8,desynch_inds,P_lf8,f,t] = locate_desynch_times_individual_v2(lf8);
        
    load ./pa_hsmm_state_seq_combined_fin_nd
    load ./pa_hsmm_state_seq7_combined_fin_nd
    hsmm_bbstate_seq8 = hsmm_bbstate_seq7;
    
    [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,Fsd,length(lf8_lf));
    dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq);
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq8);
    
    for i = 1:7
        lf8_utrig_lfp_lf(d,i,:) = get_event_trig_avg(lfp_lf(i,:),up_trans_inds8,forwardlag,backlag);
        lf8_dtrig_lfp_lf(d,i,:) = get_event_trig_avg(lfp_lf(i,:),down_trans_inds8,forwardlag,backlag);
        lf8_utrig_lfp_hf(d,i,:) = get_event_trig_avg(lfp_hf(i,:),up_trans_inds8,forwardlag,backlag);
        lf8_dtrig_lfp_hf(d,i,:) = get_event_trig_avg(lfp_hf(i,:),down_trans_inds8,forwardlag,backlag);
    end
    for i = 1:5
        lf8_utrig_lfp_csd(d,i,:) = get_event_trig_avg(csd(i,:),up_trans_inds8,forwardlag,backlag);
        lf8_dtrig_lfp_csd(d,i,:) = get_event_trig_avg(csd(i,:),down_trans_inds8,forwardlag,backlag);
    end
    lf8_utrig_wcv(d,:) = get_event_trig_avg(wcv_lf,up_trans_inds8,forwardlag,backlag);
    lf8_dtrig_wcv(d,:) = get_event_trig_avg(wcv_lf,down_trans_inds8,forwardlag,backlag);
    clear lfp_lf lfp_hf
end

cd C:\WC_Germany\sven_thomas_combined\
save combined_trig_avgs_v2 *_lf8 *_wcv maxlag Fsd *dtrig* *utrig*

%%
lags = (-backlag:forwardlag)/Fsd;
l3mec_m = find(~isnan(hpc_mua));
l3mec_Nm = l3mec(isnan(hpc_mua(l3mec)));

cmap = jet(7);
figure
for i = 1:7
    plot(lags,squeeze(mean(lf8_utrig_lfp_lf(l3mec_m,i,:))),'color',cmap(i,:))
    hold on
end
cmap = jet(7);
figure
for i = 1:7
    plot(lags,squeeze(mean(lf8_utrig_lfp_lf(l3mec_Nm,i,:))),'color',cmap(i,:))
    hold on
end

figure
for i = 1:7
    plot(lags,squeeze(mean(lf8_utrig_lfp_hf(l3mec_m,i,:))),'color',cmap(i,:))
    hold on
end
cmap = jet(7);
figure
for i = 1:7
    plot(lags,squeeze(mean(lf8_utrig_lfp_hf(l3mec_Nm,i,:))),'color',cmap(i,:))
    hold on
end

cmap = jet(5);
figure
for i = 1:5
    plot(lags,squeeze(mean(lf8_utrig_lfp_csd(l3mec_m,i,:))),'color',cmap(i,:))
    hold on
end
cmap = jet(5);
figure
for i = 1:5
    plot(lags,squeeze(nanmean(lf8_utrig_lfp_csd(l3mec_Nm,i,:))),'color',cmap(i,:))
    hold on
end