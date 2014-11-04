
clear all
% close all

%%
addpath('C:\WC_Germany\persistent_9_27_2010\')
addpath('C:\WC_Germany\new_mec\')
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd

raw_Fs = 2016;
dsf = 32;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 4;
lcf = 0.05;
lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.025;

maxLag = round(2*Fsd);
backlag = 2*Fsd;
forwardlag = 2*Fsd;
lags = -maxLag:maxLag;
NLags = length(lags);

rate_sm = round(Fsd*0.075);

min_durs = 60; %minimum segment duration
%%
d = 61; 
% d = 38;
cd(combined_dir{d})
    pwd
    load ./used_data
    if ctx_lfp(d) == 7
        lf8 = lf7;
    end
    [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    [desynch_times_lf8,desynch_inds,P_lf8,f,t] = locate_desynch_times_individual(lf8);
    
    if ~isnan(hpc_mua(d))
        load ./mua_data
        load ./sync_times.mat
        synct_d = downsample(synct,dsf);
        hpc_mua_times = mua_times{hpc_mua(d)};
        hpc_mua_times(hpc_mua_times > synct_d(end) | hpc_mua_times < synct_d(1)) = [];
        mua_rate =hist(hpc_mua_times,synct_d)*Fsd;
        mua_rate = jmm_smooth_1d_cor(mua_rate,rate_sm);
        %         mua_rate = zscore(mua_rate(:));
        
        [log_mua,offset(d)] = log_transform_sig(mua_rate);
        lmua_rate = zscore(log_mua(:));
        mua_rate = zscore(mua_rate);
        if length(mua_rate) > length(t_axis)
            mua_rate = mua_rate(1:length(t_axis));
            lmua_rate = lmua_rate(1:length(t_axis));
        end
    else
        mua_rate = nan(size(lf8_lf));
    end
    
    load ./pa_hsmm_state_seq_combined_fin.mat
    load ./pa_hsmm_state_seq7_combined_fin.mat
    hsmm_bbstate_seq8 = hsmm_bbstate_seq7;
    
    bb_Fs = 252;
    temp_lf8 = get_lf_features(lf8,raw_Fs,252,[0.05 10]);
    [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,bb_Fs,length(temp_lf8));
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq);
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq8);
    
    
    up_trans_inds8 = round(up_trans_inds8/(dsf/8));
    up_trans_inds = round(up_trans_inds/(dsf/8));
    down_trans_inds8 = round(down_trans_inds8/(dsf/8));
    down_trans_inds = round(down_trans_inds/(dsf/8));
    
%%
xl1 = [261 277];
xl2 = [234 249];
xl3 = [282 304];

hpc_mua_times_s = (hpc_mua_times-synct(1))/1e6;
winsize = 20;
t1 = [-2.2 -1.7];
t_axis_offset = t_axis - xl2(1); 
figure; set(gca,'fontname','arial','fontsize',14)
hold on
plot(t_axis_offset,lf8_lf,'k','linewidth',1.5)
plot(t_axis_offset,wcv_lf,'r','linewidth',1.5)
plot(t_axis_offset,mua_rate,'color',[0.2 0.8 0.2],'linewidth',1.5)
% plot(t_axis,lmua_rate,'color',[0.2 0.8 0.2])
% for i = 1:length(hpc_mua_times)
%     line([hpc_mua_times_s(i) hpc_mua_times_s(i)],t1,'color','k','linewidth',0.01)
% end
% for i = 1:500
%     sp = (i-1)*winsize;ep = sp+winsize;
%     xlim([sp ep])
%     pause
% end
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude (z)','fontsize',16)
% legend('Ctx LFP','MEC MP','Hpc MUA')
xlim(xl2-xl2(1))