clear all
close all

% cd C:\WC_Germany\sven_thomas_combined
% load ./combined_dir.mat
% l3mec_m = find(~isnan(hpc_mua));
% d = 60
% cd(combined_dir{d})

% load ./distal_dir
% d=3
% cd(distal_dir{d})

cd C:\WC_Germany\persistent_downs\
load ./new_pdown_dir.mat
d = 23

cd(new_pdown_dir{d})

%%
amp_threshold = 30;
max_overlap = 0.5;
% [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);
load ./mua_data3

figure
plot(avg_waveform(2:end,:)')
load ./sync_times.mat
synct_d = downsample(synct,8);
counts = cellfun(@(x) length(x),mua_amps)/range(synct)*1e6
figure
plot(2:8,counts(2:8),'o-')
 
%%
%           load ./used_data
%         clear S
%         params.Fs = 2016;
%         params.tapers = [2 3];
%         params.fpass = [0 40];
%         win = 50;
%         for i = 2:8
%             eval(['[S(i,:),f]=mtspectrumsegc(' sprintf('lf%d',i) ',win,params,1);']);
%         end
%         uds_freqs = find(f > 0.1 & f < 1);
%         peak_uds_pow = max(S(:,uds_freqs),[],2);
% plot(peak_uds_pow)
 
%%
close all
raw_Fs = 2016;
dsf = 4;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.1;
% lcf_hf = 15;
% hcf_hf = 80;
% hcf_sm = 0.025;

to_zscore = 1;
eps = 3;

% hpc_ch = 3;
% synct_d = downsample(synct,dsf);
% mua_binned = hist(mua_times{hpc_ch},synct_d);
% mua_binned([1 end]) = 0;
% mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.06));
% hpc_mua_rate = zscore(mua_rate);

load ./used_data lf8 lf5 lf7 lf6 lf4 lf3 lf2 wcv_minus_spike
[lf8_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf],to_zscore);
[lf6_lf,t_axis] = get_lf_features(lf6,raw_Fs,Fsd,[lcf hcf],to_zscore);
[lf5_lf,t_axis] = get_lf_features(lf5,raw_Fs,Fsd,[lcf hcf],to_zscore);
[lf4_lf,t_axis] = get_lf_features(lf4,raw_Fs,Fsd,[lcf hcf],to_zscore);
[lf3_lf,t_axis] = get_lf_features(lf3,raw_Fs,Fsd,[lcf hcf],to_zscore);
[lf2_lf,t_axis] = get_lf_features(lf2,raw_Fs,Fsd,[lcf hcf],to_zscore);
wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);

lcf = 40;
hcf = 200;
offset = -5*eps;
% eps = 1e-4;
% [lf8_hf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf],to_zscore);
% [lf6_hf,t_axis] = get_lf_features(lf6,raw_Fs,Fsd,[lcf hcf],to_zscore);
% [lf5_hf,t_axis] = get_lf_features(lf5,raw_Fs,Fsd,[lcf hcf],to_zscore);
[lf4_hf,t_axis] = get_lf_features(lf4,raw_Fs,Fsd,[lcf hcf],to_zscore);
[lf3_hf,t_axis] = get_lf_features(lf3,raw_Fs,Fsd,[lcf hcf],to_zscore);
[lf2_hf,t_axis] = get_lf_features(lf2,raw_Fs,Fsd,[lcf hcf],to_zscore);

window_size = 20;
window_size_n = round(Fsd*window_size);
n_windows = ceil(max(t_axis)/window_size);
for i = 1:n_windows
    cur_set = (i-1)*window_size_n + (1:window_size_n);
    cur_set(cur_set > length(t_axis)) = [];
% subplot(2,1,1)
plot(t_axis(cur_set),wcv_lf(cur_set)+2*eps)
hold on
plot(t_axis(cur_set),lf8_lf(cur_set)+eps,'r','linewidth',1)
% plot(t_axis(cur_set),lf6_lf(cur_set),'y','linewidth',1)
% plot(t_axis(cur_set),lf5_lf(cur_set)-0*eps,'g')
% plot(t_axis(cur_set),lf4_lf(cur_set)-2*eps,'c')
plot(t_axis(cur_set),lf3_lf(cur_set)-0*eps,'c')
plot(t_axis(cur_set),lf2_lf(cur_set)-1*eps,'k')
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude (V)','fontsize',16)
% subplot(2,1,2)
% % plot(t_axis(cur_set),lf4_hf(cur_set),'c')
% hold on
% plot(t_axis(cur_set),lf3_hf(cur_set)-eps,'c')
% plot(t_axis(cur_set),lf2_hf(cur_set)-2*eps,'k')
% plot(t_axis(cur_set),hpc_mua_rate(cur_set)-eps,'r','linewidth',2)
pause
clf
end
    %%
close all
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 4;
lcf = 0.1;
to_zscore = 1;
eps = 3;
    load ./all_combined_mp_uds.mat
    load ./all_combined_lf7_uds.mat
        lfp_state_seq = hmm_bbstate_seq7;
        mp_state_seq = hmm_bbstate_seq;        
        [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,Fsd,mp_state_seq);
        [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
        [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);

load ./used_data lf8 lf5 lf7 lf6 lf4 lf3 lf2 wcv_minus_spike
[lf7_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf],to_zscore);
[lf3_lf,t_axis] = get_lf_features(lf3,raw_Fs,Fsd,[lcf hcf],to_zscore);
[lf2_lf,t_axis] = get_lf_features(lf2,raw_Fs,Fsd,[lcf hcf],to_zscore);
wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);

wcv_up_log = nan(size(t_axis));
lfp_up_log = nan(size(t_axis));

for ns = 1:hmm.Nsegs
    cur_seg = new_seg_inds(ns,1):new_seg_inds(ns,2);
    wcv_up_log(cur_seg) = logical(mp_state_seq{ns}-1);
    lfp_up_log(cur_seg) = logical(lfp_state_seq{ns}-1);
end

window_size = 30;
window_size_n = round(Fsd*window_size);
n_windows = ceil(max(t_axis)/window_size);
for i = 1:n_windows
    cur_set = (i-1)*window_size_n + (1:window_size_n);
    cur_set(cur_set > length(t_axis)) = [];
plot(t_axis(cur_set),wcv_lf(cur_set)+eps)
hold on
plot(t_axis(cur_set),lf7_lf(cur_set),'r','linewidth',1)
plot(t_axis(cur_set),lf3_lf(cur_set)-eps,'c','linewidth',1)
plot(t_axis(cur_set),lf2_lf(cur_set)-2*eps,'k')
plot(t_axis(cur_set),wcv_up_log(cur_set)+eps,'k')
plot(t_axis(cur_set),lfp_up_log(cur_set),'k')
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude (V)','fontsize',16)

pause
clf
end

%%

% plot(t_axis,lf8_hf+eps,'r','linewidth',1)
hold on
% plot(t_axis,lf6_hf,'y','linewidth',1)
% plot(t_axis,wcv_hf+3)
plot(t_axis,2*lf5_hf-eps+offset,'g')
plot(t_axis,2*lf4_hf-2*eps+offset,'c')
plot(t_axis,2*lf3_hf-3*eps+offset,'b')
plot(t_axis,2*lf2_hf-4*eps+offset,'k')
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude (V)','fontsize',16)
%%
ctx_ch = 7;
hpc_ch = 4;
hpc_ch2 = 3;
synct_d = downsample(synct,dsf);

mua_binned = hist(mua_times{ctx_ch},synct_d);
mua_binned([1 end]) = 0;
mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.06));
ctx_mua_rate = zscore(mua_rate);

mua_binned = hist(mua_times{hpc_ch},synct_d);
mua_binned([1 end]) = 0;
mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.06));
hpc_mua_rate = zscore(mua_rate);

mua_binned = hist(mua_times{hpc_ch2},synct_d);
mua_binned([1 end]) = 0;
mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.06));
hpc_mua_rate2 = zscore(mua_rate);

hold on
plot(t_axis,ctx_mua_rate+4,'k')
plot(t_axis,hpc_mua_rate+4,'g')
plot(t_axis,hpc_mua_rate2+4,'c')
shg

%%
maxlag = round(Fsd*2);
for ch = 2:8
    mua_binned = hist(mua_times{ch},synct_d);
    mua_binned([1 end]) = 0;
    mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.05));
    [mua_lfp_xc(ch,:),lags] = xcov(lf8_lf,mua_rate,maxlag,'coeff');
end