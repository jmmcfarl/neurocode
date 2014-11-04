clear all
close all
addpath('C:\Code\general_functions\');

cd C:\WC_Germany\sven_thomas_combined
% load ./combined_dir.mat
% l3mec_m = find(~isnan(hpc_mua));
% d = 56
% cd(combined_dir{d})

% load ./distal_dir
% d=3
% cd(distal_dir{d})

cd C:\WC_Germany\persistent_downs\
load ./new_pdown_dir.mat
d = 6
cd(new_pdown_dir{d})



%%
amp_threshold = 30;
max_overlap = 0.5;
% [mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);
load ./mua_data3

figure
plot(avg_waveform(2:end-1,:)')
load ./sync_times.mat
synct_d = downsample(synct,8);
counts = cellfun(@(x) length(x),mua_amps)/range(synct)*1e6
figure
plot(1:7,counts(1:7),'o-')
 
        %%
close all
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 4;
lcf = 0.05;
to_zscore = 1;
load ./all_combined_mp_uds.mat
load ./all_combined_lf7_uds.mat
lfp_state_seq = hmm_bbstate_seq7;
mp_state_seq = hmm_bbstate_seq;
[new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,Fsd,mp_state_seq);
[up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
[up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);

load ./used_data lf7 wcv_minus_spike
[lf7_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf],to_zscore);
wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);

load ./aligned_heka.mat
ac_time_d = downsample(ac_time,dsf);
dc_fs = 1/median(diff(dc_time));
dc_dsf = 20;
dc_fsd = dc_fs/dc_dsf;
dc_data = decimate(dc_data,dc_dsf);
dc_time = downsample(dc_time,dc_dsf);

wcv_up_log = nan(size(t_axis));
lfp_up_log = nan(size(t_axis));

for ns = 1:hmm.Nsegs
    cur_seg = new_seg_inds(ns,1):new_seg_inds(ns,2);
    wcv_up_log(cur_seg) = logical(mp_state_seq{ns}-1);
    lfp_up_log(cur_seg) = logical(lfp_state_seq{ns}-1);
end

rate_sm = round(Fsd*0.05);

load ./mua_data3
load ./sync_times.mat
synct_d = downsample(synct,dsf);
uset = find(mua_amps{ctx_ch} > 50);
if ~isnan(new_pdown_hpcmua(d))
    hpc_mua_times = mua_times{new_pdown_hpcmua(d)};
    hpc_mua_times(hpc_mua_times > synct_d(end) | hpc_mua_times < synct_d(1)) = [];
    hpc_mua_rate =hist(hpc_mua_times,synct_d)*Fsd;
    hpc_mua_rate = jmm_smooth_1d_cor(hpc_mua_rate,rate_sm);
else
    hpc_mua_rate = nan(size(synct_d));
end
ctx_mua_times = sort([mua_times{6} mua_times{7} mua_times{8}]);
ctx_mua_rate = hist(ctx_mua_times,synct_d)*Fsd;
ctx_mua_rate = jmm_smooth_1d_cor(ctx_mua_rate,rate_sm);
if length(hpc_mua_rate) > length(synct_d)
    hpc_mua_rate = hpc_mua_rate(1:length(synct_d));
    ctx_mua_rate = ctx_mua_rate(1:length(synct_d));
end

hpc_mua_rate = zscore(hpc_mua_rate);
ctx_mua_rate = zscore(ctx_mua_rate);

% load ./used_data lf3
% lcf = 35;
% hcf = 80;
% [lf3_hf,t_axis] = get_lf_features(lf3,raw_Fs,Fsd,[lcf hcf],to_zscore);
% lf3_rpow  = zscore(abs(hilbert(lf3_hf)));
% lcf = 80;
% hcf = 200;
% [lf3_hf,t_axis] = get_lf_features(lf3,raw_Fs,Fsd,[lcf hcf],to_zscore);
% lf3_rpow2  = zscore(abs(hilbert(lf3_hf)));

%%
window_size = 30;
window_size_n = round(Fsd*window_size);
window_size_dc = round(dc_fsd*window_size);
n_windows = ceil(max(t_axis)/window_size);
for i = 1:n_windows
    cur_set = (i-1)*window_size_n + (1:window_size_n);
    cur_set(cur_set > length(t_axis)) = [];
    cur_set_dc = find(dc_time >= ac_time_d(cur_set(1)) & dc_time < ac_time_d(cur_set(end)));
    subplot(3,1,1)
    plot(dc_time(cur_set_dc),dc_data(cur_set_dc));hold on
    plot(ac_time_d(cur_set),wcv_up_log(cur_set)*20-50,'k')
    
    subplot(3,1,2)
%     plot(ac_time_d(cur_set),wcv_lf(cur_set)+5)
    hold on
    plot(ac_time_d(cur_set),lf7_lf(cur_set),'r','linewidth',1)
%     plot(ac_time_d(cur_set),wcv_up_log(cur_set)+5,'k')
    plot(ac_time_d(cur_set),lfp_up_log(cur_set),'k')
    
    subplot(3,1,3)
    hold on
    plot(ac_time_d(cur_set),ctx_mua_rate(cur_set)-5,'r')
    plot(ac_time_d(cur_set),hpc_mua_rate(cur_set)-8,'k')
%     plot(ac_time_d(cur_set),jmm_smooth_1d_cor(lf3_rpow2(cur_set),20)*2-5,'g')
    xlabel('Time (s)','fontsize',16)
    ylabel('Amplitude (V)','fontsize',16)
 
    pause
    clf
end

%%
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