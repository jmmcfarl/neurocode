clear all
close all
cd C:\WC_Germany\sven_thomas_combined
load ./distal_dir

d = 1; %U LF8. Maybe can use Ch2, but no clear peak
d = 2; %LF5, no hpc peak, good mua signals only at 50uv
d = 3; %Not usable LFP (L4, but NOT others) peak on ch 3, wierd signals even with higher threshold, but looks maybe OK in terms of UDS modulation
% d = 4; %LF8 (maybe) U? Peak on Ch2 NOPE
% d = 5; %U LF8 peak on ch 2
% d = 6; %no LFPs have usable UDS;
d = 7; %no usable LFPs noisy sigs, even up to 75 uv, probably not usable.
% d = 8; % LF5, no hpc peak. noisy sigs, but clear uds mod.
% d = 9; %U (M2, LF5) NOPE
% d = 10; %no usable LFPs good mua sigs
% d = 11; %U MAYBE LF5 (M3) NOPE
% d = 12; %? 
% d = 13; %noisy, MAYBE peak on ch 4, after increasing to 50uv. Pretty reasonable UDS mod though.
% d = 14; %noisey, but descent uds mod on some ch at least. peaks on 3 and 5, but no delay on 3

cd(distal_dir{d})
pwd

% cd C:\wc_data\2012_06_21_B\2012-6-21_16-1-44 %LEC neuron, good LFP, ?MUA, but cell wierd (depolarized MP)
% % cd C:\wc_data\2012_06_23_A\2012-6-23_18-39-42 %not usable LFP
% cd C:\wc_data\2012_06_24_B\2012-6-24_15-46-53 %good MUA, questionable LFP (at times), NEED to use LF4
% % cd C:\wc_data\2012_06_24_C\2012-6-24_16-24-8 %maybe OK MUA, questionabel LFP, (again, need LF4)
% cd C:\wc_data\2012_06_25_A\2012-6-25_12-58-45 %not usable LFP

%%
amp_threshold = 25;
max_overlap = 0.5;
[mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);

figure
plot(avg_waveform(2:end,:)')
load ./sync_times.mat
synct_d = downsample(synct,8);
counts = cellfun(@(x) length(x),mua_amps)/range(synct)*1e6
figure
plot(2:8,counts(2:8),'o-')
 
%%
          load ./used_data
        clear S
        params.Fs = 2016;
        params.tapers = [2 3];
        params.fpass = [0 40];
        win = 50;
        for i = 2:8
            eval(['[S(i,:),f]=mtspectrumsegc(' sprintf('lf%d',i) ',win,params,1);']);
        end
        uds_freqs = find(f > 0.1 & f < 1);
        peak_uds_pow = max(S(:,uds_freqs),[],2);
plot(peak_uds_pow)
 
%%
close all
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;
lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.025;
load ./used_data lf8 lf5 lf7 lf6 lf4 wcv_minus_spike
% if ctx_lfp(d) == 5
%     lf8 = lf5;
% end
[lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
% [lf5_lf,t_axis] = get_lf_features(lf5,raw_Fs,Fsd,[lcf hcf]);
wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
% plot(t_axis,lf8_lf,'r')
hold on
plot(t_axis,wcv_lf+3)
hold on
% plot(t_axis,lf5_lf,'k')
plot(t_axis,lf8_lf,'r')

%     load ./pa_hsmm_state_seq_combined_fin_nd.mat
%     load ./pa_hsmm_state_seq7_combined_fin_nd.mat
%     hsmm_bbstate_seq8 = hsmm_bbstate_seq7;
%     hsmm8 = hsmm7;
%     [new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs, Fsd,hsmm_bbstate_seq);
%     dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
%     [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq);
%     [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq8);

%     plot(t_axis(up_trans_inds8),lf8_lf(up_trans_inds8),'ro')
%%
ctx_ch = 7;
hpc_ch = 2;
hpc_ch2 = 3;

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
% plot(t_axis,hpc_mua_rate2+4,'c')
shg

%%
maxlag = round(Fsd*2);
for ch = 2:8
    mua_binned = hist(mua_times{ch},synct_d);
    mua_binned([1 end]) = 0;
    mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.05));
    [mua_lfp_xc(ch,:),lags] = xcov(lf8_lf,mua_rate,maxlag,'coeff');
end