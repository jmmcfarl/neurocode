clear all
distal_lec{1} = 'C:\wc_data\2012_06_21_A\2012-6-21_14-33-23';
% distal_lec{2} = 'C:\wc_data\2012_06_21_B\2012-6-21_16-1-44'; %LEC neuron, good LFP, ?MUA, but cell wierd (depolarized MP)
% distal_lec{3} = 'C:\wc_data\2012_06_23_A\2012-6-23_18-39-42'; %not usable LFP (MAYBE with LF5) no DEPOL
% distal_lec{4} = 'C:\wc_data\2012_06_23_B\2012-6-23_19-5-21'; %no DEPOL
distal_lec{5} = 'C:\wc_data\2012_06_23_C\2012-6-23_19-39-37'; %usable LF5
% distal_lec{6}  = 'C:\wc_data\2012_06_24_A\2012-6-24_14-49-46'; %questionable LFP
% distal_lec{7} = 'C:\wc_data\2012_06_24_B\2012-6-24_15-46-53'; %good MUA, questionable LFP (at times), NEED to use LF4
% distal_lec{8} = 'C:\wc_data\2012_06_24_C\2012-6-24_16-24-8'; %maybe OK MUA, questionabel LFP, (again, need LF4)
% distal_lec{9} = 'C:\wc_data\2012_06_25_A\2012-6-25_12-58-45'; %not usable LFP

% load ./distal_dir.mat

% for d = 1:14
d=5
cd(distal_lec{d})
% cd(distal_dir{d})
pwd


%%
amp_threshold = 25;
max_overlap = 0.5;
[mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);
% 
load ./sync_times.mat
synct_d = downsample(synct,8);
rate(d,:) = cellfun(@(x) length(x),mua_amps)/range(synct)*1e6;
figure
plot(2:8,rate(d,2:8),'o-')
 figure
plot(avg_waveform(2:end,:)')

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
        peak_uds_pow(d,:) = max(S(:,uds_freqs),[],2);
        peak_uds_pow(d,1) = nan;
% end
        figure
        plot(2:8,peak_uds_pow(d,2:8),'o-')
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
load ./used_data lf8 lf5 lf4 lf7 lf6 wcv_minus_spike
[lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
[lf5_lf,t_axis] = get_lf_features(lf5,raw_Fs,Fsd,[lcf hcf]);
wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
% plot(t_axis,lf8_lf,'r')
hold on
plot(t_axis,wcv_lf+3)
hold on
plot(t_axis,lf8_lf,'r')
plot(t_axis,lf5_lf,'k')

%%
ctx_ch = 8;
hpc_ch = 4;

mua_binned = hist(mua_times{ctx_ch},synct_d);
mua_binned([1 end]) = 0;
mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.05));
ctx_mua_rate = zscore(mua_rate);

mua_binned = hist(mua_times{hpc_ch},synct_d);
mua_binned([1 end]) = 0;
mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.05));
hpc_mua_rate = zscore(mua_rate);

hold on
plot(t_axis,ctx_mua_rate+4,'k')
plot(t_axis,hpc_mua_rate+4,'r')
shg

%%
maxlag = round(Fsd*0.5);
for ch = 2:8
    mua_binned = hist(mua_times{ch},synct_d);
    mua_binned([1 end]) = 0;
    mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.05));
    [mua_lfp_xc(ch,:),lags] = xcov(lf8_lf,mua_rate,maxlag,'coeff');
end