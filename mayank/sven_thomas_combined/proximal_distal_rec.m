clear all
close all
cd C:\wc_data\2012_07_14\2012-7-14_20-24-51

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
lcf = 0.05;
hcf = 10;
params.Fs = raw_Fs;
params.tapers = [2 3];
win = 50;

load ./used_data 

%%
[S8,f] = mtspectrumsegc(lf8,win,params);
[S4,f] = mtspectrumsegc(lf4,win,params);
[S6,f] = mtspectrumsegc(lf6,win,params);
[S2,f] = mtspectrumsegc(lf2,win,params);
[S5,f] = mtspectrumsegc(lf5,win,params);
[S3,f] = mtspectrumsegc(lf3,win,params);
[S7,f] = mtspectrumsegc(lf7,win,params);
[S1,f] = mtspectrumsegc(wcv,win,params);

%%
cmap = jet(4);
figure
subplot(2,1,1)
plot(f,S8,'color',cmap(4,:))
xlim([0 2])
ylim([0 1e-7])
hold on
plot(f,S7,'color',cmap(3,:))
plot(f,S6,'color',cmap(2,:))
plot(f,S5,'color',cmap(1,:))
subplot(2,1,2)
plot(f,S4,'color',cmap(4,:))
hold on
xlim([0 2])
ylim([0 1e-7])
plot(f,S3,'color',cmap(3,:))
plot(f,S2,'color',cmap(2,:))
plot(f,S1,'color',cmap(1,:))

figure
subplot(2,1,1)
plot(f,S8,'color',cmap(4,:))
set(gca,'xscale','log','yscale','log')
hold on
plot(f,S7,'color',cmap(3,:))
plot(f,S6,'color',cmap(2,:))
plot(f,S5,'color',cmap(1,:))
subplot(2,1,2)
plot(f,S4,'color',cmap(4,:))
hold on
set(gca,'xscale','log','yscale','log')
plot(f,S3,'color',cmap(3,:))
plot(f,S2,'color',cmap(2,:))
plot(f,S1,'color',cmap(1,:))

%%
[lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf],0);
lf7_lf = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf],0);
lf6_lf = get_lf_features(lf6,raw_Fs,Fsd,[lcf hcf],0);
lf5_lf = get_lf_features(lf5,raw_Fs,Fsd,[lcf hcf],0);
lf4_lf = get_lf_features(lf4,raw_Fs,Fsd,[lcf hcf],0);
lf3_lf = get_lf_features(lf3,raw_Fs,Fsd,[lcf hcf],0);
lf2_lf = get_lf_features(lf2,raw_Fs,Fsd,[lcf hcf],0);
lf1_lf = get_lf_features(wcv,raw_Fs,Fsd,[lcf hcf],0);

%%
eps = 4e-4;
figure
plot(t_axis,lf8_lf+4*eps,'r')
hold on
plot(t_axis,lf7_lf+3*eps,'k')
plot(t_axis,lf6_lf+2*eps,'g')
plot(t_axis,lf5_lf+1*eps,'b')
xlim([228 248])
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude (V)','fontsize',16)

figure
plot(t_axis,lf7_lf+5*eps,'k-','linewidth',2)
hold on
plot(t_axis,lf4_lf+4*eps,'r')
plot(t_axis,lf3_lf+3*eps,'k')
plot(t_axis,lf2_lf+2*eps,'g')
plot(t_axis,lf1_lf+1*eps,'b')
xlim([228 248])
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude (V)','fontsize',16)

%%
amp_threshold = 25;
max_overlap = 0.5;
[mua_times,mua_amps,mua_widths,avg_waveform,std_waveform] = extract_MiP_mua(amp_threshold,max_overlap);

figure
plot(avg_waveform(1:end,:)')
load ./sync_times.mat
synct_d = downsample(synct,8);
counts = cellfun(@(x) length(x),mua_amps)/range(synct)*1e6
figure
plot(1:8,counts(1:8),'o-')

%%
plot(t_axis,zscore(lf7_lf),'r')
hold on
plot(t_axis,zscore(lf4_lf),'b')

ctx_ch1 = 8;
ctx_ch2 = 4;

mua_binned = hist(mua_times{ctx_ch1},synct_d);
mua_binned([1 end]) = 0;
mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.05));
ctx_mua_rate = zscore(mua_rate);

mua_binned = hist(mua_times{ctx_ch2},synct_d);
mua_binned([1 end]) = 0;
mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.05));
ctx_mua_rate2 = zscore(mua_rate);

hold on
plot(t_axis,ctx_mua_rate+4,'k')
plot(t_axis,ctx_mua_rate2+4,'r')
shg

hpc_ch1 = 6;
hpc_ch2 = 2;

mua_binned = hist(mua_times{hpc_ch1},synct_d);
mua_binned([1 end]) = 0;
mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.05));
hpc_mua_rate1 = zscore(mua_rate);

mua_binned = hist(mua_times{hpc_ch2},synct_d);
mua_binned([1 end]) = 0;
mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.05));
hpc_mua_rate2 = zscore(mua_rate);

hold on
plot(t_axis,hpc_mua_rate1-4,'k')
plot(t_axis,hpc_mua_rate2-4,'r')
shg

