%% data viewer

clear all
close all

%% for hc_cort
% load('C:\WC_Germany\overall_calcs\HC_Cort\hc_cor_dir.mat')
% load C:\WC_Germany\overall_calcs\HC_Cort\UDS_dur_raw\UDS_raw_data
% load C:\WC_Germany\overall_calcs\HC_Cort\UDS_synch_state_dur\UDS_synch_state_dur_data

%% for EC
load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data
load C:\WC_Germany\overall_calcs\overall_info_file_data
% 
%%

d =62

cd(over_dir{d})
pwd

winSize = 20;

Fs = 2016;
dsf = 8;
niqf = Fs/2;
Fsd = Fs/dsf;
lcf = 5/niqf;
hcf = 40/niqf;
lcf2 = 0.05/niqf;
hcf2 = 2/niqf;
[b2,a2] = butter(2,150/niqf,'high');
[b3,a3] = butter(2,[4/niqf 20/niqf]);
[b4,a4] = butter(1,5/Fsd*2,'high');
[b,a] = butter(2,[0.05/niqf 5/niqf]);


load used_data wcv_minus_spike lf3 lf8 lf5 lf7

wcv_f = filter(b3,a3,wcv_minus_spike);
wcv_f = downsample(wcv_f,dsf);
wcv_f = zscore(wcv_f);
% wcv_f = jmm_smooth_1d(wcv_f.^2,5);
% rmsval = sqrt(mean(wcv_f.^2));
% wcv_f(wcv_f <= 0) = rmsval*0.01;
% wcv_f = log(wcv_f);
% wcv_f = zscore(wcv_f);
% wcv_der = jmm_smooth_1d([0;diff(wcv_minus_spike)],4);
% wcv_orf = filtfilt(b,a,wcv_minus_spike);
wcv_d = downsample(wcv_minus_spike,dsf);
% wcv_der = downsample(wcv_der,dsf);
wcv_d = zscore(wcv_d);
% wcv_der = zscore(wcv_der);
% 
lf8_f = filter(b,a,lf8);
lf8_d = downsample(lf8_f,dsf);
lf8_d = zscore(lf8_d);
% lf8_f2 = filter(b2,a2,lf8);
% 
% lf8_d2 = zscore(lf8_f2);
% lf8_f3 = filter(b2,a2,lf7);
% 
% lf8_d3 = zscore(lf8_f3);
% lf8_f4 = filter(b2,a2,lf5);
% 
% lf8_d4 = zscore(lf8_f4);

% lf3_f = filtfilt(b3,a3,lf3);
% lf3_d = downsample(lf3_f,dsf);
% lf3_d = zscore(lf3_d);

t = (1:length(wcv_d))/Fsd;

numWins = floor(max(t)/winSize);


figure
% plot(t,lf8_d2/8,'r')
hold on
% plot(t,lf8_d3/4-2,'g')
% plot(t,lf8_d4/8+2,'m')
hold on
plot(t,wcv_d)
% hold on
% plot(t,wcv_der,'r')
hold on
plot(t,lf8_d,'k','linewidth',2)

% plot(t,wcv_f,'k')
% plot(t(cross_points),wcv_d(cross_points),'ro')
% plot(t,wcv_d3,'g')
% plot(t,wcv_d2,'k')
% plot(t,wcv_d2,'r')
% plot(t(spike_ids),ones(size(spike_ids))*2,'k.')
plot(t(up_trans{d}),wcv_d(up_trans{d}),'ro')
plot(t(down_trans{d}),wcv_d(down_trans{d}),'go')
% plot(t(up_trans8{d}),lf8_d(up_trans8{d}),'ro')
% plot(t(down_trans8{d}),lf8_d(down_trans8{d}),'ro')

for i = 1:numWins
    xlim([(i-1)*winSize i*winSize])
%     line([(i-1)*winSize i*winSize],[-10.5 -10.5],'Color','k')
    ylim([-3 3])
    pause
end
% 

