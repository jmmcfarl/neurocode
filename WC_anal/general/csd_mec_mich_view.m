clear all
% close all

% cd C:\WC_Germany\2007-11-02_Par_EC_LFP_A\2007-11-2_13-13-6
cd C:\WC_Germany\2007-11-02_Par_EC_LFP_B\2007-11-2_14-36-23

load C:\WC_Germany\JMM_analysis_pyr\sim_record_UDS_dur_2

load used_data lf15 lf16 lf13 lf8 lf3

dsf = 8;
winSize = 15;

Fsd = 2016/dsf;
niqf = 2016/2;

params.Fs = Fsd;
params.tapers = [2 3];
params.fpass = [0 10];
params.err = [2 .05];
movingwin = [25 1];

maxlag = 2*Fsd;
lags = -maxlag:maxlag;

% [b,a] = butter(2,[0.1/niqf 15/niqf]);
% [b2,a2] = butter(2,[5/niqf 60/niqf]);
[b2,a2] = butter(2,[0.1/niqf 15/niqf]);
lf8_f = filtfilt(b2,a2,lf8);
lf3_f = filtfilt(b2,a2,lf3);
lf16_f = filtfilt(b2,a2,lf16);
% lf16_hf = filtfilt(b2,a2,lf16);
lf15_f = filtfilt(b2,a2,lf15);
lf13_f = filtfilt(b2,a2,lf13);
csd = lf16_f + lf13_f - 2*lf15_f;

lf8_d = downsample(lf8_f,dsf);
lf3_d = downsample(lf3_f,dsf);
lf16_d = downsample(lf16_f,dsf);
lf8_d = zscore(lf8_d);
lf3_d = zscore(lf3_d);
lf16_d = zscore(lf16_d);
csd_d = downsample(csd,dsf);
csd_d = zscore(csd_d);

csd_response_mat = zeros(length(up_trans8),length(lags));
for i = 1:length(up_trans8)
    
    if up_trans8(i) > maxlag & (length(csd_d)-up_trans8(i)) > maxlag
        csd_response_mat(i,:) = csd_d(up_trans8(i)-maxlag:up_trans8(i)+maxlag);
    else
        csd_response_mat(i,:) = nan;
    end
  
end

% [C16,phi,S12,S1,S2,Ct,Cf]=cohgramc(lf8_d,lf16_d,movingwin,params);
% 
% [C8,phi,S12,S1,S2,Ct,Cf]=cohgramc(lf8_d,csd_d,movingwin,params);
% [C3,phi,S12,S1,S2,Ct,Cf]=cohgramc(lf3_d,csd_d,movingwin,params);

t = (1:length(lf8_d))/Fsd;

% numWins = floor(max(t)/winSize);
% for i = 1:numWins
%     seg_beg = (i-1)*winSize*Fsd+1;
%     seg_end = i*winSize*Fsd;
%     cur_lf8_ups = up_trans8(find(up_trans8 > seg_beg & up_trans8 < seg_end));
%     cur_lf8_downs = down_trans8(find(down_trans8 > seg_beg & down_trans8 < seg_end));
%     
%     plot(t(seg_beg:seg_end),-csd_d(seg_beg:seg_end),'k','linewidth',2)
% hold on
% plot(t(seg_beg:seg_end),lf8_d(seg_beg:seg_end),'r','linewidth',2)
% plot(t(cur_lf8_ups),lf8_d(cur_lf8_ups),'ro')
% plot(t(cur_lf8_downs),lf8_d(cur_lf8_downs),'go')
% % plot(t(seg_beg:seg_end),lf16_d(seg_beg:seg_end),'k','linewidth',2)
% ylim([-4 4])
% pause
% clf
% end
