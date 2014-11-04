clear all

FieldSelectionFlags = [1 0 1 0 1];
ExtractMode = 1;
HeaderExtractionFlag = 1;


Filename = 'CSC1.Ncs';
[CSC1_Timestamps, CSC1_SampleFrequencies,CSC1_Samples, Header] = ...
    Nlx2MatCSC_v3( Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
CSC1_Timestamps = interp1(1:512:numel(CSC1_Samples),CSC1_Timestamps,1:numel(CSC1_Samples));
bad = find(isnan(CSC1_Timestamps));
CSC1_Timestamps(bad) = [];
Fs1 = median(CSC1_SampleFrequencies);
csc1 = CSC1_Samples(:)*str2num(Header{15}(13:end));
csc1(bad) = [];
clear CSC1_Samples

% Filename = 'CSC2.Ncs';
% [CSC2_Timestamps, CSC2_SampleFrequencies,CSC2_Samples, Header] = ...
%     Nlx2MatCSC_v3( Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
% CSC2_Timestamps = interp1(1:512:numel(CSC2_Samples),CSC2_Timestamps,1:numel(CSC2_Samples));
% bad = find(isnan(CSC2_Timestamps));
% CSC2_Timestamps(bad) = [];
% Fs2 = median(CSC2_SampleFrequencies);
% csc2 = CSC2_Samples(:);
% csc2(bad) = [];
% clear CSC2_Samples

% Filename = 'CSC5.Ncs';
% [CSC5_Timestamps, CSC5_SampleFrequencies,CSC5_Samples, Header] = ...
%     Nlx2MatCSC_v3( Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
% CSC5_Timestamps = interp1(1:512:numel(CSC5_Samples),CSC5_Timestamps,1:numel(CSC5_Samples));
% bad = find(isnan(CSC5_Timestamps));
% CSC5_Timestamps(bad) = [];
% Fs5 = median(CSC5_SampleFrequencies);
% csc5 = CSC5_Samples(:);
% csc5(bad) = [];
% clear CSC5_Samples

Filename = 'CSC4.Ncs';
[CSC8_Timestamps, CSC8_SampleFrequencies,CSC8_Samples, Header] = ...
    Nlx2MatCSC_v3( Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
CSC8_Timestamps = interp1(1:512:numel(CSC8_Samples),CSC8_Timestamps,1:numel(CSC8_Samples));
bad = find(isnan(CSC8_Timestamps));
CSC8_Timestamps(bad) = [];
Fs8 = median(CSC8_SampleFrequencies);
csc8 = CSC8_Samples(:)*str2num(Header{15}(13:end));
csc8(bad) = [];
clear CSC8_Samples

%%
Fs = Fs8;
dsf = 8;
% Fsd = Fs5/dsf;
% csc1_d = zscore(decimate(csc1,dsf));
% csc2_d = zscore(decimate(csc5,dsf));
% csc8_d = zscore(decimate(csc8,dsf));
% 
t8_d = downsample(CSC8_Timestamps,dsf);
% t1_d = downsample(CSC1_Timestamps,dsf*16);
% interp_t = linspace(t8_d(1),t8_d(end),length(t8_d));
% csc1_d = -interp1(t1_d',csc1_d,interp_t);
% csc2_d = interp1(t8_d',csc2_d,interp_t);
% csc8_d = interp1(t8_d',csc8_d,interp_t);
% bad = find(isnan(csc1_d));
% csc1_d(bad) = 0;
% 
% t = (1:length(csc8_d))/Fsd;
dsf = 8;
Fsd = Fs/dsf;
[b,a] = butter(2,[0.2 10]/(Fs/2));
% [b,a] = butter(2,[0.25 10]/(Fs/2));
csc1_d = downsample(filtfilt(b,a,csc1),dsf);
% % csc2_d = downsample(filtfilt(b,a,csc5),dsf);
csc8_d = downsample(filtfilt(b,a,csc8),dsf);
t = (1:length(csc8_d))/Fsd;

figure
plot(t,zscore(csc1_d))
hold on
plot(t,zscore(csc8_d)-6,'k')
xlabel('Time (s)','fontsize',16)
ylabel('Amplitude (z)','fontsize',16)

csc1_d = zscore(csc1_d);
csc8_d = zscore(csc8_d);

% [x,lags] = xcov(csc1_d,csc8_d,round(Fsd*1),'coeff');
% figure
% plot(lags/Fsd,x)
%%
% params.Fs = Fsd;
% params.tapers = [2 3];
% win = [10 5];
% [S,tt,ff] = mtspecgramc(csc8_d,win,params);
% pcolor(tt,ff,log(S'));shading flat
% caxis([-10 5])
% set(gca,'yscale','log')
% 
% win = [20 10];
% params.tapers = [5 9];
% [S,phi,S12,S1,S2,tt,ff] = cohgramc(csc1_d,csc8_d,win,params);
% pcolor(tt,ff,S');shading flat
% set(gca,'yscale','log')
% 
win = [20 10];
params.Fs = Fsd;
params.tapers = [5 9];
[S,phi,S12,S1,S2,ff] = coherencysegc(csc1_d,csc8_d,win,params);
plot(ff,S)
set(gca,'xscale','log')

figure
plot(ff,log(S1),ff,log(S2),'r')
set(gca,'xscale','log')

%%
% [b,a] = butter(2,[0.1 40]/(Fsd/2));
% csc1_d = zscore(filtfilt(b,a,csc1_d));
% % csc2_d = zscore(filtfilt(b,a,csc2_d));
% csc8_d = zscore(filtfilt(b,a,csc8_d));

plot(t,csc1_d*1e3)
hold on
plot(t,csc8_d*1e3,'k')
% plot(t,csc2_d-3,'r')

%%
% Fs2 = Fs5;
% [b,a] = butter(2,[20 80]/(Fs2/2));
% csc1_hf = filtfilt(b,a,csc1);
% csc8_hf = filtfilt(b,a,csc5);
% csc1_hf = smooth(abs(csc1_hf),round(Fs2*0.15));
% csc8_hf = smooth(abs(csc8_hf),round(Fs2*0.15));
% 
% csc1_hf = zscore(decimate(csc1_hf,dsf));
% csc8_hf = zscore(decimate(csc8_hf,dsf));
% csc1_hf = interp1(t8_d',csc1_hf,interp_t);
% csc8_hf = interp1(t8_d',csc8_hf,interp_t);
% plot(t,csc1_hf,'r')
% plot(t,csc8_hf,'g')
% 
%%
% clear all
FieldSelectionFlags = [1 0 0 1 0];
HeaderExtractionFlag = 1;
ExtractMode = 1;
Filename = 'Sc1.ntt';
[Timestamps, Samples, Header] = ...
    Nlx2MatSpike_v3( Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
% peakvals = squeeze(max(Samples)); 
conv_factor = str2num(Header{15}(13:end));
peakvals = bsxfun(@times,Samples(1:4,:),conv_factor'); 
clear Samples
[ov_peak,peak_ch] = max(peakvals);
for i = 1:4
    cur_mua = find(peak_ch == i);
%     cur_mua(ov_peak(cur_mua) < 20/1e6) = [];
    mua{i} = Timestamps(cur_mua);
end
Filename = 'Sc2.ntt';
[Timestamps, Samples, Header] = ...
    Nlx2MatSpike_v3( Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
% peakvals = squeeze(max(Samples)); 
conv_factor = str2num(Header{15}(13:end));
peakvals = bsxfun(@times,Samples(1:4,:),conv_factor'); 
clear Samples
[ov_peak,peak_ch] = max(peakvals);
for i = 1:4
    cur_mua = find(peak_ch == i);
%     cur_mua(ov_peak(cur_mua) < 20/1e6) = [];
    mua{i+4} = Timestamps(cur_mua);
end

%%
mua_sm = round(Fsd*0.1);
% for used_ch = 2:8;
used_ch = 2;
hpc_mua = mua{used_ch};
hpc_mua(hpc_mua > t8_d(end) | hpc_mua < t8_d(1)) = [];
binned_spks = hist(hpc_mua,t8_d);
mua_rate = smooth(binned_spks,mua_sm)*Fsd;
mua_rate = zscore(mua_rate);
% 
hold on
plot(t,mua_rate-16,'c')
%%
mua_sm = round(Fsd*0.01);
ctx_ch = 1;
hpc_mua = mua{ctx_ch};
hpc_mua(hpc_mua > t8_d(end) | hpc_mua < t8_d(1)) = [];
binned_spks = hist(hpc_mua,t8_d);
mua_rate = smooth(binned_spks,mua_sm)*Fsd;
mua_rate = zscore(mua_rate);

up = 1e5:8e5;
[S,phi,S12,S1,S2,ff] = coherencysegc(csc8_d(up),mua_rate(up),win,params);

%%
params.Fs = Fsd;
params.tapers = [4 7];
win = 20;
mua_sm = round(Fsd*0.1);
for i = 1:8
    i
used_ch = i;
hpc_mua = mua{used_ch};
hpc_mua(hpc_mua > t8_d(end) | hpc_mua < t8_d(1)) = [];
binned_spks = hist(hpc_mua,t8_d);
mua_rate = smooth(binned_spks,mua_sm)*Fsd;
mua_rate = zscore(mua_rate);
[S(i,:),f] = mtspectrumsegc(mua_rate,win,params);
end

%%
FieldSelectionFlags = [1 0 0 0 1];
HeaderExtractionFlag = 1;
ExtractMode = 1;
Filename = 'Sc1.ntt';
[Timestamps, Samples, Header] = ...
    Nlx2MatSpike_v3( Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
peakvals = squeeze(max(Samples)); 
conv_factor = str2num(Header{15}(13:end));
[ov_peak,peak_ch] = max(peakvals);

best_wvfrms = zeros(32,size(Samples,3));
for i = 1:4
    cur_spks = find(peak_ch == i);
    best_wvfrms(:,cur_spks) = Samples(:,i,cur_spks)*conv_factor(i);
end
best_wvfrms = bsxfun(@rdivide,best_wvfrms,max(best_wvfrms));
clear Samples

Filename = 'Sc2.ntt';
[Timestamps, Samples, Header] = ...
    Nlx2MatSpike_v3( Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
peakvals = squeeze(max(Samples)); 
conv_factor = str2num(Header{15}(13:end));
[ov_peak,peak_ch2] = max(peakvals);

best_wvfrms2 = zeros(32,size(Samples,3));
for i = 1:4
    cur_spks = find(peak_ch2 == i);
    best_wvfrms2(:,cur_spks) = Samples(:,i,cur_spks)*conv_factor(i);
end
best_wvfrms2 = bsxfun(@rdivide,best_wvfrms2,max(best_wvfrms2));
clear Samples
