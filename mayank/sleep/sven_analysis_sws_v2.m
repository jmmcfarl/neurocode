clear all
% close all
%
% cd /Users/James/Data/Sven/2012_06_19/2012-6-19_12-21-9
%
% cd /Users/James/Data/Sven/2012_06-22/2012-6-22_14-31-1
% cd /Users/James/Data/Sven/2012_06-22/2012-6-22_18-52-27
%
% cd /Users/James/Data/Sven/2012-6-23_13-45-41
% cd /Users/James/Data/Sven/2012-6-23_14-44-12
%
% cd /Users/James/Data/Sven/2012_06_27/2012-6-27_14-35-25
% cd /Users/James/Data/Sven/2012_06_27/2012-6-27_15-36-56
% cd /Users/James/Data/Sven/2012_06_27/2012-6-27_16-42-48
cd /Users/James/Data/Sven/2012_06_27/2012-6-27_17-42-24

% cd ~/Data/Sven/2012-6-30_16-6-48/
% cd ~/Data/Sven/2012-6-30_17-29-16/
% cd ~/Data/Sven/2012-6-30_18-44-18/


% cd ~/Data/Sven/2012-7-2_15-5-21/
% cd ~/Data/Sven/2012-7-2_16-19-30/
% cd ~/Data/Sven/2012-7-2_17-17-14/

%% LOAD RAW DATA
% FieldSelectionFlags = [1 0 1 0 1];
% ExtractMode = 1;
% HeaderExtractionFlag = 1;
% 
% for i = 1:8
%     fprintf('Loading CSC%d\n',i)
%     Filename = sprintf('CSC%d.Ncs',i);
%     [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
%         Nlx2MatCSC_v3(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
%     CSC_Timestamps = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
%     Fs(i) = median(CSC_SampleFrequencies);
%     csc(i,:) = CSC_Samples(:)*str2num(Header{15}(13:end));
%     clear CSC_Samples
% end
% if length(unique(Fs)) > 1
%     error('Unequal sample rates')
% end
% Fs = unique(Fs);
% 
% %% DOWN-SAMPLE
% dsf = 8;
% Fsd = Fs/dsf;
% for i = 1:8
%     cscd(:,i) = decimate(csc(i,:),dsf);
% end
% t = (1:size(cscd,1))/Fsd;
%%
temp = pwd;
sl = find(temp == '/',1,'last');
dname = strcat(temp(sl+1:end),'.mat');
load(dname)

%% COMPUTE SPECTROGRAMS
movingwin = [10 10];
params.Fs = Fsd;
params.tapers = [2 3];
clear Sg tt f
for i = 1:8
    [Sg(i,:,:),tt,f]=mtspecgramc(cscd(:,i),movingwin,params);
end

[~,uds_freq] = min(abs(f - 1.5)); 
uds_pow = squeeze(Sg(:,:,uds_freq));
low_freqs = find(f < 0.5);
lf_pow = squeeze(sum(Sg(:,:,low_freqs),3));

%% IDENTIFY UDS EPOCHS
min_uds_pow = 3e-6;
ctx_ch = 1;
uds_epochs = find(uds_pow(ctx_ch,:) >= min_uds_pow);

figure
plot(tt,uds_pow(ctx_ch,:));hold on
plot(tt(uds_epochs),uds_pow(ctx_ch,uds_epochs),'r.')

uds_inds = [];
for i = 1:length(uds_epochs)
    cur_set = find(t > tt(uds_epochs(i)) - movingwin(1)/2 & t < tt(uds_epochs(i))+movingwin(1)/2);
   uds_inds = [uds_inds cur_set];
end
uds_inds = unique(uds_inds);

%% COMPUTE POWER SPECTRA
% win = [20 10];
win = 10;
params.Fs = Fsd;
params.tapers = [4 7];
params.trialave = 0;
for i = 2:8
    fprintf('Channel %d\n',i);
%     [C(i,:),~,~,S(1,:),S(i,:),f] = coherencysegc(cscd(:,1),cscd(:,i),win,params);
    [C_temp,~,~,S1_temp,S_temp,f] = coherencysegc(cscd(:,1),cscd(:,i),win,params);
    C(i,:) = mean(C_temp,2);
    S(1,:) = mean(S1_temp,2);
    S(i,:) = mean(S_temp,2);
    if ~isempty(uds_epochs)
       C_uds(i,:) = mean(C_temp(:,uds_epochs),2);
       S_uds(1,:) = mean(S1_temp(:,uds_epochs),2);
       S_uds(i,:) = mean(S_temp(:,uds_epochs),2);
    end
end
normS = bsxfun(@rdivide,S,sum(S,2));
normS_uds = bsxfun(@rdivide,S_uds,sum(S_uds,2));
% figure
% plot(f,log(normS));set(gca,'xscale','log')



%% LOAD MUA
% FieldSelectionFlags = [1 0 0 1 0];
% HeaderExtractionFlag = 1;
% ExtractMode = 1;
% Filename = 'Sc1.ntt';
% [Timestamps, Samples, Header] = ...
%     Nlx2MatSpike_v3( Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
% conv_factor = str2num(Header{15}(13:end));
% peakvals = bsxfun(@times,Samples(1:4,:),conv_factor');
% clear Samples
% [ov_peak,peak_ch] = max(peakvals);
% for i = 1:4
%     cur_mua = find(peak_ch == i);
%     mua{i} = Timestamps(cur_mua);
% end
% Filename = 'Sc2.ntt';
% [Timestamps, Samples, Header] = ...
%     Nlx2MatSpike_v3( Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
% conv_factor = str2num(Header{15}(13:end));
% peakvals = bsxfun(@times,Samples(1:4,:),conv_factor');
% clear Samples
% [ov_peak,peak_ch] = max(peakvals);
% for i = 1:4
%     cur_mua = find(peak_ch == i);
%     mua{i+4} = Timestamps(cur_mua);
% end
% 
%% BIN MUA
% down_ts = downsample(CSC_Timestamps,dsf);
% for ch = 1:8;
%     cur_mua = mua{ch};
%     cur_mua(cur_mua > down_ts(end) | cur_mua < down_ts(1)) = [];
%     binned_spks(:,ch) = hist(cur_mua,down_ts)';
% end
% binned_spks(end,:) = mean(binned_spks);

%% MUA SPECTRA
% win = [20 10];
win = [10 5];
params.Fs = Fsd;
params.tapers = [4 7];

mec_ch = 6;
ctx_ch = 3;
for ch = 1:8
    ch
%     [Cmua_mec(ch,:),~,~,~,Smua(ch,:),f,zerosp]=coherencysegcpb(cscd(:,mec_ch),binned_spks(:,ch),win,params,0);
    [Cmua_temp,~,~,~,Smua_temp,f,zerosp]=coherencysegcpb(cscd(:,mec_ch),binned_spks(:,ch),win,params,0);
    Cmua_mec(ch,:) = mean(Cmua_temp,2);
    Smua(ch,:) = mean(Smua_temp,2);
    if ~isempty(uds_inds)
       Cmua_mec_uds(ch,:) = mean(Cmua_temp(:,uds_epochs),2);
       Smua_uds(ch,:) = mean(Smua_temp(:,uds_epochs),2);
    end
end
Smua_norm = bsxfun(@rdivide,Smua,sum(Smua,2));
Smua_uds_norm = bsxfun(@rdivide,Smua_uds,sum(Smua_uds,2));

for ch = 1:8
    ch
    [Cmua_ctx(ch,:),~,~,~,Smua(ch,:),f,zerosp]=coherencysegcpb(cscd(:,ctx_ch),binned_spks(:,ch),win,params);
end
Smua_norm = bsxfun(@rdivide,Smua,sum(Smua,2));
%% COMPUTE FILTERED LFPS AND SMOOTHED MUAS
[b,a] = butter(2,[0.4 5]/(Fsd/2));
mua_sm = round(Fsd*0.1);
for i = 1:8
    csc_f(:,i) = zscore(filtfilt(b,a,cscd(:,i)));
    sm_rate(:,i) = zscore(smooth(binned_spks(:,i),mua_sm)*Fsd);
end


%%
ctx_ch = 3;
mec_ch = 6;

% csc_f = -csc_f;
min_amp = 3;
[down_blip_amps,down_blip_locs] = findpeaks(-csc_f(uds_inds,ctx_ch),'minpeakheight',min_amp);
down_blip_locs = uds_inds(down_blip_locs);
[down_blip_amps_mec,down_blip_locs_mec] = findpeaks(-csc_f(uds_inds,mec_ch),'minpeakheight',min_amp);
down_blip_locs_mec = uds_inds(down_blip_locs_mec);
% csc_f = -csc_f;

% figure
% plot(t,csc_f(:,1))
% hold on
% plot(t(down_blip_locs),csc_f(down_blip_locs,1),'r.')
% 
maxlag = round(Fsd*0.5);
down_blip_locs(down_blip_locs < maxlag | down_blip_locs > size(csc_f,1)-maxlag);
trig_avgs = zeros(2*maxlag+1,8);
mtrig_avgs = zeros(2*maxlag+1,8);
for i = 1:length(down_blip_locs)
    cur_inds = (down_blip_locs(i)-maxlag):(down_blip_locs(i)+maxlag);
    trig_avgs = trig_avgs + csc_f(cur_inds,:);
    mtrig_avgs = mtrig_avgs + sm_rate(cur_inds,:);
end
trig_avgs = trig_avgs/length(down_blip_locs);
mtrig_avgs = mtrig_avgs/length(down_blip_locs);

down_blip_locs_mec(down_blip_locs_mec < maxlag | down_blip_locs_mec > size(csc_f,1)-maxlag);
trig_avgs_mec = zeros(2*maxlag+1,8);
mtrig_avgs_mec = zeros(2*maxlag+1,8);
for i = 1:length(down_blip_locs_mec)
    cur_inds = (down_blip_locs_mec(i)-maxlag):(down_blip_locs_mec(i)+maxlag);
    trig_avgs_mec = trig_avgs_mec + csc_f(cur_inds,:);
    mtrig_avgs_mec = mtrig_avgs_mec + sm_rate(cur_inds,:);
end
trig_avgs_mec = trig_avgs_mec/length(down_blip_locs_mec);
mtrig_avgs_mec = mtrig_avgs_mec/length(down_blip_locs_mec);

lags = (-maxlag:maxlag)/Fsd;

% figure
% subplot(2,2,1)
% plot(lags,trig_avgs);hold on;plot(lags,trig_avgs(:,1),'k','linewidth',2)
% title('Ctx DS Trig LFP','fontsize',18)
% xlabel('Time (s)','fontsize',16)
% ylabel('Relative Amp (z)','fontsize',16)
% subplot(2,2,2)
% plot(lags,mtrig_avgs);hold on;plot(lags,mtrig_avgs(:,1),'k','linewidth',2)
% title('Ctx DS Trig MUA','fontsize',18)
% xlabel('Time (s)','fontsize',16)
% ylabel('Relative Amp (z)','fontsize',16)
% subplot(2,2,3)
% plot(lags,trig_avgs_mec);hold on;plot(lags,trig_avgs_mec(:,1),'k','linewidth',2)
% title('MEC DS Trig LFP','fontsize',18)
% xlabel('Time (s)','fontsize',16)
% ylabel('Relative Amp (z)','fontsize',16)
% subplot(2,2,4)
% plot(lags,mtrig_avgs_mec);hold on;plot(lags,mtrig_avgs_mec(:,1),'k','linewidth',2)
% title('MEC DS Trig MUA','fontsize',18)
% xlabel('Time (s)','fontsize',16)
% ylabel('Relative Amp (z)','fontsize',16)
figure
subplot(2,2,1)
plot(lags,trig_avgs(:,1),'k')
hold on
for i = 2:4
    plot(lags,trig_avgs(:,i),'b')
    hold on
end
for i = 5:8
    plot(lags,trig_avgs(:,i),'r')
end
title('Ctx DS Trig LFP','fontsize',18)
xlabel('Time (s)','fontsize',16)
ylabel('Relative Amp (z)','fontsize',16)
subplot(2,2,2)
plot(lags,mtrig_avgs(:,1),'k')
hold on
for i = 2:4
    plot(lags,mtrig_avgs(:,i),'b')
    hold on
end
for i = 5:8
    plot(lags,mtrig_avgs(:,i),'r')
end
title('Ctx DS Trig MUA','fontsize',18)
xlabel('Time (s)','fontsize',16)
ylabel('Relative Amp (z)','fontsize',16)
subplot(2,2,3)
plot(lags,trig_avgs_mec(:,1),'k')
hold on
for i = 2:4
    plot(lags,trig_avgs_mec(:,i),'b')
    hold on
end
for i = 5:8
    plot(lags,trig_avgs_mec(:,i),'r')
end
title('MEC DS Trig LFP','fontsize',18)
xlabel('Time (s)','fontsize',16)
ylabel('Relative Amp (z)','fontsize',16)
subplot(2,2,4)
plot(lags,mtrig_avgs_mec(:,1),'k')
hold on
for i = 2:4
    plot(lags,mtrig_avgs_mec(:,i),'b')
    hold on
end
for i = 5:8
    plot(lags,mtrig_avgs_mec(:,i),'r')
end
title('MEC DS Trig MUA','fontsize',18)
xlabel('Time (s)','fontsize',16)
ylabel('Relative Amp (z)','fontsize',16)



