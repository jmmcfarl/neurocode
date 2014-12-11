clear all
close all

%%
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
l3mec_np = find(ismember(all_cells(uset),l3mec_np));
l3lec_np = find(ismember(all_cells(uset),l3lec_np));
combined_dir = combined_dir(uset);
hpc_mua = hpc_mua(uset);
hpc_lfp = hpc_lfp(uset);
ctx_lfp = ctx_lfp(uset);

%%
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
[b,a] = butter(2,40/(raw_Fs/2),'low');

params.Fs = Fsd;
params.tapers = [3 5];
% movingwin = [10 2.5];
movingwin = [40 2.5];
% rate_sm = round(Fsd*0.05);
%
% params2.Fs = 1/movingwin(2);
% params2.tapers = [4 7];
%

spec_Fs = 1/movingwin(2);
maxlag = round(1500*spec_Fs);

% params2.Fs = spec_Fs;
% params2.tapers = [2 3];
%%
% d = 6;
for d = 1:length(combined_dir)
    cd(combined_dir{d})
    pwd
    load ./used_data wcv_minus_spike lf8 lf7 lf6
    
    if ctx_lfp(d) == 7
        lf8 = lf7;
    elseif ctx_lfp(d) == 6
        lf8 = lf6;
    end
    
    lf8_lf = zscore(downsample(filtfilt(b,a,lf8),dsf));
    wcv_lf = zscore(downsample(filtfilt(b,a,wcv_minus_spike),dsf));
    %     if ~isnan(hpc_mua(d))
    %         load ./mua_data2
    %         load ./sync_times.mat
    %         synctt = downsample(synct,dsf);
    %         hpc_mua_times = mua_times{hpc_mua(d)};
    %         hpc_mua_times(hpc_mua_times < synctt(1) | hpc_mua_times > synctt(end)) = [];
    %         temp_mua_rate =hist(hpc_mua_times,synct)*raw_Fs;
    %         mua_rate = jmm_smooth_1d_cor(temp_mua_rate,rate_sm);
    %         mua_rate = zscore(mua_rate);
    %         if length(mua_rate) > length(lf8_lf)
    %             mua_rate = mua_rate(1:length(lf8_lf));
    %         end
    %     end
    
    [S8,t,f] = mtspecgramc(lf8_lf,movingwin,params);
    [Sw,t,f] = mtspecgramc(wcv_lf,movingwin,params);
    %     if ~isnan(hpc_mua(d))
    %         [Sm,t,f] = mtspecgramc(mua_rate,movingwin,params);
    %     end
    uds_range = find(f >= 0.2 & f <= 1);
    
    lS8 = log(S8);
    lSw = log(Sw);
    %  lS8 = S8;
    %  lSw = Sw;
    %
    [uds_pow8,uds_freq8] = max(lS8(:,uds_range),[],2);
    [uds_poww,uds_freqw] = max(lSw(:,uds_range),[],2);
    
%     uds_ipow8 = sum(log(S8(:,uds_range)),2);
%     uds_ipoww = sum(log(Sw(:,uds_range)),2);
    uds_ipow8 = trapz(lS8(:,uds_range),2)*(f(2)-f(1));
    uds_ipoww = trapz(lSw(:,uds_range),2)*(f(2)-f(1));
    
%     period_8 = uds_pow8./uds_ipow8;
%     period_w = uds_poww./uds_ipoww;
%     avg_period_w(d) = mean(period_w);
%     avg_period_8(d) = mean(period_8);
    
    [aa,bb] = corrcoef(uds_ipow8,uds_ipoww);
    instant_pow_corr(d) = aa(2,1);
    instant_pow_p(d) = bb(2,1);
    
    [aa,bb] = corrcoef(uds_freq8,uds_freqw);
    instant_freq_corr(d) = aa(2,1);
    instant_freq_p(d) = bb(2,1);
    
%     [P88,f] = mtspectrumc(uds_ipow8,params2);
    
    [Cw8(d,:),lags] = xcov(uds_pow8,uds_poww,maxlag,'coeff');
    [Cww(d,:),lags] = xcov(uds_poww,maxlag,'coeff');
    [C88(d,:),lags] = xcov(uds_pow8,maxlag,'coeff');
    [Ciw8(d,:),lags] = xcov(uds_ipow8,uds_ipoww,maxlag,'coeff');
    [Ciww(d,:),lags] = xcov(uds_ipoww,maxlag,'coeff');
    [Ci88(d,:),lags] = xcov(uds_ipow8,maxlag,'coeff');
    [Cfw8(d,:),lags] = xcov(uds_freq8,uds_freqw,maxlag,'coeff');
    [Cfww(d,:),lags] = xcov(uds_freqw,maxlag,'coeff');
    [Cf88(d,:),lags] = xcov(uds_freq8,maxlag,'coeff');
    
%     [Cw8(d,:),lags] = xcov(uds_pow8,uds_poww,maxlag);
%     [Cww(d,:),lags] = xcov(uds_poww,maxlag);
%     [C88(d,:),lags] = xcov(uds_pow8,maxlag);
%     [Ciw8(d,:),lags] = xcov(uds_ipow8,uds_ipoww,maxlag);
%     [Ciww(d,:),lags] = xcov(uds_ipoww,maxlag);
%     [Ci88(d,:),lags] = xcov(uds_ipow8,maxlag);
%     [Cfw8(d,:),lags] = xcov(uds_freq8,uds_freqw,maxlag);
%     [Cfww(d,:),lags] = xcov(uds_freqw,maxlag);
%     [Cf88(d,:),lags] = xcov(uds_freq8,maxlag);
    
    
    %     if ~isnan(hpc_mua(d))
    % %         uds_powm = sum(Sm(:,uds_range),2);
    %         [uds_powm,uds_freqm] = max(Sm(:,uds_range),[],2);
    %        [C8m(d,:),lags] = xcov(uds_pow8,uds_powm,maxlag,'coeff');
    %         [Cwm(d,:),lags] = xcov(uds_poww,uds_powm,maxlag,'coeff');
    %         [Cmm(d,:),lags] = xcov(uds_powm,maxlag,'coeff');
    %        [Cf8m(d,:),lags] = xcov(uds_freq8,uds_freqm,maxlag,'coeff');
    %         [Cfwm(d,:),lags] = xcov(uds_freqw,uds_freqm,maxlag,'coeff');
    %         [Cfmm(d,:),lags] = xcov(uds_freqm,maxlag,'coeff');
    %     end
   
    %%
%     figure
%     set(gca,'fontsize',14,'fontname','arial')
%     pcolor(t,f,log(S8'));shading flat
%     ylim([0.1 3])
%     set(gca,'yscale','log')
%     xlabel('Time (s)','fontsize',16,'fontname','arial')
%     ylabel('Frequency (Hz)','fontsize',16,'fontname','arial')
%     caxis([-5 1])
%     
%     figure
%     set(gca,'fontsize',14,'fontname','arial')
%     pcolor(t,f,log(Sw'));shading flat
%     ylim([0.1 3])
%     set(gca,'yscale','log')
%     caxis([-5 1])
%     xlabel('Time (s)','fontsize',16,'fontname','arial')
%     ylabel('Frequency (Hz)','fontsize',16,'fontname','arial')
%     
%     figure
%     set(gca,'fontsize',14,'fontname','arial')
%     plot(lags/spec_Fs,Ci88(d,:),'k')
%     hold on
%     plot(lags/spec_Fs,Ciww(d,:),'r')
%     xlim([-1000 1000])
%     xlabel('Time (s)','fontsize',16,'fontname','arial')
%     ylabel('Correlation','fontsize',16,'fontname','arial')
%     
%     figure
%     set(gca,'fontsize',14,'fontname','arial')
%     plot(lags/spec_Fs,Cf88(d,:),'k')
%     hold on
%     plot(lags/spec_Fs,Cfww(d,:),'r')
%     xlim([-1000 1000])
%      xlabel('Time (s)','fontsize',16,'fontname','arial')
%     ylabel('Correlation','fontsize',16,'fontname','arial')
%    
%     pause
%     close all

zero_lag = find(lags == 0);

end

%%
% l3mec_m = l3mec(~isnan(hpc_mua(l3mec)));

% Cw8_z = Cw8(:,zero_lag);
% 
% figure
% set(gca,'fontsize',14,'fontname','arial')
% hold on
% shadedErrorBar(lags/spec_Fs,nanmean(Cww(l3lec,:)),std(Cww(l3lec,:))/sqrt(length(l3lec)),{'b'});
% shadedErrorBar(lags/spec_Fs,nanmean(C88),std(C88)/sqrt(length(combined_dir)),{'k'});
% shadedErrorBar(lags/spec_Fs,nanmean(Cww(l3mec,:)),std(Cww(l3mec,:))/sqrt(length(l3mec)),{'r'});
% % shadedErrorBar(lags/spec_Fs,nanmean(Cmm(l3mec_m,:)),std(Cmm(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});
% xlabel('Time (s)','fontsize',16)
% ylabel('Correlation','fontsize',16)
% % ylim([-0.1 0.15])
% % xlim([-300 300])
% 
% % xlim([-300 300])
% 
% figure
% set(gca,'fontsize',14,'fontname','arial')
% hold on
% shadedErrorBar(lags/spec_Fs,nanmean(Cpww(l3lec,:)),std(Cpww(l3lec,:))/sqrt(length(l3lec)),{'b'});
% shadedErrorBar(lags/spec_Fs,nanmean(Cp88),std(C88)/sqrt(length(combined_dir)),{'k'});
% shadedErrorBar(lags/spec_Fs,nanmean(Cpww(l3mec,:)),std(Cpww(l3mec,:))/sqrt(length(l3mec)),{'r'});
% % shadedErrorBar(lags/spec_Fs,nanmean(Cmm(l3mec_m,:)),std(Cmm(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});
% xlabel('Time (s)','fontsize',16)
% ylabel('Correlation','fontsize',16)
% 
% figure
% set(gca,'fontsize',14,'fontname','arial')
% hold on
% shadedErrorBar(lags/spec_Fs,nanmean(Ciww(l3lec,:)),std(Ciww(l3lec,:))/sqrt(length(l3lec)),{'b'});
% shadedErrorBar(lags/spec_Fs,nanmean(Ci88),std(C88)/sqrt(length(combined_dir)),{'k'});
% shadedErrorBar(lags/spec_Fs,nanmean(Ciww(l3mec,:)),std(Ciww(l3mec,:))/sqrt(length(l3mec)),{'r'});
% % shadedErrorBar(lags/spec_Fs,nanmean(Cmm(l3mec_m,:)),std(Cmm(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});
% xlabel('Time (s)','fontsize',16)
% ylabel('Correlation','fontsize',16)
% 
% figure
% set(gca,'fontsize',14,'fontname','arial')
% hold on
% shadedErrorBar(lags/spec_Fs,nanmean(Cfww(l3lec,:)),std(Cfww(l3lec,:))/sqrt(length(l3lec)),{'b'});
% shadedErrorBar(lags/spec_Fs,nanmean(Cf88),std(Cf88)/sqrt(length(combined_dir)),{'k'});
% shadedErrorBar(lags/spec_Fs,nanmean(Cfww(l3mec,:)),std(Cfww(l3mec,:))/sqrt(length(l3mec)),{'r'});
% % shadedErrorBar(lags/spec_Fs,nanmean(Cfmm(l3mec_m,:)),std(Cfmm(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});
% xlabel('Time (s)','fontsize',16)
% ylabel('Correlation','fontsize',16)
% % ylim([-0.15 0.15])
% % xlim([-300 300])
% 
% 
% figure
% set(gca,'fontsize',14,'fontname','arial')
% hold on
% shadedErrorBar(lags/spec_Fs,nanmean(Cw8(l3lec,:)),std(Cw8(l3lec,:))/sqrt(length(l3lec)),{'b'});
% shadedErrorBar(lags/spec_Fs,nanmean(Cw8(l3mec,:)),std(Cw8(l3mec,:))/sqrt(length(l3mec)),{'r'});
% % shadedErrorBar(lags/spec_Fs,nanmedian(C8m(l3mec_m,:)),std(C8m(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});
% xlabel('Time (s)','fontsize',16)
% ylabel('Correlation','fontsize',16)
% % ylim([-0.2 0.25])
% xlim([-500 500])
% 
% figure
% set(gca,'fontsize',14,'fontname','arial')
% hold on
% shadedErrorBar(lags/spec_Fs,nanmean(Ciw8(l3lec,:)),std(Ciw8(l3lec,:))/sqrt(length(l3lec)),{'b'});
% shadedErrorBar(lags/spec_Fs,nanmean(Ciw8(l3mec,:)),std(Ciw8(l3mec,:))/sqrt(length(l3mec)),{'r'});
% % shadedErrorBar(lags/spec_Fs,nanmedian(C8m(l3mec_m,:)),std(C8m(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});
% xlabel('Time (s)','fontsize',16)
% ylabel('Correlation','fontsize',16)
% % ylim([-0.2 0.25])
% xlim([-500 500])
% 
% figure
% set(gca,'fontsize',14,'fontname','arial')
% hold on
% shadedErrorBar(lags/spec_Fs,nanmean(Cpw8(l3lec,:)),std(Cpw8(l3lec,:))/sqrt(length(l3lec)),{'b'});
% shadedErrorBar(lags/spec_Fs,nanmean(Cpw8(l3mec,:)),std(Cpw8(l3mec,:))/sqrt(length(l3mec)),{'r'});
% % shadedErrorBar(lags/spec_Fs,nanmedian(C8m(l3mec_m,:)),std(C8m(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});
% xlabel('Time (s)','fontsize',16)
% ylabel('Correlation','fontsize',16)
% % ylim([-0.2 0.25])
% xlim([-500 500])
% 
% figure
% set(gca,'fontsize',14,'fontname','arial')
% hold on
% shadedErrorBar(lags/spec_Fs,nanmean(Cfw8(l3lec,:)),std(Cfw8(l3lec,:))/sqrt(length(l3lec)),{'b'});
% shadedErrorBar(lags/spec_Fs,nanmean(Cfw8(l3mec,:)),std(Cfw8(l3mec,:))/sqrt(length(l3mec)),{'r'});
% % shadedErrorBar(lags/spec_Fs,nanmedian(Cf8m(l3mec_m,:)),std(Cf8m(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});
% xlabel('Time (s)','fontsize',16)
% ylabel('Correlation','fontsize',16)
% xlim([-500 500])
% 
% 
% % figure
% % hold on
% % shadedErrorBar(lags/spec_Fs,nanmedian(Cwm(l3mec_m,:)),std(Cww(l3mec_m,:))/sqrt(length(l3mec_m)),{'b'});
% % % shadedErrorBar(lags/spec_Fs,nanmedian(C8m(l3mec_m,:)),std(C8m(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});

%%
% to_disp = 1;
% close all
% smth = 20;
% n_cells = size(Cw8,1);
% Ci88_period = nan(n_cells,1);
% Ci88_mod = nan(n_cells,1);
% C88_period = nan(n_cells,1);
% C88_mod = nan(n_cells,1);
% Cf88_period = nan(n_cells,1);
% Cf88_mod =nan(n_cells,1);
% Ciww_period = nan(n_cells,1);
% Ciww_mod = nan(n_cells,1);
% Cww_period = nan(n_cells,1);
% Cww_mod = nan(n_cells,1);
% Cfww_period = nan(n_cells,1);
% Cfww_mod = nan(n_cells,1);
% 
% for i = 1:size(Cw8,1)
%     [min_amp,min_loc,max_amp,max_loc,sm_Ci88] = get_acorr_peaks(Ci88(i,:),lags,smth,1);
%     if ~isempty(max_loc) & ~isempty(min_loc)
%         Ci88_period(i) = lags(max_loc);
%         Ci88_mod(i) = max_amp+min_amp;
%     end
%     [min_amp,min_loc,max_amp,max_loc,sm_C88] = get_acorr_peaks(C88(i,:),lags,smth);
%     if ~isempty(max_loc) & ~isempty(min_loc)
%         C88_period(i) = lags(max_loc);
%         C88_mod(i) = max_amp+min_amp;
%     end
%     [min_amp,min_loc,max_amp,max_loc,sm_Cf88] = get_acorr_peaks(Cf88(i,:),lags,smth);
%     if ~isempty(max_loc) & ~isempty(min_loc)
%         Cf88_period(i) = lags(max_loc);
%         Cf88_mod(i) = max_amp+min_amp;
%     end
%     [min_amp,min_loc,max_amp,max_loc,sm_Ciww] = get_acorr_peaks(Ciww(i,:),lags,smth);
%     if ~isempty(max_loc) & ~isempty(min_loc)
%         Ciww_period(i) = lags(max_loc);
%         Ciww_mod(i) = max_amp+min_amp;
%     end
%     [min_amp,min_loc,max_amp,max_loc,sm_Cww] = get_acorr_peaks(Cww(i,:),lags,smth);
%     if ~isempty(max_loc) & ~isempty(min_loc)
%         Cww_period(i) = lags(max_loc);
%         Cww_mod(i) = max_amp+min_amp;
%     end
%     [min_amp,min_loc,max_amp,max_loc,sm_Cfww] = get_acorr_peaks(Cfww(i,:),lags,smth);
%     if ~isempty(max_loc) & ~isempty(min_loc)
%         Cfww_period(i) = lags(max_loc);
%         Cfww_mod(i) = max_amp+min_amp;
%     end
%     [min_amp,min_loc,max_amp,max_loc,sm_Ciw8] = get_acorr_peaks(Ciw8(i,:),lags,smth);
%     if ~isempty(max_loc) & ~isempty(min_loc)
%         Ciw8_period(i) = lags(max_loc);
%         Ciw8_mod(i) = max_amp+min_amp;
%     end
%     [min_amp,min_loc,max_amp,max_loc,sm_Cw8] = get_acorr_peaks(Cw8(i,:),lags,smth);
%     if ~isempty(max_loc) & ~isempty(min_loc)
%         Cw8_period(i) = lags(max_loc);
%         Cw8_mod(i) = max_amp+min_amp;
%     end
%     [min_amp,min_loc,max_amp,max_loc,sm_Cfw8] = get_acorr_peaks(Cfw8(i,:),lags,smth);
%     if ~isempty(max_loc) & ~isempty(min_loc)
%         Cfw8_period(i) = lags(max_loc);
%         Cfw8_mod(i) = max_amp+min_amp;
%     end
%     Ci88_period(i)/spec_Fs
%     Ciww_period(i)/spec_Fs
%     if to_disp == 1
%         if ismember(i,l3mec)
%             sprintf('L3MEC')
%         elseif ismember(i,l3lec)
%             sprintf('L3LEC')
%         end
%         subplot(3,2,1)
%         hold on
%         plot(lags/spec_Fs,sm_Ci88,'r','linewidth',2)
%         plot(lags/spec_Fs,Ci88(i,:))
%         subplot(3,2,2)
%         plot(lags/spec_Fs,sm_Cf88,'r','linewidth',2)
%         hold on
%         plot(lags/spec_Fs,Cf88(i,:))
%         subplot(3,2,3)
%         plot(lags/spec_Fs,sm_Ciww,'r','linewidth',2)
%         hold on
%         plot(lags/spec_Fs,Ciww(i,:))
%         subplot(3,2,4)
%         plot(lags/spec_Fs,sm_Cfww,'r','linewidth',2)
%         hold on
%         plot(lags/spec_Fs,Cfww(i,:))
%         subplot(3,2,5)
%         plot(lags/spec_Fs,sm_Ciw8,'r','linewidth',2)
%         hold on
%         plot(lags/spec_Fs,Ciw8(i,:))
%         subplot(3,2,6)
%         plot(lags/spec_Fs,sm_Cfw8,'r','linewidth',2)
%         hold on
%         plot(lags/spec_Fs,Cfw8(i,:))
%         
%         if Ci88_period(i)/spec_Fs < 100
%             disp('SHORT')
%         end
%         
%         
%         pause
%         clf
%     end
%     
% end
% 
% Ci88_period = Ci88_period/spec_Fs;
% Cf88_period = Cf88_period/spec_Fs;
% Ciww_period = Ciww_period/spec_Fs;
% Cfww_period = Cfww_period/spec_Fs;
% Ciw8_period = Ciw8_period/spec_Fs;
% Cfw8_period = Cfw8_period/spec_Fs;
% 
%%
close all
L = length(lags);
NFFT = 2^nextpow2(L);
f = spec_Fs/2*linspace(0,1,NFFT/2+1);
lim = find(1./f > 750,1,'last');
lim2 = find(1./f < 50,1,'first');
for i = 1:size(Ciw8,1)
    Y = fft(Ci88(i,:),NFFT)/length(lags);
    Y = 2*abs(Y(1:NFFT/2+1));
    [peakamps,peaklocs] = findpeaks(Y(lim:lim2));
        peaklocs = peaklocs + lim-1;
[a,b] = sort(peakamps);
    Ci88_peakpow_s(i,:) = Y;
    Ci88_peakpow(i) = a(end);
    Ci88_peakf(i) = f(peaklocs(b(end)));

    Y = fft(C88(i,:),NFFT)/length(lags);
    Y = 2*abs(Y(1:NFFT/2+1));
    [peakamps,peaklocs] = findpeaks(Y(lim:end));
        peaklocs = peaklocs + lim-1;
[a,b] = sort(peakamps);
    C88_peakpow_s(i,:) = Y;
    C88_peakpow(i) = a(end);
    C88_peakf(i) = f(peaklocs(b(end)));
    
    Y = fft(Cf88(i,:),NFFT)/length(lags);
    Y = 2*abs(Y(1:NFFT/2+1));
    [peakamps,peaklocs] = findpeaks(Y(lim:lim2));
        peaklocs = peaklocs + lim-1;
[a,b] = sort(peakamps);
    Cf88_peakpow_s(i,:) = Y;
    Cf88_peakpow(i) = a(end);
    Cf88_peakf(i) = f(peaklocs(b(end)));
    
    Y = fft(Ciww(i,:),NFFT)/length(lags);
    Y = 2*abs(Y(1:NFFT/2+1));
    [peakamps,peaklocs] = findpeaks(Y(lim:lim2));
    peaklocs = peaklocs + lim-1;
    [a,b] = sort(peakamps);
    Ciww_peakpow_s(i,:) = Y;
    Ciww_peakpow(i) = a(end);
    Ciww_peakf(i) = f(peaklocs(b(end)));
    
%             subplot(2,1,1)
%     plot(lags/spec_Fs,Ciww(i,:))
%     xlim([0 2000])
%     subplot(2,1,2)
%     plot(1./f,Y)
%     hold on
%     plot(1./f(peaklocs(b(end))),Y(peaklocs(b(end))),'ro')
%     xlim([0 2000])
%     pause
%     clf

    Y = fft(Cww(i,:),NFFT)/length(lags);
    Y = 2*abs(Y(1:NFFT/2+1));
    [peakamps,peaklocs] = findpeaks(Y(lim:end));
    peaklocs = peaklocs + lim-1;
    [a,b] = sort(peakamps);
    Cww_peakpow_s(i,:) = Y;
    Cww_peakpow(i) = a(end);
    Cww_peakf(i) = f(peaklocs(b(end)));
    
    Y = fft(Cfww(i,:),NFFT)/length(lags);
    Y = 2*abs(Y(1:NFFT/2+1));
    [peakamps,peaklocs] = findpeaks(Y(lim:lim2));
    peaklocs = peaklocs + lim-1;
    [a,b] = sort(peakamps);
    Cfww_peakpow_s(i,:) = Y;
    Cfww_peakpow(i) = a(end);
    Cfww_peakf(i) = f(peaklocs(b(end)));

    Y = fft(Ciw8(i,:),NFFT)/length(lags);
    Y = 2*abs(Y(1:NFFT/2+1));
    [peakamps,peaklocs] = findpeaks(Y(lim:lim2));
    peaklocs = peaklocs + lim-1;
    [a,b] = sort(peakamps);
    Ciw8_peakpow_s(i,:) = Y;
    Ciw8_peakpow(i) = a(end);
    Ciw8_peakf(i) = f(peaklocs(b(end)));

        Y = fft(Cw8(i,:),NFFT)/length(lags);
    Y = 2*abs(Y(1:NFFT/2+1));
    [peakamps,peaklocs] = findpeaks(Y(lim:end));
    peaklocs = peaklocs + lim-1;
    [a,b] = sort(peakamps);
    Cw8_peakpow_s(i,:) = Y;
    Cw8_peakpow(i) = a(end);
    Cw8_peakf(i) = f(peaklocs(b(end)));

        Y = fft(Cfw8(i,:),NFFT)/length(lags);
    Y = 2*abs(Y(1:NFFT/2+1));
    [peakamps,peaklocs] = findpeaks(Y(lim:lim2));
    peaklocs = peaklocs + lim-1;
    [a,b] = sort(peakamps);
    Cfw8_peakpow_s(i,:) = Y;
    Cfw8_peakpow(i) = a(end);
    Cfw8_peakf(i) = f(peaklocs(b(end)));

    cur_Cf88_period = 1/Cf88_peakf(i);
    [~,close_lag] = min(abs(cur_Cf88_period-lags/spec_Fs));
    Cfw8_peakxcorr(i) = Cfw8(i,close_lag);
     cur_Ci88_period = 1/Cf88_peakf(i);
    [~,close_lag] = min(abs(cur_Ci88_period-lags/spec_Fs));
    Ciw8_peakxcorr(i) = Ciw8(i,close_lag);
%     Cw8_peakxcorr(i) = Cw8(i,close_lag);
    
end

%%
close all
figure
set(gca,'fontsize',14,'fontname','arial')
plot(Cf88_peakpow(l3mec),Cfww_peakpow(l3mec),'ro','linewidth',1,'markersize',8);
hold on
plot(Cf88_peakpow(l3lec),Cfww_peakpow(l3lec),'bo','linewidth',1,'markersize',8);
xlabel('Ncx Infraslow Power','fontsize',16,'fontname','arial')
ylabel('EC Infraslow Power','fontsize',16,'fontname','arial')
xlim([0 0.3]); ylim([0 0.3])

figure
set(gca,'fontsize',14,'fontname','arial')
plot(1./Cf88_peakf(l3mec),1./Cfww_peakf(l3mec),'ro','linewidth',1,'markersize',8);
hold on
plot(1./Cf88_peakf(l3lec),1./Cfww_peakf(l3lec),'bo','linewidth',1,'markersize',8);
xlabel('Ncx Infraslow Period (s)','fontsize',16,'fontname','arial')
ylabel('EC Infraslow  Period (s)','fontsize',16,'fontname','arial')
xlim([0 700]); ylim([0 700])
%%
close all
figure
set(gca,'fontsize',14,'fontname','arial')
plot(C88_peakpow(l3mec),Cww_peakpow(l3mec),'ro','linewidth',1,'markersize',8);
hold on
plot(C88_peakpow(l3lec),Cww_peakpow(l3lec),'bo','linewidth',1,'markersize',8);
xlabel('Ncx Infraslow Power','fontsize',16,'fontname','arial')
ylabel('EC Infraslow Power','fontsize',16,'fontname','arial')
xlim([0 0.3]); ylim([0 0.3])

figure
set(gca,'fontsize',14,'fontname','arial')
plot(1./C88_peakf(l3mec),1./Cww_peakf(l3mec),'ro','linewidth',1,'markersize',8);
hold on
plot(1./C88_peakf(l3lec),1./Cww_peakf(l3lec),'bo','linewidth',1,'markersize',8);
xlabel('Ncx Infraslow Period (s)','fontsize',16,'fontname','arial')
ylabel('EC Infraslow  Period (s)','fontsize',16,'fontname','arial')
xlim([0 700]); ylim([0 700])

%%
cd C:\WC_Germany\sven_thomas_combined
load ./combined_core_analysis_fin_nd.mat

maxlag = round(1500*spec_Fs);
lags = -maxlag:maxlag;
%%
figure
subplot(2,2,1)
plot(Cf88_peakpow(l3mec),fract_rt2_ups(l3mec),'o')
xlabel('Cortical ISLOW-f Power','fontsize',14,'fontname','arial')
ylabel('Persistence','fontsize',14,'fontname','arial')
subplot(2,2,2)
plot(Ci88_peakpow(l3mec),fract_rt2_ups(l3mec),'o')
xlabel('Cortical ISLOW-p Power','fontsize',14,'fontname','arial')
ylabel('Persistence','fontsize',14,'fontname','arial')
subplot(2,2,3)
plot(Cfww_peakpow(l3mec),fract_rt2_ups(l3mec),'o')
xlabel('MEC ISLOW-f Power','fontsize',14,'fontname','arial')
ylabel('Persistence','fontsize',14,'fontname','arial')
subplot(2,2,4)
plot(Ciww_peakpow(l3mec),fract_rt2_ups(l3mec),'o')
xlabel('MEC ISLOW-p Power','fontsize',14,'fontname','arial')
ylabel('Persistence','fontsize',14,'fontname','arial')