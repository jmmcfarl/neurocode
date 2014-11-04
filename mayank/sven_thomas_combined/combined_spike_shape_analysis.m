clear all
close all

cur_dir_letter = 'C';
addpath('C:\WC_Germany\persistent_9_27_2010\')
addpath('C:\WC_Germany\new_mec\')
cd C:\WC_Germany\sven_thomas_combined\
% load ./combined_dir.mat
load ./combined_dir_nd.mat

%%
% l3mec_np(ismember(l3mec_np,[64 69])) = [];
uset = sort([l3mec l3lec l3mec_np l3lec_np]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
l3mec_np = find(ismember(all_cells(uset),l3mec_np));
l3lec_np = find(ismember(all_cells(uset),l3lec_np));
l3mec_np(l3mec_np == 59) = [];
l3mec_np(l3mec_np == 60) = [];
l3mec_np(l3mec_np == 62) = []; %int
int = 62;
combined_dir = combined_dir(uset);
% hpc_mua = hpc_mua(uset);
% hpc_lfp = hpc_lfp(uset);
ctx_lfp = ctx_lfp(uset);

Fs = 32258;
maxlag = round(Fs*0.1);
[b,a] = butter(2,60/(Fs/2),'high');

acorr_dt = 2e-3;
acorr_maxlag = round(0.5/acorr_dt);

isi_range = logspace(log10(2e-3),log10(5e-2),51);  

%%
for d = 1:length(combined_dir)
    cur_dir = combined_dir{d};
    cd(cur_dir)
    pwd
    
    load ./spike_time_jmm.mat
    load ./sync_times.mat
    if exist('./all_eeg_data2.mat','file')
        load ./all_eeg_data2 CSC1_Samples
    elseif exist('./all_eeg_data.mat','file')
        load ./all_eeg_data CSC1_Samples
    else
        disp('ERROR DATA DOES NOT EXIST')
        return
    end
    
    wcv = -CSC1_Samples(:);
    wcv = filtfilt(b,a,wcv);
%     wcv_d = [0; diff(wcv)];
    
    spkid(synct1id(spkid) < maxlag | synct1id(spkid) > length(wcv)-maxlag) = [];
    
    beg_spkids = find(synct1id(spkid) < length(wcv)/20);%first 20% of spikes
    
    n_spks(d) = length(spkid);
    if length(spkid) >= 10
        [spk_trg_avg(d,:)] = get_event_trig_avg(wcv,synct1id(spkid),maxlag,maxlag);
%         [spk_trg_avg(d,:),spk_trg_mat] = get_event_trig_avg(wcv,synct1id(spkid),maxlag,maxlag);
%         [spk_trg_davg,spk_trg_dmat] = get_event_trig_avg(wcv_d,synct1id(spkid),maxlag,maxlag);
    else
        spk_trg_avg(d,:) = nan;
    end
    if length(spkid) >= 10
        spk_trg_avg_beg(d,:) = get_event_trig_avg(wcv,synct1id(spkid(beg_spkids)),maxlag,maxlag);
    else
        spk_trg_avg_beg(d,:) = nan;
    end
    
%     %compute spike widths
%     spk_trg_mat_norm = bsxfun(@rdivide,spk_trg_mat,max(spk_trg_mat,[],2));
    spk_beg = nan(length(spkid),1);
    spk_end = nan(length(spkid),1);
    spk_peak = nan(length(spkid),1);
%     for i = 1:length(spkid)
% %         spk_peak(i) = maxlag - 1 + find(spk_trg_mat_norm(i,maxlag:end) > 0.5,1,'first');
% %         spk_end(i) = spk_peak(i) - 1 + find(spk_trg_mat_norm(i,spk_peak(i):end) < 0.5,1,'first');
% %         spk_beg(i) = find(spk_trg_mat_norm(i,1:spk_end(i)-1) < 0.5,1,'last');
%         spk_peak(i) = maxlag-1+find(spk_trg_dmat(i,maxlag:end) < 0,1,'first');
%         cur_seg = spk_trg_mat(i,:)/abs(spk_trg_mat(i,spk_peak(i)));
%         cur = spk_peak(i) - 1 + find(cur_seg(spk_peak(i):end) < 0.5,1,'first');
%         if ~isempty(cur)
%             spk_end(i) = cur;
%             spk_beg(i) = find(cur_seg(1:spk_end(i)-1) < 0.5,1,'last');
%         end
%     end
    spk_widths = (spk_end - spk_beg)/Fs;
    med_spk_widths(d) = nanmedian(spk_widths);
    mean_spk_widths(d) = nanmean(spk_widths);
    med_spk_width_beg(d) = nanmedian(spk_widths(beg_spkids));
    mean_spk_widths_beg(d) = nanmean(spk_widths(beg_spkids));
    
    spk_times = synct(spkid)/1e6;
    spk_taxis = (synct(1)/1e6):acorr_dt:(synct(end)/1e6); 
    spk_bin = hist(spk_times,spk_taxis);
    spk_bin(end) = 0;
    
%     [spk_acorr(d,:),lags] = xcov(spk_bin,acorr_maxlag,'coeff');
%     zp = find(lags == 0);
%     spk_acorr(d,zp) = 0;
%     
%     isis = diff(spk_times);
%     isi_dist(d,:) = histc(isis,isi_range);
%     isi_dist(d,:) = isi_dist(d,:)/sum(isi_dist(d,:));
    
    mean_rate(d) = length(spkid)/range(spk_taxis);
    
    
end

cd C:\WC_Germany\sven_thomas_combined
save spk_trg_avgs_60hz_v2 spk_trg_avg* *spk_width*
%%
% figure
% hold on
% shadedErrorBar(lags*acorr_dt,mean(spk_acorr(l3mec,:)),std(spk_acorr(l3mec,:))/sqrt(length(l3mec)),{'r'});
% shadedErrorBar(lags*acorr_dt,mean(spk_acorr(l3lec,:)),std(spk_acorr(l3lec,:))/sqrt(length(l3lec)),{'b'});
% shadedErrorBar(lags*acorr_dt,mean(spk_acorr(l3mec_np,:)),std(spk_acorr(l3mec_np,:))/sqrt(length(l3mec_np)),{'k'});
% shadedErrorBar(lags*acorr_dt,mean(spk_acorr(l3lec_np,:)),std(spk_acorr(l3lec_np,:))/sqrt(length(l3lec_np)),{'g'});
% 
% figure
% hold on
% shadedErrorBar(isi_range,nanmean(isi_dist(l3mec,:)),std(isi_dist(l3mec,:))/sqrt(length(l3mec)),{'r'});
% shadedErrorBar(isi_range,nanmean(isi_dist(l3lec,:)),std(isi_dist(l3lec,:))/sqrt(length(l3lec)),{'b'});
% shadedErrorBar(isi_range,nanmean(isi_dist(l3mec_np,:)),nanstd(isi_dist(l3mec_np,:))/sqrt(length(l3mec_np)),{'k'});
% shadedErrorBar(isi_range,nanmean(isi_dist(l3lec_np,:)),nanstd(isi_dist(l3lec_np,:))/sqrt(length(l3lec_np)),{'g'});

spk_trg_avg_norm = bsxfun(@rdivide,spk_trg_avg,max(spk_trg_avg,[],2));
spk_trg_avg_norm_beg = bsxfun(@rdivide,spk_trg_avg_beg,max(spk_trg_avg_beg,[],2));
lags = (-maxlag:maxlag)/Fs;


figure
subplot(2,1,1)
hold on
plot(lags,nanmean(spk_trg_avg_norm(l3mec,:)),'r')
plot(lags,nanmean(spk_trg_avg_norm(l3mec_np,:)),'k')
legend('Pyr','Non-pyr')
shadedErrorBar(lags,nanmean(spk_trg_avg_norm(l3mec,:)),std(spk_trg_avg_norm(l3mec,:))/sqrt(length(l3mec)),{'r'});
shadedErrorBar(lags,nanmean(spk_trg_avg_norm(l3mec_np,:)),std(spk_trg_avg_norm(l3mec_np,:))/sqrt(length(l3mec_np)),{'k'});
plot(lags,spk_trg_avg_norm(int,:),'g','linewidth',2)
title('MEC','fontsize',16)
xlim([-0.02 0.3])
xlabel('Time (s)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
ylim([-0.1 1])
% figure
subplot(2,1,2)
hold on
plot(lags,nanmean(spk_trg_avg_norm(l3lec,:)),'b')
plot(lags,nanmean(spk_trg_avg_norm(l3lec_np,:)),'y')
legend('Pyr','Non-pyr')
% shadedErrorBar(lags,nanmean(spk_trg_avg_norm(l3lec,:)),std(spk_trg_avg_norm(l3lec,:))/sqrt(length(l3lec)),{'b'});
% shadedErrorBar(lags,nanmean(spk_trg_avg_norm(l3lec_np,:)),std(spk_trg_avg_norm(l3lec_np,:))/sqrt(length(l3lec_np)),{'g'});
title('LEC','fontsize',16)
xlim([-0.02 0.3])
xlabel('Time (s)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
ylim([-0.1 1])

%%
spk_trg_avg_norm = bsxfun(@rdivide,spk_trg_avg,max(spk_trg_avg,[],2));
spk_trg_avg_norm_beg = bsxfun(@rdivide,spk_trg_avg_beg,max(spk_trg_avg_beg,[],2));
lags = (-maxlag:maxlag)/Fs;
xl = [-0.05 0.1];
ulags = find(lags >= -0.05 & lags <= 0.35);
figure
subplot(2,1,1)
hold on
plot(lags(ulags),nanmean(spk_trg_avg_norm_beg(l3mec,ulags)),'r')
plot(lags(ulags),nanmean(spk_trg_avg_norm_beg(l3mec_np,ulags)),'k')
legend('Pyr','Non-pyr')
% shadedErrorBar(lags(ulags),nanmean(spk_trg_avg_norm_beg(l3mec,ulags)),std(spk_trg_avg_norm_beg(l3mec,ulags))/sqrt(length(l3mec)),{'r'});
% shadedErrorBar(lags(ulags),nanmean(spk_trg_avg_norm_beg(l3mec_np,ulags)),std(spk_trg_avg_norm_beg(l3mec_np,ulags))/sqrt(length(l3mec_np)),{'k'});
% plot(lags,spk_trg_avg_norm(int,:),'g','linewidth',2)
title('MEC','fontsize',16)
xlim(xl)
xlabel('Time (s)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
ylim([-0.1 1])
% figure
subplot(2,1,2)
hold on
plot(lags(ulags),nanmean(spk_trg_avg_norm_beg(l3lec,ulags)),'b')
plot(lags(ulags),nanmean(spk_trg_avg_norm_beg(l3lec_np,ulags)),'g')
legend('Pyr','Non-pyr')
% shadedErrorBar(lags(ulags),nanmean(spk_trg_avg_norm_beg(l3lec,ulags)),std(spk_trg_avg_norm_beg(l3lec,ulags))/sqrt(length(l3lec)),{'b'});
% shadedErrorBar(lags(ulags),nanmean(spk_trg_avg_norm_beg(l3lec_np,ulags)),std(spk_trg_avg_norm_beg(l3lec_np,ulags))/sqrt(length(l3lec_np)),{'y'});
title('LEC','fontsize',16)
xlim(xl)
xlabel('Time (s)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
ylim([-0.1 1])

%%
for i = 1:size(spk_trg_avg)
    if ~isnan(spk_trg_avg_norm(i,1));
        start = find(spk_trg_avg_norm(i,:) > 0.5,1,'first');
        last = find(spk_trg_avg_norm(i,:) > 0.5,1,'last');
        spk_width(i) = (last - start)/Fs;
    else
        spk_width(i) = nan;
    end
    if ~isnan(spk_trg_avg_norm_beg(i,1));
        start = find(spk_trg_avg_norm_beg(i,:) > 0.5,1,'first');
        last = find(spk_trg_avg_norm_beg(i,:) > 0.5,1,'last');
        spk_width_beg(i) = (last - start)/Fs;
    else
        spk_width_beg(i) = nan;
    end
end
spk_amps = max(spk_trg_avg,[],2);
spk_amps_beg = max(spk_trg_avg_beg,[],2);

[use,ord] = sort([l3mec l3mec_np l3lec l3lec_np]);

% Y = spk_width(use)*1e3;
Y = spk_width_beg(use)*1e3;
% Y = mean_spk_widths_beg(use)*1e3;
% Y = med_spk_width_beg(use)*1e3;
G = [ones(length(l3mec),1); 2*ones(length(l3mec_np),1); 4*ones(length(l3lec),1); 5*ones(length(l3lec_np),1)];
G = G(ord);

figure
boxplot(Y,G)
set(gca,'fontsize',14,'fontname','arial')
ylabel('Spike width (ms)','fontsize',16)