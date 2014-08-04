
%% Example MEC MP trace
clear all
close all
load C:\WC_Germany\persistent_revised\pers_revised_dir

%specify session number
d = 14;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 1;
Fsd = 2016/dsf;

cd(dir_array{d})
pwd

load used_data lf8 wcv wcv_minus_spike

lf8_f = filtfilt(b,a,lf8);
% wcv_f = filtfilt(b,a,wcv_minus_spike);

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv,dsf);

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);

wcv = zscore(wcv);
lf8 = zscore(lf8);
t = (1:length(wcv_d))/2016*dsf;

plot(t-66,wcv_d,'linewidth',2)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('MEC Membrane Potential','FontSize',16)
xlim([66-66 82-66])
ylim([-3 8])

figure
plot(t-66,lf8_d,'r','linewidth',3)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('Cortical LFP','FontSize',16)
legend('MP','LFP')
xlim([66-66 82-66])
ylim([-2 3])

% used_spikes = find(spike_ids/Fsd > 66 & spike_ids/Fsd < 82);
% % plot(t(spike_ids)-66,ones(size(spike_ids))*3,'ro')
% line_bottom = 2.8;
% line_top = 3.2;
% for i = 1:length(used_spikes)
%    line([t(spike_ids(used_spikes(i)))-66 t(spike_ids(used_spikes(i)))-66],...
%        [line_bottom line_top],'Color','k','linewidth',2)
% end
shg

%% Example LEC MP trace
clear all
close all
load C:\WC_Germany\persistent_revised\pers_revised_dir

%specify session number
d = 27;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 1;
Fsd = 2016/dsf;

cd(dir_array{d})
pwd

load used_data lf8 wcv wcv_minus_spike

lf8_f = filtfilt(b,a,lf8);

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv,dsf);

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);

wcv = zscore(wcv);
lf8 = zscore(lf8);
t = (1:length(wcv_d))/2016*dsf;

plot(t,wcv_d,'linewidth',2)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('MEC Membrane Potential','FontSize',16)
% xlim([214 230])
% ylim([-3 8])
% xlim([970 986])
% xlim([430 446])
xlim([633 649])

figure
plot(t,lf8_d,'r','linewidth',3)
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (zscore)','FontSize',14)
title('Cortical LFP','FontSize',16)
legend('MP','LFP')
% xlim([214 230])
% ylim([-2 3])
% xlim([430 446])
xlim([633 649])

%% All Distributions MP

clear all
close all
load C:\WC_Germany\persistent_revised\amp_dist
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

wcv_dist = amp_dist_f;
wcv_range = amp_range_f;
lf8_dist = lf8_dist_f;
lf8_range = amp_range_f;

%for mean and std error
m_dist_mec = mean(wcv_dist(mec_cells,:));
u_dist_mec = m_dist_mec+1*std(wcv_dist(mec_cells,:))/sqrt(length(mec_cells));
l_dist_mec = m_dist_mec-1*std(wcv_dist(mec_cells,:))/sqrt(length(mec_cells));

m_dist_lec = mean(wcv_dist(lec_cells,:));
u_dist_lec = m_dist_lec+1*std(wcv_dist(lec_cells,:))/sqrt(length(lec_cells));
l_dist_lec = m_dist_lec-1*std(wcv_dist(lec_cells,:))/sqrt(length(lec_cells));

m_dist_lf8 = mean(lf8_dist([mec_cells lec_cells],:));
u_dist_lf8 = m_dist_lf8+1*std(lf8_dist([mec_cells lec_cells],:))/sqrt((length(lec_cells)+length(mec_cells)));
l_dist_lf8 = m_dist_lf8-1*std(lf8_dist([mec_cells lec_cells],:))/sqrt((length(lec_cells)+length(mec_cells)));

figure
plot(wcv_range,m_dist_mec,'linewidth',2)
hold on
plot(wcv_range,m_dist_lec,'g','linewidth',2)
plot(lf8_range,m_dist_lf8,'r','linewidth',2)
xlabel('Amplitude (z)','FontSize',14)
ylabel('Probability Density','FontSize',14)
xlim([-3 3])
ylim([0 1.0])
legend('MEC','LEC','LF8')

begp = find(wcv_range < -3,1,'last');
endp = find(wcv_range > 3,1,'first');

X = [wcv_range(begp:endp) fliplr(wcv_range(begp:endp))];
Y = [u_dist_mec(begp:endp) fliplr(l_dist_mec(begp:endp))];
fill(X,Y,'b')
X = [wcv_range(begp:endp) fliplr(wcv_range(begp:endp))];
Y = [u_dist_lec(begp:endp) fliplr(l_dist_lec(begp:endp))];
fill(X,Y,'g')
X = [lf8_range(begp:endp) fliplr(lf8_range(begp:endp))];
Y = [u_dist_lf8(begp:endp) fliplr(l_dist_lf8(begp:endp))];
fill(X,Y,'r')
xlabel('Amplitude (z)')
ylabel('Probability')

%% Overall Power Spectra
clear all
close all
load C:\WC_Germany\persistent_revised\power_spec_data
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

max_fr = 1;
max_f = find(f > max_fr,1,'first');

%convert to log power
Sw = 10*log10(Sw);
S8 = 10*log10(S8);

mean_S8 = mean(S8([mec_cells lec_cells],:));
u_S8 = mean_S8 + 1*std(S8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));
l_S8 = mean_S8 - 1*std(S8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));

mean_Sw_mec = mean(Sw(mec_cells,:));
u_Sw_mec = mean_Sw_mec + 1*std(Sw(mec_cells,:))/sqrt(length(mec_cells));
l_Sw_mec = mean_Sw_mec - 1*std(Sw(mec_cells,:))/sqrt(length(mec_cells));

mean_Sw_lec = mean(Sw(lec_cells,:));
u_Sw_lec = mean_Sw_lec + 1*std(Sw(lec_cells,:))/sqrt(length(lec_cells));
l_Sw_lec = mean_Sw_lec - 1*std(Sw(lec_cells,:))/sqrt(length(lec_cells));

plot(f(1:max_f),mean_S8(1:max_f),'r','linewidth',2)
hold on
plot(f(1:max_f),mean_Sw_mec(1:max_f),'linewidth',2)
plot(f(1:max_f),mean_Sw_lec(1:max_f),'g','linewidth',2)

xlim([0 max_fr])
X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_S8(1:max_f) fliplr(l_S8(1:max_f))];
fill(X,Y,'r')
X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_Sw_mec(1:max_f) fliplr(l_Sw_mec(1:max_f))];
fill(X,Y,'b')
X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_Sw_lec(1:max_f) fliplr(l_Sw_lec(1:max_f))];
fill(X,Y,'g')

xlabel('Frequency (Hz)')
ylabel('Power (dB)')

%% Overall Coherence
clear all
close all
load C:\WC_Germany\persistent_revised\coherence\coherence_data
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

max_fr = 1; %maximum frequency of interest in Hz
max_f = find(f > max_fr,1,'first');

mean_coh_mec = mean(Cmn(mec_cells,:));
u_coh_mec = mean_coh_mec+1*std(Cmn(mec_cells,:))/sqrt(length(mec_cells));
l_coh_mec = mean_coh_mec-1*std(Cmn(mec_cells,:))/sqrt(length(mec_cells));

mean_coh_lec = mean(Cmn(lec_cells,:));
u_coh_lec = mean_coh_lec+1*std(Cmn(lec_cells,:))/sqrt(length(lec_cells));
l_coh_lec = mean_coh_lec-1*std(Cmn(lec_cells,:))/sqrt(length(lec_cells));

% for i = 1:17
%    sig_check(i,:) = Cerr(i,1,:) > ConfC(i);
% end
% frac_sig = sum(sig_check)/17;

plot(f(1:max_f),mean_coh_mec(1:max_f),'linewidth',2)
hold on
plot(f(1:max_f),mean_coh_lec(1:max_f),'g','linewidth',2)
legend('MEC','LEC')

X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_coh_mec(1:max_f) fliplr(l_coh_mec(1:max_f))];
fill(X,Y,'b')
X = [f(1:max_f) fliplr(f(1:max_f))];
Y = [u_coh_lec(1:max_f) fliplr(l_coh_lec(1:max_f))];
fill(X,Y,'g')

xlim([0 max_fr])
xlabel('Frequency (Hz)','FontSize',14)
ylim([0 1])

%% Up State Distribution Compare
clear all
close all
load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\state_dur_hist_data

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

mean_up_mec = mean(up_hist(mec_cells,:));
u_up_mec = mean_up_mec + 1*std(up_hist(mec_cells,:))/sqrt(length(mec_cells));
l_up_mec = mean_up_mec - 1*std(up_hist(mec_cells,:))/sqrt(length(mec_cells));

mean_up_lfp = mean(up_hist8([mec_cells lec_cells],:));
u_up_lfp = mean_up_lfp + 1*std(up_hist8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));
l_up_lfp = mean_up_lfp - 1*std(up_hist8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));

mean_up_lec = mean(up_hist(lec_cells,:));
u_up_lec = mean_up_lec + 1*std(up_hist(lec_cells,:))/sqrt(length(lec_cells));
l_up_lec = mean_up_lec - 1*std(up_hist(lec_cells,:))/sqrt(length(lec_cells));

stairs(up_grid,mean_up_lfp,'r','linewidth',2)
hold on
stairs(up_grid,mean_up_mec,'linewidth',2)
stairs(up_grid,mean_up_lec,'g','linewidth',2)

legend('LFP','MEC','LEC')

epsilon = 1e-7;
u_up_lfp(u_up_lfp < epsilon) = epsilon;
l_up_lfp(l_up_lfp < epsilon) = epsilon;
mean_up_lfp(mean_up_lfp < epsilon) = epsilon;
u_up_mec(u_up_mec < epsilon) = epsilon;
l_up_mec(l_up_mec < epsilon) = epsilon;
mean_up_mec(mean_up_mec < epsilon) = epsilon;
u_up_lec(u_up_lec < epsilon) = epsilon;
l_up_lec(l_up_lec < epsilon) = epsilon;
mean_up_lec(mean_up_lec < epsilon) = epsilon;

for i = 1:length(up_grid)-1
    clear X Y
    bin_width = (up_grid(i+1)-up_grid(i));
    X = [up_grid(i) up_grid(i)+bin_width];
    Y = [u_up_lec(i) u_up_lec(i)];
    X = [X fliplr(X)];
    Y = [Y l_up_lec(i) l_up_lec(i)];
    fill(X,Y,'g','EdgeColor','none');
    hold on
end

for i = 1:length(up_grid)-1
    clear X Y
    bin_width = (up_grid(i+1)-up_grid(i));
    X = [up_grid(i) up_grid(i)+bin_width];
    Y = [u_up_lfp(i) u_up_lfp(i)];
    X = [X fliplr(X)];
    Y = [Y l_up_lfp(i) l_up_lfp(i)];
    fill(X,Y,'r','EdgeColor','none');
    hold on
end

for i = 1:length(up_grid)-1
    clear X Y
    bin_width = (up_grid(i+1)-up_grid(i));
    X = [up_grid(i) up_grid(i)+bin_width];
    Y = [u_up_mec(i) u_up_mec(i)];
    X = [X fliplr(X)];
    Y = [Y l_up_mec(i) l_up_mec(i)];
    fill(X,Y,'b','EdgeColor','none');
    hold on
end


set(gca,'yscale','log')
set(gca,'YMinorGrid','off')
set(gca,'YGrid','off')
% set(gca,'XMinorGrid','off')
ylim([1e-4 0.2])
xlim([0 10])
shg
xlabel('Up State Duration (s)','FontSize',14)
ylabel('Percent of Data','FontSize',14)

%% Down State Distribution Compare
clear all
close all
load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\state_dur_hist_data

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

mean_down_mec = mean(down_hist(mec_cells,:));
u_down_mec = mean_down_mec + 1*std(down_hist(mec_cells,:))/sqrt(length(mec_cells));
l_down_mec = mean_down_mec - 1*std(down_hist(mec_cells,:))/sqrt(length(mec_cells));

mean_down_lfp = mean(down_hist8([mec_cells lec_cells],:));
u_down_lfp = mean_down_lfp + 1*std(down_hist8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));
l_down_lfp = mean_down_lfp - 1*std(down_hist8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));

mean_down_lec = mean(down_hist(lec_cells,:));
u_down_lec = mean_down_lec + 1*std(down_hist(lec_cells,:))/sqrt(length(lec_cells));
l_down_lec = mean_down_lec - 1*std(down_hist(lec_cells,:))/sqrt(length(lec_cells));

stairs(down_grid,mean_down_lfp,'r','linewidth',2)
hold on
stairs(down_grid,mean_down_mec,'linewidth',2)
stairs(down_grid,mean_down_lec,'g','linewidth',2)

legend('LFP','MEC','LEC')

epsilon = 1e-7;
u_down_lfp(u_down_lfp < epsilon) = epsilon;
l_down_lfp(l_down_lfp < epsilon) = epsilon;
mean_down_lfp(mean_down_lfp < epsilon) = epsilon;
u_down_mec(u_down_mec < epsilon) = epsilon;
l_down_mec(l_down_mec < epsilon) = epsilon;
mean_down_mec(mean_down_mec < epsilon) = epsilon;
u_down_lec(u_down_lec < epsilon) = epsilon;
l_down_lec(l_down_lec < epsilon) = epsilon;
mean_down_lec(mean_down_lec < epsilon) = epsilon;

for i = 1:length(down_grid)-1
    clear X Y
    bin_width = (down_grid(i+1)-down_grid(i));
    X = [down_grid(i) down_grid(i)+bin_width];
    Y = [u_down_lec(i) u_down_lec(i)];
    X = [X fliplr(X)];
    Y = [Y l_down_lec(i) l_down_lec(i)];
    fill(X,Y,'g','EdgeColor','none');
    hold on
end

for i = 1:length(down_grid)-1
    clear X Y
    bin_width = (down_grid(i+1)-down_grid(i));
    X = [down_grid(i) down_grid(i)+bin_width];
    Y = [u_down_lfp(i) u_down_lfp(i)];
    X = [X fliplr(X)];
    Y = [Y l_down_lfp(i) l_down_lfp(i)];
    fill(X,Y,'r','EdgeColor','none');
    hold on
end

for i = 1:length(down_grid)-1
    clear X Y
    bin_width = (down_grid(i+1)-down_grid(i));
    X = [down_grid(i) down_grid(i)+bin_width];
    Y = [u_down_mec(i) u_down_mec(i)];
    X = [X fliplr(X)];
    Y = [Y l_down_mec(i) l_down_mec(i)];
    fill(X,Y,'b','EdgeColor','none');
    hold on
end


set(gca,'yscale','log')
set(gca,'YMinorGrid','off')
set(gca,'YGrid','off')
% set(gca,'XMinorGrid','off')
ylim([1e-4 0.2])
xlim([0 10])
shg
xlabel('Down State Duration (s)','FontSize',14)
ylabel('Percent of Data','FontSize',14)


%% MP Triggered UP Trans Overall
clear all
close all
load C:\WC_Germany\persistent_revised\trig_avgs\trig_avg_data_wideband
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

first_look_point = find(lags > -1,1,'first');

mec_avg_utrig = nanmean(mp_utrig_lf8(mec_cells,:));
mec_uc_utrig = mec_avg_utrig+nanstd(mp_utrig_lf8(mec_cells,:))/sqrt(length(mec_cells));
mec_lc_utrig = mec_avg_utrig-nanstd(mp_utrig_lf8(mec_cells,:))/sqrt(length(mec_cells));

[mval,mloc] = max(mec_avg_utrig);
[minval,minloc] = min(mec_avg_utrig(1:mloc));
mid_pt = (mval+minval)/2;
mid_loc_mec = lags(first_look_point+find(mec_avg_utrig(first_look_point:end) > mid_pt,1,'first'));
mid_loc_mec_u = lags(first_look_point+find(mec_uc_utrig(first_look_point:end) > mid_pt,1,'first'));
mid_loc_mec_l = lags(first_look_point+find(mec_lc_utrig(first_look_point:end) > mid_pt,1,'first'));

lec_avg_utrig = nanmean(mp_utrig_lf8(lec_cells,:));
lec_uc_utrig = lec_avg_utrig+nanstd(mp_utrig_lf8(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_utrig = lec_avg_utrig-nanstd(mp_utrig_lf8(lec_cells,:))/sqrt(length(lec_cells));

[mval,mloc] = max(lec_avg_utrig);
[minval,minloc] = min(lec_avg_utrig(1:mloc));
mid_pt = (mval+minval)/2;
mid_loc_lec = lags(first_look_point+find(lec_avg_utrig(first_look_point:end) > mid_pt,1,'first'));
mid_loc_lec_u = lags(first_look_point+find(lec_uc_utrig(first_look_point:end) > mid_pt,1,'first'));
mid_loc_lec_l = lags(first_look_point+find(lec_lc_utrig(first_look_point:end) > mid_pt,1,'first'));

figure
plot(lags,mec_avg_utrig,'linewidth',2)
hold on
plot(lags,lec_avg_utrig,'g','linewidth',2)
legend('MEC','LEC')

X = [lags fliplr(lags)];
Y = [mec_uc_utrig fliplr(mec_lc_utrig)];
fill(X,Y,'b')

Y = [lec_uc_utrig fliplr(lec_lc_utrig)];
fill(X,Y,'g')

xlim([-2 2])
% xlabel('Time Lag (s)')
ylabel('Amplitude (z)')
line([0 0],[-1 2],'Color','k')
line([-2 2],[0 0],'Color','k')


%% LF8 Triggered UP Trans Overall
clear all
close all
load C:\WC_Germany\persistent_revised\trig_avgs\lf8_trig_avg_data_wideband
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

first_look_point = find(lags > 0,1,'first');

lfp_avg_utrig = nanmean(lf8_utrig_lf8([mec_cells lec_cells],:));
lfp_uc_utrig = lfp_avg_utrig+nanstd(lf8_utrig_lf8([mec_cells lec_cells],:))/sqrt(28);
lfp_lc_utrig = lfp_avg_utrig-nanstd(lf8_utrig_lf8([mec_cells lec_cells],:))/sqrt(28);


mec_avg_utrig = nanmean(lf8_utrig_mp(mec_cells,:));
mec_uc_utrig = mec_avg_utrig+nanstd(lf8_utrig_mp(mec_cells,:))/sqrt(length(mec_cells));
mec_lc_utrig = mec_avg_utrig-nanstd(lf8_utrig_mp(mec_cells,:))/sqrt(length(mec_cells));

[mval,mloc] = max(mec_avg_utrig);
[minval,minloc] = min(mec_avg_utrig(1:mloc));
mid_pt = (mval+minval)/2;
mid_loc_mec = lags(first_look_point+find(mec_avg_utrig(first_look_point:end) > mid_pt,1,'first'));
mid_loc_mec_u = lags(first_look_point+find(mec_uc_utrig(first_look_point:end) > mid_pt,1,'first'));
mid_loc_mec_l = lags(first_look_point+find(mec_lc_utrig(first_look_point:end) > mid_pt,1,'first'));

lec_avg_utrig = nanmean(lf8_utrig_mp(lec_cells,:));
lec_uc_utrig = lec_avg_utrig+nanstd(lf8_utrig_mp(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_utrig = lec_avg_utrig-nanstd(lf8_utrig_mp(lec_cells,:))/sqrt(length(lec_cells));

[mval,mloc] = max(lec_avg_utrig);
[minval,minloc] = min(lec_avg_utrig(1:mloc));
mid_pt = (mval+minval)/2;
mid_loc_lec = lags(first_look_point+find(lec_avg_utrig(first_look_point:end) > mid_pt,1,'first'));
mid_loc_lec_u = lags(first_look_point+find(lec_uc_utrig(first_look_point:end) > mid_pt,1,'first'));
mid_loc_lec_l = lags(first_look_point+find(lec_lc_utrig(first_look_point:end) > mid_pt,1,'first'));

figure
plot(lags,mec_avg_utrig,'linewidth',2)
hold on
plot(lags,lec_avg_utrig,'g','linewidth',2)
plot(lags,lfp_avg_utrig,'r','linewidth',2)
legend('MEC','LEC','LFP')

X = [lags fliplr(lags)];
Y = [mec_uc_utrig fliplr(mec_lc_utrig)];
fill(X,Y,'b')

Y = [lec_uc_utrig fliplr(lec_lc_utrig)];
fill(X,Y,'g')
Y = [lfp_uc_utrig fliplr(lfp_lc_utrig)];
fill(X,Y,'r')

xlim([-2 2])
% xlabel('Time Lag (s)')
ylabel('Amplitude (z)')
line([0 0],[-1 1.5],'Color','k')
line([-2 2],[0 0],'Color','k')

%% MP Triggered DOWN Trans Overall
clear all
close all
load C:\WC_Germany\persistent_revised\trig_avgs\trig_avg_data_wideband
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

first_look_point = find(lags > -1,1,'first');

mec_avg_dtrig = nanmean(mp_dtrig_lf8(mec_cells,:));
mec_uc_dtrig = mec_avg_dtrig+nanstd(mp_dtrig_lf8(mec_cells,:))/sqrt(length(mec_cells));
mec_lc_dtrig = mec_avg_dtrig-nanstd(mp_dtrig_lf8(mec_cells,:))/sqrt(length(mec_cells));

[mval,mloc] = min(mec_avg_dtrig);
[maxval,maxloc] = max(mec_avg_dtrig(1:mloc));
mid_pt = (mval+maxval)/2;
mid_loc_mec = lags(first_look_point+find(mec_avg_dtrig(first_look_point:end) < mid_pt,1,'first'));
mid_loc_mec_u = lags(first_look_point+find(mec_uc_dtrig(first_look_point:end) < mid_pt,1,'first'));
mid_loc_mec_l = lags(first_look_point+find(mec_lc_dtrig(first_look_point:end) < mid_pt,1,'first'));

lec_avg_dtrig = nanmean(mp_dtrig_lf8(lec_cells,:));
lec_uc_dtrig = lec_avg_dtrig+nanstd(mp_dtrig_lf8(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_dtrig = lec_avg_dtrig-nanstd(mp_dtrig_lf8(lec_cells,:))/sqrt(length(lec_cells));

[mval,mloc] = min(lec_avg_dtrig);
[maxval,maxloc] = max(lec_avg_dtrig(1:mloc));
mid_pt = (mval+maxval)/2;
mid_loc_lec = lags(first_look_point+find(lec_avg_dtrig(first_look_point:end) < mid_pt,1,'first'));
mid_loc_lec_u = lags(first_look_point+find(lec_uc_dtrig(first_look_point:end) < mid_pt,1,'first'));
mid_loc_lec_l = lags(first_look_point+find(lec_lc_dtrig(first_look_point:end) < mid_pt,1,'first'));

figure
plot(lags,mec_avg_dtrig,'linewidth',2)
hold on
plot(lags,lec_avg_dtrig,'g','linewidth',2)
legend('MEC','LEC')

X = [lags fliplr(lags)];
Y = [mec_uc_dtrig fliplr(mec_lc_dtrig)];
fill(X,Y,'b')

Y = [lec_uc_dtrig fliplr(lec_lc_dtrig)];
fill(X,Y,'g')

ylim([-1 1.5])
xlim([-2 2])
xlabel('Time Lag (s)')
ylabel('Amplitude (z)')
line([0 0],[-1 2],'Color','k')
line([-2 2],[0 0],'Color','k')

%% LF8 Triggered DOWN Trans Overall
clear all
close all
load C:\WC_Germany\persistent_revised\trig_avgs\lf8_trig_avg_data_wideband
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

first_look_point = find(lags > 0,1,'first');

lfp_avg_dtrig = nanmean(lf8_dtrig_lf8([lec_cells mec_cells],:));
lfp_uc_dtrig = lfp_avg_dtrig+nanstd(lf8_dtrig_lf8([lec_cells mec_cells],:))/sqrt(28);
lfp_lc_dtrig = lfp_avg_dtrig-nanstd(lf8_dtrig_lf8([lec_cells mec_cells],:))/sqrt(28);


mec_avg_dtrig = nanmean(lf8_dtrig_mp(mec_cells,:));
mec_uc_dtrig = mec_avg_dtrig+nanstd(lf8_dtrig_mp(mec_cells,:))/sqrt(length(mec_cells));
mec_lc_dtrig = mec_avg_dtrig-nanstd(lf8_dtrig_mp(mec_cells,:))/sqrt(length(mec_cells));

[mval,mloc] = min(mec_avg_dtrig);
[maxval,maxloc] = max(mec_avg_dtrig(1:mloc));
mid_pt = (mval+maxval)/2;
mid_loc_mec = lags(first_look_point+find(mec_avg_dtrig(first_look_point:end) < mid_pt,1,'first'));
mid_loc_mec_u = lags(first_look_point+find(mec_uc_dtrig(first_look_point:end) < mid_pt,1,'first'));
mid_loc_mec_l = lags(first_look_point+find(mec_lc_dtrig(first_look_point:end) < mid_pt,1,'first'));

lec_avg_dtrig = nanmean(lf8_dtrig_mp(lec_cells,:));
lec_uc_dtrig = lec_avg_dtrig+nanstd(lf8_dtrig_mp(lec_cells,:))/sqrt(length(lec_cells));
lec_lc_dtrig = lec_avg_dtrig-nanstd(lf8_dtrig_mp(lec_cells,:))/sqrt(length(lec_cells));

[mval,mloc] = min(lec_avg_dtrig);
[maxval,maxloc] = max(lec_avg_dtrig(1:mloc));
mid_pt = (mval+maxval)/2;
mid_loc_lec = lags(first_look_point+find(lec_avg_dtrig(first_look_point:end) < mid_pt,1,'first'));
mid_loc_lec_u = lags(first_look_point+find(lec_uc_dtrig(first_look_point:end) < mid_pt,1,'first'));
mid_loc_lec_l = lags(first_look_point+find(lec_lc_dtrig(first_look_point:end) < mid_pt,1,'first'));

figure
plot(lags,mec_avg_dtrig,'linewidth',2)
hold on
plot(lags,lec_avg_dtrig,'g','linewidth',2)
plot(lags,lfp_avg_dtrig,'r','linewidth',2)
legend('MEC','LEC')

X = [lags fliplr(lags)];
Y = [mec_uc_dtrig fliplr(mec_lc_dtrig)];
fill(X,Y,'b')

Y = [lec_uc_dtrig fliplr(lec_lc_dtrig)];
fill(X,Y,'g')
Y = [lfp_uc_dtrig fliplr(lfp_lc_dtrig)];
fill(X,Y,'r')

ylim([-1 1.5])
xlim([-2 2])
xlabel('Time Lag (s)')
ylabel('Amplitude (z)')
line([0 0],[-1 2],'Color','k')
line([-2 2],[0 0],'Color','k')


%% Quantization Figure Overall
clear all
close all
load C:\WC_Germany\persistent_revised\corresponding_lfp_state_data
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
mec_cells(21) = [];


hist_range = linspace(0,7,1000);
binsize = hist_range(2)-hist_range(1);
for i = 1:n_mec+n_lec
    cell_hist(i,:) = hist(mp_updur_corresp_lfp{i},hist_range);
    cell_hist(i,:) = cell_hist(i,:)/sum(cell_hist(i,:))/binsize;
    cell_hist(i,:) = jmm_smooth_1d_cor(cell_hist(i,:),25);
end

mean_hist_mec = mean(cell_hist(mec_cells,:));
u_mec = mean_hist_mec+1*std(cell_hist(mec_cells,:))/sqrt(length(mec_cells));
l_mec = mean_hist_mec-1*std(cell_hist(mec_cells,:))/sqrt(length(mec_cells));

mean_hist_lec = mean(cell_hist(lec_cells,:));
u_lec = mean_hist_lec+1*std(cell_hist(lec_cells,:))/sqrt(length(lec_cells));
l_lec = mean_hist_lec-1*std(cell_hist(lec_cells,:))/sqrt(length(lec_cells));

plot(hist_range,mean_hist_mec,'linewidth',2)
hold on
plot(hist_range,mean_hist_lec,'g','linewidth',2)
legend('MEC','LEC')

epsilon = 1e-10;
cell_hist(cell_hist < epsilon) = epsilon;
u_mec(u_mec < epsilon) = epsilon;
l_mec(l_mec < epsilon) = epsilon;
u_lec(u_lec < epsilon) = epsilon;
l_lec(l_lec < epsilon) = epsilon;

% X = [hist_range fliplr(hist_range)];
% Y = [u_mec fliplr(l_mec)];
% fill(X,Y,'b')
% Y = [u_lec fliplr(l_lec)];
% fill(X,Y,'g')
% 
set(gca,'yscale','log')
xlim([0 6.2])
ylim([1e-3 3])

[a,b] = findpeaks(mean_hist_mec);
peak_locs = hist_range(b);

% for i = 1:length(peak_locs)
%     line([peak_locs(i) peak_locs(i)],[5e-4 2],'Color','k')
% end
xlabel('MP Up state duration (LFP cycles)')
ylabel('Probability')
%% Quantization Figure Single Example
clear all
close all
load C:\WC_Germany\persistent_revised\corresponding_lfp_state_data
load C:\WC_Germany\persistent_revised\pers_revised_dir

%good examples are 1,2,5,8,14,15,19
hist_range = linspace(0,7,1000);
binsize = hist_range(2)-hist_range(1);
for i = 1:n_mec+n_lec
    cell_hist(i,:) = hist(mp_updur_corresp_lfp{i},hist_range);
    cell_hist(i,:) = cell_hist(i,:)/sum(cell_hist(i,:))/binsize;
    cell_hist(i,:) = jmm_smooth_1d_cor(cell_hist(i,:),25);
end

examp_ind = 5;
x = linspace(0,7,100);
binsize = x(2)-x(1);

y = hist(mp_updur_corresp_lfp{examp_ind},x);
sy = sum(y);

figure
bar(x,y/sy/binsize,1,'k')
hold on
temp = jmm_smooth_1d_cor(y,2);
temp = temp/sy;
% fine_axis = linspace(0,7,1000);
% smooth_temp = spline(x,temp,fine_axis);
% plot(x,jmm_smooth_1d_cor(y,2),'r','linewidth',2)
plot(hist_range,cell_hist(examp_ind,:),'r','linewidth',2)
xlim([0 4.2])
% ylim([0 0.15])
shg


%% Histogram of Persistent Activity
clear all
close all
load C:\WC_Germany\persistent_revised\corresponding_lfp_state_data
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

numBins = 20;
% pers_range = linspace(0,0.6,numBins);
pers_range = linspace(0,1.,numBins);

pers_mec = histc(pers_fract_within_cycle(mec_cells),pers_range);
pers_lec = histc(pers_fract_within_cycle(lec_cells),pers_range);
pers_mec_np = histc(pers_fract_within_np(mec_cells),pers_range);

stairs(pers_range*100,pers_mec,'linewidth',2)
hold on
stairs(pers_range*100,pers_lec,'g','linewidth',2)
stairs(pers_range*100,pers_mec_np,'c','linewidth',2)
% xlim([0 65])
ylim([0 5])
xlabel('Fraction of Up States Persisting','FontSize',14)
ylabel('Number of Cells','FontSize',14)
legend('MEC','LEC')

%% WCV Up Trig Mat example
clear all
close all

%%for new method
load C:\WC_Germany\persistent_revised\trig_avgs\trig_avg_single_examp
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data
load C:\WC_Germany\persistent_revised\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\persistent_revised\pers_revised_dir

d = 14;
dsf = 8;
Fsd = 2016/dsf;

cur_mp_mat = mp_utrig_mp_mat;
cur_lfp_mat = mp_utrig_lf8_mat;

bad_rows = 1;
cur_mp_mat(bad_rows,:) = [];
cur_lfp_mat(bad_rows,:) = [];
synch_up_dur{d}(bad_rows) = [];

[dummy,up_order] = sort(synch_up_dur{d});

Fig = figure(1)
clf
set(Fig,'PaperUnits','centimeters');
set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
imagesc(lags,(1:length(synch_ups{d})),(cur_mp_mat(up_order,:)));shading flat;
caxis([-3 3]);colorbar
xlim([-2 10])
xlabel('Time (s)','FontSize',14)
ylabel('MP Up Duration Percentile','FontSize',14)

Fig2 = figure(2)
clf
set(Fig2,'PaperUnits','centimeters');
set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
set(Fig2,'PaperPosition',[0,0,(get(Fig2,'PaperSize'))])
imagesc(lags,(1:length(synch_ups{d})),(cur_lfp_mat(up_order,:)));shading flat;
caxis([-2 4]);colorbar
xlim([-2 10])

figure(3)
plot(synch_up_dur{d}(up_order),1-(1:length(synch_up_dur{d}))/length(synch_up_dur{d}),'r','linewidth',2)
line([0 0],[0 1],'Color','k')
xlim([-2 10])
xlabel('Time (s)','FontSize',14)
ylabel('MP Up Duration Percentile','FontSize',14)

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 8;
Fsd = 2016/dsf;

cd(dir_array{d})
pwd

load used_data lf8 wcv_minus_spike
load spike_time_jmm
spike_ids = round(spkid/dsf);

lf8_f = filtfilt(b,a,lf8);
wcv_f = filtfilt(b,a,wcv_minus_spike);

lf8_d = downsample(lf8_f,dsf);
wcv_d = downsample(wcv_f,dsf);

lf8_d = zscore(lf8_d);
wcv_d = zscore(wcv_d);
t = (1:length(lf8_d))/Fsd;

% wcv_d(wcv_d > 2) = 2;

%at 17%
trans_num = 24;
prev_lfp_tran = find(up_trans8{d}<up_trans{d}(trans_num),1,'last');
% cur_lag = up_trans{d}(trans_num)-up_trans8{d}(prev_lfp_tran);
st = up_trans{d}(trans_num)/Fsd-4;
et = st+9;
begpt = find(t > 100.8,1,'first');
endpt = find(t > 103.5,1,'first');
figure
plot(t(begpt:endpt),wcv_d(begpt:endpt),'linewidth',3)
hold on
plot(t(begpt:endpt),lf8_d(begpt:endpt)-1,'r','linewidth',3)
hold on
xlim([100.8 103.5])
ylim([-3 3])
used_spikes = find(spike_ids/Fsd > 100.8 & spike_ids/Fsd < 103.5);
line_bottom = 2.2;
line_top = 2.6;
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
       [line_bottom line_top],'Color','k','linewidth',2)
end


%at 60%
trans_num = 94;
prev_lfp_tran = find(up_trans8{d}<up_trans{d}(trans_num),1,'last');
cur_lag = up_trans{d}(trans_num)-up_trans8{d}(prev_lfp_tran);
st = up_trans{d}(trans_num)/Fsd-3;
et = st+10;
begpt = find(t > 495,1,'first');
endpt = find(t > 499.3,1,'first');
figure
plot(t(begpt:endpt),wcv_d(begpt:endpt),'linewidth',2)
hold on
plot(t(begpt:endpt),lf8_d(begpt:endpt)-1,'r','linewidth',3)
xlim([495 499.3])
ylim([-3 3])
used_spikes = find(spike_ids/Fsd > 495 & spike_ids/Fsd < 499.3);
line_bottom = 2.2;
line_top = 2.6;
for i = 1:length(used_spikes)
    line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
        [line_bottom line_top],'Color','k','linewidth',2)
end


%at 80%
trans_num = 31;
prev_lfp_tran = find(up_trans8{d}<up_trans{d}(trans_num),1,'last');
cur_lag = up_trans{d}(trans_num)-up_trans8{d}(prev_lfp_tran);
st = up_trans{d}(trans_num)/Fsd-6;
et = st+15;
figure
begpt = find(t > 125.3,1,'first');
endpt = find(t > 131.1,1,'first');
plot(t(begpt:endpt),wcv_d(begpt:endpt),'linewidth',2)
hold on
plot(t(begpt:endpt),lf8_d(begpt:endpt)-1,'r','linewidth',3)
xlim([125.3 131.1])
ylim([-3 3])
used_spikes = find(spike_ids/Fsd > 125.3 & spike_ids/Fsd < 131.1);
line_bottom = 2.2;
line_top = 2.6;
for i = 1:length(used_spikes)
    line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
        [line_bottom line_top],'Color','k','linewidth',2)
end



%at 95%
trans_num = 116;
prev_lfp_tran = find(up_trans8{d}<up_trans{d}(trans_num),1,'last');
cur_lag = up_trans{d}(trans_num)-up_trans8{d}(prev_lfp_tran);
st = up_trans{d}(trans_num)/Fsd-6;
et = st+15;
begpt = find(t > 609.1,1,'first');
endpt = find(t > 609.1+8,1,'first');
figure
plot(t(begpt:endpt)-609.1,wcv_d(begpt:endpt),'linewidth',2)
hold on
plot(t(begpt:endpt)-609.1,lf8_d(begpt:endpt)-1,'r','linewidth',3)
xlim([0 7.5])
ylim([-3 3])
used_spikes = find(spike_ids/Fsd > 609.1 & spike_ids/Fsd < 616.4);
line_bottom = 2.2;
line_top = 2.6;
for i = 1:length(used_spikes)
   line([t(spike_ids(used_spikes(i))) t(spike_ids(used_spikes(i)))],...
       [line_bottom line_top],'Color','k','linewidth',2)
end


%% UP TRAN PHASE
clear all
close all
load C:\WC_Germany\persistent_revised\lag_phase_data
load C:\WC_Germany\persistent_revised\pers_revised_dir

lf8_down_phase = 144;

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

phase_range = linspace(0,360,1000);
dphase = phase_range(2)-phase_range(1);
for d = 1:28
    phase_hist(d,:) = histc(mp_up_phase{d},phase_range);
    phase_hist(d,:) = phase_hist(d,:)/sum(phase_hist(d,:))/dphase;
    sm_phase_hist(d,:) = [phase_hist(d,:) phase_hist(d,:)];
    sm_phase_hist(d,:) = jmm_smooth_1d_cor(sm_phase_hist(d,:),30);
%     phase_hist(d,:) = fgsmooth(phase_hist(d,:),4);
end

mean_phase_mec = mean(sm_phase_hist(mec_cells,:));
u_phase_mec = mean_phase_mec + 1*std(sm_phase_hist(mec_cells,:))/sqrt(length(mec_cells));
l_phase_mec = mean_phase_mec - 1*std(sm_phase_hist(mec_cells,:))/sqrt(length(mec_cells));
mean_phase_lec = mean(sm_phase_hist(lec_cells,:));
u_phase_lec = mean_phase_lec + 1*std(sm_phase_hist(lec_cells,:))/sqrt(length(lec_cells));
l_phase_lec = mean_phase_lec - 1*std(sm_phase_hist(lec_cells,:))/sqrt(length(lec_cells));

phase_range = [phase_range (phase_range+360)];
% mean_phase_mec = [mean_phase_mec mean_phase_mec(2:end)];
% mean_phase_lec = [mean_phase_lec mean_phase_lec(2:end)];
% u_phase_mec = [u_phase_mec u_phase_mec];
% u_phase_lec = [u_phase_lec u_phase_lec];
% l_phase_mec = [l_phase_mec l_phase_mec];
% l_phase_lec = [l_phase_lec l_phase_lec];

plot(phase_range,mean_phase_mec,'linewidth',2)
hold on
plot(phase_range,mean_phase_lec,'g','linewidth',2)
legend('MEC','LEC')

X = [phase_range fliplr(phase_range)];
Y = [u_phase_lec fliplr(l_phase_lec)];
fill(X,Y,'g')
X = [phase_range fliplr(phase_range)];
Y = [u_phase_mec fliplr(l_phase_mec)];
fill(X,Y,'b')
xlim([250 610])
line([lf8_down_phase lf8_down_phase],[0 0.015],'Color','k')
line([360 360],[0 0.015],'Color','k')
line([360+lf8_down_phase 360+lf8_down_phase],[0 0.015],'Color','k')

xlabel('LFP Phase (degrees)')
ylabel('Probability')
%% DOWN TRAN PHASE
load C:\WC_Germany\persistent_revised\lag_phase_data
load C:\WC_Germany\persistent_revised\pers_revised_dir

lf8_down_phase = 144;

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

phase_range = linspace(0,360,1000);
dphase = phase_range(2)-phase_range(1);
for d = 1:28
    phase_hist(d,:) = histc(mp_down_phase{d},phase_range);
    phase_hist(d,:) = phase_hist(d,:)/sum(phase_hist(d,:))/dphase;
    sm_phase_hist(d,:) = [phase_hist(d,:) phase_hist(d,:)];
    sm_phase_hist(d,:) = jmm_smooth_1d_cor(sm_phase_hist(d,:),30);
end


mean_phase_mec = mean(sm_phase_hist(mec_cells,:));
u_phase_mec = mean_phase_mec + 1*std(sm_phase_hist(mec_cells,:))/sqrt(length(mec_cells));
l_phase_mec = mean_phase_mec - 1*std(sm_phase_hist(mec_cells,:))/sqrt(length(mec_cells));
mean_phase_lec = mean(sm_phase_hist(lec_cells,:));
u_phase_lec = mean_phase_lec + 1*std(sm_phase_hist(lec_cells,:))/sqrt(length(lec_cells));
l_phase_lec = mean_phase_lec - 1*std(sm_phase_hist(lec_cells,:))/sqrt(length(lec_cells));

phase_range = [phase_range (phase_range+360)];

plot(phase_range,mean_phase_mec,'linewidth',2)
hold on
plot(phase_range,mean_phase_lec,'g','linewidth',2)
legend('MEC','LEC')
X = [phase_range fliplr(phase_range)];
Y = [u_phase_lec fliplr(l_phase_lec)];
fill(X,Y,'g')
X = [phase_range fliplr(phase_range)];
Y = [u_phase_mec fliplr(l_phase_mec)];
fill(X,Y,'b')
xlim([30 390])
line([lf8_down_phase lf8_down_phase],[0 0.014],'Color','k')
line([360 360],[0 0.014],'Color','k')
line([360+lf8_down_phase 360+lf8_down_phase],[0 0.014],'Color','k')
xlabel('LFP Phase (degrees)')
ylabel('Probability')

%% UP TRAN LAG
clear all
close all
load C:\WC_Germany\persistent_revised\corresponding_lfp_state_data
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

lag_range = linspace(-3,3,1000);
dlag = lag_range(2)-lag_range(1);
for d = 1:28
    lag_hist(d,:) = histc(lfp_up_lag{d},lag_range);
    lag_hist(d,:) = lag_hist(d,:)/sum(lag_hist(d,:))/dlag;
        non_pers_states = find(mp_updur_corresp_lfp{d} < 1);
    lag_hist_np(d,:) = histc(lfp_up_lag{d}(non_pers_states),lag_range);
    lag_hist_np(d,:) = lag_hist_np(d,:)/sum(lag_hist_np(d,:))/dlag;
    sm_lag_hist(d,:) = jmm_smooth_1d_cor(lag_hist(d,:),10);
    sm_lag_hist_np(d,:) = jmm_smooth_1d_cor(lag_hist_np(d,:),10);
%     phase_hist(d,:) = fgsmooth(phase_hist(d,:),4);
end

mean_lag_mec = mean(sm_lag_hist(mec_cells,:));
u_lag_mec = mean_lag_mec + 1*std(sm_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
l_lag_mec = mean_lag_mec - 1*std(sm_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
mean_lag_lec = mean(sm_lag_hist(lec_cells,:));
u_lag_lec = mean_lag_lec + 1*std(sm_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
l_lag_lec = mean_lag_lec - 1*std(sm_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
mean_lag_np_mec = mean(sm_lag_hist_np(mec_cells,:));
u_lag_np_mec = mean_lag_np_mec + 1*std(sm_lag_hist_np(mec_cells,:))/sqrt(length(mec_cells));
l_lag_np_mec = mean_lag_np_mec - 1*std(sm_lag_hist_np(mec_cells,:))/sqrt(length(mec_cells));


plot(lag_range,mean_lag_mec,'linewidth',2)
hold on
plot(lag_range,mean_lag_lec,'g','linewidth',2)
plot(lag_range,mean_lag_np_mec,'c','linewidth',2)
% legend('MEC','LEC')

X = [lag_range fliplr(lag_range)];
Y = [u_lag_lec fliplr(l_lag_lec)];
fill(X,Y,'g')
X = [lag_range fliplr(lag_range)];
Y = [u_lag_mec fliplr(l_lag_mec)];
fill(X,Y,'b')
X = [lag_range fliplr(lag_range)];
Y = [u_lag_np_mec fliplr(l_lag_np_mec)];
fill(X,Y,'c')
xlim([-0.5 1.5])
ylim([0 3])
line([0 0],[0 3],'Color','k')
xlabel('Lag (s)')
ylabel('Probability')
%% DOWN TRAN LAG
clear all
close all
load C:\WC_Germany\persistent_revised\corresponding_lfp_state_data
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

lag_range = linspace(-1,5,1000);
dlag = lag_range(2)-lag_range(1);
for d = 1:28
    lag_hist(d,:) = histc(lfp_down_lag{d},lag_range);
    lag_hist(d,:) = lag_hist(d,:)/sum(lag_hist(d,:))/dlag;
    non_pers_states = find(mp_updur_corresp_lfp{d} < 1);
    lag_hist_np(d,:) = histc(lfp_down_lag{d}(non_pers_states),lag_range);
    lag_hist_np(d,:) = lag_hist_np(d,:)/sum(lag_hist_np(d,:))/dlag;
    up_lag_hist(d,:) = histc(lfp_up_lag{d},lag_range);
    up_lag_hist(d,:) = up_lag_hist(d,:)/sum(up_lag_hist(d,:))/dlag;
    sm_lag_hist(d,:) = jmm_smooth_1d_cor(lag_hist(d,:),20);
    sm_lag_hist_np(d,:) = jmm_smooth_1d_cor(lag_hist_np(d,:),20);
    sm_up_lag_hist(d,:) = jmm_smooth_1d_cor(up_lag_hist(d,:),20);
%     phase_hist(d,:) = fgsmooth(phase_hist(d,:),4);
end

mean_lag_mec = mean(sm_lag_hist(mec_cells,:));
u_lag_mec = mean_lag_mec + 1*std(sm_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
l_lag_mec = mean_lag_mec - 1*std(sm_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
mean_lag_lec = mean(sm_lag_hist(lec_cells,:));
u_lag_lec = mean_lag_lec + 1*std(sm_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
l_lag_lec = mean_lag_lec - 1*std(sm_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
mean_uplag_mec = mean(sm_up_lag_hist(mec_cells,:));
u_uplag_mec = mean_uplag_mec + 1*std(sm_up_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
l_uplag_mec = mean_uplag_mec - 1*std(sm_up_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
mean_uplag_lec = mean(sm_up_lag_hist(lec_cells,:));
u_uplag_lec = mean_uplag_lec + 1*std(sm_up_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
l_uplag_lec = mean_uplag_lec - 1*std(sm_up_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
mean_lag_np_mec = mean(sm_lag_hist_np(mec_cells,:));
u_lag_np_mec = mean_lag_np_mec + 1*std(sm_lag_hist_np(mec_cells,:))/sqrt(length(mec_cells));
l_lag_np_mec = mean_lag_np_mec - 1*std(sm_lag_hist_np(mec_cells,:))/sqrt(length(mec_cells));


plot(lag_range,mean_lag_mec,'linewidth',2)
hold on
plot(lag_range,mean_lag_lec,'g','linewidth',2)
plot(lag_range,mean_uplag_mec,'r','linewidth',2)
plot(lag_range,mean_lag_np_mec,'c','linewidth',2)

X = [lag_range fliplr(lag_range)];
Y = [u_lag_lec fliplr(l_lag_lec)];
fill(X,Y,'g')
X = [lag_range fliplr(lag_range)];
Y = [u_lag_mec fliplr(l_lag_mec)];
fill(X,Y,'b')
X = [lag_range fliplr(lag_range)];
Y = [u_uplag_mec fliplr(l_uplag_mec)];
fill(X,Y,'r')
X = [lag_range fliplr(lag_range)];
Y = [u_lag_np_mec fliplr(l_lag_np_mec)];
fill(X,Y,'c')
xlim([-0.5 2.5])
line([0 0],[0 2],'Color','k')
xlabel('Lag (s)')
ylabel('Probability')

%% Supplementary figure all MP autocorr
clear all
close all
cd C:\WC_Germany\persistent_revised
load C:\WC_Germany\persistent_revised\time_domain\time_domain_data
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

zero_lag = find(lags/Fsd >= 0,1,'first');

m_corr_mec = mean(tot_wcv_acorr(mec_cells,:));
u_corr_mec = m_corr_mec + 1*std(tot_wcv_acorr(mec_cells,:))/sqrt(length(mec_cells));
l_corr_mec = m_corr_mec - 1*std(tot_wcv_acorr(mec_cells,:))/sqrt(length(mec_cells));

[min_mec,min_loc] = min(m_corr_mec(zero_lag:end));
[mec_peaks,mec_peak_locs] = findpeaks(m_corr_mec(zero_lag+min_loc:end));
for i = 1:length(mec_peaks)
    peak_test_mec(i) = ttest(tot_wcv_acorr(mec_cells,zero_lag+min_loc+mec_peak_locs(i)));
end

m_corr_lec = mean(tot_wcv_acorr(lec_cells,:));
u_corr_lec = m_corr_lec + 1*std(tot_wcv_acorr(lec_cells,:))/sqrt(length(lec_cells));
l_corr_lec = m_corr_lec - 1*std(tot_wcv_acorr(lec_cells,:))/sqrt(length(lec_cells));

[min_lec,min_lec_loc] = min(m_corr_lec(zero_lag:end));
[peak_lec,peak_loc_lec] = max(m_corr_lec(zero_lag+min_lec_loc:end));
peak_loc_lec = peak_loc_lec + zero_lag+min_lec_loc;

m_corr_lfp = mean(tot_lf8_acorr([mec_cells lec_cells],:));
u_corr_lfp = m_corr_lfp + 1*std(tot_lf8_acorr([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));
l_corr_lfp = m_corr_lfp - 1*std(tot_lf8_acorr([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));

[min_lfp,min_lfp_loc] = min(m_corr_lfp(zero_lag:end));
[peak_lfp,peak_loc_lfp] = max(m_corr_lfp(zero_lag+min_lfp_loc:end));
peak_loc_lfp = peak_loc_lfp + zero_lag+min_lfp_loc;



plot(lags/Fsd,m_corr_mec,'linewidth',2)
hold on
plot(lags/Fsd,m_corr_lec,'g','linewidth',2)
plot(lags/Fsd,m_corr_lfp,'r','linewidth',2)

% X = [lags/Fsd fliplr(lags)/Fsd];
% Y = [l_corr_mec fliplr(u_corr_mec)];
% fill(X,Y,'b')
% Y = [l_corr_lec fliplr(u_corr_lec)];
% fill(X,Y,'g')
% Y = [l_corr_lfp fliplr(u_corr_lfp)];
% fill(X,Y,'r')

line([0 0],[-0.4 1],'Color','k')
line([-8 8],[0 0],'Color','k')
xlim([-8 8])
ylim([-0.4 1])
ylabel('Correlation Coefficient')
xlabel('Time (s)')

%% Supplementary figure for MP/LFP xcorr
clear all
close all
load C:\WC_Germany\persistent_revised\time_domain\time_domain_data
load C:\WC_Germany\persistent_revised\pers_revised_dir

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

first_point = find(lags/Fsd > -1,1,'first');

[all_peaks,all_peak_locs] = max(tot_w8_x(:,first_point:end),[],2);
peak_lags = lags(all_peak_locs+first_point)/Fsd;

m_corr_mec = mean(tot_w8_x(mec_cells,:));
u_corr_mec = m_corr_mec + 1*std(tot_w8_x(mec_cells,:))/sqrt(length(mec_cells));
l_corr_mec = m_corr_mec - 1*std(tot_w8_x(mec_cells,:))/sqrt(length(mec_cells));

[pk_mec,pk_mec_loc] = max(m_corr_mec);
corr_at_peak_mec = tot_w8_x(1:21,pk_mec_loc);

pk_mec_u = u_corr_mec(pk_mec_loc);
pk_mec_l = l_corr_mec(pk_mec_loc);
pk_shift_mec = lags(pk_mec_loc)/Fsd;

m_corr_lec = mean(tot_w8_x(lec_cells,:));
u_corr_lec = m_corr_lec + 1*std(tot_w8_x(lec_cells,:))/sqrt(length(lec_cells));
l_corr_lec = m_corr_lec - 1*std(tot_w8_x(lec_cells,:))/sqrt(length(lec_cells));

[pk_lec,pk_lec_loc] = max(m_corr_lec);
pk_lec_u = u_corr_lec(pk_lec_loc);
pk_lec_l = l_corr_lec(pk_lec_loc);
pk_shift_lec = lags(pk_lec_loc)/Fsd;
corr_at_peak_lec = tot_w8_x(lec_cells,pk_lec_loc);

plot(lags/Fsd,m_corr_mec,'linewidth',2)
hold on
plot(lags/Fsd,m_corr_lec,'g','linewidth',2)

X = [lags/Fsd fliplr(lags)/Fsd];
Y = [l_corr_mec fliplr(u_corr_mec)];
fill(X,Y,'b')
Y = [l_corr_lec fliplr(u_corr_lec)];
fill(X,Y,'g')

line([0 0],[-0.4 1],'Color','k')
line([-8 8],[0 0],'Color','k')
xlim([-4 4])
ylim([-0.4 0.8])
ylabel('Correlation Coefficient')
xlabel('Time (s)')


%% Up State Distribution Compare NO CUT
clear all
close all
load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data_no_cut

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

mean_up_mec = mean(up_hist(mec_cells,:));
u_up_mec = mean_up_mec + 1*std(up_hist(mec_cells,:))/sqrt(length(mec_cells));
l_up_mec = mean_up_mec - 1*std(up_hist(mec_cells,:))/sqrt(length(mec_cells));

mean_up_lfp = mean(up_hist8([mec_cells lec_cells],:));
u_up_lfp = mean_up_lfp + 1*std(up_hist8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));
l_up_lfp = mean_up_lfp - 1*std(up_hist8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));

mean_up_lec = mean(up_hist(lec_cells,:));
u_up_lec = mean_up_lec + 1*std(up_hist(lec_cells,:))/sqrt(length(lec_cells));
l_up_lec = mean_up_lec - 1*std(up_hist(lec_cells,:))/sqrt(length(lec_cells));

stairs(up_grid,mean_up_lfp,'r','linewidth',2)
hold on
stairs(up_grid,mean_up_mec,'linewidth',2)
stairs(up_grid,mean_up_lec,'g','linewidth',2)

legend('LFP','MEC','LEC')

epsilon = 1e-7;
u_up_lfp(u_up_lfp < epsilon) = epsilon;
l_up_lfp(l_up_lfp < epsilon) = epsilon;
mean_up_lfp(mean_up_lfp < epsilon) = epsilon;
u_up_mec(u_up_mec < epsilon) = epsilon;
l_up_mec(l_up_mec < epsilon) = epsilon;
mean_up_mec(mean_up_mec < epsilon) = epsilon;
u_up_lec(u_up_lec < epsilon) = epsilon;
l_up_lec(l_up_lec < epsilon) = epsilon;
mean_up_lec(mean_up_lec < epsilon) = epsilon;

% for i = 1:length(up_grid)-1
%     clear X Y
%     bin_width = (up_grid(i+1)-up_grid(i));
%     X = [up_grid(i) up_grid(i)+bin_width];
%     Y = [u_up_lec(i) u_up_lec(i)];
%     X = [X fliplr(X)];
%     Y = [Y l_up_lec(i) l_up_lec(i)];
%     fill(X,Y,'g','EdgeColor','none');
%     hold on
% end
% 
% for i = 1:length(up_grid)-1
%     clear X Y
%     bin_width = (up_grid(i+1)-up_grid(i));
%     X = [up_grid(i) up_grid(i)+bin_width];
%     Y = [u_up_lfp(i) u_up_lfp(i)];
%     X = [X fliplr(X)];
%     Y = [Y l_up_lfp(i) l_up_lfp(i)];
%     fill(X,Y,'r','EdgeColor','none');
%     hold on
% end
% 
% for i = 1:length(up_grid)-1
%     clear X Y
%     bin_width = (up_grid(i+1)-up_grid(i));
%     X = [up_grid(i) up_grid(i)+bin_width];
%     Y = [u_up_mec(i) u_up_mec(i)];
%     X = [X fliplr(X)];
%     Y = [Y l_up_mec(i) l_up_mec(i)];
%     fill(X,Y,'b','EdgeColor','none');
%     hold on
% end
% 

set(gca,'yscale','log')
set(gca,'YMinorGrid','off')
set(gca,'YGrid','off')
% set(gca,'XMinorGrid','off')
ylim([1e-4 0.2])
xlim([0 10])
line([0.3 0.3],[1e-4 0.2],'Color','k')
shg
xlabel('Up State Duration (s)','FontSize',14)
ylabel('Percent of Data','FontSize',14)

%% Down State Distribution Compare NO CUT
clear all
close all
load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data_no_cut

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

mean_down_mec = mean(down_hist(mec_cells,:));
u_down_mec = mean_down_mec + 1*std(down_hist(mec_cells,:))/sqrt(length(mec_cells));
l_down_mec = mean_down_mec - 1*std(down_hist(mec_cells,:))/sqrt(length(mec_cells));

mean_down_lfp = mean(down_hist8([mec_cells lec_cells],:));
u_down_lfp = mean_down_lfp + 1*std(down_hist8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));
l_down_lfp = mean_down_lfp - 1*std(down_hist8([mec_cells lec_cells],:))/sqrt(length(mec_cells)+length(lec_cells));

mean_down_lec = mean(down_hist(lec_cells,:));
u_down_lec = mean_down_lec + 1*std(down_hist(lec_cells,:))/sqrt(length(lec_cells));
l_down_lec = mean_down_lec - 1*std(down_hist(lec_cells,:))/sqrt(length(lec_cells));

stairs(down_grid,mean_down_lfp,'r','linewidth',2)
hold on
stairs(down_grid,mean_down_mec,'linewidth',2)
stairs(down_grid,mean_down_lec,'g','linewidth',2)

legend('LFP','MEC','LEC')

epsilon = 1e-7;
u_down_lfp(u_down_lfp < epsilon) = epsilon;
l_down_lfp(l_down_lfp < epsilon) = epsilon;
mean_down_lfp(mean_down_lfp < epsilon) = epsilon;
u_down_mec(u_down_mec < epsilon) = epsilon;
l_down_mec(l_down_mec < epsilon) = epsilon;
mean_down_mec(mean_down_mec < epsilon) = epsilon;
u_down_lec(u_down_lec < epsilon) = epsilon;
l_down_lec(l_down_lec < epsilon) = epsilon;
mean_down_lec(mean_down_lec < epsilon) = epsilon;

% for i = 1:length(down_grid)-1
%     clear X Y
%     bin_width = (down_grid(i+1)-down_grid(i));
%     X = [down_grid(i) down_grid(i)+bin_width];
%     Y = [u_down_lec(i) u_down_lec(i)];
%     X = [X fliplr(X)];
%     Y = [Y l_down_lec(i) l_down_lec(i)];
%     fill(X,Y,'g','EdgeColor','none');
%     hold on
% end
% 
% for i = 1:length(down_grid)-1
%     clear X Y
%     bin_width = (down_grid(i+1)-down_grid(i));
%     X = [down_grid(i) down_grid(i)+bin_width];
%     Y = [u_down_lfp(i) u_down_lfp(i)];
%     X = [X fliplr(X)];
%     Y = [Y l_down_lfp(i) l_down_lfp(i)];
%     fill(X,Y,'r','EdgeColor','none');
%     hold on
% end
% 
% for i = 1:length(down_grid)-1
%     clear X Y
%     bin_width = (down_grid(i+1)-down_grid(i));
%     X = [down_grid(i) down_grid(i)+bin_width];
%     Y = [u_down_mec(i) u_down_mec(i)];
%     X = [X fliplr(X)];
%     Y = [Y l_down_mec(i) l_down_mec(i)];
%     fill(X,Y,'b','EdgeColor','none');
%     hold on
% end


set(gca,'yscale','log')
set(gca,'YMinorGrid','off')
set(gca,'YGrid','off')
% set(gca,'XMinorGrid','off')
ylim([1e-4 0.2])
line([0.3 0.3],[1e-4 0.2],'Color','k')
xlim([0 10])
shg
xlabel('Down State Duration (s)','FontSize',14)
ylabel('Percent of Data','FontSize',14)

